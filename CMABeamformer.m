%Author: Rob Baummer
%Algorithm Source: Vishwanath Venkataraman et al.  Adaptive Algorithms for 
%OFDMA Wireless Ad Hoc Networks with Multiple Antennas.  
%
%CMABeamformer uses the Constant Modulus Algorithm to beamform an OFDM
%signal arriving at an Antenna Array
%
%properties
%   F: FFT Matrix Fi,j = exp(-1i*2*pi*(i-1)*(j-1)/N)
%   G: FFT Matrix Gi,j = exp(1i*2*pi*(i-1)*(j-1)/N)/N
%default constructor: CMABeamformer(OFDM_inst) takes OFDM object and
%   creates matrices F&G
classdef CMABeamformer
    properties
        %FFT Matrix
        F = 0;
        %IFFT Matrix
        G = 0;
    end
    
    methods
        %constructor
        function obj = CMABeamformer(OFDM)
            N = OFDM.N;
            
            %Create FFT matrix of size NxN
            obj.F = exp(-1i*2*pi*(0:N-1)'*(0:N-1)/N);
            %Create IFFT matrix of size NxN
            obj.G = exp(1i*2*pi*(0:N-1)'*(0:N-1)/N)/N;
        end
        
        %Time-Domain CMA Beamformer
        %AntennaConfiguration contains the antenna array object
        %OFDM contains the OFDM object
        %array waveform is in the time domain at passband and contains 
        %the time domain samples at each antenna NxSamples
        %mu is the step size of the CMA algorithm
        function output = CMA_TimeDomain(obj, AntennaConfiguration, OFDM, array_waveform, mu)
            %check upsampled FFT size size matches F/G size
            assert(OFDM.N == length(obj.F), 'F&G must be NxN');
            
            %Define permutation matrix
            %a permutation matrix contains a 1 in every row and column,
            %zeros elsewhere.  Round is necessary due to numerical
            %representation leaving very small numbers (1e-17) but not exactly 0
            P = round(obj.G*conj(obj.F));
            
            %Determine OFDM symbol locations, organize into time domain
            %OFDM symbols at each antenna and remove cyclic prefix
            r = zeros(OFDM.N, OFDM.num_symbols - 3, AntennaConfiguration.N);
            for i = 1:AntennaConfiguration.N
                %time domain passband OFDM symbols
                received_symbol = OFDM.serialtoparallel(array_waveform(i,:));
                %save number of OFDM symbols - 4 (minimum length after
                %synchronization)
                r(:,:,i) = received_symbol(:, 1:OFDM.num_symbols - 4);
            end
            
            %number of OFDM symbols
            [~,l,~] = size(r);
            
            %adaptively beamform the array input
            %identity matrix the same size as upampled OFDM symbol
            I = diag(ones(OFDM.N,1));
            v = [ones(1,AntennaConfiguration.N); zeros(OFDM.N-1,AntennaConfiguration.N)];
            for k = 1:l
                %multiply ofdm symbol by weights at each antenna
                for m = 1:AntennaConfiguration.N
                    y_bar(:,m) = conj(circ(v(:,m)))*r(:,k,m);
                end
                
                %sum OFDM symbol from each antenna
                y_bar = 1/AntennaConfiguration.N*sum(y_bar,2);
                
                %Calculate weight update equation
                A = circ(P*conj(y_bar));
                D = I - A*circ(y_bar)*A;
                for m = 1:AntennaConfiguration.N
                    %r is the time domain ofdm symbol at antenna m
                    %update the beamformer weights v from previous
                    %beamformer weights
                    v(:,m) = v(:,m) + mu*D*r(:,k,m);
                end
            end
        end
       
        %Frequency Domain Constant Modulus Algorithm performed at baseband
        function w = oCMA_FrequencyDomain(obj, AntennaConfiguration, OFDM, array_waveform, mu)
            %Determine OFDM symbol locations, organize into time domain
            %OFDM symbols at each antenna and remove cyclic prefix
            s = zeros(OFDM.N, OFDM.num_symbols - 4, AntennaConfiguration.N);
            for i = 1:AntennaConfiguration.N
                %time domain passband OFDM symbols
                r = OFDM.serialtoparallel(array_waveform(i,:));
                %frequency domain baseband OFDM symbols
                received_symbol = 1/OFDM.N*fft(r);
                %save number of OFDM symbols - 4 (minimum length after
                %synchronization)
                s(:,:,i) = received_symbol(:, 1:OFDM.num_symbols - 4);
            end
            
            %normalize the power across time and antennas
            p = mean(mean(abs(s),2),3);
            for i = 1:OFDM.N
                s(i,:,:) = s(i,:,:)/p(i);
            end
            
            %number of OFDM symbols
            [~,l,~] = size(s);
            
            %adaptively beamform the array input
            %initial weights
            W = [];
            %W is a tall block matrix [diag(w1); diag(w2)...] where w1,
            %w2... are the weights for antennas 1,2...
            for m = 1:AntennaConfiguration.N
                %Use sparse matrices to significantly reduce processing
                %time.  Sparse matrix is equivalent to 
                %[W; diag(1/AntennaConfiguration.N*ones(OFDM.N,1))]
                W = [W; sparse(1:OFDM.N, 1:OFDM.N, 1/AntennaConfiguration.N*ones(OFDM.N,1))];
            end
            
            %initial inverse correlation matrix
            %each channel of the OFDM signal has its own correlation matrix
            %correlation matrix is a square matrix sized to the number of
            %elements in the array
            R = repmat(diag(ones(AntennaConfiguration.N,1)), [1 1 OFDM.N]);
            %create a stacked square matrix
            R = square_to_stacked(R);
            
            %test R and weight gen
            R_test = diag(ones(AntennaConfiguration.N,1));
            w_test = 1/AntennaConfiguration.N*ones(AntennaConfiguration.N,1);
            
            %set forgetting factor
            a = 1-0.985;
            
            %Weight update loop
            for k = 1:l
                %Construct the tall block matrix of the frequency domain
                %symbol.  [diag(S1); diag(S2)...] where S1,S2 are the
                %frequency domain OFDM sybmols for antennas 1,2...
                S = [];
                for m = 1:AntennaConfiguration.N
                    %Use sparse matrices to significantly reduce processing
                    %time.  Sparse matrix is equivalent to [S; diag(s(:,k,m))]
                    S = [S; sparse(1:OFDM.N,1:OFDM.N,s(:,k,m))];
                end
                %output Y is the sum of the weights x OFDM symbol at each
                %antenna in column vector form
                y(:,k) = diag(W'*S);
                
                %conj operator moved from weight update equation to improve
                %memory usage
                d = conj((diag(ones(OFDM.N,1)) - diag(y(:,k).*conj(y(:,k))))*y(:,k));
                %Construct the tall block matrix D = diag(d.'; d.'...)
                D = [];
                for m = 1:AntennaConfiguration.N
                    D = [D; d];
                end
                %Use sparse matrix to significantly reduce processing time
                %sparse matrix is equivalent to diag(D);
                D = sparse(1:AntennaConfiguration.N*OFDM.N,1:AntennaConfiguration.N*OFDM.N,D);
                
                %Calculate the inverse correlation matrix
                %calculate scalar divisor, due to the stacked matrix format
                %this results in a OFDM.NxOFDM.N diagonal matrix
                %scalar_to_stacked puts it in a format that the schur
                %product of the inverse yields the same result as matrix/scalar
                scalar_term = diag((1-a)*ones(OFDM.N,1)) + a*S'*R*S;
                r = scalar_to_stacked(inv(scalar_term), AntennaConfiguration.N);
                %Inverse correlation matrix estimate
                %R = R/(1-a) - 1/(1-a)*(a*R*(S*S')*R)/r; move division to r
                R = R/(1-a) - 1/(1-a)*(a*R*(S*S')*R).*r;
                
                %Test inverse correlation matrix stacked format by
                %calculating for a single channel
                x = squeeze(s(64,k,:));
                c = 1-a + a*x'*R_test*x;
                R_test = R_test/(1-a) - 1/(1-a)*(a*R_test*(x*x')*R_test)/c;
                
                %Test weight update for a single channel
                y_test = w_test'*x;
                err(k) = y_test - y_test*abs(y_test)^2;
                w_test = w_test + mu*R_test*x*conj(err(k));
                
                %weight update equation
                %W = W + mu*conj(D)*S; move conj to d to improve memory usage
                W = W + mu*R*D*S;
                
                %weights extracted from matrix form
                w = obj.ExtractWeights(AntennaConfiguration, OFDM, W);
                %save one set of weights
                v(:,k) = w(:,1);
                
                if sum(isnan(abs(v(:,k)))) > 0
                    stop = 1;
                end
                %Plot every 5th beam pattern
                if mod(k,500) == 0
                    figure;
                    f = OFDM.CalculateCenterFrequencies;
                    subplot(2,2,1);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,64), 3e8/f(1));
                    title('64');
                    subplot(2,2,2);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,60), 3e8/f(4));
                    title('60');
                    subplot(2,2,3);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,59), 3e8/f(5));
                    title('59');
                    subplot(2,2,4);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,55), 3e8/f(7));
                    title('55');
                end
            end
        end
        
        %Frequency Domain Constant Modulus Algorithm performed at baseband
        %different stacked array format to support stacked square matrices
        function w = CMA_FrequencyDomain(obj, AntennaConfiguration, OFDM, array_waveform, mu)
            %Determine OFDM symbol locations, organize into time domain
            %OFDM symbols at each antenna and remove cyclic prefix
            s = zeros(OFDM.N, OFDM.num_symbols - 4, AntennaConfiguration.N);
            for i = 1:AntennaConfiguration.N
                %time domain passband OFDM symbols
                r = OFDM.serialtoparallel(array_waveform(i,:));
                %frequency domain passband OFDM symbols
                received_symbol = 1/OFDM.N*fft(r);
                %save number of OFDM symbols - 4 (minimum length after
                %synchronization)
                s(:,:,i) = received_symbol(:, 1:OFDM.num_symbols - 4);
            end
            
            %normalize the power across time and antennas
            p = mean(mean(abs(s),2),3);
            for i = 1:OFDM.N
                s(i,:,:) = s(i,:,:)/p(i);
            end
            
            %number of OFDM symbols
            [~,l,~] = size(s);
            
            %adaptively beamform the array input
            %initial weights
            W = [];
            %W is a tall block matrix [diag(w1); diag(w2)...] where w1,
            %w2... are the weights for antennas 1,2...
            for m = 1:AntennaConfiguration.N
                %Use sparse matrices to significantly reduce processing
                %time.  Sparse matrix is equivalent to 
                %[W; diag(1/AntennaConfiguration.N*ones(OFDM.N,1))]
                W = [W; sparse(1:OFDM.N, 1:OFDM.N, 1/AntennaConfiguration.N*ones(OFDM.N,1))];
            end
            %Weight update loop
            for k = 1:l
                %Construct the tall block matrix of the frequency domain
                %symbol.  [diag(S1); diag(S2)...] where S1,S2 are the
                %frequency domain OFDM sybmols for antennas 1,2...
                S = [];
                for m = 1:AntennaConfiguration.N
                    %Use sparse matrices to significantly reduce processing
                    %time.  Sparse matrix is equivalent to [S; diag(s(:,k,m))]
                    S = [S; sparse(1:OFDM.N,1:OFDM.N,s(:,k,m))];
                end
                %output Y is the sum of the weights x OFDM symbol at each
                %antenna in column vector form
                y(:,k) = diag(W'*S);
                
                %conj operator moved from weight update equation to improve
                %memory usage
                d = conj((diag(ones(OFDM.N,1)) - diag(y(:,k).*conj(y(:,k))))*y(:,k));
                %Construct the tall block matrix D = diag(d.'; d.'...)
                D = [];
                for m = 1:AntennaConfiguration.N
                    D = [D; d];
                end
                %Use sparse matrix to significantly reduce processing time
                %sparse matrix is equivalent to diag(D);
                D = sparse(1:AntennaConfiguration.N*OFDM.N,1:AntennaConfiguration.N*OFDM.N,D);
                
                %weight update equation
                %W = W + mu*conj(D)*S; move conj to d to improve memory usage
                W = W + mu*D*S;
                
                %weights extracted from matrix form
                w = obj.ExtractWeights(AntennaConfiguration, OFDM, W);
                %save one set of weights
                v(:,k) = w(:,5);
                
                if sum(isnan(abs(v(:,k)))) > 0
                    stop = 1;
                end
                %Plot every 500th beam pattern
                if mod(k,500) == 0
                    figure;
                    f = OFDM.CalculateCenterFrequencies;
                    subplot(2,2,1);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,1), 3e8/f(1));
                    title('1');
                    subplot(2,2,2);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,4), 3e8/f(4));
                    title('4');
                    subplot(2,2,3);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,5), 3e8/f(5));
                    title('5');
                    subplot(2,2,4);
                    AntennaConfiguration.PlotBeampattern2D(0, w(:,7), 3e8/f(7));
                    title('7');
                end
            end
        end
        
        %The weight matrix W is a sparse
        %OFDM.N*AntennaConfiguration.NxOFDM.N matrix
        function w = ExtractWeights(obj, AntennaConfiguration, OFDM, W)
            %there is a set of weights for each sample of the OFDM symbol
            for i = 1:OFDM.N
                w(:,i) = W(i:OFDM.N:OFDM.N*AntennaConfiguration.N,i);
            end
        end
    end
end