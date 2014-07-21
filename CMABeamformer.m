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
        
        %Frequency Domain Constant Modulus Algorithm performed at baseband
        function [waveform, w] = oCMA_FrequencyDomain(obj, AntennaConfiguration, OFDM, array_waveform, mu)
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
            
            %remove zero+noise channels from symbol
            %these are present from interpolation factor of 4 in IFFT size
            %midpoint of channels
            mid = floor(OFDM.channels/2);
            %remove obj.N-obj.channels zeros(+noise) from symbol
            s(mid+1:OFDM.N-mid,:,:) = [];
            
            %normalize the power across time and antennas
            p = mean(mean(abs(s),2),3);
            for i = 1:OFDM.channels
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
                %[W; diag(1/AntennaConfiguration.N*ones(OFDM.channels,1))]
                W = [W; sparse(1:OFDM.channels, 1:OFDM.channels, 1/AntennaConfiguration.N*ones(OFDM.channels,1))];
            end
            
            %initial inverse correlation matrix
            %each channel of the OFDM signal has its own correlation matrix
            %correlation matrix is a square matrix sized to the number of
            %elements in the array
            R = repmat(diag(ones(AntennaConfiguration.N,1)), [1 1 OFDM.channels]);
            %create a stacked square matrix
            R = square_to_stacked(R);
            
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
                    S = [S; sparse(1:OFDM.channels,1:OFDM.channels,s(:,k,m))];
                end
                %output Y is the sum of the weights x OFDM symbol at each
                %antenna in column vector form
                y(:,k) = diag(W'*S);
                
                %save error to display learning curve
                e(:,k) = (ones(OFDM.channels,1) - y(:,k).*conj(y(:,k))).*y(:,k);
                
                %conj operator moved from weight update equation to improve
                %memory usage
                d = conj((diag(ones(OFDM.channels,1)) - diag(y(:,k).*conj(y(:,k))))*y(:,k));
                %Construct the tall block matrix D = diag(d.'; d.'...)
                D = [];
                for m = 1:AntennaConfiguration.N
                    D = [D; d];
                end
                %Use sparse matrix to significantly reduce processing time
                %sparse matrix is equivalent to diag(D);
                D = sparse(1:AntennaConfiguration.N*OFDM.channels,1:AntennaConfiguration.N*OFDM.channels,D);
                
                %Calculate the inverse correlation matrix
                %calculate scalar divisor, due to the stacked matrix format
                %this results in a OFDM.channelsxOFDM.channels diagonal matrix
                %scalar_to_stacked puts it in a format that the schur
                %product of the inverse yields the same result as matrix/scalar
                scalar_term = diag((1-a)*ones(OFDM.channels,1)) + a*S'*R*S;
                r = scalar_to_stacked(inv(scalar_term), AntennaConfiguration.N);
                %Inverse correlation matrix estimate
                %R = R/(1-a) - 1/(1-a)*(a*R*(S*S')*R)/r; move division to r
                R = R/(1-a) - 1/(1-a)*(a*R*(S*S')*R).*r;

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
%                 if mod(k,500) == 0
%                     figure;
%                     f = OFDM.CalculateCenterFrequencies;
%                     subplot(2,2,1);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,1), 3e8/f(1));
%                     title(sprintf('oCMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w(:,1)'*w(:,1))));
%                     subplot(2,2,2);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,16), 3e8/f(4));
%                     title(sprintf('oCMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w(:,16)'*w(:,16))));
%                     subplot(2,2,3);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,32), 3e8/f(5));
%                     title(sprintf('oCMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w(:,32)'*w(:,32))));
%                     subplot(2,2,4);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,96), 3e8/f(7));
%                     title(sprintf('oCMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w(:,96)'*w(:,96))));
%                 end
            end
            
            figure;
            f = OFDM.CalculateCenterFrequencies;
            subplot(2,2,1);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,1), 3e8/f(1));
            title(sprintf('oCMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w(:,1)'*w(:,1))));
            subplot(2,2,2);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,16), 3e8/f(4));
            title(sprintf('oCMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w(:,16)'*w(:,16))));
            subplot(2,2,3);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,32), 3e8/f(5));
            title(sprintf('oCMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w(:,32)'*w(:,32))));
            subplot(2,2,4);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,96), 3e8/f(7));
            title(sprintf('oCMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w(:,96)'*w(:,96))));
            
            %apply final set of weights to OFDM frequency domain symbols
            waveform = zeros(64,l);
            for i = 1:OFDM.channels
                waveform(i,:) = w(:,i)'*squeeze(s(i,:,:)).';
            end
            
            figure;
            plot(mean(abs(e),1));
            title('oCMA Learning Curve');
            xlabel('OFDM Symbols')
            ylabel('mean |error|');
        end
        
        function [waveform, w] = LCoCMA_FrequencyDomain(obj, AntennaConfiguration, OFDM, array_waveform, mu, angle)
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
            
            %remove zero+noise channels from symbol
            %these are present from interpolation factor of 4 in IFFT size
            %midpoint of channels
            mid = floor(OFDM.channels/2);
            %remove obj.N-obj.channels zeros(+noise) from symbol
            s(mid+1:OFDM.N-mid,:,:) = [];
            
            %normalize the power across time and antennas
            p = mean(mean(abs(s),2),3);
            for i = 1:OFDM.channels
                s(i,:,:) = s(i,:,:)/p(i);
            end
            
            %generate the constraint vector for the specified angle
            C = AntennaConfiguration.CalculateNormalizedReplicaVector(angle);
            
            %generate the blocking matrix for a generalized sidelobe canceller using
            %the constraint C
            %projection matrix onto constraint subspace C
            Pc = C*(C'*C)^-1*C';
            Pc_orth = diag(ones(AntennaConfiguration.N,1)) - Pc;

            %find the orthonormalization of Pc_orth
            [Q,~] = qr(Pc_orth);
            %Blocking matrix B is the first N-1 columns of Q
            b = Q(:,1:AntennaConfiguration.N-1);
            %create stacked blocking matrix
            B = blocking_to_stacked(repmat(b, [1 1 OFDM.channels]));

            %weights for constraint vector
            w_c = C*(C'*C)^-1;
            
            %number of OFDM symbols
            [~,l,~] = size(s);
            
            %adaptively beamform the array input
            %constraint weights
            W_c = [];
            %initial adaptive weights
            W = [];
            %W_c is a tall block matrix 
            for m = 1:AntennaConfiguration.N
                %Use sparse matrices to significantly reduce processing time
                W_c = [W_c; sparse(1:OFDM.channels, 1:OFDM.channels, w_c(m)*ones(OFDM.channels,1))];
            end
            %W is a tall block matrix [diag(w1); diag(w2)...] where w1,
            %w2... are the weights for antennas 1,2...
            for m = 1:AntennaConfiguration.N-1
                %Use sparse matrices to significantly reduce processing
                %time.  Sparse matrix is equivalent to 
                %[W; diag(1/AntennaConfiguration.N*ones(OFDM.channels,1))]
                W = [W; sparse(1:OFDM.channels, 1:OFDM.channels, 1/(AntennaConfiguration.N-1)*ones(OFDM.channels,1))];
            end
            
            %initial inverse correlation matrix
            %each channel of the OFDM signal has its own correlation matrix
            %correlation matrix is a square matrix sized to the number of
            %elements in the array
            R = repmat(diag(ones(AntennaConfiguration.N-1,1)), [1 1 OFDM.channels]);
            %create a stacked square matrix
            R = square_to_stacked(R);
            
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
                    S = [S; sparse(1:OFDM.channels,1:OFDM.channels,s(:,k,m))];
                end
                
                %Use the blocking matrix to remove the constraint subspace
                %from the input vectors
                S_reduced = B'*S;
                
                %output Y is the sum of the weights x OFDM symbol at each
                %antenna in column vector form
                y(:,k) = diag(W_c'*S - W'*S_reduced);
                
                %save error to display learning curve
                e(:,k) = (ones(OFDM.channels,1) - y(:,k).*conj(y(:,k))).*y(:,k);

                %conj operator moved from weight update equation to improve
                %memory usage
                d = conj((diag(ones(OFDM.channels,1)) - diag(y(:,k).*conj(y(:,k))))*y(:,k));
                %d = ((diag(y(:,k).*conj(y(:,k))) - diag(ones(OFDM.channels,1)))*conj(y(:,k)));
                %Construct the tall block matrix D = diag(d.'; d.'...)
                D = [];
                for m = 1:AntennaConfiguration.N-1
                    D = [D; d];
                end
                %Use sparse matrix to significantly reduce processing time
                %sparse matrix is equivalent to diag(D);
                D = sparse(1:(AntennaConfiguration.N-1)*OFDM.channels,1:(AntennaConfiguration.N-1)*OFDM.channels,D);
                
                %Calculate the inverse correlation matrix
                %calculate scalar divisor, due to the stacked matrix format
                %this results in a OFDM.channelsxOFDM.channels diagonal matrix
                %scalar_to_stacked puts it in a format that the schur
                %product of the inverse yields the same result as matrix/scalar
                scalar_term = diag((1-a)*ones(OFDM.channels,1)) + a*S_reduced'*R*S_reduced;
                r = scalar_to_stacked(inv(scalar_term), AntennaConfiguration.N-1);
                %Inverse correlation matrix estimate
                %R = R/(1-a) - 1/(1-a)*(a*R*(S*S')*R)/r; move division to r
                R = R/(1-a) - 1/(1-a)*(a*R*(S_reduced*S_reduced')*R).*r;
                
                %weight update equation
                %W = W + mu*R*conj(D)*S_reduced; move conj to d to improve memory usage
                W = W + mu*R*D*S_reduced;
                
                %weights extracted from matrix form
                w = obj.ExtractReducedWeights(AntennaConfiguration, OFDM, W);
                
                %save one set of weights
                v(:,k) = w_c - b*w(:,5);

                if sum(isnan(abs(v(:,k)))) > 0
                    stop = 1;
                end
                %Plot every 500th beam pattern
%                 if mod(k,500) == 0
%                     figure;
%                     f = OFDM.CalculateCenterFrequencies;
%                     subplot(2,2,1);
%                     w_gsc =  w_c - b*w(:,1)/sum(abs(w(:,1)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(1));
%                     title(sprintf('LCoCMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                     subplot(2,2,2);
%                     w_gsc =  w_c - b*w(:,16)/sum(abs(w(:,16)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(4));
%                     title(sprintf('LCoCMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                     subplot(2,2,3);
%                     w_gsc =  w_c - b*w(:,32)/sum(abs(w(:,32)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(5));
%                     title(sprintf('LCoCMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                     subplot(2,2,4);
%                     w_gsc =  w_c - b*w(:,96)/sum(abs(w(:,96)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(7));
%                     title(sprintf('LCoCMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                 end
            end
            
            figure;
            f = OFDM.CalculateCenterFrequencies;
            subplot(2,2,1);
            w_gsc =  w_c - b*w(:,1)/sum(abs(w(:,1)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(1));
            title(sprintf('LCoCMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            subplot(2,2,2);
            w_gsc =  w_c - b*w(:,16)/sum(abs(w(:,16)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(4));
            title(sprintf('LCoCMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            subplot(2,2,3);
            w_gsc =  w_c - b*w(:,32)/sum(abs(w(:,32)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(5));
            title(sprintf('LCoCMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            subplot(2,2,4);
            w_gsc =  w_c - b*w(:,96)/sum(abs(w(:,96)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(7));
            title(sprintf('LCoCMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            
            %apply final set of weights to OFDM frequency domain symbols
            waveform = zeros(64,l);
            for i = 1:OFDM.channels
                w_gsc = w_c - b*w(:,i);
                waveform(i,:) = w_gsc'*squeeze(s(i,:,:)).';
            end
            
            figure;
            plot(mean(abs(e),1));
            title('LCoCMA Learning Curve');
            xlabel('OFDM Symbols')
            ylabel('mean |error|');
        end
        
        %Frequency Domain Constant Modulus Algorithm performed at baseband
        function [waveform, w] = CMA_FrequencyDomain(obj, AntennaConfiguration, OFDM, array_waveform, mu)
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
            
            %remove zero+noise channels from symbol
            %these are present from interpolation factor of 4 in IFFT size
            %midpoint of channels
            mid = floor(OFDM.channels/2);
            %remove obj.N-obj.channels zeros(+noise) from symbol
            s(mid+1:OFDM.N-mid,:,:) = [];
            
            %normalize the power across time and antennas
            p = mean(mean(abs(s),2),3);
            for i = 1:OFDM.channels
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
                %[W; diag(1/AntennaConfiguration.N*ones(OFDM.channels,1))]
                W = [W; sparse(1:OFDM.channels, 1:OFDM.channels, 1/AntennaConfiguration.N*ones(OFDM.channels,1))];
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
                    S = [S; sparse(1:OFDM.channels,1:OFDM.channels,s(:,k,m))];
                end
                %output Y is the sum of the weights x OFDM symbol at each
                %antenna in column vector form
                y(:,k) = diag(W'*S);
                
                %save error to display learning curve
                e(:,k) = (ones(OFDM.channels,1) - y(:,k).*conj(y(:,k))).*y(:,k);
                
                %conj operator moved from weight update equation to improve
                %memory usage
                d = conj((diag(ones(OFDM.channels,1)) - diag(y(:,k).*conj(y(:,k))))*y(:,k));
                %Construct the tall block matrix D = diag(d.'; d.'...)
                D = [];
                for m = 1:AntennaConfiguration.N
                    D = [D; d];
                end
                %Use sparse matrix to significantly reduce processing time
                %sparse matrix is equivalent to diag(D);
                D = sparse(1:AntennaConfiguration.N*OFDM.channels,1:AntennaConfiguration.N*OFDM.channels,D);
                
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
%                 %Plot every 500th beam pattern
%                 if mod(k,500) == 0
%                     figure;
%                     f = OFDM.CalculateCenterFrequencies;
%                     subplot(2,2,1);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,1), 3e8/f(1));
%                     title(sprintf('CMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w(:,1)'*w(:,1))));
%                     subplot(2,2,2);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,16), 3e8/f(4));
%                     title(sprintf('CMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w(:,16)'*w(:,16))));
%                     subplot(2,2,3);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,32), 3e8/f(5));
%                     title(sprintf('CMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w(:,32)'*w(:,32))));
%                     subplot(2,2,4);
%                     AntennaConfiguration.PlotBeampattern2D(0, w(:,96), 3e8/f(7));
%                     title(sprintf('CMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w(:,96)'*w(:,96))));
%                 end
            end
            
            figure;
            f = OFDM.CalculateCenterFrequencies;
            subplot(2,2,1);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,1), 3e8/f(1));
            title(sprintf('CMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w(:,1)'*w(:,1))));
            subplot(2,2,2);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,16), 3e8/f(4));
            title(sprintf('CMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w(:,16)'*w(:,16))));
            subplot(2,2,3);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,32), 3e8/f(5));
            title(sprintf('CMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w(:,32)'*w(:,32))));
            subplot(2,2,4);
            AntennaConfiguration.PlotBeampattern2D(0, w(:,96), 3e8/f(7));
            title(sprintf('CMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w(:,96)'*w(:,96))));
            
            %apply final set of weights to OFDM frequency domain symbols
            waveform = zeros(OFDM.channels,l);
            for i = 1:OFDM.channels
                waveform(i,:) = w(:,i)'*squeeze(s(i,:,:)).';
            end
            
            figure;
            plot(mean(abs(e),1));
            title('CMA Learning Curve');
            xlabel('OFDM Symbols')
            ylabel('mean |error|');
        end
        
        function [waveform, w] = LCCMA_FrequencyDomain(obj, AntennaConfiguration, OFDM, array_waveform, mu, angle)
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
            
            %remove zero+noise channels from symbol
            %these are present from interpolation factor of 4 in IFFT size
            %midpoint of channels
            mid = floor(OFDM.channels/2);
            %remove obj.N-obj.channels zeros(+noise) from symbol
            s(mid+1:OFDM.N-mid,:,:) = [];
            
            %normalize the power across time and antennas
            p = mean(mean(abs(s),2),3);
            for i = 1:OFDM.channels
                s(i,:,:) = s(i,:,:)/p(i);
            end
            
            %generate the constraint vector for the specified angle
            C = AntennaConfiguration.CalculateNormalizedReplicaVector(angle);
            
            %generate the blocking matrix for a generalized sidelobe canceller using
            %the constraint C
            %projection matrix onto constraint subspace C
            Pc = C*(C'*C)^-1*C';
            Pc_orth = diag(ones(AntennaConfiguration.N,1)) - Pc;

            %find the orthonormalization of Pc_orth
            [Q,~] = qr(Pc_orth);
            %Blocking matrix B is the first N-1 columns of Q
            b = Q(:,1:AntennaConfiguration.N-1);
            %create stacked blocking matrix
            B = blocking_to_stacked(repmat(b, [1 1 OFDM.channels]));

            %weights for constraint vector
            w_c = C*(C'*C)^-1;
            
            %number of OFDM symbols
            [~,l,~] = size(s);
            
            %adaptively beamform the array input
            %constraint weights
            W_c = [];
            %initial adaptive weights
            W = [];
            %W_c is a tall block matrix 
            for m = 1:AntennaConfiguration.N
                %Use sparse matrices to significantly reduce processing time
                W_c = [W_c; sparse(1:OFDM.channels, 1:OFDM.channels, w_c(m)*ones(OFDM.channels,1))];
            end
            %W is a tall block matrix [diag(w1); diag(w2)...] where w1,
            %w2... are the weights for antennas 1,2...
            for m = 1:AntennaConfiguration.N-1
                %Use sparse matrices to significantly reduce processing
                %time.  Sparse matrix is equivalent to 
                %[W; diag(1/AntennaConfiguration.N*ones(OFDM.channels,1))]
                W = [W; sparse(1:OFDM.channels, 1:OFDM.channels, 1/(AntennaConfiguration.N-1)*ones(OFDM.channels,1))];
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
                    S = [S; sparse(1:OFDM.channels,1:OFDM.channels,s(:,k,m))];
                end
                
                %Use the blocking matrix to remove the constraint subspace
                %from the input vectors
                S_reduced = B'*S;
                
                %output Y is the sum of the weights x OFDM symbol at each
                %antenna in column vector form
                y(:,k) = diag(W_c'*S - W'*S_reduced);
                
                %save error to display learning curve
                e(:,k) = (ones(OFDM.channels,1) - y(:,k).*conj(y(:,k))).*y(:,k);
                
                %conj operator moved from weight update equation to improve
                %memory usage
                d = conj((diag(ones(OFDM.channels,1)) - diag(y(:,k).*conj(y(:,k))))*y(:,k));
                %d = ((diag(y(:,k).*conj(y(:,k))) - diag(ones(OFDM.channels,1)))*conj(y(:,k)));
                %Construct the tall block matrix D = diag(d.'; d.'...)
                D = [];
                for m = 1:AntennaConfiguration.N-1
                    D = [D; d];
                end
                %Use sparse matrix to significantly reduce processing time
                %sparse matrix is equivalent to diag(D);
                D = sparse(1:(AntennaConfiguration.N-1)*OFDM.channels,1:(AntennaConfiguration.N-1)*OFDM.channels,D);
                
                %weight update equation
                %W = W + mu*conj(D)*S_reduced; move conj to d to improve memory usage
                W = W + mu*D*S_reduced;
                
                %weights extracted from matrix form
                w = obj.ExtractReducedWeights(AntennaConfiguration, OFDM, W);
                %save one set of weights
                v(:,k) = w_c - b*w(:,5);
                
                if sum(isnan(abs(v(:,k)))) > 0
                    stop = 1;
                end
%                 %Plot every 500th beam pattern
%                 if mod(k,500) == 0
%                     figure;
%                     f = OFDM.CalculateCenterFrequencies;
%                     subplot(2,2,1);
%                     w_gsc =  w_c - b*w(:,1)/sum(abs(w(:,1)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(1));
%                     title(sprintf('LCCMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                     subplot(2,2,2);
%                     w_gsc =  w_c - b*w(:,16)/sum(abs(w(:,16)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(4));
%                     title(sprintf('LCCMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                     subplot(2,2,3);
%                     w_gsc =  w_c - b*w(:,32)/sum(abs(w(:,32)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(5));
%                     title(sprintf('LCCMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                     subplot(2,2,4);
%                     w_gsc =  w_c - b*w(:,96)/sum(abs(w(:,96)));
%                     AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(7));
%                     title(sprintf('LCCMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
%                 end
            end
            
            figure;
            f = OFDM.CalculateCenterFrequencies;
            subplot(2,2,1);
            w_gsc =  w_c - b*w(:,1)/sum(abs(w(:,1)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(1));
            title(sprintf('LCCMA Channel 1 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            subplot(2,2,2);
            w_gsc =  w_c - b*w(:,16)/sum(abs(w(:,16)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(4));
            title(sprintf('LCCMA Channel 16 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            subplot(2,2,3);
            w_gsc =  w_c - b*w(:,32)/sum(abs(w(:,32)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(5));
            title(sprintf('LCCMA Channel 32 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            subplot(2,2,4);
            w_gsc =  w_c - b*w(:,96)/sum(abs(w(:,96)));
            AntennaConfiguration.PlotBeampattern2D(0, w_gsc, 3e8/f(7));
            title(sprintf('LCCMA Channel 96 Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
            
            %apply final set of weights to OFDM frequency domain symbols
            waveform = zeros(64,l);
            for i = 1:OFDM.channels
                w_gsc = w_c - b*w(:,i);
                waveform(i,:) = w_gsc'*squeeze(s(i,:,:)).';
            end
            
            figure;
            plot(mean(abs(e),1));
            title('LCCMA Learning Curve');
            xlabel('OFDM Symbols')
            ylabel('mean |error|');
        end
        
        %The weight matrix W is a sparse
        %OFDM.channels*AntennaConfiguration.NxOFDM.channels matrix
        function w = ExtractWeights(obj, AntennaConfiguration, OFDM, W)
            %there is a set of weights for each sample of the OFDM symbol
            for i = 1:OFDM.channels
                w(:,i) = W(i:OFDM.channels:OFDM.channels*AntennaConfiguration.N,i);
            end
        end
        
        %The weight matrix W is a sparse
        %OFDM.channels*AntennaConfiguration.NxOFDM.channels matrix
        function w = ExtractReducedWeights(obj, AntennaConfiguration, OFDM, W)
            %there is a set of weights for each sample of the OFDM symbol
            for i = 1:OFDM.channels
                w(:,i) = W(i:OFDM.channels:OFDM.channels*(AntennaConfiguration.N-1),i);
            end
        end
    end
end