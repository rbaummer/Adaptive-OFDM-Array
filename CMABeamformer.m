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
            N = OFDM.N*OFDM.L;
            
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
        function output = TimeDomainCMA(obj, AntennaConfiguration, OFDM, array_waveform, mu)
            %check upsampled FFT size size matches F/G size
            assert(OFDM.N*OFDM.L == length(obj.F), 'F&G must be NxN');
            
            %Define permutation matrix
            P = obj.G*conj(obj.F);
            
            %Determine OFDM symbol locations, organize into time domain
            %OFDM symbols at each antenna and remove cyclic prefix
            r = zeros(OFDM.N*OFDM.L, OFDM.num_symbols - 3, AntennaConfiguration.N);
            for i = 1:AntennaConfiguration.N
                %time domain passband OFDM symbols
                received_symbol = OFDM.serialtoparallel(array_waveform(i,:));
                %save number of OFDM symbols - 3 (minimum length after
                %synchronization)
                r(:,:,i) = received_symbol(:, 1:OFDM.num_symbols - 3);
            end
            
            %number of OFDM symbols
            [~,~,l] = size(r);
            
            %adaptively beamform the array input
            %identity matrix the same size as upampled OFDM symbol
            I = diag(ones(OFDM.N*OFDM.L,1));
            v = [ones(1,AntennaConfiguration.N); zeros(OFDM.N*OFDM.L-1,AntennaConfiguration.N)];
            for k = 1:l
                %multiply ofdm symbol by weights at each antenna
                for m = 1:AntennaConfiguration.N
                    R(:,m) = conj(circ(v(:,m)))*r(:,k,m);
                end
                
                %sum OFDM symbol from each antenna
                y_bar = sum(R,2);
                
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

        function FrequencyDomainCMA(obj, AntennaConfiguration, OFDM, array_waveform, mu)
            %Determine OFDM symbol locations, organize into time domain
            %OFDM symbols at each antenna and remove cyclic prefix
            s = zeros(OFDM.N*OFDM.L, OFDM.num_symbols - 3, AntennaConfiguration.N);
            for i = 1:AntennaConfiguration.N
                %time domain passband OFDM symbols
                r = OFDM.serialtoparallel(array_waveform(i,:));
                %frequency domain passband OFDM symbols
                received_symbol = fft(r);
                %save number of OFDM symbols - 3 (minimum length after
                %synchronization)
                s(:,:,i) = received_symbol(:, 1:OFDM.num_symbols - 3);
            end
            
            %number of OFDM symbols
            [~,~,l] = size(r);
            
            %adaptively beamform the array input
            %identity matrix the same size as upampled OFDM symbol
            I = diag(ones(OFDM.N*OFDM.L,1));
            v = [ones(1,AntennaConfiguration.N); zeros(OFDM.N*OFDM.L-1,AntennaConfiguration.N)];
            for k = 1:l
                %multiply ofdm symbol by weights at each antenna
                for m = 1:AntennaConfiguration.N
                    R(:,m) = conj(circ(v(:,m)))*r(:,k,m);
                end
                
                %sum OFDM symbol from each antenna
                y_bar = sum(R,2);
                
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
    end
end