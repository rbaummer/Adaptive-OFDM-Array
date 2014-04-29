%Author: Rob Baummer
%
classdef AntennaConfiguration
    
    properties
        Configuration
        N;
        M;
        fd;
        G;
        
    end
    
    properties(Constant = true, Hidden = true)
        c = 3e8;
    end
    
    methods
        %constructor
        function obj = AntennaConfiguration(Configuration, N, M, fd)
            obj.Configuration = Configuration;
            obj.N = N;
            obj.M = M;
            %designed frequency of array; lambda = c/fd
            obj.fd = fd;
            if strcmp(Configuration, 'Linear') == true
                %set distance between antennas in x plane to lambda/2
                dx = obj.c/(2*fd);
                dy = 0;
                dz = 0;
                
                %matrix containing the position of each antenna
                obj.G = ((0:N-1) - (N-1)/2)'*[dx dy dz];
            elseif strcmp(Configuration, 'Rectangular') == true
                %set distance between antennas in x plane to lambda/2
                dx = obj.c/(2*fd);
                %set distance between antennas in y plane to lambda/2
                dy = obj.c/(2*fd);
                dz = 0;
                
                %matrix containing the position of each antenna
                for i = 0:N-1
                    for j = 0:M-1
                        obj.G(i*M+j+1, :) = [(i-(N-1)/2)*dx (j-(M-1)/2)*dy dz];
                    end
                end
            elseif strcmp(Configuration, 'Circular') == true
                %calculate the circumference N*d where d equals lambda/2
                d = obj.c/(2*fd);
                cir = N*d;
                r = cir/(2*pi);
                theta = 2*pi/N*(0:N-1)';
                
                %matrix containing the position of each antenna
                obj.G = [r*cos(theta) r*sin(theta) zeros(N,1)];
            end
            
        end
        
        %Calculate the replica vectors for the array
        %V = CalculateReplicaVectors(az, el, lambda)
        function V = CalculateReplicaVectors(obj, az, el, lambda)
            %unit length vector in direction of propagation
            p = [cos(az).*cos(el); sin(az).*cos(el); sin(el)];
            
            %replica vector (array response) for wavelength lambda
            V = exp(-1i*2*pi*obj.G*p*(1./lambda));
        end
        
        %Perform broadband simulation of waveform arriving at
        %array inputs.  Performs time shifts in frequency domain
        %Returns frequency domain so multiple signals can be summed
        function waveform_fft = BroadbandSimulation(obj, az, el, OFDM_inst)
            %size of broadband waveform
            len_fft = length(OFDM_inst.waveform);
            %sampling period of broadband waveform
            Ts = OFDM_inst.Ts/(OFDM_inst.L*OFDM_inst.N);
            %sampling frequency
            Fs = 1/Ts;
            
            %Calculate the frequency bins for FFT of length len
            %in wavelength
            f = Fs*linspace(0,1-1/len_fft,len_fft);
            mid = round(len_fft/2)+1;
            %subtract Fs to set range to +/-Fs
            f(mid:end) = f(mid:end) - Fs;
            lambda = 3e8./f;
            
            %Calculate timeshifts in frequency domain
            %timeshift in frequency domain is exp(j*omega*t) which is 
            %equivalent to the replica vector for narrowband omega
            %Calculate replica vectors for each frequency bin of FFT
            V = obj.CalculateReplicaVectors(az, el, lambda);
            %V is array [Array size by len]
            [len_array,~] = size(V);
            
            %Frequency domain of broadband waveform from OFDM object
            s = fft(OFDM_inst.waveform);
            
            %each fft bin is multiplied by corresponding replica vector
            waveform_fft = V.*(ones(len_array,1)*s.');            
        end
        
        %K-Omega Plot (frequency vs wavenumber)
        %x = array input Number of Sensors x Number of Samples
        %L = block length
        %M = block overlap
        %twin, swin = temporal and spatial windows
        %ntFFT, nsFFT = size of temporal and spatial FFTs
        function PlotKOmega(obj, x, Fs, L, M, twin, swin, ntFFT, nsFFT)
            %Change array size to Samples x Sensors
            x = x.';
            
            %Number of samples and sensors
            [N_samples, N_sensors] = size(x);
                        
            %Temporal and Spatial FFTs
            %Number of new samples per block
            new_samples = L - M;

            %determine number of blocks
            %round up to next integer
            num_blocks = ceil(N_samples/new_samples);

            %append L zeros to x to cover last block with partial data
            x = [x; zeros(L,N_sensors)];

            %make sure twin is a column vector and normalize to 1
            twin = twin(:)/sum(twin);
            %expand twin to a matrix L x number of sensors
            twin = twin*ones(1,N_sensors);

            %make sure swin is a row vector and normalize to 1
            swin = (swin(:)/sum(swin)).';
            swin = ones(ntFFT/2,1)*swin;

            %seperate x into overlapping blocks
            for i = 1:num_blocks
                %create blocks of x of lenght L with overlap M and apply temporal
                %window
                %Samples x Sensors x Blocks
                blocks(:,:,i) = x(new_samples*(i-1)+1 : L+new_samples*(i-1),:).*twin;

                %perform temporal FFT of length ntFFT on block at each sensor
                temporal_fft(:,:,i) = fft(blocks(:,:,i), ntFFT);

                %apply spatial window to rows
                %since data is real you can disgard half of the temporal FFT
                temporal_fft_real(:,:,i) = temporal_fft(1:ntFFT/2,:,i).*swin;

                %perform spatial FFT of length nsFFT on block at each frequency
                %takes FFT of rows, apply fftshift only on rows
                spatial_fft(:,:,i) = fftshift(fft(temporal_fft_real(:,:,i), nsFFT, 2),2);
            end

            %average the results over the blocks
            %sum |spatial_fft|^2
            freq_kz_spectra = 1/num_blocks*sum(abs(spatial_fft).^2,3);
            
            %Plotting
            %generate frequency index
            f = Fs/2*linspace(0,1-1/(ntFFT/2),(ntFFT/2));
            
            %distance between antennas
            d = obj.c/(2*obj.fd);
            
            %generate spatial index psi
            psi = 2*pi*linspace(0,1-1/nsFFT,nsFFT);
            mid = round(nsFFT/2)+1;
            %subtract 2*pi to set axis to -pi/ to pi
            psi(mid:end) = psi(mid:end) - 2*pi;
            psi = fftshift(psi);
            %calculate the wavenumber kz
            kz = -1/d*psi;

            figure;
            %plot the frequency-wavenumber spectra on log scale
            %transpose necessary to put frequency on x axis
            imagesc(f,kz,10*log10(freq_kz_spectra.'));
            xlabel('Frequency (Hz)');
            ylabel('kz');
            title('Frequency Wavenumber Beamformer');
            colorbar;
        end
        
        %Plot the 2D Beampattern at specified elevation
        function PlotBeampattern2D(obj, el, w, lambda)
            %if weights not specified use conventional beamformer
            if nargin < 3
                if obj.M > 0
                    NumElements = obj.N*obj.M;
                else
                    NumElements = obj.N;
                end
                w = 1/NumElements*ones(NumElements,1);
                lambda = obj.c/obj.fd;
            end
            
            %azimuth angles from 0 to PI
            az = (0:pi/1024:pi);
            el = el*ones(1,1025);
            
            %Replica vectors for array
            V = obj.CalculateReplicaVectors(az, el, lambda);
            %Beampattern = weights'*ReplicaVectors
            BP = w'*V;
            
            %plot beam pattern
            %figure;
            plot(180*az/pi,10*log10(abs(BP).^2));
            title('Beampattern');
            xlabel('Azimuth');
        end
        
        %Plot the 3D Beampattern
        function PlotBeampattern3D(obj, w, lambda)
            %if weights not specified use conventional beamformer
            if nargin < 2
                if obj.M > 0
                    NumElements = obj.N*obj.M;
                else
                    NumElements = obj.N;
                end
                w = 1/NumElements*ones(NumElements,1);
                lambda = obj.c/obj.fd;
            end
            
            %
            az = (0:pi/256:pi);
            el = (0:pi/256:pi);
            %calculate for all evelations
            for i = 1:257
                %calculate replica vectors for constant elevation
                V = obj.CalculateReplicaVectors(az,el(i)*ones(1,257),lambda);
                BP(i,:) = w'*V;                
            end
            
            %convert azimuth, elevation and vector length to cartesian
            %coordinates
            [x,y,z] = sph2cart(ones(257,1)*az,el'*ones(1,257),abs(BP).^2);
            
            %3D beampattern plot
            figure;
            mesh(x,y,z);
            title('3D Beampattern Plot');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            axis([-1 1 -1 1 0 1]);          
        end        
    end    
end