%Author: Rob Baummer
%
%Create OFDM object
%
%properties
%   channels: number of channels in system
%   assigned_channels: channels assigned to this node
%   pilot_channels: define pilot channels
%   pilot_symbols: define pilot symbols
%   num_symbols: number of OFDM symbols generated
%   N: IFFT/FFT size, N > channels
%   delta_f: spacing of individual channels, delta_f = 1/Ts
%   Ts: symbol duration
%   Tg: guard duration
%   Fs: sampling frequency of waveform, Fs = L*N*delta_f
%   Mn: bits per channel
%   waveform: OFDM waveform with random data
%default constructor: OFDM(channels, carrier_spacing, symbol_duration, 
%     sampling_frequency, ModulationType, len)
%method set: obj.assigned_channels = new_channels
%   new_channels must be a subset of obj.channels
%method set: obj.waveform = len
%   creates a random bit string len OFDM symbols 
classdef OFDM
    %OFDM properties
    properties
        channels
        assigned_channels
        pilot_channels
        pilot_symbols
        num_symbols
        N
        delta_f
        Ts
        Tg
        Fc
        Mn
        symbol_set
        waveform
    end
    
    properties (Constant = true, Hidden = true)
        %Samples per period
        L = 8;
        %Carrier to OFDM sample period ratio
        k = 2;
        %QPSK signal constellation from DVB-S standard
        Sn_2bits = [1+1i; 1-1i; -1+1i; -1-1i];
        Sn_3bits = [1; 1+1i; 1i; -1+1i; -1; -1-1i; -1i; 1-1i];
        %16 QAM signal constellation from DVB-S standard
        Sn_4bits = [3+3i; 3+1i; 1+3i; 1+1i; 3-3i; 3-1i; 1-3i; 1-1i; -3+3i; -3+1i; -1+3i; -1+1i; -3-3i; -3-1i; -1-3i; -1-1i];
    end

    methods
        %constructor
        function obj = OFDM(channels, pilots, FFT_size, Ts, Tg, Mn)
            obj.channels = channels;
            obj.assigned_channels = 0;
            obj.pilot_channels = pilots;
            obj.N = FFT_size;
            obj.delta_f = 1/Ts;
            obj.Ts = Ts;
            obj.Tg = Tg;
            obj.Fc = obj.k*FFT_size*(1/Ts);
            %obj.Fc = 2e6;
            obj.Mn = Mn;
            %determine symbol set
            if Mn == 2
                obj.symbol_set = obj.Sn_2bits;
            elseif Mn == 3
                obj.symbol_set = obj.Sn_3bits;
            elseif Mn == 4
                obj.symbol_set = obj.Sn_4bits;
            else
                obj.symbol_set = obj.Sn_2bits;
            end
            obj.waveform = 0;
            obj.pilot_symbols = 0;
        end
    
        %assign new channels to OFDM object
        function obj = set.assigned_channels(obj, new_channels)
            if max(new_channels) > obj.channels %#ok<*MCSUP>
                error('Total System Channels must contain new channel assigments');
            else
                channel = new_channels;
                %remove pilot channels from channel assigment
                channel(obj.pilot_channels) = [];
                obj.assigned_channels = channel;
            end
        end
        
        %generate a waveform random bits with length equal to len
        function obj = set.waveform(obj, len)
            %if channels aren't assigned don't create a waveform
            if obj.assigned_channels == 0
                obj.waveform = 0;
            %call modulate method to generate a waveform of len bits
            else
                obj.num_symbols = len;
                obj.pilot_symbols = obj.symbol_set(randi(2^obj.Mn,len,1));
                obj.waveform = obj.modulate(len);
            end 
        end
        
        %generate a OFDM waveform for the object
        function waveform = modulate(obj, len)
            %size of OFDM symbols in bits
            size_OFDM_sym = length(obj.assigned_channels)*obj.Mn;
            %Number of OFDM symbols
            num_OFDM_syms = len;
            
            %generate binary bit stream
            bits = randi([0 1], len*size_OFDM_sym, 1);

            %empty waveform;
            wf = [];
            
            %build OFDM symbols
            for i = 1:num_OFDM_syms-1
                start = (i-1)*size_OFDM_sym;
                %generate the symbols for each channel
                %preallocate complex symbol array
                Sn = zeros(length(obj.assigned_channels),1);
                for j = 1:length(obj.assigned_channels)
                    %grab the next Mn bits
                    b = bits(start + (j-1)*obj.Mn + 1: start + j*obj.Mn);
                    %convert to integer
                    d= bi2de(b');
                    %map integer to complex symbol
                    Sn(j,1) = obj.symbol_set(d+1);
                end

                %Generate zero vector number channels x 1
                OFDM_sym = zeros(obj.channels,1);
                %Set assigned channels to their respective symbols
                %This is the frequency domain OFDM symbol
                OFDM_sym(obj.assigned_channels) = Sn;
                
                %Add in pilot signals
                OFDM_sym(obj.pilot_channels) = obj.pilot_symbols(i);

                %zero stuff C to FFT size leaving OFDM_sym at center
                z = obj.N - obj.channels;
                mid = floor(obj.channels/2);
                %construct complex symbol, interpolates by N/channels
                C = obj.N/obj.channels*[OFDM_sym(1:mid); zeros(z,1); OFDM_sym(mid+1:end)];

                %Generate time domain OFDM symbol and append it to
                %previously generated symbols
                ofdm_sym = obj.N*ifft(C,obj.N);

                %generate the cyclic prefix
                %CP_size = Tg/OFDM_Ts (ofdm sample rate is Ts/N pre CP)
                CP_size = round(obj.Tg/(obj.Ts/obj.N));
                ofdm_sym_cp = [ofdm_sym(end - CP_size + 1:end); ofdm_sym];

                %append time domain ofdm symbol to previously generated
                %symbols
                wf = [wf; ofdm_sym_cp];
            end
            
            %upsample and filter OFDM time domain
            wf_baseband = interp(wf,obj.L);
% Note AWGN is added in broadband array simulation so not necessary here           
%             %add AWGN
%             n = randn(length(wf)*obj.L,1) + 1i*randn(length(wf)*obj.L,1);
%             wf_baseband = wf_baseband + n;

            %upconvert to carrier frequency
            %create time array for num_OFDM symbols sampled at Fs
            t = (0:obj.Ts/(obj.L*obj.N):num_OFDM_syms*(obj.Ts + obj.Tg))';
            carrier = exp(1i*2*pi*obj.Fc*t);

            %assign generated waveform to class property
            waveform = real(wf_baseband.*carrier(1:length(wf_baseband)));

            %debug plots
            figure;
            subplot(3,1,1)
            pwelch(wf, [], [], [], 1/(obj.Ts/(obj.N)));
            subplot(3,1,2)
            pwelch(wf_baseband, [], [], [], 1/(obj.Ts/(obj.L*obj.N)));
            subplot(3,1,3)
            pwelch(waveform, [], [], [], 1/(obj.Ts/(obj.L*obj.N)));
        end
        
        %Synchronize the OFDM symbols at the higher sample rate
        %and return the parallelized form after removing the CP.
        %One symbol per column
        function ofdm_sym = serialtoparallel(obj, waveform)
            %make sure waveform is a column vector
            waveform = waveform(:);
            
            %Downconvert the passpand signal to baseband and downsample to
            %OFDM symbol size
            baseband_waveform = obj.down_conversion(waveform);
            
            %time synchronization
            time_offset = obj.time_synchronizaiton(baseband_waveform);
            
            %serial to parallel operation forming time domain OFDM symbols
            ofdm_sym = obj.framing(baseband_waveform, time_offset);
        end

        %Generate scatter plots of waveform
        %note waveform does not need to be the obj.waveform but it should
        %be the same length
        function OFDM_sym = demodulate(obj, waveform)
            %make sure waveform is a column vector
            waveform = waveform(:);
            
            %Downconvert the passpand signal to baseband and downsample to
            %OFDM symbol size
            baseband_waveform = obj.down_conversion(waveform);
            
            %debug plots
            figure;
            pwelch(baseband_waveform, [], [], [], 1/(obj.Ts/obj.N));
            
            %time synchronization
            time_offset = obj.time_synchronizaiton(baseband_waveform);
            
            %serial to parallel operation forming time domain OFDM symbols
            ofdm_sym = obj.framing(baseband_waveform, time_offset);
            
            %generate frequency domain symbols from time domain using FFT
            OFDM_sym = fft(ofdm_sym,obj.N);
            
            %estimate the channel from pilots
            %interpolate across all obj.channels
            H = obj.channel_equalization(OFDM_sym(:,1));            

            %midpoint of channels
            mid = floor(obj.channels/2);
            %remove obj.N-obj.channels zeros(+noise) from symbol
            OFDM_sym(mid+1:obj.N-mid,:) = [];
            
            %remove channel based on estimate
            for i = 1:size(OFDM_sym,2)
                OFDM_sym(:,i) = OFDM_sym(:,i).*H./abs(H).^2;
            end
            
            %remove pilots
            OFDM_sym(obj.pilot_channels,:) = [];
            
            %generate scatter plots
            scatterplot(OFDM_sym(:,1));
        end
        
        %Perform downconversion of waveform: Heterodyne to baseband and
        %downsample to OFDM symbol size
        function baseband_waveform = down_conversion(obj, waveform)
            %find length of passband waveform
            l = length(waveform);
            
            %generate carrier
            t = obj.Ts/(obj.L*obj.N)*(0:l-1)';
            carrier = exp(-1i*2*pi*obj.Fc*t);
            
            %downconvert back to baseband
            wf_baseband = carrier.*waveform;
            
            %decimate by interpolation factor
            baseband_waveform = decimate(wf_baseband, obj.L);
        end
        
        %Use the cyclic prefix to synchronize the OFDM frame
        %Returns time offset of start of OFDM symbol
        function time_offset = time_synchronizaiton(obj, waveform)
            %Cyclic prefix size
            CP_size = round(obj.Tg/(obj.Ts/obj.N));
            frame_size = obj.N + CP_size;
            
            %time synchronization
            %since the cylic prefix is a repetition of the end of the OFDM
            %frame a correlation separated by N peaks at the end of the CP
            for i = 1:2*obj.N + CP_size
            %for i = 1:length(wf)-obj.N-CP_size%2*obj.N + CP_size
               %correlation value corrected by the average power
               corr(i) = sum(waveform(i:i+CP_size-1).*conj(waveform(i+obj.N:i+frame_size-1)));
               pavg(i) = 0.5*sum(abs(waveform(i:i+CP_size-1)).^2 + abs(waveform(i+obj.N:i+frame_size-1)).^2);
            end
            
            %time offset to beginning of OFDM frame is the maximum of the
            %correlation function
            [~,time_offset] = max(abs(corr) - pavg);
%             time_offset
%             figure;
%             plot(abs(corr)-pavg);
        end
        
        %Perform serial to parallel operation to build OFDM symbols
        %Returns OFDM symbols in time domain
        function ofdm_sym = framing(obj, waveform, time_offset)
            %Cyclic prefix size
            CP_size = round(obj.Tg/(obj.Ts/obj.N));
            frame_size = obj.N + CP_size;
            
            %find length of baseband waveform            
            len = length(waveform);
            
             %loop until no more frames to process
             %loop until no more frames to process
            i = 0;  %loop counter
            %subtract skipped samples from beginning
            len = len - time_offset;   
            while len > frame_size
                index = i*frame_size + time_offset;
                %grab frame size baseband time domain bits
                ofdm_sym_cp = waveform(index:index + frame_size-1);
                %remove cyclic prefix
                ofdm_sym(:,i+1) = ofdm_sym_cp(CP_size+1:end);
                
                len = len - frame_size;
                i = i + 1;
            end
        end
        
        %Estimate the channel from the pilots
        %Returns channel estimate
        function H = channel_equalization(obj, OFDM_sym)
            %extract pilots
            %channels 1:mid are at beginning of symbol
            %channels mid+1:end are at end of symbol
            mid = floor(obj.channels/2);
            %upper index of lower pilot channels
            m = length(obj.pilot_channels(obj.pilot_channels <= mid));
            %calculate indices of pilot channels in size N symbol
            pilot_index = [obj.pilot_channels(1:m) (obj.N - obj.channels) + obj.pilot_channels(m+1:end)];
            
            %extract pilot symbols
            pilots = OFDM_sym(pilot_index);
            
            %channel estimate is 
            %transmitted pilot symbol * conj(received pilot symbol)
            H_hat = obj.pilot_symbols(1).*conj(pilots);
            
%             %angle rotation error
%             theta = unwrap(angle(chan)).';
%             %2nd order polynomial interpolation for first half of symbol
%             coef1 = polyfit(obj.pilot_channels(1:m), theta(1:m), 2);
%             %2nd order polynomial interpolation for second half of symbol
%             coef2 = polyfit(obj.pilot_channels(m+1:end), theta(m+1:end), 2);
%             
%             %first half correction
%             H1 = exp(1i*(coef1(1)*(1:mid).^2 + coef1(2)*(1:mid) + coef1(3))).';
%             %second half correction
%             H2 = exp(1i*(coef2(1)*(mid+1:obj.channels).^2 + coef2(2)*(mid+1:obj.channels) + coef2(3))).';
% 
%             %correct channel phase rotation
%             OFDM_sym(1:mid) = OFDM_sym(1:mid).*H1;
%             OFDM_sym(obj.N - mid + 1:end) = OFDM_sym(obj.N - mid + 1:end).*H2;

            %2nd order polynomial interpolation for first half of symbol
            coef1 = polyfit(obj.pilot_channels(1:m), H_hat(1:m).', 2);
            %2nd order polynomial interpolation for second half of symbol
            coef2 = polyfit(obj.pilot_channels(m+1:end), H_hat(m+1:end).', 2);
            
            %first half correction
            H1 = (coef1(1)*(1:mid).^2 + coef1(2)*(1:mid) + coef1(3)).';
            %second half correction
            H2 = (coef2(1)*(mid+1:obj.channels).^2 + coef2(2)*(mid+1:obj.channels) + coef2(3)).';
            
            H = [H1; H2];
        end
        
        %Calculate the center frequencies for each channel in the OFDM
        %waveform.  
        function f = CalculateCenterFrequencies(obj)
            %half the channels are to the left of the carrier and half are
            %to the right
            lower_frequency = obj.Fc - ceil(obj.channels/2)*obj.delta_f;
            
            %lower_frequency is the center frequency of the first channel
            %so subtract 1 from assigned_channels
            %f = (obj.assigned_channels-1)*obj.delta_f + lower_frequency;            
            f = (1:obj.channels)*obj.delta_f + lower_frequency;            
        end
    end
end