%Author: Rob Baummer
%
%AccessPoint describes a user in an adhoc network
%
%properties
%   location: [x y z]
%   AntennaConfiguration: describes the antenna array configuration for the
%     user
%   Channel: place holder
%   OFDM_inst: describes the orthogonal frequency division multiplexing
%     setup for the user
%
%default constructor: AccessPoint(location, OFDM_inst)
%   if OFDM_inst is not provided a default is used
classdef AccessPoint
    properties
        %location x,y,z
        location = [0 0 0];
        AntennaConfiguration
        Channel = 1;
        OFDM_inst
    end
    
    properties(Constant = true, Hidden = true)
        %16 Channels, FFT size 64, Ts 224us, Tg 224 us, Mn 2
        %OFDM_default = OFDM(16, [1 8 9 16],64,224e-6,224e-6,2);
        %OFDM_default = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
        %128 channels, FFT size 512, Ts = 224us, Tg = 28us, Mn = 2
        %pilot channels [1 8 32 48 64 65 86 128]
        OFDM_default = OFDM(128, [1 8 32 48 64 65 86 128],512,224e-6,28e-6,2);
        %OFDM_default = OFDM(512, [1 128 256 257 384 512],2048,224e-6,28e-6,2);
        AntennaConfiguration_default = AntennaConfiguration('Linear', 10, 0, 5e6);
    end
    
    methods
        %constructor
        function obj = AccessPoint(location, AntennaConfiguration_inst, OFDM_inst)
            obj.location = location;
            if nargin < 2
                obj.AntennaConfiguration = obj.AntennaConfiguration_default;
            else
                obj.AntennaConfiguration = AntennaConfiguration_inst;
            end
            if nargin < 3
                obj.OFDM_inst = obj.OFDM_default;
                %obj.OFDM_inst.assigned_channels = 1:128;
                obj.OFDM_inst.assigned_channels = 1:obj.OFDM_inst.channels;
                obj.OFDM_inst.waveform = 2000;
            else
                obj.OFDM_inst = OFDM_inst;
            end
        end
        
    end
end