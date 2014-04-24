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
        %128 channels, FFT size 512, Ts = 224us, Tg = 28us, Mn = 2
        %pilot channels [1 32 64 65 96 128]
        OFDM_default = OFDM(128, [1 32 64 65 96 128],512,224e-6,28e-6,2);
        %OFDM_default = OFDM(512, [1 128 256 257 384 512],2048,224e-6,28e-6,2);
        AntennaConfiguration_default = AntennaConfiguration('Linear', 10, 0, 20e6);
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
                obj.OFDM_inst.assigned_channels = 1:128;
                obj.OFDM_inst.waveform = 10;
            else
                obj.OFDM_inst = OFDM_inst;
            end
        end
        
    end
end