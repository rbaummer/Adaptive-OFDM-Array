%Author: Rob Baummer
%
function ofdm_test_script

close all;

% OFDM_inst = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
% OFDM_inst.assigned_channels = 1:16;
OFDM_inst = OFDM(128, [1 8 32 48 64 65 86 128],512,224e-6,28e-6,2);
OFDM_inst.assigned_channels = 1:128;
OFDM_inst.waveform = 1000;

f = OFDM_inst.CalculateCenterFrequencies;

OFDM_inst.demodulate(OFDM_inst.waveform);

%test serial to parallel at non-baseband
parallel_waveform = OFDM_inst.serialtoparallel(OFDM_inst.waveform);

OFDM_inst2 = OFDM(128, [1 8 32 48 64 65 86 128],512,224e-6,28e-6,4);
OFDM_inst2.assigned_channels = 1:128;
OFDM_inst2.waveform = 1000;

OFDM_inst.demodulate(OFDM_inst2.waveform);

figure;
subplot(2,1,1)
pwelch(OFDM_inst.waveform, [], [], [], 1/(OFDM_inst.Ts/(8*OFDM_inst.N)));
subplot(2,1,2)
pwelch(OFDM_inst2.waveform, [], [], [], 1/(OFDM_inst.Ts/(8*OFDM_inst.N)));