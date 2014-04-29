%Author: Rob Baummer
%
function ofdm_test_script

close all;

%OFDM_inst = OFDM(64, [1 8 32 33 48 64],256,224e-6,28e-6,2);
OFDM_inst = OFDM(128, [1 8 32 48 64 65 86 128],512,224e-6,28e-6,2);
OFDM_inst.assigned_channels = 1:128;
OFDM_inst.waveform = 100;

f = OFDM_inst.CalculateCenterFrequencies;

% CMA = CMABeamformer(OFDM_inst);
% 
% t = 0:1/2048:1-1/2048;
% x = exp(1i*2*pi*125*t);
% 
% figure;
% plot(abs(fft(x.')),'g');
% hold on;
% plot(abs(CMA.F*x.'),'--r');
% 
% X = CMA.F*x.';
% 
% figure;
% plot(real(x),'g');
% hold on;
% plot(real(CMA.G*X),'--r');

OFDM_inst.demodulate(OFDM_inst.waveform);

%test serial to parallel at non-baseband
parallel_waveform = OFDM_inst.serialtoparallel(OFDM_inst.waveform);

OFDM_inst2 = OFDM(128, [1 8 32 48 64 65 86 128],512,224e-6,28e-6,4);
OFDM_inst2.assigned_channels = 1:128;
OFDM_inst2.waveform = 100;

OFDM_inst.demodulate(OFDM_inst2.waveform);

figure;
subplot(2,1,1)
pwelch(OFDM_inst.waveform, [], [], [], 1/(OFDM_inst.Ts/(64*OFDM_inst.N)));
subplot(2,1,2)
pwelch(OFDM_inst2.waveform, [], [], [], 1/(OFDM_inst.Ts/(64*OFDM_inst.N)));