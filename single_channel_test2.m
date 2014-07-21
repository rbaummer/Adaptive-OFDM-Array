function single_channel_test2

close all

%Angles of Arrival for QPSK signals
angles = pi/180*[20 45 74 107 137];
%Signal Power
sigma_s = [1 1 1 1 1];

%Run length (Number OFDM symbols)
len = 1000;
%Array length
N = 10;
%Noise power
sigma_n = 0.05; %26 dB SNR

%% generate test OFDM sources arriving at N element array
%s1
OFDM_inst1 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
OFDM_inst1.assigned_channels = 1:16;
OFDM_inst1.waveform = len;

%s2
OFDM_inst2 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
OFDM_inst2.assigned_channels = 1:16;
OFDM_inst2.waveform = len;

%s3
OFDM_inst3 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
OFDM_inst3.assigned_channels = 1:16;
OFDM_inst3.waveform = len;

%s4
OFDM_inst4 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
OFDM_inst4.assigned_channels = 1:16;
OFDM_inst4.waveform = len;

%s5
OFDM_inst5 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
OFDM_inst5.assigned_channels = 1:16;
OFDM_inst5.waveform = len;

%% broadband simulation of array inputs
ant = AntennaConfiguration('Linear', 10, 0, 600e3);
%Ts = OFDM_inst1.Ts/(OFDM_inst1.L*OFDM_inst1.N);
%Fs = 1/Ts;
%ant = AntennaConfiguration('Linear', 10, 0, Fs);

w1_fft = ant.BroadbandSimulation(angles(1),0,OFDM_inst1);
w2_fft = ant.BroadbandSimulation(angles(2),0,OFDM_inst2);
w3_fft = ant.BroadbandSimulation(angles(3),0,OFDM_inst3);
w4_fft = ant.BroadbandSimulation(angles(4),0,OFDM_inst4);
w5_fft = ant.BroadbandSimulation(angles(5),0,OFDM_inst5);

w_fft = w1_fft + w2_fft + w3_fft + w4_fft + w5_fft;

wf = (ifft(w_fft,[],2));

%signal power mean at each antenna
p = mean(abs(w_fft),2);
%siganl power mean
p_avg = mean(p);

SNR = -20;
n = sqrt(10^(SNR/20)*p_avg)*randn(size(w_fft));

waveform = wf + n;

%%
for i = 1:ant.N
    %time domain passband OFDM symbols
    r = OFDM_inst1.serialtoparallel(waveform(i,:));
    %frequency domain baseband OFDM symbols
    received_symbol = 1/OFDM_inst1.N*fft(r);
    %save number of OFDM symbols - 4 (minimum length after
    %synchronization)
    s(:,:,i) = received_symbol(:, 1:OFDM_inst1.num_symbols - 4);
end

s(9:56,:,:) = [];

%normalize the power across time and antennas
p = mean(mean(abs(s),2),3)/2 %#ok<NOPRT>
for i = 1:OFDM_inst1.channels
    s(i,:,:) = s(i,:,:)/p(i);
end

x = squeeze(s(16,:,:)).';

%%

% %% test CMA algorithm with QPSK
% %set intitial weights to CBF
% w_init = 1/N*ones(N,1);
% %w_init = [0 0 0 0 0 1 0 0 0 0]';
% 
% %Calculate step size
% %mu from "Convergence Behavior of the Constant Modulus Algorithm"
% R = V_qpsk*V_qpsk';
% [~,D] = eig(R);
% mu = 0.25/(6*real(max(max(D))^2));
% 
% %Run CMA
% [w, err] = CMA(w_init, mu, x);
% w_end = w(:,end);
% 
% %check to make sure algorithm converged
% assert(isfinite(w_end(1)), sprintf('Error: CMA: QPSK did not converge mu = %f',mu)); 
% 
% %plot before and after scatter plot
% h = scatterplot(w_init'*x,1,0,'r.');
% hold on;
% scatterplot(w_end'*x,1,0,'b.',h)
% title('CMA Signal Constellation');
% legend('Before','After');
% 
% %plot beam pattern after convergence
% figure;
% bf_plot(w_end, angles);
% title(sprintf('QPSK CMA Beam Pattern\n WNG: %2.1f', 1/(w_end'*w_end)));
% xlabel('Degrees')
% ylabel('Power (dB)');
% figure;
% plot(abs(err));
% title('CMA Learning Curve');
% xlabel('Iteration');
% ylabel('|Error|');


%% test orthogonal CMA algorithm with QPSK
%set initial weights to CBF
w_init = 1/N*ones(N,1);
%w_init = [0 0 0 0 0 1 0 0 0 0]';
%set step size
mu = 0.01;
%set forgetting factor
alpha = 1-0.985;
%set initial inverse correlation matrix
R = diag(ones(N,1));

%Run o-CMA
[w, err, R_inv] = oCMA(w_init, R, mu, alpha, x);
w_end = w(:,end);

%check to make sure algorithm converged
assert(isfinite(w_end(1)), 'Error: o-CMA: QPSK did not converge'); 

%plot before and after scatter plot
h = scatterplot(w_init'*x,1,0,'r.');
hold on;
scatterplot(w_end'*x,1,0,'b.',h)
title('o-CMA Signal Constellation');
legend('Before','After');

%plot beam pattern after convergence
figure;
bf_plot(w_end, angles);
title(sprintf('QPSK O-CMA Beam Pattern\n WNG: %2.1f', 1/(w_end'*w_end)));
xlabel('Degrees')
ylabel('Power (dB)');
figure;
plot(abs(err));
title('O-CMA Learning Curve');
xlabel('Iteration');
ylabel('|Error|');

