function single_channel_test

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

% %s2
% OFDM_inst2 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
% OFDM_inst2.assigned_channels = 1:16;
% OFDM_inst2.waveform = len;
% 
% %s3
% OFDM_inst3 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
% OFDM_inst3.assigned_channels = 1:16;
% OFDM_inst3.waveform = len;
% 
% %s4
% OFDM_inst4 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
% OFDM_inst4.assigned_channels = 1:16;
% OFDM_inst4.waveform = len;
% 
% %s5
% OFDM_inst5 = OFDM(16, [1 8 9 16],64,224e-6,28e-6,2);
% OFDM_inst5.assigned_channels = 1:16;
% OFDM_inst5.waveform = len;

%pick single channel for each
r = OFDM_inst1.serialtoparallel(OFDM_inst1.waveform);
received_symbol = 1/OFDM_inst1.N*fft(r);
s(:,1) = received_symbol(1, 1:len - 4);

%r = OFDM_inst2.serialtoparallel(OFDM_inst2.waveform);
%received_symbol = 1/OFDM_inst2.N*fft(r);
s(:,2) = received_symbol(2, 1:len - 4);

%r = OFDM_inst3.serialtoparallel(OFDM_inst3.waveform);
%received_symbol = 1/OFDM_inst3.N*fft(r);
s(:,3) = received_symbol(3, 1:len - 4);

%r = OFDM_inst4.serialtoparallel(OFDM_inst4.waveform);
%received_symbol = 1/OFDM_inst4.N*fft(r);
s(:,4) = received_symbol(4, 1:len - 4);

%r = OFDM_inst5.serialtoparallel(OFDM_inst5.waveform);
%received_symbol = 1/OFDM_inst5.N*fft(r);
s(:,5) = received_symbol(5, 1:len - 4);

%normalize the power across time and antennas
p = mean(abs(s),1);
for i = 1:5
    s(:,i) = s(:,i)./p(i);
end

s = s.';

%% Array
%sensor spacing d = lambda/2
pz = ((0:N-1) - (N-1)/2) * 0.5;
%angles in kz space, lambda normalized
lambda = 1;
kz = -2*pi./lambda*cos(angles);
%replica vectors
V = (ones(N,1)*sigma_s.^2).*exp(-1i*pz'*kz);

% Output
%AWGN
x = 0.05*(randn(N,len-4)+1i*randn(N,len-4));
%sources at each sensor
x = x + V*s;

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

