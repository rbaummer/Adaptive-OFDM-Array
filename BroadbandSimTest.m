function BroadbandSimTest

close all;

network = adhoc_network(AccessPoint([0 0 0], AntennaConfiguration('Linear', 30, 0, 20e6)));

%% Test Data
N = 30;
c = 3e8;
d = c/(2*network.nodes(1).AntennaConfiguration.fd);
Ts = network.nodes(1).OFDM_inst.Ts/(network.nodes(1).OFDM_inst.L*network.nodes(1).OFDM_inst.N);
Fs = 1/Ts;
%create 10 second time vector sampled at Fs
%t = 0:1/fs:10;
t = 0:Ts:Ts*50000;

%generate signal at 45 degrees in non aliasing region (f<fd)
% f = 18e6;
% az = 100*pi/180;
% sig1 = cos(2*pi*f*t);

f0 = 18e6;    %start frequency
f1 = 20e6;   %stop frequency
k = (f1-f0)/t(end);

az = 25*pi/180;
sig1 = cos(2*pi*(f0*t + k/2*t.^2));

%% Broadband BF test
%size of broadband waveform
len_fft = length(sig1);

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
V = network.nodes(1).AntennaConfiguration.CalculateReplicaVectors(az, 0, lambda);
%V is array [Array size by len]
[len_array,~] = size(V);

%Frequency domain of broadband waveform from OFDM object
s = fft(sig1.');

%each fft bin is multiplied by corresponding replica vector
waveform_fft = V.*(ones(len_array,1)*s.');   

%% Calc array input test
n=sqrt(0.00001)*randn(size(waveform_fft));

%cumulative_fft is a Nxlength matrix.  Take IFFT across rows
%array_waveform = real(ifft(cumulative_fft,[],2)) + n;
array_waveform = (ifft(waveform_fft,[],2)) + n;

%K-Omega plot of input at node 1 array
L = 1024;
M = L/2;
twin = ones(L,1);
%swin = ones(size(array_waveform,1),1);
swin = hamming(size(array_waveform,1));
% twin = hamming(L);
% swin = hamming(nrcv);
ntFFT = 2048;
nsFFT = 2048;
Ts = network.nodes(1).OFDM_inst.Ts/(network.nodes(1).OFDM_inst.L*network.nodes(1).OFDM_inst.N);
network.nodes(1).AntennaConfiguration.PlotKOmega(array_waveform, 1/Ts, L, M, twin, swin, ntFFT, nsFFT)
