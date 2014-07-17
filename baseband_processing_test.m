function baseband_processing_test

close all;

%verify that beamforming can be done at baseband since time delays are
%preserved when downconverting
array = AntennaConfiguration('Linear',8,0,10e6);

Fs = 64e6;

%1 ms sampled at 64 MHz
t = 0:1/Fs:1e-3-1/Fs;

sig = [zeros(1,10) ones(1,200) zeros(1,length(t)-210)];

%2 MHz carrier
%car = exp(1i*2*pi*2e+06*t);
car = cos(2*pi*2e+06*t);
passband = car.*sig;
passband = passband(:);

%% broadband array setup
az = pi/4;
el = 0;
len_fft = length(passband);
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
V = array.CalculateReplicaVectors(az, el, lambda);
%V is array [Array size by len]
[len_array,~] = size(V);

%Frequency domain of broadband waveform from OFDM object
%sig =sig.*feval('tukeywin',length(sig),0.1);
s = fft(passband,len_fft);

%each fft bin is multiplied by corresponding replica vector
waveform_fft = V.*(ones(len_array,1)*s.');

%waveform_fft is a Nxlength matrix.  Take IFFT across rows
array_waveform = real(ifft(waveform_fft,[],2));

%% plot delays at passband
range = 50:150;
figure;
hold on;
plot(array_waveform(1,range),'b');
plot(array_waveform(2,range),'g');
plot(array_waveform(3,range),'c');
plot(array_waveform(4,range),'r');
plot(array_waveform(5,range),'k');
plot(array_waveform(6,range),'-.b');
plot(array_waveform(7,range),'-.g');
plot(array_waveform(8,range),'-.c');

%% downconvert to near baseband
sig2 = exp(1i*2*pi*1e+06*t);

Hf = fdesign.lowpass('N,Fc',200,0.08);
H = design(Hf);

array_waveform2 = real(array_waveform.*(ones(8,1)*sig2));
for i = 1:8
    array_waveform2(i,:) = filter(H,array_waveform2(i,:));
end

figure;
subplot(2,1,1);
plot(f,abs(fft(array_waveform(1,:))));
subplot(2,1,2);
plot(f,abs(fft(array_waveform2(1,:))));

range = 150:250;
figure;
hold on;
plot(array_waveform2(1,range),'b');
plot(array_waveform2(2,range),'g');
plot(array_waveform2(3,range),'c');
plot(array_waveform2(4,range),'r');
plot(array_waveform2(5,range),'k');
plot(array_waveform2(6,range),'-.b');
plot(array_waveform2(7,range),'-.g');
plot(array_waveform2(8,range),'-.c');

%phase is preserved when changing carrier frequency as seen in the two
%plots hence you can heterdyne the waveform at each antenna and then
%beamform



