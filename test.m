%DVB-T 2K Transmission
%The available bandwidth is 8 MHz
%2K is intended for mobile services
clear all;
close all;
%DVB-T Parameters
Tu=224e-6; %useful OFDM symbol period
T=Tu/2048; %baseband elementary period
G=0; %choice of 1/4, 1/8, 1/16, and 1/32
delta=G*Tu; %guard band duration
Ts=delta+Tu; %total OFDM symbol period
Kmax=1705; %number of subcarriers
Kmin=0;
FS=4096; %IFFT/FFT length
q=10; %carrier period to elementary period ratio
fc=q*1/T; %carrier frequency
Rs=4*fc; %simulation period
t=0:1/Rs:Tu;
%Data generator (A)
M=Kmax+1;
rand('state',0);
a=-1+2*round(rand(M,1)).'+1i*(-1+2*round(rand(M,1))).';
A=length(a);
info=zeros(FS,1);
info(1:(A/2)) = [ a(1:(A/2)).']; %Zero padding
info((FS-((A/2)-1)):FS) = [ a(((A/2)+1):A).'];
%Subcarriers generation (B)
carriers=FS.*ifft(info,FS);
tt=0:T/2:Tu;
figure(1);
subplot(211);
stem(tt(1:20),real(carriers(1:20)));
subplot(212);
stem(tt(1:20),imag(carriers(1:20)));
figure(2);
f=(2/T)*(1:(FS))/(FS);
subplot(211);
plot(f,abs(fft(carriers,FS))/FS);
subplot(212);
pwelch(carriers,[],[],[],2/T);


% D/A simulation
L = length(carriers);
chips = [ carriers.';zeros((2*q)-1,L)];
p=1/Rs:1/Rs:T/2;
g=ones(length(p),1); %pulse shape
figure(3);
stem(p,g);
dummy=conv(g,chips(:));
u=[dummy(1:length(t))]; % (C)
figure(4);
subplot(211);
plot(t(1:400),real(u(1:400)));
subplot(212);
plot(t(1:400),imag(u(1:400)));
figure(5);
ff=(Rs)*(1:(q*FS))/(q*FS);
subplot(211);
plot(ff,abs(fft(u,q*FS))/FS);
subplot(212);
pwelch(u,[],[],[],Rs);
[b,a] = butter(13,1/20); %reconstruction filter
[H,F] = FREQZ(b,a,FS,Rs);
figure(6);
plot(F,20*log10(abs(H)));
uoft = filter(b,a,u); %baseband signal (D)
figure(7);
subplot(211);
plot(t(80:480),real(uoft(80:480)));
subplot(212);
plot(t(80:480),imag(uoft(80:480)));
figure(8);
subplot(211);
plot(ff,abs(fft(uoft,q*FS))/FS);
subplot(212);
pwelch(uoft,[],[],[],Rs);
%Upconverter
s_tilde=(uoft.').*exp(1i*2*pi*fc*t);
s=real(s_tilde); %passband signal (E)
figure(9);
plot(t(80:480),s(80:480));
figure(10);
subplot(211);
%plot(ff,abs(fft(((real(uoft).').*cos(2*pi*fc*t)),q*FS))/FS);
%plot(ff,abs(fft(((imag(uoft).').*sin(2*pi*fc*t)),q*FS))/FS);
plot(ff,abs(fft(s,q*FS))/FS);
subplot(212);
%pwelch(((real(uoft).').*cos(2*pi*fc*t)),[],[],[],Rs);
%pwelch(((imag(uoft).').*sin(2*pi*fc*t)),[],[],[],Rs);
pwelch(s,[],[],[],Rs);