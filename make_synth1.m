%%% Output file
outfile='synth1.mat';
writenow=false;

%%% Simulation parameters
c0=1500;				% sound speed m/s
fs=150;					% sampling frequency (Hz)
Tw=20;
Tpad=1;
win='hanning';

%%% Signal parameters
src(1).type='cw';			% Specify the signals
src(1).Tw=Tw;
src(1).fc=20;
src(1).thetadeg=90+25;
src(1).amp=2;

src(2).type='bbnoise';			% Specify the signals
src(2).Tw=Tw;
src(2).thetadeg=60;
src(2).sigsq=3;

src(3).type='lfm';
src(3).Tw=20;
src(3).fc=50;
src(3).bw=30;
src(3).thetadeg=95;
src(3).amp=5;

%%% Noise parameters
noise.sigwsq=.5;

%%% Receiver parameters
N=40;					% Specify the receiver
rcvr.d=30;
rcvr.z=(0:N-1)'*rcvr.d;
rcvr.c=c0;

x=make_synth_data(src,rcvr,noise,Tw,Tpad,fs,win);

dz=rcvr.d;
[nsig,nrcv]=size(x);

if writenow
  save(outfile,'fs','dz','nrcv','nsig','x','-V6');
end

params.fs=150;
%params.blocklength=1024;
params.blocklength=512;
params.skip=params.blocklength/2;
params.timewindow=blackman(params.blocklength);
params.spacewindow=hanning(40);
params.Ntempfft=4096;
params.Nspatfft=1024;
params.detrend='removemean';
params.dz=30;
params.c0=1500;

goodchan=ones(40,1);
fmaxplot=params.fs/2;

%%% Compute komega spectrum
clf;
makeplots=1;
if makeplots; figure(2); end;
[komega_spectra,temp_spectra,fvec,kzvec,nblk]=...
    calc_komega_v2(x,goodchan,params);
if makeplots
  imagesc(fvec,kzvec,10*log10(abs(komega_spectra)))
  hold on
  plot(fvec,2*pi*fvec/params.c0,'w-','Linewidth',1.1)
  plot(fvec,-2*pi*fvec/params.c0,'w-','LineWidth',1.1)
  plot(fvec,-2*pi*fvec/params.c0*cos(pi/180*115),'k--','Linewidth',1.1)
  plot(fvec,-2*pi*fvec/params.c0*cos(pi/180*80),'k--','Linewidth',1.1)
  hold off
  axis([0 fmaxplot -pi/params.dz pi/params.dz])
  set(gca,'YDir','normal')
  caxis([max(caxis)-60 max(caxis)])
  colorbar
  xlabel('Frequency (Hz)')
  ylabel('Wavenumber k_z')
%  tstr=sprintf('k-\\Omega for %s',datafile);
%  title(tstr)
end  
