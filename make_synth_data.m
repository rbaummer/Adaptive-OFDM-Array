function x=make_synth_data(src,rcvr,noise,Tw,Tpad,fs,win)
%%% File to make simulated data for Project I

nsrc=length(src);
taxis=0:1/fs:(Tw+2*Tpad);
startind=find(taxis==Tpad);
nrcv=length(rcvr.z);

%%% Loop through all sources and add them up
for srcind=1:nsrc

  %%% Get the time domain response of each source and transform it
  sig=zeros(length(taxis),1);

  if strcmp(src(srcind).type,'cw')
    tim=0:1/fs:src(srcind).Tw;
%    s=sin(2*pi*src(srcind).fc*tim).*feval('tukeywin',length(tim),0.1)';
    s=sin(2*pi*src(srcind).fc*tim);
    s=s*src(srcind).amp;
  elseif strcmp(src(srcind).type,'lfm')
    [senv,s,tim]=lfmpulse(src(srcind).fc,src(srcind).bw,...
	src(srcind).Tw,fs);
    s=s*src(srcind).amp;
  elseif strcmp(src(srcind).type,'bbnoise')
    tim=0:1/fs:src(srcind).Tw;
    s=randn(size(tim))*sqrt(src(srcind).sigsq);
  end

  s=s.*feval('tukeywin',length(tim),0.1)';
  subplot(211)
  plot(s)
  subplot(212)
  plot(abs(fft(s)))
  sig(startind:startind+length(tim)-Tpad,1)=s.';

  sigfft=fft(sig);
  if srcind==1
    K=length(sigfft);
    dtomega=(2*pi/K)*[0:(K-1)]';			% DT frequency vector
    mid=ceil(K/2)+1;
    dtomega(mid:K)=dtomega(mid:K)-2*pi;
    Xdft=zeros(K,nrcv);
  end
  
  for k=1:K
    rcvr.freq=dtomega(k)*fs/(2*pi);
    vsrc=replicaz(rcvr,'theta',src(srcind).thetadeg*pi/180);
    Xdft(k,:)=sigfft(k)*vsrc.'+Xdft(k,:);
  end
end

n=sqrt(noise.sigwsq)*randn(size(Xdft));
x=real(ifft(Xdft))+n;
