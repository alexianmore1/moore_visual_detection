function [Pxx, Pyy, Pxy, XYphi, Cxy, F, nTapers]= multitaperSpectrum(X,Y,Fs,bandWidth,NFFT,removeTemporalMean,RemoveEnsembleMean,nTapers)
% [Pxx, Pyy, Pxy, XYphi, Cxy, F, nTapers]= multitaperSpectrum(X,Y,Fs,bandWidth,NFFT,RemoveEnsembleMean,nTapers)
% Single window multitaper spectra for bivariate sequences: 
% X(timepoints,trials) and Y(timepoints,trials)
% Fs := sampling rate in Hertz
% BandWidth := full bandwidth, e.g. 4 Hz
% RemoveTemporalMean = true -> subtract the temporal mean from each trial
% RemoveEnsembleMean = true -> subtract the ensemble mean from each trial
% nTapers := number of tapers or simply nTapers = []
% Output:
% Pxx:= X Power spectral density (one-sided)
% Pyy:= Y power spectrum density
% Pxy:= cross-spectral density
% XYphi:= spectral phase
% Cxy:= SQUARED spectral coherence, 
% F:= frequency axis (one-sided) 
% nTapers:=2*N*W-1, where W=BandWidth/(2*Fs), if the input nTapers = []

% Note: Units in (.)^2/Hz 2 * dt^2 * 1/T * Pxx;
% If X is a stationary (binned) point process, for high frequency the above corresponds to
% 2 * P[x=1 in dt] * dt = 2 * lambda*dt * dt where lambda is the mean rate of the process.
% To obtain the power in units of a rate (events/second or Hz), 
% (useful given that the power spectrum spectrum of a stationary point process asymptotes the mean rate as the frequency \to \infty)
% one simply renormalizes as power = 1/2 * 1/dt^2 * Pxx

% Example: Power spectrum of a 50Hz sinusoid signal plus Gaussian white noise estimated via multitaper
% clear all,close all
% T=1; Fs = 1000; dt=1/Fs; N=T*Fs;
% t=0:dt:T-dt;
% f=50;
% x=sin(2*pi*f*t)' + 0.1*randn(N,1);
% % multitaper
% R = 4;%full bandwidth
% NFFT=N;%no zero padding
% % [Pxx, Pyy, Pxy, XYphi, Cxy, F, nTapers]= multitaperSpectrum(X,Y,Fs,bandWidth,NFFT,removeTemporalMean,RemoveEnsembleMean,nTapers)
% cd C:\Dropbox\workCodes\spectral_analysis\multitaper\
% [Pxx, ~, ~, ~, ~,F, nTapers]= multitaperSpectrum(x,x,Fs,R,NFFT,1,1,[]);
% figure(1),set(gcf,'color',[1 1 1])
% plot(F,10*log10(Pxx),'k.-')
% axis([0 100 -100 0])
% xlabel('Frequency (Hz)')
% ylabel('dB')

% This code was written for use in the course NEUR 2110 Statistical Neuroscience
% This code prioritizes transparency. For efficiency/general applications one might use the Chronux toolbox.
% WT, Brown University

if nargin<8,display('incorrect number of arguments; exiting ...');return,end
if size(X)~=size(Y), display('X and Y have different dimensions, exiting now ...');return; end
if size(X,1)==1, X=X(:);end
if size(Y,1)==1, Y=Y(:);end
[N,ntrials]=size(X);

T=N/Fs;%Duration of each realization or trial; in seconds
dt=1/Fs;
df=1/T;%frequency resolution, in per second
F0=df;%lowest nonzero frequency
Nyquist=Fs/2;
% Compute the true frequency resolution which is also the lowest estimated frequency F0=df=1/T=Fs/N; 
% This is computed before zero-padding. Zero padding will lead to a finer
% frequency bin, but the true resolution remains the same
F0=Fs/N; 

% The half-bandwidth W is expressed in terms of frequencies in [-0.5, 0.5]. So we
% divide bandWidth by the sampling rate Fs. ("bandWidth" is the full bandwidth)
W = 0.5 * bandWidth/Fs;

% Get the dpss
[E,V] = dpss(N,round(N*W));

% We need tapers such that the integral of the square of each taper equals 1; 
% dpss computes orthogonal tapers, i.e. the SUM of squares equals 1 - Thus we need
% to multiply the tapers by sqrt(Fs) to get the right normalization
% E = E*sqrt(Fs);
E = E*sqrt(N);

if isempty(nTapers)
    k = min(round(2*N*W),N); % By convention, the first 2*NW eigenvalues/vectors are stored
    k = max(k-1,1);
    % nTapers=round(2*N*W-1);
    nTapers=k;    
elseif nTapers > size(E,2)
    Pxx=[]; Pyy=[]; Pxy=[]; XYphi=[]; Cxy=[]; F=[];
    display('Too many tapers requested.'); return
end

V = V(1:nTapers);
E = E(:,1:nTapers);

if isempty(NFFT) NFFT = 2^nextpow2(N);end %zero padding 
if RemoveEnsembleMean
    MX=mean(X')';MY=mean(Y')';
    X=X-repmat(MX,1,ntrials);Y=Y-repmat(MY,1,ntrials);
end

Pxx=zeros(NFFT,1);Pyy=Pxx;Pxy=Pxx;Cxy=Pxx;

for trial=1:ntrials
    
    if removeTemporalMean
        x=X(:,trial)-mean(X(:,trial));
        y=Y(:,trial)-mean(Y(:,trial));
    else
        x=X(:,trial);
        y=Y(:,trial);
    end
    % Compute the windowed DFTs.
    if N<=NFFT
        xf=fft(E.*x(:,ones(1,nTapers)),NFFT);
        yf=fft(E.*y(:,ones(1,nTapers)),NFFT);
    else  % Wrap the data modulo nfft if N > nfft
        % use CZT to compute DFT on nfft evenly spaced samples around the 
        % unit circle:
        xf=czt(E.*x(:,ones(1,nTapers)),NFFT);
        yf=czt(E.*y(:,ones(1,nTapers)),NFFT);
    end 
    wt = ones(nTapers,1);
    % Park estimate
    %wt=V(:);
    
    % Xx2 = (xf.*conj(xf))*wt/nTapers
    Xx2 = (abs(xf).^2)*wt/nTapers;%average across tappers
    Yy2 = (abs(yf).^2)*wt/nTapers;
    Xy2 = (yf.*conj(xf))*wt/nTapers;
    
    Pxx = Pxx + Xx2;
    Pyy = Pyy + Yy2;
    Pxy = Pxy + Xy2;
end

Pxx=Pxx/ntrials;
Pyy=Pyy/ntrials;
Pxy=Pxy/ntrials;

XYphi=angle(Pxy);

Cxy = abs(Pxy).^2./(Pxx.*Pyy);% squared coherence function estimate 

% Select first half
% factor 2 to get one-sided
% In the periodogram, a factor 2 is not needed for dc. 
% However, in the case of the multitaper, the dc is within a larger
% bandwidth, so we use a factor 2 in this case too.

if imag(X)==0 & imag(Y)==0,   % Real X and Y
    if rem(NFFT,2),    % nfft odd
        select = [1:(NFFT+1)/2];nf=length(select);
        df=Nyquist/nf;        
          
%         Pxx_unscaled = Pxx(select,:); % Take only [0,pi] or [0,pi)
%         Pxx = [Pxx_unscaled(1,:); 2*Pxx_unscaled(2:end,:)];  % Only DC is a unique point and doesn't get doubled
% 
%         Pyy_unscaled = Pyy(select,:); % Take only [0,pi] or [0,pi)
%         Pyy = [Pyy_unscaled(1,:); 2*Pyy_unscaled(2:end,:)];  % Only DC is a unique point and doesn't get doubled        
%         
%         Pxy_unscaled = Pxy(select,:); % Take only [0,pi] or [0,pi)
%         Pxy = [Pxy_unscaled(1,:); 2*Pxy_unscaled(2:end,:)];  % Only DC is a unique point and doesn't get doubled        
    else % Even
        select = [1:NFFT/2+1];   % include DC AND Nyquist
        nf=length(select);     
        
%         Pxx_unscaled = Pxx(select,:); % Take only [0,pi] or [0,pi)
%         Pxx = [Pxx_unscaled(1,:); 2*Pxx_unscaled(2:end-1,:); Pxx_unscaled(end,:)]; % Don't double unique Nyquist point        
%         
%         Pyy_unscaled = Pyy(select,:); % Take only [0,pi] or [0,pi)
%         Pyy = [Pyy_unscaled(1,:); 2*Pyy_unscaled(2:end-1,:); Pyy_unscaled(end,:)]; % Don't double unique Nyquist point        
%         
%         Pxy_unscaled = Pxy(select,:); % Take only [0,pi] or [0,pi)
%         Pxy = [Pxy_unscaled(1,:); 2*Pxy_unscaled(2:end-1,:); Pxy_unscaled(end,:)]; % Don't double unique Nyquist point        
    end
    % DC component, F0=lowest frequency given frequency resolution (depends on
    % the length of the realization), up to the Nyquist frequency
    F=[0 linspace(F0,Nyquist,length(select)-1)];
else % Complex X and Y
    select = 1:NFFT;nf=length(select);
    F=[0 linspace(F0,2*Nyquist,nf-1)];
end

Pxx=Pxx(select,:);
Pyy=Pyy(select,:);
Pxy=Pxy(select,:);
Cxy=Cxy(select,:);

% Compute the PSD [Power/freq]
% Power = energy in a give time period (length of a realization or trial), 
% thus divide by T; The power spectrum density has units of (signal units)^/Hz;
% This is equivalent to psd = 2 * dt^2 * 1/T * Pxx;
% The multiplication by 2 gives the one-sided spectrum
% The multiplication by dt^2 is done because Matlab's definition of DFT
% does not include the factor dt. By computing the power, that becomes dt^2
%
%
Pxx = 2 * dt^2 * 1/T * Pxx; % Divide by the window length
Pyy = 2 * dt^2 * 1/T * Pyy;
Pxy = 2 * dt^2 * 1/T * Pxy;
      



return