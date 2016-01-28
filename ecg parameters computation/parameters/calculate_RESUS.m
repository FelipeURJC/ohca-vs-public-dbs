function [bCP, bWT, bW] = calculate_RESUS(signal,fs)

% Slope and Frequency domain Features, from
% "A high-temporal resolution algorithm to discriminate shockable
% from nonshockable rhythms in adults and children", Irusta U. et al.
%
% INPUT:
% - signal: filtered ecg [mV]
% - fs: sampling frequency [sample/s]
%
% OUTPUT:
% - bCP, bWT, bW: features 
%
% by Unai Irusta and Unai Ayala (unai.irusta@ehu.es)

signal = signal(:);

L     = length(signal);
L_w   = 2*fs;   %%% Usamos subventanas de 2s en bCP y bWT
valC  = 0.0055; %%% Para bCP
valPT = 47;  %%% Para bWT

%frequency domain feature
Nfft = 4096;
Sfft = fft(signal.*hamming(L),Nfft);

ff   = fs/2*linspace(0,1,Nfft/2+1);
Sfft = Sfft(1:length(ff));
%ff   = linspace(0,fs/2,Nfft);
Pss = abs(Sfft).^2;

%%% Suppress components and normalize
f_max = 40; f_min = 0.75;
Pss(ff > f_max | ff < f_min)=0;
Pss = Pss/sum(Pss);

%%% Calculamos ancho banda entre 20-80 potencia
P_sum = cumsum(Pss);
f_min = ff(find(P_sum < 0.2,1,'last'));
f_max = ff(find(P_sum > 0.8,1,'first'));

bW = f_max - f_min;

%%% Slope domain feature
slope = diff(signal);
slope = slope.^2; slope(end+1)=slope(end);

n_w = floor(length(slope)/L_w);
for i=1:n_w
    slp_w = slope((i-1)*L_w+1:i*L_w);
    bCP(i) = sum(slp_w < valC*max(slp_w))/L_w;
end
bCP = min(bCP);

%%% Time domain feature (baseline)
tr_sec = 2;
signal   = [ones(1,tr_sec*fs)*signal(1) signal(:)' ones(1,tr_sec*fs)*signal(end)];
[b, a]   = butter(5,2*[6.5 30]/fs,'bandpass');
signal   = filter(b,a,signal);
signal   = signal(tr_sec*fs+1:end-tr_sec*fs);

n_w = floor(length(signal)/L_w);
for i=1:n_w
    sig_w  = signal((i-1)*L_w+1:i*L_w);
    bWT(i) = (prctile(sig_w,50+valPT/2)-prctile(sig_w,50-valPT/2))/max(abs(sig_w));
end
bWT = max(bWT);