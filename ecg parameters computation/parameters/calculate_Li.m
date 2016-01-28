function Li = calculate_Li(signal,fs,verbose)

% Feature from
% "An Algorithm Used for Ventricular Fibrillation
% Detection Without Interrupting Chest Compression", Li et al.
%
% INPUT:
% - signal: filtered ecg [mV]
% - fs: sampling frequency [sample/s]
% - verbose: flag variable, for debugging purposes
%
% OUTPUT:
% - Li: mean residues (modified from original) 
%
% by Unai Ayala and Unai Irusta (unai.irusta@ehu.es)

signal = signal(:);

% Parameters
nMax=6; t_win=0.3; scale=8; wave='db5';

% Transformada wavelet
cwt_ecg = cwt(signal,scale,wave);
tfm_ecg = abs(cwt_ecg).^2;

% Calculamos picos y seleccionamos nMax mayores
% Debe poder seleccionarse intervalo (+,-)t_win*fs alrededor
[pks, locs] = findpeaks(tfm_ecg(t_win*fs+1:end-t_win*fs),'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',round(0.15*fs));
[~  , idx ] = sort(pks,'descend');
if length(pks) < nMax; nMax=length(pks); disp('hello world');end
locs        = locs(idx(1:nMax))+t_win*fs;

% Calculamos el template
for i = 1:nMax
    interval  = (-t_win*fs+1:t_win*fs)+locs(i);
    cw_j(i,:) = cwt_ecg(interval);
end
tw_j = mean(cw_j,1);

% Calculamos las correlaciones y residuos
xx_tw = xcorr(tw_j); xx_tw = xx_tw/max(xx_tw);
for i = 1:nMax
    xy_tw(i,:) = xcorr(tw_j, cw_j(i,:)); 
    xy_tw(i,:)=xy_tw(i,:)/max(xy_tw(i,:));
    res(i) = sum(abs(xx_tw-xy_tw(i,:)));
end

%%% Modified from original
Li = prctile(res,50);


if verbose
    figure;
    t = (0:length(signal)-1)/fs;
    subplot(211); plot(t, signal);
    subplot(212); plot(t, tfm_ecg); hold on;
    for i=1:length(locs)
        plot([1 1]*locs(i)/fs, get(gca,'Ylim'),':r');
    end
    
    keyboard;
end