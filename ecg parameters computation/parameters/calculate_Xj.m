function [x1,x2] = calculate_Xj(signal,fs)

% Slope and Frequency domain Features, from
% "A Reliable Method for Rhythm Analysis During Cardiopulmonary 
% Resuscitation", Ayala U. et al.
%
% INPUT:
% - signal: filtered ecg [mV]
% - fs: sampling frequency [sample/s]
%
% OUTPUT:
% - x1, x2: slope domain features 
%
% by Unai Ayala and Unai Irusts (unai.irusta@ehu.es)

signal = signal(:);
% Promediador y cálculo de la pendiente
aMs=100;    % Promediamos 100ms
bP=ones(1,round(aMs*fs/1000));

%%% Calculo pendiente al cuadrado
slope = diff(signal); slope(end+1)=slope(end);
slope = slope.^2;

%%% Promediado de la pendiente (eliminando transitorio)
tr_sec  = 2; %%% 2 segundos de transitorio
slope   = [ones(1,tr_sec*fs)*slope(1) slope(:)' ones(1,tr_sec*fs)*slope(end)]; 
slope   = filter(bP,1,slope);
slope   = slope(tr_sec*fs+1:end-tr_sec*fs);

%%% ensbCP es el percentil 10
prC=10;
x1 = prctile(slope,prC)/max(slope);

% nPA Modificado
minAmp     = 0.1*max(slope);  % 10% del normalizado
minDist    = fs/10; 
[pks, locs]= findpeaks(slope,'MINPEAKHEIGHT',minAmp,'MINPEAKDISTANCE',minDist);
x2         = length(pks);