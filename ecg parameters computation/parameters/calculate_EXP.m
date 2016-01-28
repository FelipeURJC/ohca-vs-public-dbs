function N = calculate_EXP(signal,fs,verbose)

% EXPONENTIAL parameter, based on:
%
% 1) "Reliability of old and new ventricular fibrillation detection 
%    algorithms for automated external defibrillators"
%    A. Amann, R. Tratning, and K. Unterkofler,
%    Biomed Eng Online, 4(60), 2005.
%
% INPUT:
% - signal: ecg signal (preprocessed)
% - fs: sampling frequency
% - wL: window length, in seconds 
%
% OUTPUT
% - exp parameter
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

signal = signal(:);

% max value
[M,nm] = max(abs(signal)); 

% tau = 3s
tau = 3;
L = length(signal);
n = 1:L;

Es = M.*exp(-abs(nm-n)./(tau*fs));

% Intersections
above   = Es>=signal';
inter   = find(diff(above)~=0);
peaks   = length(inter);

N = peaks*60/(L/fs); %crossing per minutes

if verbose
    figure(1);
    plot(n/fs,signal); hold on;
    plot(n/fs,Es,'r');
    plot(inter/fs,signal(inter),'ro');
    xlabel('time (sec)');
    ylabel('ECG (a.u.)')
    keyboard;
    clf
end


