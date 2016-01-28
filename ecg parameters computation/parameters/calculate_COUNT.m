function [c1,c2,c3b] = calculate_COUNT(signal,fs,L,verbose)

% Band-pass filtering and auxiliary counts, based on 
% "Real time detection of ventricular fibrillation and tachycardia"
% I. Jekova, V. Krasteva, Physiol. Meas. 25, 1167-1178, 2004.
%
% INPUT:
% - signal: filtered ecg
% - fs: sampling frequency
% - L: window length 
% - verbose: flag variable, for debugging purposes
%
% OUTPUT
% - c1, c2 and c3: calculated count1, count2 and count3 parameters
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

% filter coefs definition, according to paper
B = [0.5 0 -0.5];
A = [8 -14 7];

% FSi signal
Fsignal = filtfilt(B,A,signal); 
Fsignal_Abs = abs(Fsignal);

Count1 = zeros(1,L);
Count2 = zeros(1,L);
Count3 = zeros(1,L);
for count = 0:L-1   
    Intervalo_Fsignal_Abs = Fsignal_Abs(count*fs+1:count*fs+fs);
    max_Fsignal = max(Intervalo_Fsignal_Abs);
    mean_Fsignal = mean(Intervalo_Fsignal_Abs);
    md_Fsignal = mean(abs(Intervalo_Fsignal_Abs - mean_Fsignal));
    
    Count1(count+1) = sum(Intervalo_Fsignal_Abs >= 0.5*max_Fsignal); 
    Count2(count+1) = sum(Intervalo_Fsignal_Abs >= mean_Fsignal);
    Count3(count+1) = sum((Intervalo_Fsignal_Abs >= (mean_Fsignal - md_Fsignal)) .* (Intervalo_Fsignal_Abs <= (mean_Fsignal + md_Fsignal)));
end
Count1 = sum(Count1);
Count2 = sum(Count2);
Count3 = sum(Count3);

% MODIFIED FROM ORIGINAL: since c1, c2 and c3 depends on the window length
% (L)they are normalized by L

c1 = Count1/L;
c2 = Count2/L;
c3 = Count3/L;

c3b = c1*c2/c3; % not interested in c3, 
                % but on c1*c2/c3 (see figure 4 of the orignal paper).

% Plot for debugging
if verbose
    
    % replicating figure 2 of the paper
    [H,f] = freqz(B,A,0:fs/1000:50,fs);
    figure;
    mg2 = abs(H).^2;
    plot(f,mg2/max(mg2));
    xlabel('frequency (Hz)');
    ylabel('Magnitude Squared');
    grid on;
    keyboard;
end





