function [x3,x4,x5] = calculate_Xi(signal,fs,L,verbose)

% Slope and Frequency domain Features, from
% "A Reliable Method for Rhythm Analysis During Cardiopulmonary 
% Resuscitation", Ayala U. et al.
%
% INPUT:
% - signal: filtered ecg [mV]
% - fs: sampling frequency [sample/s]
% - L: analysis window length [s] 
% - verbose: flag variable, for debugging purposes
%
% OUTPUT:
% - x3, x4, x5: frequency domain features 
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

signal = signal(:);

%frequency domain features
Nfft = 2048;
Sfft = fftshift(fft(signal.*hamming(L*fs),Nfft));

ff   = linspace(-fs/2,fs/2,Nfft);
Pss = abs(Sfft).^2;

df   = ff(2)-ff(1);
area = sum(Pss*df)/2;   %f>0
Pss  = Pss./area;       %normalization to unit Pss

ind = find( (ff>=0) & (ff<=30) ); 
f = ff(ind);            % just consider 0<f<30 Hz
Pss = Pss(ind);

% frequency ranges
f_1_10  = find(f>=1 & f<=10);
f_2_7   = (f>=2.5 & f<=7.5);
f_12_30 = (f>=12 & f<=30);

[val,pos] = max(Pss(f_1_10));
x3 = f(f_1_10(pos));
x4 = sum(Pss(f_2_7)*df);
x5 = sum(Pss(f_12_30)*df);

if verbose
    
    t = (0:L*fs-1)/fs;
    
    plot(f,Pss,x3,val,'ro');
    msg = sprintf('x_4 = %2.2f\nx_5 = %2.2f',x4,x5);
    text(20,0.5,msg,'FontSize',14)
    xlabel('f (Hz)','Interpreter','Latex'); ylabel('$P_{ss}$','Interpreter','Latex')
    
    keyboard;
end




