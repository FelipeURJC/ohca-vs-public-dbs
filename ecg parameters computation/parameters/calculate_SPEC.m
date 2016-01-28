function spec = calculate_SPEC(senal,fs)

% SPECTRAL ALGORITHM, based on:
%
% 1) "Algorithmic sequential decision making in the frequency domain for
%    life threatening ventricular arrythmias and imitative artifacts: a
%    diagnostic system",
%    S. Barro, R. Ruiz, D. Cabello, and J. Mira, 
%    Journal of Biomed Eng., 11(4): 320-8, 1989
%
% 2) "Reliability of old and new ventricular fibrillation detection 
%    algorithms for automated external defibrillators"
%    A. Amann, R. Tratning, and K. Unterkofler,
%    Biomed Eng Online, 4(60), 2005.
%
% INPUT:
% - senal: ecg signal (preprocessed)
% - fs: sampling frequency
%
% OUTPUT
% - M, A1, A2, A3 parameters
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

senal = senal(:);

L       = length(senal);
ventana = hamming(L); 

% Multiply by hamming window and calculate FFT
NFFT = 2^13; 
Y    = fft(senal.*ventana,NFFT);  %senal enventanada
Y    = Y(1:NFFT/2);
f    = linspace(0,fs/2,NFFT/2);

%definition of amplitude spectrum according to the original paper
Amp = abs(real(Y)) + abs(imag(Y));

% insignificant components reduction: those below 5% of the max. ampl. 
iomega = find( (f>=0.5) & (f<=9) );
romega = f(iomega);

[Amax,pos]  = max(Amp(iomega));
th          = 0.05*Amax;
Amp(Amp<th) = 0;

% Parameters
omega   = romega(pos);
jmax    = min(20*omega,100);

% M: 0-min(20*omega,100)
js = find(f<=jmax);
ws = f(js);
aj = Amp(js)';

M = (1/omega)*aj*ws'/sum(aj);

% A1
jsnum = find((f>=0.5)& (f<=omega/2));
ws    = f(jsnum);
ajnum = Amp(jsnum)';

jsden = (f>=0.5)& (f<=jmax);
ajden = Amp(jsden);

A1 = (1/omega)*ajnum*ws'/sum(ajden);

% A2
jsnum = find((f>=0.7*omega)& (f<=1.4*omega));
ws    = f(jsnum);
ajnum = Amp(jsnum)';

A2 = (1/omega)*ajnum*ws'/sum(ajden);

% A3
jsnum = find((f>=2*omega)& (f<=8*omega));
ws    = f(jsnum);
ajnum = Amp(jsnum)';

A3 = (1/omega)*ajnum*ws'/sum(ajden);

% all
spec = [M, A1, A2, A3];




