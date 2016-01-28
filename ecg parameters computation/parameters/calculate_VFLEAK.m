function vf = calculate_VFLEAK(senal)

% VF LEAK parameter, based on:
%
% 1) "Computer detection of ventricular fibrillation"
%    S. Kuo, D. Dillman, Computers in Cardiology, 1978, 2747-2750.
%
% 2) "Reliability of old and new ventricular fibrillation detection 
%    algorithms for automated external defibrillators"
%    A. Amann, R. Tratning, and K. Unterkofler,
%    Biomed Eng Online, 4(60), 2005.
%
% INPUT:
% - senal: ecg signal (preprocessed)
%
% OUTPUT
% - vfleak parameter
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso


num = sum (abs(senal(2:end)));
den = sum (abs( senal(2:end) - senal(1:end-1) ));

N = floor ( ( pi*(num)/(den) ) + 1/2 );

num = sum( abs( senal(N+1:end) + senal(1:end-N)) );
den = sum( abs(senal(N+1:end)) + abs(senal(1:end-N)));

vf = num/den;


