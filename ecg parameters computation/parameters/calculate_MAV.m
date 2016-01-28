function mav = calculate_MAV(signal,fs)

% MEAN ABSOLUTE VALUE parameter, based on:
%
% 1) "Sequential algorithm for life threatening cardiac pathologies 
%    detection based on mean signal strength and EMD functions"
%    E. Anas, S. Lee, and M. Hasan,
%    Biomed Eng Online, 9(1), 2010.
%
% INPUT:
% - signal: ecg signal (preprocessed)
% - fs: sampling frequency
%
% OUTPUT
% - mav parameter
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso


L = length(signal);
T = L/fs;           
N = T-1;

Le = 2;
Ne = (Le-1)*fs;

mavi = zeros(1,N);

for j=0:N-1
    xa = signal(j*Ne+1:(j+2)*Ne);
    xa = xa./max(abs(xa));
    mavi(j+1) = (1/Ne)*sum( abs(xa) );
end
mav = mean(mavi);  

