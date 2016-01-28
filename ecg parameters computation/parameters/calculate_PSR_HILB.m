function [psr,hilb] = calculate_PSR_HILB(signal,fs)

% PHASE SPACE RECONSTRUCTION and HILBERT TRANSFORM parameter, based on:
%
% 1) "Detecting ventricular fibrillation by time-delay methods"
%    A. Amann, R. Tratning, and K. Unterkofler,
%    IEEE Trans Biomed Eng, 54(1): 174-77, 2007.
%
% 2) "A new ventricular fibrillation algorithm for automated external 
%    defibrillators"
%    A. Amann, R. Tratning, and K. Unterkofler,
%    Computers in Cardiology, 559-562, 2005.
%
% INPUT:
% - signal: ecg signal (preprocessed)
% - fs: sampling frequency
%
% OUTPUT
% - psr, hilb parameters
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso


%first downsampling ECG to 50 Hz;
fd = 50;
M = round(fs/fd);

signal = signal(1:M:end);
fs = fd;

%--- HILB
X = hilbert(signal);

xr = real(X); % Señal x(t)
xi = imag(X); % Señal TH[x(t)]

xr = xr - min(xr); % Restamos porque la primera casilla es (0,0)
xi = xi - min(xi);

deltaxr = max(xr)/39; % El "paso" para normalizar
deltaxi = max(xi)/39;

xr = floor(xr/deltaxr) + 1; % Sale entre 1 y 40 => OK
xi = floor(xi/deltaxi) + 1; % Sale entre 1 y 40 => OK
%--- 

%--- PSR
signal = xr;

tau  = 0.5; % delay (seconds)
ntau = tau*fs; % delay (samples)
N    = length(signal);

nt     = 1:(N-ntau);
ntplus = ntau+1:N; 

xt  = signal(nt);
xtp = signal(ntplus);
%---


a = find(xt>40);
b = find(xtp>40);
c = find(xi>40);
d = find(xr>40);

if ~isempty(a) 
    xt(a)  = 40;
end
if ~isempty(b)
    xtp(b) = 40;
end
if ~isempty(c)
    xi(c) = 40;
end
if ~isempty(d)
    xr(d)  = 40;
end

A = zeros(40,40);
I = sub2ind(size(A),xt,xtp);
A(I) = 1;
 
B = zeros(40,40);
I = sub2ind(size(B),xr,xi);
B(I) = 1;

psr  = sum(sum(A))/1600;
hilb = sum(sum(B))/1600;
