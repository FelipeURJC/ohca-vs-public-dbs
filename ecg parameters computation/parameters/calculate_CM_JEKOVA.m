function parameters = calculate_CM_JEKOVA(signal,fs,wL)

% COMPLEXITY MEASURE parameters, based on: 
%
% 1) "Detecting ventricular tachycardia and fibrillation by complexity
%    parameter"
%    X.S. Zhang, Y.S. Zhu, N.V. Thakor and Z.Z. Wang
%    IEEE Trans. Biomed. Eng. 46(5): 548-55, 1999    
%
% 2) "Reliability of old and new ventricular fibrillation detection 
%    algorithms for automated external defibrillators"
%    A. Amann, R. Tratning, and K. Unterkofler,
%    Biomed Eng Online, 4(60), 2005.
%
% COVARIANCE, FREQUENCY, AREA, and KURTOSIS parameters, based on:
%
% 3) "Shock advisory tool: detection of life-threatening cardiac
%    arrhyrhmias and shock success prediction by means of a common
%    parameter set"
%    I. Jekova, Biomed. Sig. Proc. Control, 2:25-33, 2007
% 
% 4) "Ventricular Fibrillation and Tachycardia Classification uning a
%    Machine Learning Approach", 
%    Q. Li, C. Rajagopalan and G.D. Clifford,
%    IEEE Trans. Biomed. Eng. (In Press)
%
% INPUT:
% - signal: ecg signal (preprocessed)
% - wL: window length, in seconds 
%
% OUTPUT
% - parameters: [1x5] vector containing: cm cvbin, frqbin, abin, 
%                                        and kurt paramters.
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es) and
%    Eduardo Morgado (eduardo.morgado@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

%--- Binary signal
n = length(signal);
s = zeros(1,n);

xi = signal - mean(signal);

xmax = max(xi); 
xmin = min(xi);

Pc = length(find( 0 < xi < 0.1*xmax ) );
Nc = length(find( 0.1*xmin < xi < 0 ) );

if (Pc+Nc)<0.4*n
    th = 0;
elseif Pc < Nc
    th = 0.2*xmax;
else
    th = 0.2*xmin;
end

s(xi>=th) = 1;

% Complexity measure
[~,cm] = kolmogorov(s);

% Convariance calculation
cvbin = var(s);

% Frequency calculation
frqbin = sum(diff(s) == 1);    
frqbin = frqbin / wL;

% Area calculation
N = sum(s);
abin = max(N,fs*wL-N);

% Kurtosis calculation
m_s = mean(signal);
s_s = std(signal);
kurt = mean((signal - m_s).^4)/s_s^4 - 3;

parameters = [cm,cvbin,frqbin,abin,kurt];

