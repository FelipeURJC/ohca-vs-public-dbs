function Na = calculate_TCSC(xf,fs,wL,verbose)

% THRESHOLD CROSSING SAMPLE COUNT parameter, based on:
%
% 1) "A simple time domain algorithm for the detection of ventricular
%    fibrillation in electrocardiogram"
%    Arafat M, Chowdhury A and Hasan M, 
%    Signal, Image and Video Processing, 5: 1-10, 2011.
%
%
% INPUT:
% - xf: ecg signal (preprocessed)
% - fs: sampling frequency
% - wL: window length, in seconds. 
% - verbose: debugging variable (1: plot; 0: default, not ploting)
%
% OUTPUT
% - tcsc parameter
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

xf = xf(:);

Ls = 3;
Le = wL;

wls = Ls*fs;
w   = zeros(1,Ls*fs);

tt = linspace(0,Le,Le*fs);
t  = linspace(0,Ls,Ls*fs);
ia = find( (t<0.25) );
ib = (t>=0.25) & (t <= (Ls-0.25)) ;
ic = find( (t>Ls-0.25) );

% 1/4 segundos
w(ia) = 0.5*( 1-cos(4*pi*t(ia)) );
w(ib) = 1;
w(ic) = 0.5*( 1-cos(4*pi*t(ic)) );


V0 = 0.2;

N = zeros(1,Le-2); 

for j=0:Le-3 
    
    wsamples = j*1*fs+1:(j+3)*1*fs;
    
    wecg = xf(wsamples)'.*w;
    
    maxv = max( abs(wecg) );
    wecg = wecg./maxv;
    
    becg = zeros(size(wecg));
    becg(abs(wecg) > V0 ) = 1;
    
    N(j+1) = sum(becg)*100/wls;
    
    %Plot data
    if verbose
        figure(1)
        
        time = tt(wsamples);
        
        plot(time,becg,'r'), hold on, 
        plot(time,wecg),
        plot(time,V0*ones(size(wecg)),'g')
        plot(time,-V0*ones(size(wecg)),'g')
        
        axis([time(1) time(end) -1.1 1.1])
        
        msg = sprintf('N = %2.2f',N(j+1));
        title(msg)
        keyboard;
        clf;
    end
    
end

Na = mean(N);

