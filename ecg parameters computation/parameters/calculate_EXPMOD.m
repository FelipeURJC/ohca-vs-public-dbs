function N = calculate_EXPMOD(signal,fs,verbose)

% MODIFIED EXPONENTIAL parameter, based on:
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
% - expmod parameter
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

signal = signal(:);

%normalization
maximo = max(signal);
signal = signal./maximo;

%peak detection
nm    = peakdet(signal,0.2); 
mxpos = nm(:,1);

if mxpos(1)==1
    mxpos = mxpos(2:end);
end

% tau = 0.2s
tau = 0.2;
L = length(signal);

if isempty(mxpos)
    N = 1*60/(L/fs);
    return; 
end

% initialization, until the first maxima
nm = mxpos;
nm1 = nm(1);
En = zeros(L,1);
En(1:nm1-1) = signal(1:nm1-1);
maxpos = nm1;


% remaining samples
n = nm1:L;
%offset = nm1-1;
nmj = nm1;
%i = 1;
N = 1;
fin = 0;

%th = fs/5;

while ~fin
    
    Mj = signal(nmj);
    En(n) = Mj*exp(-(n-nmj)./(tau*fs));    
       
    ncj = find(signal > En );
    val = find (ncj-nmj > 10);
    ncj = ncj(val);        
    
    if isempty(ncj)
        fin=1;
    else
        ncj = ncj(1);
        aux = signal; aux(1:ncj-1) = zeros(1,ncj-1);
        nm = peakdet(aux,0.3);
           
        if isempty(nm)
           fin = 1;
           En(ncj:L) = signal(ncj:L);
        else
            
            nmj = nm(nm>ncj,1);
            
            if isempty(nmj)
                fin = 1;
                En(ncj:L) = signal(ncj:L);
            else
                
                
                nmj = nm(nm>ncj,1);
                nmj = nmj(1);
                En(ncj:nmj-1) = signal(ncj:nmj-1);
                n = nmj:L;
                %offset = nmj-1;
                
                N = N +1;
                
                maxpos = [maxpos nmj];
                               
            end
            
            
        end
    end
    
end


N = N*60/(L/fs);

%plot data
if verbose
    figure(1);
    plot(signal); hold on;
    plot(maxpos,signal(maxpos),'rx');
    plot(En,'g');
    msg = sprintf('MEA value: %2.2f',N);
    title(msg)
    keyboard;
    clf
end


