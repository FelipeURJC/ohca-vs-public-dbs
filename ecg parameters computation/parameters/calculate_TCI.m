function tci = calculate_TCI(xf,fs,wL,verbose)

% THRESHOLD CROSSING INTERVAL parameter, based on:
%
% 1) "Ventricular tachycardia and fibrillation detection by a sequential 
%    hypothesis testing algorithm"
%    N.V Thakor, Y.S. Zhu, and K.Y. Pan, 
%    IEEE Trans Biomed Eng, 37(9): 837-843, 1990.
%
% 2) "Reliability of old and new ventricular fibrillation detection 
%    algorithms for automated external defibrillators"
%    A. Amann, R. Tratning, and K. Unterkofler,
%    Biomed Eng Online, 4(60), 2005.
%
% INPUT:
% - xf: ecg signal (preprocessed)
% - fs: sampling frequency
% - wL: window length, in seconds 
% - verbose: debugging variable (1: plot; 0: default, not ploting)
%
% OUTPUT
% - tci parameter
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

wl = 1*fs;  % 1-sec window samples
wa = 3*fs;  % 3-sec window samples


L = wL-3+1; % number of 3-sec windows in wl segment 

%becg = zeros(1,wa);
becg1 = zeros(1,wl);
becg2 = zeros(1,wl);
becg3 = zeros(1,wl);

tci6 = zeros(1,L);

for j=0:L-1
    
    wsamples1 = j*wl+1:(j+1)*wl;
    wsamples2 = (j+1)*wl+1:(j+2)*wl;
    wsamples3 = (j+2)*wl+1:(j+3)*wl;
    
    stage1 = xf(wsamples1)-mean(xf(wsamples1)); maxv = max(stage1); 
    th1 = 0.2*maxv; becg1(stage1>th1) = 1;
    stage2 = xf(wsamples2)-mean(xf(wsamples2)); maxv = max(stage2); 
    th2 = 0.2*maxv; becg2(stage2>th2) = 1;
    stage3 = xf(wsamples3)-mean(xf(wsamples3)); maxv = max(stage3); 
    th3 = 0.2*maxv; becg3(stage3>th3) = 1;
    
    becg = [becg1 becg2 becg3];
    
    aux = [0 diff(becg)];
    
    s1 = find(aux(1:wl)==-1);
    if isempty(s1)
        t1 = 1;
    else
        t1 = (wl-s1(end))/fs;
    end
    
    index = find(aux(wl+1:2*wl));
    s2 = aux(wl+1:2*wl);
    pulses = s2(index);
    
    if pulses(1) == -1 && pulses(end) == -1
        t2 = 0;
        t3 = (wl-index(end))/fs;
        N = (length(pulses)+1)/2;
        %disp('caso 1')
        %flag = 1;
    elseif pulses(1) == 1 && pulses(end) == 1
        t2 = index(1)/fs;
        t3 = 0;
        N = (length(pulses)+1)/2;
        %disp('caso 2')
        %flag = 1;
    elseif pulses(1) == -1 && pulses(end) == 1
        t2 = 0;
        t3 = 0;
        N = (length(pulses)+2)/2;
        %disp('caso 3')
        %flag = 1;
    elseif pulses(1) == 1 && pulses(end) == -1
        t2 = index(1)/fs;
        t3 = (wl-index(end))/fs;
        N = (length(pulses))/2;
        %disp('caso normal')
    else
        disp('This should not be happening!')
        keyboard; %better to debug
    end
       
    s4 = find(aux(2*wl+1:3*wl)==1);
    if isempty(s4)
        t4 = 1;
    else
        t4 = s4(1)/fs;
    end
        
    tci6(j+1) = 1000/((N-1)+(t2/(t1+t2))+(t3/(t3+t4)));
    
    
    %Plot data
    if verbose
        
        f = figure;
        t = [wsamples1 wsamples2 wsamples3]/fs;
        stage = [stage1' stage2' stage3'];
        
        subplot(211)        
        plot(t,stage); hold on;
        plot(t,[th1*ones(1,wl) th2*ones(1,wl) th3*ones(1,wl)],'r');
        xlabel('time');
        ylabel('ECG and threshold');
        
        subplot(212);
        plot(t,becg,'k'); hold on; stem(t,aux,'r'); hold on;
        line([j+1 j+1],[-1.2 1.2]);line([j+2 j+2],[-1.2 1.2]);
        axis([t(1) t(end) -1.2 1.2])
        xlabel('times');
        ylabel('pulses');
        msg = sprintf('t_1=%2.2f\t\t t_2=%2.2f\t\t t_3=%2.2f\t\t t_4=%2.2f',...
            t1,t2,t3,t4);
        title(msg)
        text(j+1.2,-0.5,['TCI = ' num2str(tci6(j+1))])
        
        
        hold off
        keyboard;
        close(f);
    end
    
    %becg = zeros(1,wa);
    becg1 = zeros(1,wa/3);
    becg2 = zeros(1,wa/3);
    becg3 = zeros(1,wa/3);
    
end

tci = mean(tci6);