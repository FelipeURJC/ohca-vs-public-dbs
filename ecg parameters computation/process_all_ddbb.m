function process_all_ddbb

% This function compute the ECG for the processed ddbbs
% contained in path_ddbb. Results are stored in path_res
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

clc, close all;

path_ddbb = '../data/preprocessed ddbbs/';
path_res  = '../data/computed parameters/';

db_names = {'vfdb','cudb', 'ahadb'};
w_length = [4, 8];

for i = 1:length(db_names)
    for j = 1:length(w_length)
        r_filename  = sprintf('%s%s_%d', path_res , db_names{i} ,w_length(j));
        db_filename = sprintf('%s%s'   , path_ddbb, db_names{i});
        calculate_parameters(db_filename, w_length(j), r_filename);
    end
end

end

function [X,y,names,Lr,ECG] = calculate_parameters(db_filename, wlength,r_filename)

% Function for calculating a number of parameters from the preprocessed
% ECG database
%
% INPUT:
% - db_filename: name of the database to be loaded {cudb,vfdb,ahadb}
% - wlength: analysis window length, in seconds.
%
% OUPUT:
% - X: [N,K] matrix of computed parameters: N number of samples, K number
% of parameters.
% - y: [Nx1] matrix of labels
% - names: names of the database records.
% - Lr: number of samples for each episode
% - ECG: for each window 
%
% LABELS have been assigned according to this code:       
%         cases = {...
%             '(AB';...       %1
%             '(AFIB';...     %2
%             '(AFL';...      %3
%             '(ASYS';...     %4
%             '(B';...        %5
%             '(BI';...       %6
%             '(BII';...      %7
%             '(HGEA';...     %8
%             '(IVR';...      %9
%             '(N';...        %10
%             '(NOD';...      %11
%             '(NOISE';...    %12
%             '(P';...        %13
%             '(PREX';...     %14
%             '(SBR';...      %15
%             '(SVTA';...     %16
%             '(T';...        %17
%             '(VER';...      %18
%             '(VF';...       %19
%             '(VFL';...      %20
%             '(VT';...       %21
%             '(sTV';...      %22
%             '(others';...   %23
%             '(fineVF'};     %24
%
% by Felipe Alonso-Atienza (felipe.alonso@urjc.es)
% www.tsc.urjc.es/~felipe.alonso

% Modified (march 2014) by Grupo GSC
% changes: (1) names calculated different (some registers are too short <4s)
%          (2) included train and ped variables to identify pediatric and train/test dbs
%          (3) save results in filename
% contact:
%               Unai Irusta (unai.irusta@ehu.es)
%               Unai Ayala  (unai.ayala@ehu.es)


% Load data
load(db_filename);
[~, db_name, ~] = fileparts(db_filename);
eval(sprintf('database=%s',db_name));

L = length(database(:));    % number of database records

names = {};                 % Not all records are included (depends on wlength).
fs = database.fs;           % sampling frequency

% variables definition and initialization
wsamples = wlength*fs;      

X  = [];         
y  = [];
Lr = [];
ECG = [];

addpath('./parameters')
disp('Processing database...')
for n=1:L
        
        disp(database(n).name)
        
        ecg_signal = database(n).ecg;   
        ecg_labels = database(n).label; % as many labels as ECG samples
      
        % For each ECG signal record, we construct Nw non-overlapping
        % segments 
        Nw     = floor(length(ecg_signal)/wsamples); % number of windows
        
        if Nw > 0  % Small change ahadb has short ECG samples (UIZ)
            
            [X_i,label,ecg_w] = calculate_parameters_for_each_record(ecg_signal,...
                ecg_labels,Nw,wsamples,wlength,fs);
            
            names{end+1} = database(n).name;
            
            X  = [X;X_i];
            y  = [y;label'];
            
            ECG = [ECG; ecg_w];
            
            Lr = [Lr;length(label)];
            
        end
       
        
end

%X = X'; % [N x K]
%y = y'; % [N x 1]

names = names';

if nargin == 3
    save(r_filename, 'X','y','names','Lr','ECG');
end

end

function [A, label,ECG] = calculate_parameters_for_each_record(ecg_signal,...
    ecg_labels,Nw,wsamples,wlength,fs)
    
k = 1;
verbose = 0;

for j=1:Nw
    
    window = (j-1)*wsamples+1:j*wsamples;
    
    ecg_in_window = ecg_signal(window); % by default the 1st channel
    labels_in_window = ecg_labels(window); 
    
    % ecg segments with more than one rhythm are discarded for
    % further analysis.
    
    %keyboard;
    if not_mixed_rhythms(labels_in_window,j,window)
        
        
        % Threshold Crossing Interval
        tci(k)    = calculate_TCI(ecg_in_window,fs,wlength,verbose);
        
        % Threshold Crossing Sample Count
        tcsc(k)   = calculate_TCSC(ecg_in_window,fs,wlength,verbose);
        
        % Exponential Parameter
        exp(k)    = calculate_EXP(ecg_in_window,fs,verbose);
        
        % Modified Exponential Parameters
        expmod(k) = calculate_EXPMOD(ecg_in_window,fs,verbose);
        
        % Complexity and Jekova Parameters (based on CM)
        param     = calculate_CM_JEKOVA(ecg_in_window,fs,wlength);
        cm(k)     = param(1); % Complexity Measure
        cvbin(k)  = param(2); % Covariance calculation
        frqbin(k) = param(3); % Frequency calculation
        abin(k)   = param(4); % Area Calculation
        kurt(k)   = param(5); % Kurtosis
        
        % VF-filter Leak
        vfleak(k) = calculate_VFLEAK(ecg_in_window);
        
        % Spectral algorithm
        spec      = calculate_SPEC(ecg_in_window,fs);
        M(k)      = spec(1);
        A1(k)     = spec(2);
        A2(k)     = spec(3);
        A3(k)     = spec(4);
        
        % Mean Absolute Value
        mav(k)    = calculate_MAV(ecg_in_window,fs);
        
        % Phase Space Reconstruction and Hilbert Transform
        [psr(k),hilb(k)] = calculate_PSR_HILB(ecg_in_window,fs);
        
        % Sample Entropy
        SamEn(k) = calculate_SampEn(ecg_in_window);
        
        % Spectral analysis
        [x3(k), x4(k), x5(k)] = ...
            calculate_Xi(ecg_in_window,fs,wlength,verbose);
        
        % New parameters
        [x1(k), x2(k)]          = calculate_Xj   (ecg_in_window,fs);
        [bCP(k), bWT(k), bW(k)] = calculate_RESUS(ecg_in_window,fs);
        Li(k)                   = calculate_Li   (ecg_in_window,fs,verbose);
        
        % Band-pass filtering and auxiliary counts
        [count1(k),count2(k),count3(k)] =...
            calculate_COUNT(ecg_in_window,fs,wlength,verbose);
        
        % LABEL ASSIGNMENT
        label(k) = mode(labels_in_window);
        
        ECG(k,:) = ecg_in_window;
        
        % Updated counter
        k = k + 1;
    end
    
end

% Matrix of features: [samples x Feats]
A = [tci', tcsc', exp', expmod', cm', cvbin', frqbin', abin', kurt',...
    vfleak', M', A1', A2', A3', mav', psr', hilb', SamEn', x1', x2', x3',...
    x4', x5', bCP', bWT', bW', Li', count1', count2', count3'];
end

function output = not_mixed_rhythms(labels,j,window)

% This function examines whether the window under analysis contains both
% Shock and NShock rhythms. In this case, the window is discarded

% Also, the window is discarded if it contains rhythms labelled as slow TV 
% (sTV, #22), fine FV (#23), noise (#12) and asystolia (#4).

blabels = binary_labels(labels,0);  % {+1 Shock, 0 NShock}.
unique_labels = unique(blabels);

if ( length(unique_labels) == 1 ) ...
        && ( sum(ismember(labels,4))== 0 ) ...
        && ( sum(ismember(labels,12))== 0 ) ...
        && ( sum(ismember(labels,22))== 0 ) ...
        && ( sum(ismember(labels,24))== 0 )
    output = true;
else
    output = false;
    %msg = sprintf('Window %d discarded, samples: (%d,%d)',...
    %    j,window(1),window(end));
    %disp(msg);
end

end

function y = binary_labels(labels,value)

% Binary label assignment: Shock (labels 19 20 21) vs NonShock rhythms
shock  = find( labels == 19 | labels == 20 | labels == 21 );
nshock = setdiff(1:length(labels),shock);

y = value*ones(size(labels));
y(shock) = 1;

end

