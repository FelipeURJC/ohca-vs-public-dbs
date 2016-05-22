function table2

% This function provides the information of table 2 of the paper
%
% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

close all; clear all; clc;

data_path = '../data/';

w_length = [4, 8];
k = 1;

% load data and calculate performance: Se, Sp and BER
Se = zeros(30,4);
Sp = zeros(30,4);
BER = zeros(30,4);
for j = 1:length(w_length)

    % load data
    
    filename  = sprintf('%sdata_%d',data_path, w_length(j));
    load(filename);
        
    %Analyze performance of individual features
    public_db = 1:samples_for_dbs.ahadb(end);
    ohca_db   = samples_for_dbs.ohcadb;
    
    Tabla = feature_transformation(Tabla);
    
    % public dbs
    [Se(:,k), Sp(:,k), BER(:,k)] = ...
        get_metrics(Tabla,public_db,'PUBLICs',w_length(j));

    % ohca db   
    [Se(:,k+1), Sp(:,k+1), BER(:,k+1)] =...
        get_metrics(Tabla,ohca_db,'OHCA',w_length(j));
    
    k=k+2;
end

% Sort features by their BER performance
MeanBER = mean(BER,2);
[~,sorted_index] = sort(MeanBER,'ascend');

VarNames  = Tabla.Properties.VariableNames;
FeatNames = VarNames(1:end-5);
Nfeats    = length(FeatNames);

Se = Se(sorted_index,:);
Sp = Sp(sorted_index,:);
FeatNames = FeatNames(sorted_index);

msg = sprintf('\n\tPUBLIC-4   PUBLIC-8   OHCA-4     OHCA-8');
disp(msg)
msg = sprintf('\tSe/sp      Se/Sp      Se/Sp      Se/Sp');
disp(msg)

for i=1:Nfeats;
    
    msg=sprintf('%s\t%2.1f/%2.1f  %2.1f/%2.1f  %2.1f/%2.1f  %2.1f/%2.1f',...
        FeatNames{i},Se(i,1),Sp(i,1),...
        Se(i,3),Sp(i,3),...
        Se(i,2),Sp(i,2),...
        Se(i,4),Sp(i,4));
    disp(msg)
    
end

% Boxplots of Se/Sp values for inidividual features
Se = [Se(:,1) Se(:,3) Se(:,2) Se(:,4)];
Sp = [Sp(:,1) Sp(:,3) Sp(:,2) Sp(:,4)];

ddbb_names = {'public-4s','public-8s','ohca-4s','ohca-8s'};

% sensitivity
figure(1)
boxplot(Se(1:end-2,:),ddbb_names);    % remove outliers
xlabel('databases');
ylabel('Se (%)');

% specificity
figure(2)
boxplot(Sp,ddbb_names);    
xlabel('databases');
ylabel('Sp (%)');


end

function [Se,Sp,BER] = get_metrics(table,samples_db,option,wlength)

VarNames  = table.Properties.VariableNames;
FeatNames = VarNames(1:end-5);
Nfeats    = length(FeatNames);

Se = zeros(Nfeats,1);
Sp = zeros(Nfeats,1);
BER = zeros(Nfeats,1);
for i=1:Nfeats;

    xi = table{samples_db,i};
    y  = table.y(samples_db);
        
    % Build a threshold detector
    [~,~,stats] = bayesdetector(xi,y,FeatNames{i},option,wlength);
        
    Se(i) = stats.sen;
    Sp(i) = stats.esp;
    BER(i) = stats.ber;
end

end

function T = feature_transformation(T)

T.kurt   = log10(T.kurt + 3);
T.M      = log10(T.M + 1);
T.A3     = sqrt(T.A3);
T.count3 = log10(T.count3);
T.x1     = log10(T.x1);
T.x3     = log10(T.x3);
T.x5     = sqrt(T.x5);

end

function [th,predicted_scores,stats] = bayesdetector(x,scores,feat,ddbb,wlength)

x = x(:);

[xmin, pos] = min(x);
xmax = max(x);

N = 100;
t = linspace(xmin,xmax,N);

val = 1;
if scores(pos) == -1;
    val = -1;
end

if strcmp(feat,'tci') || strcmp(feat,'abin')...
            || strcmp(feat,'kurt') || strcmp(feat,'vfleak')...
            || strcmp(feat,'m') || strcmp(feat,'a3')...
            || strcmp(feat,'bcp') || strcmp(feat,'x5'),
    val = 1;
end

if ( strcmp(feat,'x4') && strcmp(ddbb,'PUBLICs') ) 
   val = -1*val;
end

if strcmp(feat,'x2') && strcmp(ddbb,'OHCA')
   val = -1*val;
end

if strcmp(feat,'x3') && strcmp(ddbb,'OHCA') && wlength==8
   val = -1*val;
end

predicted_scores = -1*ones(size(scores)).*val;
ber = zeros(1,N);
for i=1:N
    
    predicted_scores(x<t(i)) = 1*val;
    metrics = compute_metrics(predicted_scores, scores,-1);    
    ber(i)  = metrics.ber;
    
    % re-initialization
    predicted_scores = -1*ones(size(scores)).*val;

end

[~,pos] = max(abs(50-ber));
th = t(pos);

predicted_scores = -1*ones(size(scores)).*val;
predicted_scores(x<th) = val;

stats = compute_metrics(predicted_scores, scores,-1,x);

end

function stats = compute_metrics(labels,scores,value,decision_values)

% positive class: Shockable rhythms
fv = find(scores == 1);
tp = sum(labels(fv)==1);
fn = sum(labels(fv)==value);
pc = tp + fn;

% negative class: Others
rs = find(scores == value);
tn = sum(labels(rs)==value);
fp = sum(labels(rs)==1);
nc = tn + fp;

% metrics
stats.sen = tp/(tp+fn) * 100; %Sensitivity, Also known as Recall
stats.esp = tn/(tn+fp) * 100; %Especificity,  Also known as False Positives Rate
stats.pp  = tp/(tp+fp) * 100; %Pos. Predictivity, Also as Precision
stats.acc = (tp + tn) / (pc + nc) * 100; % Accuracy
stats.err = 100-stats.acc;               % Error rate
stats.fsc = 2*stats.pp*stats.sen / (stats.pp + stats.sen);
stats.ber = 0.5* (fn/pc + fp/nc ) * 100 ;  
stats.gme = sqrt(stats.sen*stats.esp);

% ROC curve
if nargin > 3
    

    labels = [ones(1,length(fv)), zeros(1,length(rs))];
    scores = [decision_values(fv)', decision_values(rs)'];
    
    if mean(decision_values(fv)) < mean(decision_values(rs))
        [pfa,pd,~,auc] = perfcurve(labels,scores,0);
    else
        [pfa,pd,~,auc] = perfcurve(labels,scores,1);
    end
    
    stats.auc = auc*100;
    stats.pfa = pfa;
    stats.pd  = pd;
end

end


