function visualize_analyze_data

% This function diplays a number of performance metrics for each 
% of the computed ECG parameters, and visualize their histograms 
% for both shockable and non-shockable rhythms.
%
% For each parameters, the optimum threshold is calculated to 
% minimize the Balanced Error Rate (map criterion) using the ad-hoc function
% bayesdetector.m
%
% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es
%
% THIS FUNCTION REQUIRES MATLAB R2015a (predefined histogram function)

close all; clear all; clc;

w_length = [4, 8];
k = 1;

for j = 1:length(w_length)

    % load data
    data_path = '../data/';
    filename  = sprintf('%sdata_%d',data_path, w_length(j));
    load(filename);
        
    %Analyze performance of individual features
    public_db = 1:samples_for_dbs.ahadb(end);
    ohca_db   = samples_for_dbs.ohcadb;
    
    Tabla = feature_transformation(Tabla);
    
    % public dbs
    msg = sprintf('\n\tPUBLIC-%ds\n\t-------------------',...
        w_length(j));
    disp(msg)
    plot_histograms(Tabla,public_db,'PUBLICs',k,w_length(j));

    % ohca db   
    msg = sprintf('\n\tOHCA-%ds\n\t-------------------',...
        w_length(j));
    disp(msg)
    plot_histograms(Tabla,ohca_db,'OHCA',k+1,w_length(j));
        
    k = k + 2;
    
end


end

function plot_histograms(table,samples_db,option,k,wlength)

VarNames = table.Properties.VariableNames;
FeatNames = VarNames(1:end-5);

data = table{:,1:length(FeatNames)};

msg = sprintf('\tSEN\tESP\tPP\tACC\tERR\tf1-sc\tBER\tGmean\tAUC\tTHRES');
disp(msg)

for i=1:length(FeatNames);
    
    xi = table{samples_db,i};
    y  = table.y(samples_db);
    
    % Build a threshold detector
    [th,~,stats] = bayesdetector(xi,y,FeatNames{i},option,wlength);
    
    
    msg=sprintf('%s\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f',...
        FeatNames{i},stats.sen,stats.esp,stats.pp,stats.acc,stats.err,...
        stats.fsc,stats.ber,stats.gme,stats.auc,th);
    disp(msg)
    
    AUC(i) = stats.auc;
    
    figure(i)
    subplot(2,2,k)
        
    h1 = histogram(xi(y==-1), 25, 'Normalization', 'probability');
    hold on;
    h2 = histogram(xi(y==+1), 25, 'Normalization', 'probability');
    
    ymax = max(max(h1.Values),max(h2.Values));
    ymax = ( ceil( ymax*100 ) ) / 100;
    legend({'NSh','Sh'},'Location','northoutside','orientation','horizontal')
    legend boxoff
    
    set(gca, 'XLim', [min(xi) max(xi)],...
        'YTick',0:0.1:ymax, 'YLim', [0 ymax],'Fontsize', 12);
    set (gca,'TicklabelInterpreter','tex')
    plot([th th],[0 ymax],'k:','linewidth',2)
    
    msg = sprintf('%s',FeatNames{i});
    xlabel(msg,'FontSize',14);
    msg = sprintf('%s, %d sec',option,wlength);
    ylabel(msg,'Interpreter','Latex','FontSize',15);
    
    
end

% plot ROC values
[val,id] = sort(AUC,'descend');

figure(i+1);
subplot(2,2,k)

barh(val);
set(gca,'YTick',1:length(FeatNames));
set(gca,'YTickLabel',FeatNames(id),'Fontsize',14);
xlabel('AUC','Interpreter','Latex','FontSize',12);
msg = sprintf('%s, %d sec',option,wlength);
title(msg,'Interpreter','Latex','FontSize',12);
axis tight;


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