function SVM_detector

% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

clc;

% display results
info.verbose = 0;

% add mysvm library (from )
addpath('./mysvm/')

% load data
w_length = [4, 8];
info.k = 1;         % counter

ber = cell(1,4);
for j = 1:length(w_length)
    
    % load data
    data_path = '../data/';
    filename  = sprintf('%sdata_%d',data_path, w_length(j));
    load(filename);
    
    %Analyze performance of individual features
    public_db = 1:samples_for_dbs.ahadb(end);
    ohca_db   = samples_for_dbs.ohcadb;
    
    T = feature_transformation(Tabla);
    T = scale_table(T);
    VarNames = T.Properties.VariableNames;
    FeatNames = VarNames(1:end-5);
    
    % public dbs
    info.ddbb = sprintf('PUBLIC-%ds',w_length(j));
    
    if info.verbose
        disp(info.ddbb)
        msg = sprintf('SEN\tESP\tPP\tACC\tERR\tf1-sc\tBER\tGmean\tAUC\tC\tGamma');
        disp(msg)
    end
    
    X = T{public_db,FeatNames};
    y = T.y(public_db);
    [C,G,ber{info.k}]= SVMBinaryClass(X,y,info);
    info.k = info.k + 1;
        
    % ohca dbs
    info.ddbb = sprintf('OHCA-%ds',w_length(j));
    
    if info.verbose
        disp(info.ddbb)
        msg = sprintf('SEN\tESP\tPP\tACC\tERR\tf1-sc\tBER\tGmean\tAUC\tC\tGamma');
        disp(msg)
    end
    
    X = T{ohca_db,FeatNames};
    y = T.y(ohca_db);
    [C,G,ber{info.k}] = SVMBinaryClass(X,y,info);
    info.k = info.k + 1;
    
    % conclusion: C = 10, G = 0.2;
    
                
end

figure(2)
color = {'b-o','r-o','b--d','r-d'};
for g=1:numel(G)
    subplot(2,2,g)
    for k=1:4
        plot(log10(C), ber{k}(g,:),color{k}); hold on;
    end
    xlabel('log(C)','Interpreter','Latex','Fontsize',14);
    ylabel('BER','Interpreter','Latex','Fontsize',14);
    title_msg = sprintf('G = %1.2f',G(g));
    title(title_msg,'Interpreter','Latex','Fontsize',14)
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

function [C,G,out_of_sample_ber] = SVMBinaryClass(X,y,info)


% finer grid
%C = 10.^(-2:0.5:1.5);
%G = 10.^(-2:0.2:0.5);

% This grid
C = 10.^(-2:0.5:1.5);
G = [0.01, 0.1, 0.2, 0.3];

%C = 10;
%G = 0.2;

out_of_sample_ber = zeros(length(G),length(C));

rng('default') % for reproducibility 
cvmodel = cvpartition(y,'k',5);

for g = 1:length(G)
    for c = 1:length(C)
            
        Parameters = ['-s 0 -t 2 -w-1 1 -w1 3 -j 1 -c ', num2str(C(c)),...
            ' -g ', num2str(G(g))];
        
        stats = zeros(cvmodel.NumTestSets,9);
        
        for n = 1:cvmodel.NumTestSets
            trIdx = cvmodel.training(n);
            teIdx = cvmodel.test(n);
                        
            svm_model = mysvmtrain(y(trIdx),X(trIdx,:),Parameters);
            [y_hat,~,decision_values]= mysvmpredict(y(teIdx),X(teIdx,:),svm_model);
                        
            stats(n,:) = compute_metrics(y_hat,y(teIdx),-1,decision_values);
            
        end
        
        out_of_sample_ber(g,c) = mean(stats(:,7)); % mean ber value
        
        if info.verbose
            disp_stats(mean(stats),C(c),G(g));
        end
    
    end
  
end

figure(1);
subplot(2,2,info.k);
imagesc(log10(C),log10(G),out_of_sample_ber);
colorbar
xlabel('log(C)','Interpreter','Latex','Fontsize',14);
ylabel('$\log(\gamma)$','Interpreter','Latex','Fontsize',14);
text_title = sprintf('BER: %s',info.ddbb);
title(text_title,'Interpreter','Latex')
grid on;

end

function T = scale_table(T)

A = table2array(T(:,1:end-5));
A = bsxfun(@minus,A,mean(A));
A = bsxfun(@times,A,1./std(A));

T(:,1:end-5) = array2table(A);

end


function disp_stats(stats,c,gamma)
% all
msg=sprintf('%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f',...
    stats(1),stats(2),stats(3),stats(4),stats(5),...
    stats(6),stats(7),stats(8),stats(9),...
    c,gamma);
disp(msg)

end


function stats_vector = compute_metrics(labels,scores,value,decision_values)

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

stats_vector = [stats.sen,...
    stats.esp,...
    stats.pp,...
    stats.acc,...
    stats.err,...
    stats.fsc,...
    stats.ber,...
    stats.gme,...
    stats.auc];
end
