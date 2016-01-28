function L1_LR_detector

% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

clc; close all;

addpath('./glmnet_matlab')

w_length = [4, 8];

% display results
msg = sprintf('\t\tSEN\tESP\tPP\tACC\tERR\tf1-sc\tBER\tGmean\tAUC');
disp(msg)

% color
color = {'b','k','b:','k:'};

% output
file_to_save = 'free_parameters_lasso';
lm.family  = 'binomial';
lm.options = [];


% Get metrics
lambda_opt  = zeros(1,4);
cv_error    = zeros(1,4);
for j = 1:length(w_length)
    
    % load data
    data_path = '../data/';
    filename  = sprintf('%sdata_%d',data_path, w_length(j));
    load(filename);
    
    %Analyze performance of individual features
    public_db = 1:samples_for_dbs.ahadb(end);
    ohca_db   = samples_for_dbs.ohcadb;
    
    T = feature_transformation(Tabla);
    VarNames = T.Properties.VariableNames;
    FeatNames = VarNames(1:end-5);
    
    % public dbs
    Xs = T{public_db,FeatNames};   
    Xs = scale_data(Xs);
    y = T.y(public_db);
    ybool = y;
    ybool(y==-1) = 0;
    
    
    
    %info = sprintf('PublicDB %ds',w_length(j));
    info = sprintf('public_db_w%d',w_length(j));
    [lambda_opt(j*2-1), cv_error(j*2-1)] =...
        lasso_function(Xs,ybool,color{j*2-1},j*2-1,...
        sprintf('PUBLIC-%ds',w_length(j)),...
        FeatNames);
    
    expression = ['lm.lambda_opt.' info '= lambda_opt;'];
    eval(expression)
        
    % ohca dbs
    Xs = T{ohca_db,FeatNames};
    Xs = scale_data(Xs);
    y = T.y(ohca_db);
    ybool = y;
    ybool(y==-1) = 0;
    
    info = sprintf('ohca_db_w%d',w_length(j));
    [lambda_opt(j*2), cv_error(j*2)] =...
        lasso_function(Xs,ybool,color{j*2},j*2,...
        sprintf('OHCA-%ds\t',w_length(j)),...
        FeatNames);
    
    expression = ['lm.lambda_opt.' info '= lambda_opt;'];
    eval(expression)
    
end

figure(1)
title('L1-Logistic Regression','Interpreter','latex','Fontsize',16)
ylabel('Missclassification Error','Interpreter','latex','Fontsize',16)
xlabel('$log(\lambda)$','Interpreter','latex','Fontsize',16);
legend('PUBLIC-4s','OHCA-4s','PUBLIC-8s','OHCA-8s','Location','northwest')
legend boxoff
grid on;

for k=1:4
    
    if mod(k,2) == 0 
    semilogx(lambda_opt(k),cv_error(k),...
        'ko','MarkerFaceColor','r',...
        'LineWidth',2,...
        'Markersize',10);
    else
     semilogx(lambda_opt(k),cv_error(k),...
        'bo','MarkerFaceColor','r',...
        'LineWidth',2,...
        'Markersize',10);
    end
end

%save(file_to_save,'lm');

end

function disp_stats(stats,info)
% all
msg=sprintf('%s\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f',...
    info, stats.sen,stats.esp,stats.pp,stats.acc,stats.err,...
    stats.fsc,stats.ber,stats.gme,stats.auc);
disp(msg)

end

function [Xs,mu,sigma] = scale_data(X,mu,sigma)

if nargin <2
    mu = mean(X);
    sigma = std(X);
end

Xs = bsxfun(@minus,X,mu);
Xs = bsxfun(@times,Xs,1./sigma);

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


function [lambda_opt, cv_error] = lasso_function(X,y,c,k,info,FeatNames)

rng('default') % for reproducibility
cvfit = cvglmnet(X, y,'binomial',[],'class');
lambda_opt = cvfit.lambda_1se;
cv_error = cvfit.cvm(cvfit.lambda==cvfit.lambda_1se);

% normalized model coefs
active_coefs = cvglmnetCoef(cvfit);
betaj = abs(active_coefs(2:end));
betaj = betaj./max(betaj);

[betaj,id] = sort(betaj,'descend');

% performance
decision_values = cvglmnetPredict(cvfit, X, 'lambda_1se');
ypred           = sign(decision_values-eps);
ypred(ypred==-1)= 0;
stats = compute_metrics(ypred,y,0,decision_values);
disp_stats(stats,info)

figure(1);
semilogx(cvfit.lambda,cvfit.cvm,c,'Linewidth',2); hold on;

figure(2)
subplot(2,2,k);
barh(betaj);
set(gca,'YTick',1:length(FeatNames));
set(gca,'YTickLabel',FeatNames(id),'Fontsize',14);
xlabel('$|\beta_j|$','Interpreter','latex')
title(info)
axis tight

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