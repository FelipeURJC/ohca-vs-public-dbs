function bootstraping
% Run five classifiers (L1-penalized logistic regression, bagging trees,
% random forests, boosted trees and SVM) to classify the data in ../data/.
%
% Comments:
%
%   - Free parameters are computed out of this code and hardcoded here.
%   - Performance is measured using bootstrap resampling.
%   - Feature selection is carried out before bootstraping, with two possible
%     methods (boosted trees and L1-LR).
%   - In the same execution results for both feature selection methods
%     and with no feature selection (all features) are computed, using
%     paired boostrap.
%   - Results are saved in ../data folder.
%   - Auxiliary functions are at the end of this file.
%   - Only one additional function is needed: feature_selection.m
%
% Options: use the first lines of this script to set the following options:
%
%   - Window length (w_length): 4 or 8 seconds.
%   - Database (ddbb): public or OHCA.
%   - Do feature selection (do_feat_sel): if true it run the feature
%     selection process. If false, it loads the feature selection results
%     from ../data/.
%
% This code, by Felipe Alonso-Atienza and Carlos Figuera
% felipe.alonso@urjc.es, carlos.figuera@urjc.es
%                                               Nov 2015

%% Initializations
clear all; clc; close all;

% Options
w_length      = 8;           % options: {4,8}
ddbb          = 'public_db'; % options: {public_db,ohca_db}
do_feat_sel   = false;       % options: {true,false}

% Add paths
addpath('./glmnet_matlab')
addpath('./mysvm');


%% Load data
data_path = './../data/';
filename  = sprintf('%sdata_%d',data_path, w_length);
load(filename);

% Samples for each ddbb
public_db = 1:samples_for_dbs.ahadb(end);
ohca_db   = samples_for_dbs.ohcadb;

%% Feature transformation
T = feature_transformation(Tabla);
VarNames = T.Properties.VariableNames;
FeatNames = VarNames(1:end-4); %v02

if strcmp(ddbb,'public_db')
    % public dbs
    X = T{public_db,FeatNames};
    y = T.y(public_db);
    id = T.patients_id(public_db);
else
    % ohca dbs
    X = T{ohca_db,FeatNames};
    y = T.y(ohca_db);
    id = T.patients_id(ohca_db);
end
     
fprintf('\n\b%s with WL = %d [s]...DONE\n',ddbb,w_length);

%% Set free parameters for L1-log regression and boosted trees;

% L1-log param for each ddbb
public_db_w4 = 2.6580e-04;
ohca_db_w4 = 0.0011;
public_db_w8 = 2.3821e-04;
ohca_db_w8 = 0.0020;
expression = sprintf('%s_w%d',ddbb,w_length);
params.lambda_opt = eval(expression);

% Boosted trees params
params.method  = 'LogitBoost';
params.lr      = 1;
params.iter    = 50;
params.minleaf = 100;
    
%% Randomize samples
[Nsamples,Nfeat] = size(X);
mix_idx = randperm(Nsamples);
X = X(mix_idx,:);
y = y(mix_idx);
id = id(mix_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature selection with L1-Logistic Regression and BOOSTED TREES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_name = ['../data/selected_features_' num2str(w_length) '_' ddbb];
if do_feat_sel,
    feats_BST = feature_selection(X,y,id,'BST', params, 'ber', FeatNames);
    feats_LLR = feature_selection(X,binary_labels(y),id,'LLR', params, 'ber', FeatNames);
    save (mat_name, 'feats_BST', 'feats_LLR');
else
    load(mat_name)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrap resampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = 500;
disp('Starting bootstrap for algorithms...');
tic
for b=1:B
    
    fprintf('b = %d/%d.', b, B);
    
    %% Extract the resample
    data = unique(id);
    Nsamples = length(data);
    id_pat_train = datasample(data,Nsamples);
    
    idx_train = [];
    for j=1:Nsamples
        idx_train = [idx_train; find(id==id_pat_train(j))]; %#ok<AGROW>
    end
    
    idx_val = setdiff(1:length(id),idx_train);
    
    Xt = X(idx_train,:);
    yt = y(idx_train);
    
    Xv = X(idx_val,:);
    yv = y(idx_val);
    
    [Xts,mu,sigma] = scale_data(Xt);
    Xvs = scale_data(Xv,mu,sigma);
    
    
    %% Machine learning algorithms
    for s = 1:3, % s = 1 --> all features
                 % s = 2 --> selected features with BST
                 % s = 3 --> selected features with L1-LR
                 % If you want to run the algorithms only with "all
                 % features" change the for statement to: 
                 %           for s = 1:1,
                 
        switch s,
            case 1
                feats = 1:Nfeat;
            case 2
                feats = feats_BST.feat_index;
            case 3
                feats = feats_LLR.feat_index;
        end
        
        % L1-logistic regression
        alg(1).output(b,s) = lasso_detector(Xts(:,feats),binary_labels(yt),Xvs(:,feats),binary_labels(yv),params.lambda_opt);
        alg(1).name = 'L_1-LR';
        
        % random forest
        alg(2).output(b,s) = random_forest(Xts(:,feats),yt,Xvs(:,feats),yv);
        alg(2).name = 'RF';
        
        % bagging
        alg(3).output(b,s) = bagged_trees(Xts(:,feats),yt,Xvs(:,feats),yv);
        alg(3).name = 'BAG';
        
        % boosted trees
        alg(4).output(b,s) = boosted_trees(Xts(:,feats),yt,Xvs(:,feats),yv, params);
        alg(4).name = 'BST';
        
        % svms
        alg(5).output(b,s) = svm_detector(Xts(:,feats),yt,Xvs(:,feats),yv);
        alg(5).name = 'SVM';
    end

    fprintf(' Tiempo restante: %2.2f minutos. \n', toc*(B-b)/(60*b));

end

nombre = ['../data/Paired_bootstrap_' num2str(w_length) '_' ddbb];
save(nombre)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = feature_transformation(T)

T.kurt   = log10(T.kurt + 3);
T.M      = log10(T.M + 1);
T.A3     = sqrt(T.A3);
T.count3 = log10(T.count3);
T.x1     = log10(T.x1);
T.x3     = log10(T.x3);
T.x5     = sqrt(T.x5);

end

function [Xs,mu,sigma] = scale_data(X,mu,sigma)

if nargin <2
    mu = mean(X);
    sigma = std(X);
end

Xs = bsxfun(@minus,X,mu);
Xs = bsxfun(@times,Xs,1./sigma);

end

function ybool = binary_labels(y)

    ybool = y;
    ybool(y==-1) = 0;
    
end

%% Machine learning algorithms
function output = lasso_detector(Xt,yt,Xv,yv,lambda_opt)

% options
family = 'binomial';
options = [];

% y must be boolean {0,1}
lm_model = glmnet(Xt,yt,family,options);

ind = find(lm_model.lambda < lambda_opt,1);
betas = lm_model.beta(:,ind);

decision_values = glmnetPredict(lm_model,Xv,lambda_opt);
predicted_label = sign(decision_values-eps);
predicted_label(predicted_label==-1)= 0;

stats = compute_metrics(predicted_label,yv,0,decision_values);

output.predicted_label = predicted_label;
output.decision_values = decision_values;
output.stats = stats;
output.beta = betas;

end

function output = random_forest(Xt,yt,Xv,yv)

Ntree = 50;
Nleaf = 10;

% random forest
rf_tree = TreeBagger(Ntree,Xt,yt,...
    'Method', 'classification',...
    'minleaf', Nleaf,...
    'prior', [0.8 0.2]);
        
[predClass,classifScore] = rf_tree.predict(Xv);
decision_values = classifScore(:,2);

% convert to double
S = sprintf('%s*', predClass{:});
predicted_label = sscanf(S, '%f*');

stats = compute_metrics(predicted_label,yv,-1,decision_values);

output.predicted_label = predicted_label;
output.decision_values = decision_values;
output.stats = stats;

end

function output = bagged_trees(Xt,yt,Xv,yv)

Ntree = 50;
Nleaf = 10;
        
% bagging
bg_tree = TreeBagger(Ntree,Xt,yt,...
    'Method', 'classification',...
    'minleaf', Nleaf,...
    'NvarToSample','all',...        
    'prior', [0.8 0.2]);

[predClass,classifScore] = bg_tree.predict(Xv);
decision_values = classifScore(:,2);

% convert to double
S = sprintf('%s*', predClass{:});
predicted_label = sscanf(S, '%f*');

stats = compute_metrics(predicted_label,yv,-1,decision_values);

output.predicted_label = predicted_label;
output.decision_values = decision_values;
output.stats = stats;

end

function output = boosted_trees(Xt,yt,Xv,yv, params)

method  = params.method; 
lr      = params.lr;               %learning rate
iter    = params.iter;
minleaf = params.minleaf;

tree  = templateTree('minleaf',minleaf); 

bt_model = fitensemble(Xt,yt,method,iter,tree,...
        'LearnRate',lr,'type','classification');
%importance = bt_model.predictor_importance;
    
[predicted_label,scores] = predict(bt_model,Xv);
decision_values = scores(:,1); %poor stimation of auc

stats = compute_metrics(predicted_label,yv,-1,decision_values);

output.predicted_label = predicted_label;
output.decision_values = decision_values;
output.stats = stats;
%output.imp = importance;

end

function output = svm_detector(Xt,yt,Xv,yv)

Parameters = '-s 0 -t 2 -w-1 1 -w1 3 -j 1 -c 10 -g 0.2';
        
svm_model = mysvmtrain(yt,Xt,Parameters);
[predicted_label,~,decision_values]=... 
    mysvmpredict(yv,Xv,svm_model);

stats = compute_metrics(predicted_label,yv,-1,decision_values);

output.predicted_label = predicted_label;
output.decision_values = decision_values;
output.stats = stats;

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
    stats.gme];%,...
    %stats.auc];
end