function bootstraping
% Run five classifiers (L1-penalized logistic regression, bagging trees,
% random forests, boosted trees and SVM) to classify the data in ../data/.
%
% Comments:
%
%   - Free parameters are computed out of this code and loaded here.
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
%                                                 Mar 2016

%% Initializations
clear all; clc; close all;
rng(2); 

% Add glmnet package
addpath('./glmnet_matlab')
addpath('./mysvm/')

% Free params path
data_path = './../data/';
free_parametets_path = [data_path 'FreeParameters/'];

% Options
w_length = 4;                       % options: {4,8}
ddbb     = 'public';                % options: {public,ohca}
do_feat_sel = false;                % options: {true,false}


%% Load data
filename = ['./../data/data_' ddbb '_' num2str(w_length)];
load(filename);

X = data.X;
y = data.y;
patient_id = data.patient_id;

%% Separate patients for train and validation sets

train_rate = 0.8;

if strcmp(ddbb,'public')
    
    % Take train patients from CUDB, VFDB and AHADB
    
    idx_cudb = strcmp(data.ddbb,'cudb');
    idx_vfdb = strcmp(data.ddbb,'vfdb');
    idx_ahdb = strcmp(data.ddbb,'ahadb');
    
    cu_pat_ids = unique(data.patient_id(idx_cudb));
    vf_pat_ids = unique(data.patient_id(idx_vfdb));
    ah_pat_ids = unique(data.patient_id(idx_ahdb));
    
    n_pats_cu = length(cu_pat_ids);
    n_pats_vf = length(vf_pat_ids);
    n_pats_ah = length(ah_pat_ids);
    
    pats_id_cudb_train = cu_pat_ids(randperm(n_pats_cu, round(n_pats_cu * train_rate)));
    pats_id_vfdb_train = vf_pat_ids(randperm(n_pats_vf, round(n_pats_vf * train_rate)));
    pats_id_ahdb_train = ah_pat_ids(randperm(n_pats_ah, round(n_pats_ah * train_rate)));
    
    pats_id_train = [pats_id_cudb_train;
                     pats_id_vfdb_train;
                     pats_id_ahdb_train];
    
else % If OHCA
    all_pat_ids = unique(patient_id);           % All the ids for all the patients
    n_pats = length(all_pat_ids);               % Number of patients
    n_pats_train = round(n_pats * train_rate);  % Number of patients for feature selection, training algorithms, etc.
    pats_id_train = all_pat_ids(randperm(n_pats, n_pats_train));
end

idx_train = [];
for p = 1:length(pats_id_train),
    idx_train = [idx_train; find(patient_id==pats_id_train(p))]; %#ok<AGROW>
end

idx_val = setdiff(1:length(patient_id),idx_train);
    
Xt = X(idx_train,:);
yt = y(idx_train);
pt = patient_id(idx_train);
dbt = data.ddbb(idx_train);

Xv = X(idx_val,:);
yv = y(idx_val);
     
[Xts,mu,sigma] = scale_data(Xt);
Xvs = scale_data(Xv,mu,sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Free parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

free_parameters_file = [free_parametets_path 'params_' num2str(w_length) '_' ddbb '.mat'];
load(free_parameters_file); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature selection with L1-Logistic Regression and BOOSTED TREES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_name = ['../RESULTS/selected_features_' num2str(w_length) '_' ddbb];
if do_feat_sel,
    
    feats_BST = feature_selection(Xts,yt,pt,dbt,'BST', 'ber', data.FeatNames);
    feats_LLR = feature_selection(Xts,binary_labels(yt),pt,dbt,'LLR', 'ber', data.FeatNames);
    save (mat_name, 'feats_BST', 'feats_LLR');
           
else
    load(mat_name) %#ok<UNRCH>
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train
%    with all feats --> s = 1
%    with feats. selected with BST --> s = 2
%    with feats. selected with LLR --> s = 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:3, % s = 1 --> all features
             % s = 2 --> selected features with BST
             % s = 3 --> selected features with LLR
             
     switch s,
         case 1
             feats = 1:size(X,2);
         case 2
             feats = feats_BST.feat_index;
         case 3
             feats = feats_LLR.feat_index;
     end
     
     lm_model{s}  = lasso_detector_train(Xts(:,feats),binary_labels(yt));
     rf_tree{s}   = random_forest_train(Xts(:,feats),yt,params.RF);
     bg_tree{s}   = bagged_trees_train (Xts(:,feats),yt,params.BG);
     bt_tree{s}   = boosted_trees_train(Xts(:,feats),yt,params.BT);
     svm_model{s} = svm_detector_train (Xts(:,feats),yt,params.SVM);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrap resampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = 500;

disp('Starting bootstrap for algorithms...');

for b = 1:B,

    fprintf('\n b = %d', b);
    
    % Resample the test set
    idx_test = datasample(1:length(idx_val), length(idx_val));
    
    Xtest = Xvs(idx_test,:);
    ytest = yv(idx_test);


    cosa(b) = length(ytest);
    
    for s = 1:3, % s = 1 --> all features
                 % s = 2 --> selected features with BST
                 % s = 3 --> selected features with LLR
                 
        switch s,
            case 1
                feats = 1:size(X,2);
            case 2
                feats = feats_BST.feat_index;
            case 3
                feats = feats_LLR.feat_index;
        end      
        
        %L1-logistic regression
        alg(1).output(b,s) = lasso_detector_test(lm_model{s},Xtest(:,feats),binary_labels(ytest),params.LR.lambda);
        alg(1).name = 'L_1-LR';
        
        %random forest
        alg(2).output(b,s) = random_forest_test(rf_tree{s},Xtest(:,feats),ytest);
        alg(2).name = 'RF';
        
        %bagging
        alg(3).output(b,s) = bagged_trees_test(bg_tree{s},Xtest(:,feats),ytest);
        alg(3).name = 'BAG';
        
        %boosted trees
        alg(4).output(b,s) = boosted_trees_test(bt_tree{s},Xtest(:,feats),ytest);
        alg(4).name = 'BST';
       
        %svms
        alg(5).output(b,s) = svm_detector_test(svm_model{s},Xtest(:,feats),ytest);
        alg(5).name = 'SVM';  

    end

end
fprintf('\n');
nombre = ['../RESULTS/Paired_bootstrap_v6_' num2str(w_length) '_' ddbb];
save(nombre)

end % of function




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

%
% LASSO
%

function lm_model = lasso_detector_train(Xt,yt)

% options
family = 'binomial';
options = [];

% y must be boolean {0,1}
lm_model = glmnet(Xt,yt,family,options);

end

function output = lasso_detector_test(lm_model,Xv,yv,lambda_opt)

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

%
% RANDOM FOREST
%

function rf_tree = random_forest_train(Xt,yt,p)

Ntree = p.Ntree;%50;
Nleaf = p.Nleaf;%10;

% random forest
rf_tree = TreeBagger(Ntree,Xt,yt,...
    'Method', 'classification',...
    'minleaf', Nleaf,...
    'prior', [0.8 0.2]);

end

function output = random_forest_test(rf_tree,Xv,yv)
     
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


%
% BAGGED TREES
%

function bg_tree = bagged_trees_train(Xt,yt,p)

Ntree = p.Ntree;%50;
Nleaf = p.Nleaf;%10;
        
bg_tree = TreeBagger(Ntree,Xt,yt,...
    'Method', 'classification',...
    'minleaf', Nleaf,...
    'NvarToSample','all',...        
    'prior', [0.8 0.2]);

end


function output = bagged_trees_test(bg_tree,Xv,yv)

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


%
% BOOSTED TREES
%

function bt_model = boosted_trees_train(Xt,yt,p)

method  = 'LogitBoost'; 
lr      = p.lr;%1;               %learning rate
iter    = p.iter;%50;
minleaf = p.minleaf;%100;

tree  = templateTree('minleaf',minleaf); 

bt_model = fitensemble(Xt,yt,method,iter,tree,...
        'LearnRate',lr,'type','classification');
end

function output = boosted_trees_test(bt_model,Xv,yv)

% method  = 'LogitBoost'; 
% lr      = 1;               %learning rate
% iter    = 50;
% minleaf = 100;
    
[predicted_label,scores] = predict(bt_model,Xv);
decision_values = scores(:,1); %poor stimation of auc

stats = compute_metrics(predicted_label,yv,-1,decision_values);

output.predicted_label = predicted_label;
output.decision_values = decision_values;
output.stats = stats;

end


%
% SVM
%

function svm_model = svm_detector_train(Xt,yt,p)

Parameters = ['-s 0 -t 2 -w-1 1 -w1 3 -j 1 -c ' num2str(p.C) ' -g ' num2str(p.G)];
%Parameters = '-s 0 -t 2 -w-1 1 -w1 3 -j 1 -c 10 -g 0.2';
svm_model = mysvmtrain(yt,Xt,Parameters);

end


function output = svm_detector_test(svm_model,Xv,yv)

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
