function bootstraping_per_ddbb
% This code, by Carlos Figuera
% carlos.figuera@urjc.es
%                     Mar 2016

%% Initializations
clear all; clc; close all; rng(2); 

% Add glmnet package
addpath('./glmnet_matlab')
addpath('./mysvm/')

% Options
w_length = 4;                     % options: {4,8}
train_rate = 0.8;

%% Load data and free parameters
load(['./../DATA/data_public_' num2str(w_length)]);
load(['./../DATA/FreeParameters/params_' num2str(w_length) '_public.mat']);

X = data.X;
y = data.y;
patient_id = data.patient_id;

%% Separate patients for train and validation sets

db_names = {'cudb', 'vfdb', 'ahadb'};
idx_val  = cell(3,1);
Xts      = cell(3,1);
yt       = cell(3,1);
Xvs      = cell(3,1);
yv       = cell(3,1);

for d = 1:length(db_names),
    
    idx = strcmp(data.ddbb,db_names{d});
    pat_ids = unique(data.patient_id(idx));
    n_pats = length(pat_ids);
    pats_id_train = pat_ids(randperm(n_pats, round(n_pats * train_rate)));
    
    idx_train = [];
    for p = 1:length(pats_id_train),
        idx_train = [idx_train; find(patient_id==pats_id_train(p))]; %#ok<AGROW>
    end
    
    idx_val{d} = setdiff(find(idx),idx_train);
    
    Xt = X(idx_train,:);
    yt{d} = y(idx_train);
    
    Xv = X(idx_val{d},:);
    yv{d} = y(idx_val{d});
    
    [Xts{d},mu,sigma] = scale_data(Xt);
    Xvs{d} = scale_data(Xv,mu,sigma);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1:3,
    lm_model{d}  = lasso_detector_train(Xts{d},binary_labels(yt{d}));
    rf_tree{d}   = random_forest_train(Xts{d},yt{d},params.RF);
    bg_tree{d}   = bagged_trees_train (Xts{d},yt{d},params.BG);
    bt_tree{d}   = boosted_trees_train(Xts{d},yt{d},params.BT);
    svm_model{d} = svm_detector_train (Xts{d},yt{d},params.SVM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrap resampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = 500;

disp('Starting bootstrap for algorithms...');

for d = 1:3   
    
    for b = 1:B,
        
        fprintf('\n d = %d, b = %d', d, b);
        
        % Resample the test set
        idx_test = datasample(1:length(idx_val{d}), length(idx_val{d}));
        Xtest = Xvs{d}(idx_test,:);
        ytest = yv{d}(idx_test);
        
        %L1-logistic regression
        alg(1).output(b,d) = lasso_detector_test(lm_model{d},Xtest,binary_labels(ytest),params.LR.lambda);
        alg(1).name = 'L_1-LR';

        %random forest
        alg(2).output(b,d) = random_forest_test(rf_tree{d},Xtest,ytest);
        alg(2).name = 'RF';

        %bagging
        alg(3).output(b,d) = bagged_trees_test(bg_tree{d},Xtest,ytest);
        alg(3).name = 'BAG';

        %boosted trees
        alg(4).output(b,d) = boosted_trees_test(bt_tree{d},Xtest,ytest);
        alg(4).name = 'BST';

        %svms
        alg(5).output(b,d) = svm_detector_test(svm_model{d},Xtest,ytest);
        alg(5).name = 'SVM';
        
    end
    
end

fprintf('\n');
nombre = ['../RESULTS/Bootstrap_per_public_bbdd_' num2str(w_length)]
save(nombre)

print_paired_bootstrap_results(nombre)

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

