function search_free_parameters()

clear all; clc; close all;
rng(2); % Use the same seed as in bootstrap function

% Add glmnet package
addpath('./glmnet_matlab')
addpath('./mysvm/')

% Options
w_length = 8;                       % options: {4,8}
ddbb     = 'ohca';                % options: {public,ohca}

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

    
Xt = X(idx_train,:);
yt = y(idx_train);
pt = patient_id(idx_train);
     
[Xts,~,~] = scale_data(Xt);

%% Free parameters search for each algorithm:

LR_params(Xts,yt,pt);
RF_params(Xts,yt,pt);
BG_params(Xts,yt,pt);
BT_params(Xts,yt,pt);
SVM_params(Xts,yt,pt);

end % OF FUNCION

%% SVM

function SVM_params(X,y,pat_id)

% add mysvm library (from )
addpath('./mysvm/')

% This grid
C = logspace(-1.5,3,8);
G = logspace(-9,-3,8);

ber = zeros(length(G),length(C));

pat_unique_ids = unique(pat_id);
cvmodel = cvpartition(length(pat_unique_ids),'k',10);

for g = 1:length(G)
    for c = 1:length(C)
        
        Parameters = ['-s 0 -t 2 -w-1 1 -w1 3 -j 1 -c ', num2str(C(c)),...
            ' -g ', num2str(G(g))];
        
        stats = zeros(cvmodel.NumTestSets,9);
        
        for n = 1:cvmodel.NumTestSets
            
            % Read the ids of training-set patients
            trn_pat_id = pat_unique_ids(cvmodel.training(n));

            % Read the registers associated to that patients
            idx_train = [];
            for j=1:length(trn_pat_id)
                idx_train = [idx_train; find(pat_id==trn_pat_id(j))]; %#ok<AGROW>
            end
    
            % The rest of registers are the test set.
            idx_val = setdiff(1:length(y),idx_train);
                        
            svm_model = svm_detector_train(X(idx_train,:),y(idx_train), Parameters);
            [y_hat, decision_values] = svm_detector_test(svm_model,X(idx_val,:),y(idx_val));
                        
            stats(n,:) = compute_metrics(y_hat,y(idx_val),-1,decision_values);
            
        end
        
        ber(g,c) = mean(stats(:,7)); % mean ber value
            
    end
  
end

figure;
imagesc(log10(C),log10(G),ber);
colorbar
xlabel('log(C)','Interpreter','Latex','Fontsize',14);
ylabel('$\log(\gamma)$','Interpreter','Latex','Fontsize',14);
grid on;


end

%% LOGISTIC REGRESSION

function LR_params(X,y,pat_id)

% This grid
LAMBDA = logspace(-6,-1,30);              

ber = zeros(length(LAMBDA),1);

pat_unique_ids = unique(pat_id);
cvmodel = cvpartition(length(pat_unique_ids),'k',10);

for l = 1:length(LAMBDA)
        
    lambda = LAMBDA(l);
    
    stats = zeros(cvmodel.NumTestSets,9);
    
    for n = 1:cvmodel.NumTestSets
        
        % Read the ids of training-set patients
        trn_pat_id = pat_unique_ids(cvmodel.training(n));
        
        % Read the registers associated to that patients
        idx_train = [];
        for j=1:length(trn_pat_id)
           idx_train = [idx_train; find(pat_id==trn_pat_id(j))]; %#ok<AGROW>
        end
        
        % The rest of registers are the test set.
        idx_val = setdiff(1:length(y),idx_train);
                
        lm_model = lasso_detector_train(X(idx_train,:),y(idx_train));
        [y_hat, decision_values] = lasso_detector_test(lm_model,X(idx_val,:),lambda);
        
        stats(n,:) = compute_metrics(y_hat,y(idx_val),-1,decision_values);
        
    end
    
    ber(l) = mean(stats(:,7)); % mean ber value
    
end


semilogx(LAMBDA,ber,'.-','MarkerSize',25);
xlabel('\lambda'); ylabel('BER with 10-fold CV');

end

%% RANDOM FOREST

function RF_params(X,y,pat_id)

% This grid
NTREE      = [25 50 75];              
NLEAF    = [10 20 30];

ber = zeros(length(NTREE),length(NLEAF));

pat_unique_ids = unique(pat_id);
cvmodel = cvpartition(length(pat_unique_ids),'k',10);

for Ntree = 1:length(NTREE)
        
    for Nleaf = 1:length(NLEAF)
            
            params.Ntree = NTREE(Ntree);
            params.Nleaf = NLEAF(Nleaf);
            
            stats = zeros(cvmodel.NumTestSets,9);
            
            for n = 1:cvmodel.NumTestSets
                
                % Read the ids of training-set patients
                trn_pat_id = pat_unique_ids(cvmodel.training(n));
                
                % Read the registers associated to that patients
                idx_train = [];
                for j=1:length(trn_pat_id)
                    idx_train = [idx_train; find(pat_id==trn_pat_id(j))]; %#ok<AGROW>
                end
                
                % The rest of registers are the test set.
                idx_val = setdiff(1:length(y),idx_train);
                
                rf_tree = random_forest_train(X(idx_train,:),y(idx_train),params);
                [y_hat, decision_values] = random_forest_test(rf_tree,X(idx_val,:));
                
                stats(n,:) = compute_metrics(y_hat,y(idx_val),-1,decision_values);
                
            end
            
            ber(Ntree,Nleaf) = mean(stats(:,7)); % mean ber value
        
    end
end

end

%% BAGGING TREES

function BG_params(X,y,pat_id)

% This grid
NTREE      = [25 50 75];              
NLEAF    = [10 20 30];

ber = zeros(length(NTREE),length(NLEAF));

pat_unique_ids = unique(pat_id);
cvmodel = cvpartition(length(pat_unique_ids),'k',10);

for Ntree = 1:length(NTREE)
    
    for Nleaf = 1:length(NLEAF)
            
            params.Ntree = NTREE(Ntree);
            params.Nleaf = NLEAF(Nleaf);
            
            stats = zeros(cvmodel.NumTestSets,9);
            
            for n = 1:cvmodel.NumTestSets
                
                % Read the ids of training-set patients
                trn_pat_id = pat_unique_ids(cvmodel.training(n));
                
                % Read the registers associated to that patients
                idx_train = [];
                for j=1:length(trn_pat_id)
                    idx_train = [idx_train; find(pat_id==trn_pat_id(j))]; %#ok<AGROW>
                end
                
                % The rest of registers are the test set.
                idx_val = setdiff(1:length(y),idx_train);
                
                rf_tree = bagged_trees_train(X(idx_train,:),y(idx_train),params);
                [y_hat, decision_values] = bagged_trees_test(rf_tree,X(idx_val,:));
                
                stats(n,:) = compute_metrics(y_hat,y(idx_val),-1,decision_values);
                
            end
            
            ber(Ntree,Nleaf) = mean(stats(:,7)); % mean ber value
                    
    end
end

end


%% BOOSTING TREES

function BT_params(X,y,pat_id)

% This grid
LR      = [0.1 0.5 1];               %learning rate
ITER    = [10 50 100];
MINLEAF = [50 100 200];

ber = zeros(length(LR),length(ITER),length(MINLEAF));

pat_unique_ids = unique(pat_id);
cvmodel = cvpartition(length(pat_unique_ids),'k',10);

for lr = 1:length(LR)
    for iter = 1:length(ITER)
        for minleaf = 1:length(MINLEAF),
            
            params.lr = LR(lr);
            params.iter = ITER(iter);
            params.minleaf = MINLEAF(minleaf);
            
            stats = zeros(cvmodel.NumTestSets,9);
            
            for n = 1:cvmodel.NumTestSets
                
                % Read the ids of training-set patients
                trn_pat_id = pat_unique_ids(cvmodel.training(n));
                
                % Read the registers associated to that patients
                idx_train = [];
                for j=1:length(trn_pat_id)
                    idx_train = [idx_train; find(pat_id==trn_pat_id(j))]; %#ok<AGROW>
                end
                
                % The rest of registers are the test set.
                idx_val = setdiff(1:length(y),idx_train);
                
                bt_model = boosted_trees_train(X(idx_train,:),y(idx_train),params);
                [y_hat, decision_values] = boosted_trees_test(bt_model,X(idx_val,:));
                
                stats(n,:) = compute_metrics(y_hat,y(idx_val),-1,decision_values);
                
            end
            
            ber(lr,iter,minleaf) = mean(stats(:,7)); % mean ber value
                        
        end
    end
end


end

%% AUXILIARY FUNCTIONS

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

function [Xs,mu,sigma] = scale_data(X,mu,sigma)

if nargin <2
    mu = mean(X);
    sigma = std(X);
end

Xs = bsxfun(@minus,X,mu);
Xs = bsxfun(@times,Xs,1./sigma);

end

function svm_model = svm_detector_train(Xt,yt, Parameters)
svm_model = mysvmtrain(yt,Xt,Parameters);
end

function [predicted_label, decision_values] = svm_detector_test(svm_model,Xv,yv)
[predicted_label,~,decision_values]= mysvmpredict(yv,Xv,svm_model);
end


function lm_model = lasso_detector_train(Xt,yt)
family = 'binomial'; options = [];
lm_model = glmnet(Xt,yt,family,options);
end

function [predicted_label, decision_values] = lasso_detector_test(lm_model,Xv,lambda_opt)
decision_values = glmnetPredict(lm_model,Xv,lambda_opt);
predicted_label = sign(decision_values-eps);
end

function rf_tree = random_forest_train(Xt,yt,params)
rf_tree = TreeBagger(params.Ntree,Xt,yt,'Method', 'classification','minleaf', params.Nleaf,'prior', [0.8 0.2]);
end

function [predicted_label,decision_values] = random_forest_test(rf_tree,Xv)     
[predClass,classifScore] = rf_tree.predict(Xv);
decision_values = classifScore(:,2);
% convert to double
S = sprintf('%s*', predClass{:});
predicted_label = sscanf(S, '%f*');
end

function bg_tree = bagged_trees_train(Xt,yt,params)
bg_tree = TreeBagger(params.Ntree,Xt,yt,'Method', 'classification','minleaf', params.Nleaf,'NvarToSample','all','prior', [0.8 0.2]);
end

function [predicted_label, decision_values] = bagged_trees_test(bg_tree,Xv)
[predClass,classifScore] = bg_tree.predict(Xv);
decision_values = classifScore(:,2);
% convert to double
S = sprintf('%s*', predClass{:});
predicted_label = sscanf(S, '%f*');
end

function bt_model = boosted_trees_train(Xt,yt,params)
method  = 'LogitBoost'; 
lr      = params.lr;               %learning rate
iter    = params.iter;
minleaf = params.minleaf;
tree  = templateTree('minleaf',minleaf); 
bt_model = fitensemble(Xt,yt,method,iter,tree,'LearnRate',lr,'type','classification');
end

function [predicted_label, decision_values] = boosted_trees_test(bt_model,Xv)
[predicted_label,scores] = predict(bt_model,Xv);
decision_values = scores(:,1);
end