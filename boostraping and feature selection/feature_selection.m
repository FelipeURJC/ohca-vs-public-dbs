function [output] = feature_selection(X,y,id,db, method, metric, FeatNames)
%
% feature_selection: performs feature selection for a binary classification
% problem using boosting trees or L1-penalized logistic regression.
%
% Input:
%   X: features matrix
%   y: labels (column vector)
%   id: patient label for each sample
%   db: database name for each sample
%   method: selection method. Options: {'BST','LLR'}:
%       BST: boosting trees.
%       LLR: L1-logistic regression.
%   metric: error metric to perform selection. Options: {'ber', 'error'}
%   FeatNames: feature names (only for drawing)
%
% Outputs:
%
%   output: a struct with the following fields:
%      · selected_features_index: index of the selected features
%      · score: measured importance for the whole set of features
%      · mean_error: mean error for different sizes features subsets (i.e. 1:all, 2: all but
%                  worst, 3: all but two worsts...)
%      · std_error: standardad deviation of error for different sizes features subsets
%      · metric: (input) error metric
%      · method: (input) method for feature selection
%      · params: (input) params for method 
%
%                                               Carlos Figuera, 2016.
%

%% (0) Set the free parameters
tree_method  = 'LogitBoost';
lr      = 1; %params.lr; %learning rate
iter    = 100; %params.iter;
minleaf = 50; %params.minleaf;
tree  = templateTree('minleaf',minleaf); % a measure of the max depth of the tree

if strcmp(db{1},'ohcadb'),
    lambda_opt = 5e-3;
else
    lambda_opt = 1e-2;
    
end
params.lr = lr;
params.iter = iter;
params.minleaf = minleaf;
params.lambda_opt = lambda_opt;

[~, Nfeat] = size(X);


%% (I) Bootstrap resampling for computing the relevance of each feature

B = 500;

imp = zeros(B, Nfeat);
sortedIndex = zeros(B, Nfeat);
bootstrap_metrics = zeros(B, Nfeat, 8);

disp('Bootstraping feature importance...');

for b=1:B
    
    fprintf('b = %d/%d\n', b, B);

    %%% 1. Build resample and scale data
    if strcmp(db{1},'ohcadb'), % if OHCA
        id_pat_train = datasample(unique(id),length(unique(id)));
    else % if Public: take samples from CUDB, VFDB and AHADB
        cu_pat_ids = unique(id(strcmp(db,'cudb')));
        vf_pat_ids = unique(id(strcmp(db,'vfdb')));
        ah_pat_ids = unique(id(strcmp(db,'ahadb')));
        
        pats_id_cudb_train = datasample(cu_pat_ids,length(cu_pat_ids));
        pats_id_vfdb_train = datasample(vf_pat_ids,length(vf_pat_ids));
        pats_id_ahdb_train = datasample(ah_pat_ids,length(ah_pat_ids));
        
        id_pat_train = [pats_id_cudb_train;
            pats_id_vfdb_train;
            pats_id_ahdb_train];
    end
    
    idx_train = [];
    for j=1:length(id_pat_train),
        idx_train = [idx_train; find(id==id_pat_train(j))]; %#ok<*AGROW>
    end
    
    idx_val = setdiff(1:length(id),idx_train);    
    
    Xt = X(idx_train,:);
    yt = y(idx_train);
    
    Xv = X(idx_val,:);
    yv = y(idx_val);
    
    [Xts,mu,sigma] = scale_data(Xt);
    Xvs = scale_data(Xv,mu,sigma);  
    
    %%% 2. Compute feature importance in training
    switch method,
        case 'BST',
            % Create the model
            bst_model = fitensemble(Xts,yt,tree_method,iter,tree,...
                'PredictorNames',FeatNames,'LearnRate',lr,...
                'type','classification');
            % Store features importance
            imp(b,:) = bst_model.predictorImportance;
        case 'LLR',
            % Create the model
            modelo = glmnet(Xts,yt,'binomial',[]);
            % Find betas for lambda_opt
            [~, idx] = min(abs(modelo.lambda - lambda_opt));
            betas = modelo.beta(:,idx);
            imp(b,:) = abs(betas);
    end
    
    [~, sortedIndex_b] = sort(imp(b,:));
    
    sortedIndex(b,:) = sortedIndex_b;
    
    %%% 3. Compute error and ber for different number of features    
    for n=1:Nfeat
        idx = sortedIndex_b(n:Nfeat);
        
        switch method,
            case 'BST'
                reduced_model = fitensemble(Xts(:,idx),yt,tree_method,iter,tree,...
                    'PredictorNames',FeatNames(idx),'LearnRate',lr,...
                    'type','classification');
                ypred = predict(reduced_model,Xvs(:,idx));
                neg_value = -1;
            case 'LLR',
                reduced_model = glmfit(Xts(:,idx),yt,'binomial'); % Fit the data with logistic regression
                ypred = double((glmval(reduced_model,Xvs(:,idx),'logit'))>0.5);            
                neg_value = 0;
        end
        
        performance_metrics = compute_metrics(ypred,yv,neg_value);
        error(b,n) = performance_metrics(5);
        ber(b,n)   = performance_metrics(7);        
        bootstrap_metrics(b,n,:) = performance_metrics;
       
   end
        
end

%% (II) Select features according to the selected metric

switch metric
    case 'ber'
        err = ber;
    case 'error'
        err = error;
end

merr = mean(err);
serr = std(err);

[minerr, idx] = min(merr);
threshold = minerr + serr(idx);
opt_n = find(merr<threshold,1,'last');
Nselected = Nfeat-opt_n+1;

% Select the relevant features for each resample
preselected = sortedIndex(:, opt_n:Nfeat);
% Score = how many times each feature is selected
mean_score = hist(preselected(:), 1:Nfeat);
mean_score = mean_score ./ max(mean_score) * 100;
[~, sortedFeatures] = sort(mean_score,'descend');
feat_index = sortedFeatures(1:Nselected);


output.feat_index = feat_index;
output.score = mean_score;
output.mean_error = merr;
output.std_error = serr;
output.metric = metric;
output.params = params;
output.method = method;
output.bootstrap_metrics = bootstrap_metrics;
output.bootstrap_sortedIndex = sortedIndex;

end %Function



function [Xs,mu,sigma] = scale_data(X,mu,sigma)

if nargin <2
    mu = mean(X);
    sigma = std(X);
end

Xs = bsxfun(@minus,X,mu);
Xs = bsxfun(@times,Xs,1./sigma);

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
