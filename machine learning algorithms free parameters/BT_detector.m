function BT

% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

clc;

% load data
w_length = [4, 8];
info.k = 1;         % counter

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
    info.ddbb = sprintf('PUBLIC-%ds',w_length(j));
    X = T{public_db,FeatNames};
    X = scale_data(X);
    y = T.y(public_db);
    BoostedTrees(X,y,info,FeatNames);
    info.k = info.k + 1;
        
    % ohca dbs
    info.ddbb = sprintf('OHCA-%ds',w_length(j));
    X = T{ohca_db,FeatNames};
    X = scale_data(X);
    y = T.y(ohca_db);
    BoostedTrees(X,y,info,FeatNames);
    info.k = info.k + 1;
                
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

function [Xs,mu,sigma] = scale_data(X,mu,sigma)

if nargin <2
    mu = mean(X);
    sigma = std(X);
end

Xs = bsxfun(@minus,X,mu);
Xs = bsxfun(@times,Xs,1./sigma);

end

function BoostedTrees(X,y,info,FeatNames)

% compare ensemble methods
method = {'GentleBoost','AdaBoostM1', 'LogitBoost'};                        

% constant
lr    = 1;               %learning rate
iter  = 50;
color = 'brgkm';

rng('default') % for reproducibility
tree  = templateTree('minleaf',100); % a measure of the max depth of the tree

for ii=1:length(method)
    
    
    rng('default') % for reproducibility       
    cv_model = fitensemble(X,y,method{ii},iter,tree,...
        'PredictorNames',FeatNames,'LearnRate',lr,...
        'type','classification','kfold',5);
       
    figure(3)
    subplot(2,2,info.k)
    Lcv = kfoldLoss(cv_model,'mode','cumulative');
    plot(Lcv,color(ii));
    hold on;
    

end
legend(method)
xlabel('Number of trees','Interpreter','Latex','Fontsize',14);
ylabel('MSE','Interpreter','Latex','Fontsize',14);
text_title = sprintf('%s',info.ddbb); 
title(text_title,'Interpreter','Latex')
grid on;


% compare minleaf values
method = 'LogitBoost';                        

% constant
lr    = 1;               %learning rate
iter  = 50;
color = 'brgkmc';

minleaf = [1,5,10,50,100,200];

for ii=1:length(minleaf)
    
    
    rng('default') % for reproducibility
    tree  = templateTree('minleaf',minleaf(ii)); % a measure of the max depth of the tree
    
    rng('default') % for reproducibility
    cv_model = fitensemble(X,y,method,iter,tree,...
        'PredictorNames',FeatNames,'LearnRate',lr,...
        'type','classification','kfold',5);
       
    figure(4)
    subplot(2,2,info.k)
    Lcv = kfoldLoss(cv_model,'mode','cumulative');
    plot(Lcv,color(ii));
    hold on;
    

end
legend('1','5','10','50','100','200');
xlabel('Number of trees','Interpreter','Latex','Fontsize',14);
ylabel('MSE','Interpreter','Latex','Fontsize',14);
text_title = sprintf('%s',info.ddbb); 
title(text_title,'Interpreter','Latex')
grid on;

% conclusion: lr = 1, Ntree = 50, LogitBoost, minleaf = 100;
end