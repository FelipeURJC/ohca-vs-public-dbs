function BAG_RF

% this code, by Felipe Alonso-Atienza
% felipe.alonso@urjc.es

close all; clc;

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
    BaggingTrees(X,y,info,FeatNames);
    info.k = info.k + 1;
        
    % ohca dbs
    info.ddbb = sprintf('OHCA-%ds',w_length(j));
    X = T{ohca_db,FeatNames};
    X = scale_data(X);
    y = T.y(ohca_db);
    BaggingTrees(X,y,info,FeatNames);
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


function BaggingTrees(X,y,info,FeatNames)


% free parameters
Ntree = 50;

min_leaf = [1,5,10,50]; % a measure of the max depth of the tree
Nleaf    = numel(min_leaf);

color = 'bkrg';

for i=1:Nleaf
     
    % random forest
    rng('default'); % for reproducibility
    rf_tree = TreeBagger(Ntree,X,y,...
        'Method', 'classification',...
        'OOBPred','On',...
        'minleaf', min_leaf(i),...
        'prior', [0.8 0.2]); %[0.8 0.2]
            
    % bagging
    rng('default'); % for reproducibility
    bg_tree = TreeBagger(Ntree,X,y,...
        'Method', 'classification',...
        'OOBPred','On',...
        'minleaf', min_leaf(i),...
        'NvarToSample','all',...        % Bagging
        'prior', [0.8 0.2]);
                           
    figure(1)
    subplot(2,2,info.k)
    plot(rf_tree.oobError,color(i),'linewidth',2); hold on;
    plot(bg_tree.oobError,[color(i) 'o-'],'linewidth',2);
        
end

xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');
legend({'rf 1', 'bg 1', 'rf 5', 'bg 5', 'rf 10', 'bg 10', 'rf 50', 'bg 50'},...
    'Location','NorthEast');
grid on; axis tight; 

text_title = sprintf('%s',info.ddbb); 
title(text_title)

% variable importance
Nleaf = 10; % from previous analysis we select Nleaf = 10;

% random forest
rng('default'); % for reproducibility
rf_tree = TreeBagger(Ntree,X,y,...
    'Method', 'classification',...
    'OOBPred','On',...
    'minleaf', Nleaf,...
    'oobvarimp','on',...
    'prior', [0.8 0.2]);
        
% bagging
rng('default'); % for reproducibility
bg_tree = TreeBagger(Ntree,X,y,...
    'Method', 'classification',...
    'OOBPred','On',...
    'minleaf', Nleaf,...
    'oobvarimp','on',...
    'NvarToSample','all',...        
    'prior', [0.8 0.2]);
        
figure(2);
subplot(2,2,info.k)
barh(rf_tree.OOBPermutedVarDeltaError);
xlabel('Out-of-bag feature importance');
set(gca,'YTick',1:length(FeatNames));
set(gca,'YTickLabel',FeatNames,'Fontsize',8);
text_title = sprintf('RF: %s',info.ddbb); 
title(text_title)
axis tight;

figure(3);
subplot(2,2,info.k)
barh(bg_tree.OOBPermutedVarDeltaError);
xlabel('Out-of-bag feature importance');
set(gca,'YTick',1:length(FeatNames));
set(gca,'YTickLabel',FeatNames,'Fontsize',8);
text_title = sprintf('BAG: %s',info.ddbb); 
title(text_title)
axis tight;
end