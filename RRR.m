function [cvLoss,ops] = RRR(X,Y,ops)
% scripts adapted from RRR from Semedo et al.

addpath(genpath('/usr/people/bichanw/SpikeSorting/Codes/witten/communication-subspace'));
if nargin < 3
	ops = struct;
end
ops.method  = getOr(ops,'method','RRR');
ops.if_plot = getOr(ops,'if_plot',false);
ops.nfolds  = getOr(ops,'nfolds',10);

SET_CONSTS
% zero-mean X, Y
X = X - mean(X,1);
Y = Y - mean(Y,1);

% make sure there is enough number of trials
min_required_tr = max(ceil([size(X,2)/(ops.nfolds-1)*ops.nfolds size(Y,2)/(ops.nfolds-1)*ops.nfolds]));
if size(X,1) < min_required_tr
	error('Trial number is not enough to do crossval');
end



%% Cross-validate Reduced Rank Regression

% choose regression method
if strcmpi(ops.method,'RRR')
	ops.dim = getOr(ops,'dim',[0:10 size(Y,2)]);
	numDimsUsedForPrediction = ops.dim;
	regressMethod = @ReducedRankRegress;
elseif strcmpi(ops.method,'ridge')
	ops.lambda = getOr(ops,'lamdba',0:0.2:2);
	numDimsUsedForPrediction = ops.lambda;  
	regressMethod = @RidgeRegress;
else
	error('Not recognizing regression method');
end

% Number of cross validation folds.
cvNumFolds = ops.nfolds;

% Initialize default options for cross-validation.
cvOptions = statset('crossval');

% If the MATLAB parallel toolbox is available, uncomment this line to
% enable parallel cross-validation.
% cvOptions.UseParallel = true;

 

% Auxiliary function to be used within the cross-validation routine (type
% 'help crossval' for more information). Briefly, it takes as input the
% the train and test sets, fits the model to the train set and uses it to
% predict the test set, reporting the model's test performance. Here we
% use NSE (Normalized Squared Error) as the performance metric. MSE (Mean
% Squared Error) is also available.
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
	numDimsUsedForPrediction, 'LossMeasure', 'NSE','Scale',false);

% Cross-validation routine.
cvl = crossval(cvFun, Y, X, 'KFold', cvNumFolds, 'Options', cvOptions);

% Stores cross-validation results: mean loss and standard error of the
% mean across folds.
cvLoss = [ mean(cvl); std(cvl)/sqrt(cvNumFolds) ];

% To compute the optimal dimensionality for the regression model, call
% ModelSelect:
ops.optDimReducedRankRegress = ModelSelect...
	(cvLoss, numDimsUsedForPrediction);

% Plot Reduced Rank Regression cross-validation results
x = numDimsUsedForPrediction;
y = 1-cvLoss(1,:);
e = cvLoss(2,:);

if ops.if_plot
	errorbar(x, y, e, 'o--', 'Color', COLOR(V2,:), ...
	    'MarkerFaceColor', COLOR(V2,:), 'MarkerSize', 10);

	xlabel('Number of predictive dimensions');
	ylabel('Predictive performance');
end



return


