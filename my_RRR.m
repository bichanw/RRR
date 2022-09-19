function [B_,coeff] = my_RRR(X,Y,rank_oi)
% scripts written by Bichan for RRR
% easier to pull out coefficient

if nargin<3
	rank_oi = 1:2;
end

Bfull = inv(X'*X) * X' *Y;
Yhat  = X*Bfull;
[coeff,score,latent,tsquared,explained,mu] = pca(Yhat);
B_ = Bfull * coeff(:,rank_oi)* coeff(:,rank_oi)';

