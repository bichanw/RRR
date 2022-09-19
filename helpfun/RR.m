function Xout = RR(X,nPC)
% reduce rank by PCA

coeff = pca(X);
Xout  = X* coeff(:,1:nPC);

end