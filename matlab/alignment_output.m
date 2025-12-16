function [outputalignmentidx, commfrac] = alignment_output(X,Y,W,if_plot)
% alignment_output computes the alignment index of the output population
% outputalignmentidx = alignment_output(X,Y,W)
% Input:
%       X: input population (n_samples x n_input_features)
%       Y: output population (n_samples x n_output_features)
%       W: communication weights (n_input_features x n_output_features)
%       if_plot: whether to plot the variance profile
% Output:
%   outputalignmentidx: raw alignment index 

% have code here for random generation

if nargin < 4
    if_plot = false; % default is not to plot
end

% ----- calculate raw alignment index -----
% Do PCA on population
[upop,spopvec] = svd(cov(Y),'vector');
spopcum = cumsum(spopvec);
spopvec_nrm = spopvec/sum(spopvec);
spopcum_nrm = cumsum(spopvec_nrm);

% Project communication covariance onto PCs
cov_predicted = cov(X*W);
scomvec = diag(upop'*cov_predicted*upop);
% scomvec_nrm = scomvec./sum(scomvec);
scomcum_nrm = cumsum(scomvec)/sum(scomvec);


a_raw = spopvec' * scomvec;

% compute communication fraction
commfrac = sum(scomvec)/sum(spopvec);


% ----- normalize alignment index -----
% compute max possible alignment
totcom = sum(scomvec);
ii = find(spopcum>totcom+1e-10,1);  % find how many dimensions we'd need for maximally aligned communication
scommax = spopvec;
scommax(ii:end) = 0;
scommax(ii) = totcom-sum(scommax); % fix up final bin
a_max = spopvec' * scommax; % max possible alignment score

% compute min possible alignment (by flipping order of pop eigenvalues)
spopvec_rev = flipud(spopvec);
spopcum_rev = cumsum(spopvec_rev);
ii = find(spopcum_rev>totcom+1e-10,1);  % find how many dimensions we'd need for minimally aligned communication
scommin = spopvec_rev;
scommin(ii:end) = 0;
scommin(ii) = totcom-sum(scommin); % fix up final bin
scommin = flipud(scommin); % flip back to normal ordering
a_min = spopvec' * scommin; % min possible alignment score

% rescale alignment scores above and below zero
outputalignmentidx = (a_raw - a_min) / (a_max - a_min);

end