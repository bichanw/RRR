function [ops,cvloss] = averaged_RRR(Xsup,Xdeep,ops)
if nargin < 3
	ops = struct('dim',0:10,'if_plot',false);
end

% RRR across multiple period
ops.twin = getOr(ops,'twin',[1 5;6 25;26 40;41 47;48 77]);
ops.win_name = getOr(ops,'win_name',{'start','cue','delay','arm','outcome'});

% RRR performance
clear cvloss ridge_loss
for iwin = 1:numel(ops.win_name)
	X = squeeze(mean(Xsup(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));
	Y = squeeze(mean(Xdeep(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));

	% RRR
	[cvloss(iwin,:,:),ops] = RRR(X,Y,ops);

	% % ridge for full model
	% lambda = 0:.1:20;
	% [lambdaOpt,tmp] = RegressModelSelect(@RidgeRegress, Y, X, lambda,'Scale', false, 'LossMeasure', 'NSE');
	% ridge_loss(iwin) = 1-tmp(1,lambda==lambdaOpt);
end


% plot loss by rank
ax = np;colors = cbrewer2('spectral',size(cvloss,3));
arrayfun(@(i) errorbar(1:numel(ops.win_name), squeeze(1-cvloss(:,1,i)), squeeze(cvloss(:,2,i)), 'o--', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 3),1:size(cvloss,3));
% figure setting
set(ax,'XLim',[0.5 numel(ops.win_name)+0.5],'XTick',1:numel(ops.win_name),'XTickLabel',ops.win_name,'XTickLabelRotation',45);
colormap(colors); c = colorbar(ax);c.Label.String = '# rank'; % set(c,'Ticks',[0 1],'TickLabels',{'0',num2str(ops.dim(end))});
set(gcf,'Position',[0 0 250 180]);


% plot loss by period
% ax = np; colors = cbrewer2('spectral',size(cvloss,1));
% arrayfun(@(i) errorbar(1:size(cvloss,3), squeeze(1-cvloss(i,1,:)), squeeze(cvloss(i,2,:)), 'o--', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 3),1:size(cvloss,1));
% legend(ops.win_name);