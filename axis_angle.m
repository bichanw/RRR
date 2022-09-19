addpath(genpath('communication-subspace'));
addpath(genpath('/Users/bichanwu/Desktop/classes/Rotation Witten lab/Codes/communication-subspace'));
% load('/jukebox/scratch/bichanw/Results.mat');
load('/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_47/11092021_T10/Results.mat');

Xdeep = catcell(Epoch.Spike_Epoch_sm(Deep==1));
Xsup  = catcell(Epoch.Spike_Epoch_sm(Deep==0));

% some manipulation to save less data
tmp = struct('left',Epoch.left,'right',Epoch.right,'correct',Epoch.correct,'incorrect',Epoch.incorrect);
Epoch = tmp;
save('tmp.mat');

% RRR across multiple period
ops = struct;
ops.twin = getOr(ops,'twin',[1 5;6 25;26 40;41 47;48 77]);
ops.win_name = getOr(ops,'win_name',{'start','cue','delay','arm','outcome'});

% calculate distribution
ind = {Epoch.left, Epoch.right};
ind = {Epoch.correct, Epoch.incorrect};
theta_dist = [];
for iwin = 3
	% iwin = 3;
	X = squeeze(mean(Xsup(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));
	Y = squeeze(mean(Xdeep(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));

	X = X - mean(X,1);
	Y = Y - mean(Y,1);
	
	% calculate communication axis
		% reduced rank regression
		[B_,coeff] = my_RRR(X,Y);

		% ridge regression
		% lambda = 0:.1:20;
		% [lambdaOpt,tmp] = RegressModelSelect(@RidgeRegress, Y, X, lambda,'Scale', false, 'LossMeasure', 'NSE');
		% ridge_loss(iwin) = 1-tmp(1,lambda==lambdaOpt);
		% B_ = RidgeRegress(Y, X, lambdaOpt, 'Scale',false);
		% B_ = B_(2:end,:); % remove the intercept

	% cue activity direction
		% cloud center location
		% choice_ax = (mean(X(ind{1},:),1)-mean(X(ind{2},:),1))';

		% d prime
		% choice_ax = my_dp(X(ind{1},:),X(ind{2},:))';

		% what if using svm to determine choice axis
		X_svm = X([ind{1} ind{2}],:);
		choice = [ones(size(ind{1})) 2*ones(size(ind{2}))];
		Mdl = fitcsvm(X_svm,choice);
		choice_ax = Mdl.Beta;

	theta_dist(end+1,:) = arrayfun(@(ideep) vec_theta(B_(:,ideep),choice_ax), 1:size(B_,2));
end

% eigenvector vs choice ax
	[U,S,V] = svd(B_); % A = U*S*V'

	% U(:,1) U(:,2) column vectors as eigenbasis
	alpha = 0:.01:1;
	theta_base_choice = arrayfun(@(a) vec_theta(U(:,1)*a + U(:,2)*(1-a),choice_ax), alpha);
		
	% min angle
	[m,I] = min(theta_base_choice);
	fprintf('min angle as %.1f when alpha = %.2f\n',m,alpha(I));

	% plot
	% ax = np;
	% plot(alpha,theta_base_choice);
	% xlabel('\alpha','Interpreter','tex'); 
	% ylabel('\theta','Interpreter','tex'); 

return
	% plot 5 distributions
		% ax = np(5);
		% arrayfun(@(i) histogram(ax(i),theta_dist(i,:)),1:5);
		% arrayfun(@(h) set(h,'XLim',[0 90]), ax);
		% arrayfun(@(i) title(ax(i), ops.win_name{i}), 1:5);

% permute and rerun
	% ind = randperm(size(X,1));
	% X_perm = X(randperm(size(X,1)),:);
	% Y_perm = Y(randperm(size(X,1)),:);
	% B_permuted  = my_RRR(X_perm,Y_perm);
	% choice_perm = (mean(X_perm(ind{1},:),1)-mean(X_perm(ind{2},:),1))';
	% theta_perm  = arrayfun(@(ideep) rem(vec_theta(B_(:,ideep),choice_perm),90), 1:size(B_,2));

% ax = np;
% scatter(arrayfun(@(ideep) vec_theta(B_(:,ideep),choice_ax), 1:size(B_,2)),...
% 		arrayfun(@(ideep) vec_theta(B_(:,ideep),choice_perm), 1:size(B_,2)))


% choice selectivity vs. inference
	% sort by influence
	[~,I] = sort(sum(abs(B_),1));


	% plot 5 to see
	% ax = np(9); 
	% for j =0:numel(ax)-1
	% 	% mdl = fitlm(abs(choice_ax),abs(B_(:,I(end-j))));
	% 	mdl = fitlm(choice_ax,B_(:,I(end-j)));
	% 	axes(ax(j+1));
	% 	h = plot(mdl);
	% 	legend('off');

	% 	xlabel('choice selectivity');
	% 	ylabel('\beta','Interpreter','tex');
	% 	title('');
	% 	set(h(1),'Color',[0 0 0],'Marker','.');
	% 	% scatter(ax(j+1),abs(choice_ax),abs(B_(:,I(end-j))));
	% end


	% plot the most important one for poster
	mdl = fitlm(choice_ax,B_(:,I(end-2)));
	% ax = np;
	close all;
	h = plot(mdl);
	% xlabel('B_{outcome decoder}','Interpreter','tex');
	xlabel('B_{choice decoder}','Interpreter','tex'); 
	ylabel('B_{RRR}','Interpreter','tex'); title('');
	set(h(1),'Color',[0 0 0],'Marker','.');
	legend('off');
	box off;
	set(gca,'Position',[0.2545    0.3077    0.6505    0.6173]);
	set(gcf,'Position',[0 0 185 141]);
	export_fig tmp.eps -painters
	saveas(gcf,'tmp1.fig');

% plot angle between subspace and choice axis, all sessions
	load mat/theta_plane_choice
	colors = [0.5 0.5 0.5;
			  0.0156 0.5117 0.9805;
			  0.7695 0.4688 0.0273;
			  0.5 0.5 0.5;
			  0.5 0.5 0.5];

	ax = np;
	h = plot(theta_plane_choice','.-','Color',[0.2 0.2 0.2],'MarkerSize',10);
	plot([0.5 2.5],[90 90],'k--','LineWidth',0.7);
	set(ax,'XTick',[1 2],'XLim',[0.5 2.5],'YLim',[35 91],'YTick',[35 60 90],'XTickLabel',{'choice','outcome'});
	ylabel('\theta_{comm plane, behavior axis}');



% angle between permutation axis and original choice axis

	rng shuffle;
	perm_dt = [];
	for i=1:2000
		% X_perm = X(randperm(size(X,1)),:);

		% permute L/R
			% choice_perm = (mean(X_perm(ind{1},:),1)-mean(X_perm(ind{2},:),1))';
			% perm_dt(end+1) = vec_theta(choice_perm,choice_ax);

		% just generate random orientation
			% perm_dt(end+1) = vec_theta(randn(size(choice_ax)),choice_ax);

		% shuffle weights
			perm_dt(end+1) = vec_theta(B_(randperm(size(B_,1)),randperm(size(B_,2),1)),choice_ax);
	end
	[N,edges] = histcounts(perm_dt,'Normalization','probability');
	edges = edges(1:end-1) + (edges(2)-edges(1))/2;
	% ax = np(2);
	% histogram(ax(1),perm_dt,'Normalization','probability');
	% axes(ax(2)); histfit(perm_dt);
	% export_fig tmp.pdf

	ax = np; 
	yyaxis left;
	histogram(theta_dist,'BinEdges',colon_right(48,5,90),'FaceColor',[0.7695 0.4688 0.0273]);
	vline(m,ax,'linespec','k--','linewidth',0.7);
	yyaxis right;
	plot(edges,N,'Color',[0.2 0.2 0.2]);
	ax.YAxis(1).Color = [0 0 0];
	ax.YAxis(2).Visible = 'off';
	% ax.Position(2) = 0.3;
	xlabel('\theta_{comm - choice}');
	set(gca,'Position',[ 0.1664    0.3182    0.6483    0.6068]);
	set(gcf,'Position',[0 0 159 124]);
	saveas(gcf,'tmp2.fig');
	export_fig tmp.eps -painters


	% test the angle difference among clustered 75 (only works for outcome of best session)
		% v_test = find(theta_dist<75);
		% theta_cluster = NaN(numel(v_test));
		% for i=1:numel(v_test)
		% 	for j=i+1:numel(v_test)
		% 		theta_cluster(i,j) = vec_theta(B_(:,v_test(i)), B_(:,v_test(j)));
		% 	end
		% end
		% c = colorbar; c.Label.String = 'angle between comm axes';


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