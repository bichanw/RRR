% inspect Ryan's data
addpath(genpath('communication-subspace'));
load('/jukebox/scratch/bichanw/Results.mat');
% load('/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_47/11092021_T10/Results.mat');

Xdeep = catcell(Epoch.Spike_Epoch_sm(Deep==1));
Xsup  = catcell(Epoch.Spike_Epoch_sm(Deep==0));

% run CCA
% [A,B,r,U,V,stats] = canoncorr(X,Y);
iepos = 26:40;
X = squeeze(mean(Xsup(:,iepos,:),2));
Y = squeeze(mean(Xdeep(:,iepos,:),2));



% run ridge
	ops = struct('method','ridge','if_plot',true,'lamdba',[0:50:500]);
	ax = np; [cvloss,ops] = RRR(X,Y,ops); xlabel('\lambda');
	export_fig tmp.pdf

	% use their own function to determine lambda in ridge?
	[lambdaOpt,cvloss] = RegressModelSelect(@RidgeRegress, Y, X, 0:.1:20,'Scale', false, 'LossMeasure', 'NSE');
	
	
% run RRR
	ops = struct('dim',0:10,'if_plot',true);
	[cvloss,ops] = RRR(X,Y,ops);

	% examine raw results
	ind = [Epoch.left Epoch.right...
		   find(~ismember(1:size(X,1),[Epoch.left Epoch.right]))]; % early ending trial
	X = X(ind,:);
	Y = Y(ind,:);
	n_left  = numel(Epoch.left);
	n_right = numel(Epoch.right);
	save tmp_XY X Y


	
return


% use rank 1, B1 to sort cell
	load 220708/tmp_XY; f = arrayfun(@(i) load(['220708/tmp_B' num2str(i) '.mat']), 1:3); n_left = 107;
	rank_oi = 2;
	B = cat(3,f(:).B);
	B_ = squeeze(B(2:end,(1:238)+238*rank_oi,1));
	B_ = my_RRR(X,Y);
	Y_pred = X*B_;

	% sort by cell influence
	% input
		X2Y = sum(abs(B_),2);
		[~,I_X] = sort(X2Y);
		sorted_X = X(:,I);

	% output
		Y2X = sum(abs(B_),1);
		[~,I_Y] = sort(Y2X);
		sorted_Y = Y(:,I);


	% output predicted
	tmp = sorted_Y;
	tmp = Y_pred(:,I);


	% plotting

	% plot d prime as hist
		ax = np(2);
		histogram(ax(1),my_dp(sorted_X(1:n_left,1:30),sorted_X(n_left+1,1:30)));
		histogram(ax(1),my_dp(sorted_X(1:n_left,56:110),sorted_X(n_left+1,56:110)));

		histogram(ax(2),my_dp(sorted_Y(1:n_left,1:119),sorted_Y(n_left+1,1:119)));
		histogram(ax(2),my_dp(sorted_Y(1:n_left,120:238),sorted_Y(n_left+1,120:238)));


	ind = {Epoch.correct, Epoch.incorrect};
	% plot d prime as line
		ax = np(2,1);
		plot(ax(1),abs(my_dp(X(ind{1},I_X),X(ind{2},I_X))));
		plot(ax(2),abs(my_dp(Y(ind{1},I_Y),Y(ind{2},I_Y))));
		plot(ax(2),abs(my_dp(Y_pred(ind{1},I_Y),Y_pred(ind{2},I_Y))));
		% plot(ax(2),abs(my_dp(sorted_Y(1:n_left,:),sorted_Y(n_left+1:end,:))));
		% plot(ax(2),abs(my_dp(Y_pred(1:n_left,I),Y_pred(n_left+1:end,I))));
		title(ax(1),"x d'");
		title(ax(2),"y d'");
		legend(ax(2),{'y','y_{pred}'});

	% plot raw activity
		ax = np;
		imagesc(tmp);
		set(ax,'XLim',[0.5 size(tmp,2)+0.5], 'YLim',[0.5 size(tmp,1)+0.5],'YTick',[50 150],'YTickLabel',{'L','R'});
		xlabel(ax,'cells with increasing B');
		hline(n_left,ax,'linespec','k--','linewidth',0.75);
		norm_color;

	% plot distribution
		ax = np(2);
		axes(ax(1));
		% violinplot(tmp(1:n_left,1:10),[],'ViolinColor',[0 0.4470 0.7410]);
		% violinplot(tmp(n_left+1:end,1:10),[],'ViolinColor',[0.8500 0.3250 0.0980]);
		violinplot(sorted_X(1:n_left,end-9:end),[],'ViolinColor',[0 0.4470 0.7410]);
		violinplot(sorted_X(n_left+1:end,end-9:end),[],'ViolinColor',[0.8500 0.3250 0.0980]);
		axes(ax(2));
		violinplot(sorted_Y(1:n_left,end-9:end),[],'ViolinColor',[0 0.4470 0.7410]);
		violinplot(sorted_Y(n_left+1:end,end-9:end),[],'ViolinColor',[0.8500 0.3250 0.0980]);

		set(ax(2),'Ytick',[],'XTickLabel',[]);
		xlabel(ax(2),{'Cells in deep ACC'});
		ax(2).Position([1 2 4]) = [0.5 0.25 0.7];

		set(ax(1),'YLim',ax(2).YLim,'XTickLabel',[]);
		ax(1).Position([2 3 4]) = ax(2).Position([2 3 4]);
		xlabel(ax(1),{'Cells in superficial ACC'});


		% figure setting
		f = gcf; f.Position([3 4]) = [2000 300];
		arrayfun(@(h) set(h,'XTick',[]), ax);
		export_fig tmp.png -m3


		M(1,:) = mean(tmp(1:n_left,:),1);
		V(1,:) = std(tmp(1:n_left,:),[],1);
		M(2,:) = mean(tmp(n_left+1:end,:),1);
		V(2,:) = std(tmp(n_left+1:end,:),[],1);
		ax = np;
		arrayfun(@(i) plot(ax,M(i,:),'-','Color',colors(i,:)),1:2)


% examine B
	f = arrayfun(@(i) load(['tmp_B' num2str(i) '.mat']), 1:3);
	B = cat(3,f(:).B);

	[ax,r,c] = np(1,3);

	arrayfun(@(i) imagesc(ax(i),squeeze(B(2:end,(1:238)+238,i))), 1:(c-1));
	imagesc(ax(3),squeeze(B(2:end,(1:238)+238,1)-B(2:end,(1:238)+238,2)));
	title(ax(1),'B1'); title(ax(2),'B2'); title(ax(3),'B1-B2');

	colormap(flip(cbrewer2('RdBu')));
	arrayfun(@(h) set(h,'CLim',[-1 1]*max(abs(h.CLim)),'YDir','reverse',...
					    'XLim',[0.5 238.5],'YLim',[0.5 size(B,1)-0.5]), ax);
	c = arrayfun(@(h) colorbar(h), ax);
		
	imagesc(B(2:end,(1:238)+238));

% examine rank of activity
	[coeff,score,latent,tsquared,explained_x,mu] = pca(X);
	[coeff,score,latent,tsquared,explained_y,mu] = pca(Y);
	ax = np;
	plot(cumsum(explained_x));
	plot(cumsum(explained_y));
	xlabel('# dim'); ylabel('% explained');
	plot([45 45;0 45]',[0 83.6; 83.6 83.6],'k--','LineWidth',0.5);
	ax.XTick = [0 45 100 200];
	ax.YTick = [0 83.6 100];
	ax.Position = [0.25    0.2500    0.65    0.600];
	l = legend({'sup','deep'},'Box','off');
	export_fig tmp.pdf


% pick correlation with choice
	iepo = 30;
	X = squeeze(Xsup (:,iepo,:));
	Y = squeeze(Xdeep(:,iepo,:));


	iepo = 30;
	X = Xsup;
	NCells = size(X,3);
	[ax,r,c] = np(NCells);
	for i=1:NCells
		histogram(ax(i),squeeze(X(Epoch.left,iepo,i)));
		histogram(ax(i),squeeze(X(Epoch.right,iepo,i)));
	end
	export_fig tmp.pdf

% randomly select 10 cells from each?
	iepo = 30;
	ops.if_plot = true;
	RRR(squeeze(Xsup(:,iepo,randperm(size(Xsup,3),10))),squeeze(Xdeep(:,iepo,randperm(size(Xdeep,3),10))),ops);
	% RRR(squeeze(Xsup(:,iepo,:)),squeeze(Xdeep(:,iepo,:)),ops);

return

% calculate correlation with choice value?




% PCA
	NPC = 10;
	Xsup  = permute(catcell(arrayfun(@(iepo) RR(squeeze(Xsup(:,iepo,:)),10), 1:77,'UniformOutput',false)), [1 3 2]);
	Xdeep = permute(catcell(arrayfun(@(iepo) RR(squeeze(Xdeep(:,iepo,:)),10), 1:77,'UniformOutput',false)), [1 3 2]);


return


% run cvloss across all tp
	ops = struct('method','RRR','dim',0:10);
	cvloss = catcell(arrayfun(@(iepo) RRR(squeeze(Xsup(:,iepo,:)),squeeze(Xdeep(:,iepo,:)),ops), 1:77,'UniformOutput',false));

% plot cvloss lines
	colors = cbrewer2('spectral',size(cvloss,2));
	ax = np;
	arrayfun(@(i) errorbar(1:77, squeeze(1-cvloss(1,i,:)), squeeze(cvloss(2,i,:)), 'o--', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6),1:size(cvloss,2));
	vline([5.5 25.5 40.5 47.5],ax,'linewidth',0.7,'linespec','k--'); % mark events
	% arrayfun(@(i) plot(squeeze(1-cvloss(1,i,:))','color',colors(i,:)), 1:size(cvloss,2));
	colormap(colors); c = colorbar; c.Label.String = '# rank'; set(c,'Ticks',[0 1],'TickLabels',{'0','full'});

% plot cvloss image
	ax = np;
	imagesc(1:77,0:10,squeeze(1-cvloss(1,:,:)));
	colormap(flip(cbrewer2('RdBu')));
	set(ax,'CLim',[-1 1]*max(abs(ax.CLim)),'YDir','reverse','YLim',[-0.5 10.5],'XLim',[0.5 77.5]);
	c = colorbar; 
	xlabel('time epoch'); ylabel('rank included');

return
% compare cvloss of 2 tp
	clear cvloss; 
	iepo = 30;cvloss(1,:,:) = RRR(squeeze(Xsup(:,iepo,:)),squeeze(Xdeep(:,iepo,:)),ops);
	iepo = 10;[cvloss(2,:,:),ops] = RRR(squeeze(Xsup(:,iepo,:)),squeeze(Xdeep(:,iepo,:)),ops);
	
	ax = np;
	errorbar(ax,ops.dim,1-squeeze(cvloss(1,1,:)),squeeze(cvloss(1,2,:)),'.--');
	errorbar(ax,ops.dim,1-squeeze(cvloss(2,1,:)),squeeze(cvloss(2,2,:)),'.--');
	export_fig tmp.pdf

% examine rank of activity
	for i = 1:77
		[coeff{i},score{i},latent{i},tsquared{i},explained{i},mu{i}] = pca(squeeze(Xdeep(:,i,:)));
	end
	explained = catcell(explained);
	plot(cumsum(explained,1));



% or process all epochs
	ops.if_plot = true;
	rank(reshape(Xsup(:,1:77,:),size(Xsup,1),[])); % 278x8470
	rank(reshape(Xdeep(:,1:77,:),size(Xdeep,1),[])) % 278x8470
	RRR(reshape(Xsup(:,1:77,:),size(Xsup,1),[]),reshape(Xdeep(:,1:77,:),size(Xdeep,1),[]),ops);

% plot input and output 
	ax = np(2); colormap(flip(cbrewer2('RdBu')));
	imagesc(ax(1),X');imagesc(ax(2),Y');
	arrayfun(@(h) set(h,'XLim',[0 size(X,1)],'CLim',[-1 1]*max(abs(h.CLim))),ax);


% reduced rank: examine predicted vs real
	ax = np;
	arrayfun(@(i) my_scatter(i,Y(:,i),ax), 1:size(Y,2));


	% variance problem?
		B_ = Bfull*V*V';
		ax = np(1,3); colormap(flip(cbrewer2('RdBu')));
		imagesc(ax(1),Bfull);
		imagesc(ax(2),B_);
		imagesc(ax(3),B_-Bfull);
		title(ax(1),'Bfull');
		title(ax(2),'B_');
		title(ax(3),'B\_-Bfull');

		arrayfun(@(h) colorbar(h), ax);
		arrayfun(@(h) set(h,'CLim',[-1 1]*2.5,'XLim',[0.5 size(Bfull,2)],'YLim',[0.5 size(Bfull,1)]), ax);

	% delta Yhat - Y
		Yhat = B; Y = B_;
		n = 1320/110;
		[ax,r,c] = np(n); colormap(flip(cbrewer2('RdBu')));
		arrayfun(@(i) imagesc(ax(i),Y-Yhat(:,(1:110)+(i-1)*110)), 1:n);
		arrayfun(@(i) colorbar(ax(i)), 1:n);
		arrayfun(@(i) set(ax(i),'CLim',[-50 50],'YLim',[0 size(Y,1)]), 1:n); %arrayfun(@(i) set(ax(i),'CLim',[-15 15],'YLim',[0 size(Y,1)]), 1:n);
		ind = sub2ind([c r],1,r);
		xlabel(ax(ind),'cell ID'); ylabel(ax(ind),'trial ID');

	% line plot
		ax = np;
		colors = cbrewer2('spectral',10);
		predicted_ymean = catcell(arrayfun(@(i) mean(Yhat(:,(1:110)+i*110),1), 0:9,'UniformOutput',false));

		colormap(colors);
		arrayfun(@(i) plot(predicted_ymean(i,:),'Color',colors(i,:)), 1:10);

		plot(ax,mean(Y,1),'k','LineWidth',2);
