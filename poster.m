% RRR explaining schematics
	n_samples = 100;
	rng(1);
	a = normrnd(0,0.5,n_samples,1);
	b = normrnd(0,0.5,n_samples,1);
	c = normrnd(0,0.5,n_samples,1);

	x = 0.5 * a + b;
	y = b + 0.1 * c;
	z = 0.1*x - 0.2 * b + c;
	x = a;
	y = b;
	z = c;
	center = mean([x y z],2);

	%% plot
	coeff = [1 1 1 ; 1 1 0;1 0.8 -1];
	sz    = [421 258;360 124;360 124];
	% sz    = [421 258;182 122;182 122];
	% cmap  = {'RdGy','RdYlGn','BrBG'};
	cmap  = {cbrewer2('RdGy')};
	cmap  = {flip(cmap{1},1), cmap{1}(:,[2 1 3]), cmap{1}(:,[3 2 1])};
	ltype = {'-','--','-.'};
	ip = 3;
	close all; h = [];

	% artificial left / right trials
		latent = [1 -1 0.5];
		L_ind = [x y z]*latent'>0;
		choice_center = [mean([x(L_ind) y(L_ind) z(L_ind)],1); mean([x(~L_ind) y(~L_ind) z(~L_ind)],1)]*3;
		scatter3(x,y,z,[],L_ind,'filled'); hold on;
		colormap([0 0 1; 1 0 0]);
		h(end+1) = plot3(choice_center(:,1),choice_center(:,2),choice_center(:,3),'k--','LineWidth',1.5); % choice axis
		grid on; hold on;
		ip = 1;

	% comm axis
	% scatter3(x,y,x-y,[],'filled'); hold on;
	scatter3(x,y,z,[],[x,y,z]*coeff(ip,:)','filled'); hold on;
	c = colorbar; c.Ticks = []; 
	% if ip==1
		c.Label.String = ['Activity of deep neuron ' num2str(ip)];
	% else
		% c.Position = [0.93    0.1100    0.04    0.8150];
	% end
	h(end+1) = plt_coeff(coeff(ip,:),center,['k' ltype{ip}],'LineWidth',1.5); % communication axis
	ax = gca;
	set(ax,'XLim',[-1.2686 1.2381],'YLim',[-1.1550 1.2027],'ZLim',[-1.0000 1.0354],...
		'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'View',[-73.4862 31.8924]);

	f = gcf; 
	f.Position([3 4]) = sz(1,:);
	colormap(cmap{ip});




	return
	f.Position([3 4]) = [165 121];
	colormap(cbrewer2('Accent',2));
	l = legend(h,{'cue axis','comm axis'},'box','off');


	% plot subspace
	hold on;
	h = arrayfun(@(i) plt_coeff(coeff(i,:),center,['k' ltype{i}]), 1:3);
	set(gca,'XLim',[-1.2686 1.2381],'YLim',[-1.1550 1.2027],'ZLim',[-1.0000 1.0354],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'View',[-73.4862 31.8924]);
	grid on;
	l = legend(h,{'deep cell 1','deep cell 2','deep cell 3'},'box','off','Position',[0.5709    0.6750    0.3801    0.2771]);
	l.Title.String = 'Communication Axes';
	f = gcf; f.Position([3 4]) = sz(1,:) + [-50 0];

% communication is low rank and restricted to cue info
	load('fig2a.mat');
	ax = np; colors = cbrewer2('spectral',size(cvloss,1));
	colors = [0.5 0.5 0.5;
			  0.0156 0.5117 0.9805;
			  0.7695 0.4688 0.0273;
			  0.5 0.5 0.5;
			  0.5 0.5 0.5];
	% h = arrayfun(@(i) errorbar((1:size(cvloss,3))-1, squeeze(1-cvloss(i,1,:)), squeeze(cvloss(i,2,:)), '.-', 'MarkerSize', 10,'LineWidth',2,'Color',colors(i,:),'MarkerFaceColor',colors(i,:)),3);
	h = arrayfun(@(i) errorbar((1:size(cvloss,3))-1, squeeze(1-cvloss(i,1,:)), squeeze(cvloss(i,2,:)), '.-', 'MarkerSize', 10,'LineWidth',2,'Color',colors(i,:),'MarkerFaceColor',colors(i,:)),1:size(cvloss,1));
	h(1).LineStyle = '--';
	h(4).LineStyle = '-.';

	arrayfun(@(h_) set(h_,'LineWidth',1),h([1 4 5]));
	plot([0 10],[0 0],'-','Color',[0.2 0.2 0.2 0.2],'LineWidth',0.7);
	set(ax,'YLim',[-0.4 0.4],'XTick',[0 5 11],'XTickLabel',{'0','5','ridge'})
	ylim([-0.4 0.4]);
	xlabel('# rank'); 
	ylabel('Performance (1-NSE)');

	arrayfun(@(i) scatter(11,ridge_loss(i),'^','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','None'), 1:5);

	l = legend(h,ops.win_name); l.Position = [0.1865    0.1619    0.3171    0.3340];
	l.Box = 'off';
	f = gcf; f.Position([3 4]) = [287 263];

% figure 2B

	% rank 2
	cv_2_plt = squeeze(cvloss(:,:,1,ops.dim==2));

	colors = [0.5 0.5 0.5;
			  0.0156 0.5117 0.9805;
			  0.7695 0.4688 0.0273;
			  0.5 0.5 0.5;
			  0.5 0.5 0.5];
	ax = np;
	plot(1-cv_2_plt','Color',[0.8 0.8 0.8 0.2]);
	arrayfun(@(i) my_scatter(i,1-cv_2_plt(:,i),ax,'filled','MarkerFaceColor',colors(i,:)), 1:size(cv_2_plt,2));

	set(ax,'XTick',1:5,'XTickLabel',ops.win_name,'XTickLabelRotation',45,'XLim',[0.5 5.5]);
	ylabel('rank 2 performance');
	f = gcf; f.Position([3 4]) = [287 270];

% figure C
	load 220708/tmp_XY; f = arrayfun(@(i) load(['220708/tmp_B' num2str(i) '.mat']), 1:3); n_left = 107; 
	rank_oi = 2;
	B = cat(3,f(:).B);
	B_ = squeeze(B(2:end,(1:238)+238*rank_oi,1));
	% B_ = my_RRR(X,Y);

	% sort by cell influence
	% input
		X2Y = sum(abs(B_),2);
		[~,I] = sort(X2Y);
		sorted_X = X(:,I);

	% output
		Y2X = sum(abs(B_),1);
		[~,I] = sort(Y2X);
		sorted_Y = Y(:,I);


	ax = np(2,1);
	axes(ax(1));
	colors = [0 0 1;1 0 0];
	n_2_plot = 10;
	dx = 0.1;
	% violinplot(tmp(1:n_left,1:10),[],'ViolinColor',[0 0.4470 0.7410]);
	% violinplot(tmp(n_left+1:end,1:10),[],'ViolinColor',[0.8500 0.3250 0.0980]);
	violinplot_poster(sorted_X(1:n_left,end-(n_2_plot-1):end),[],-dx,'ViolinColor',colors(1,:));
	violinplot_poster(sorted_X(n_left+1:end,end-(n_2_plot-1):end),[],dx,'ViolinColor',colors(2,:));
	axes(ax(2));
	violinplot_poster(sorted_Y(1:n_left,end-(n_2_plot-1):end),[],-dx,'ViolinColor',colors(1,:));
	violinplot_poster(sorted_Y(n_left+1:end,end-(n_2_plot-1):end),[],dx,'ViolinColor',colors(2,:));

	set(ax(2),'XLim',[0.5 n_2_plot+0.5],'Ytick',[],'XTickLabel',[],'YLim',[0 3]);
	xlabel(ax(2),{'Cells in deep ACC'});
	ylabel(ax(2),'Delay Activity');
	ax(2).Position(2) = 0.15;
	% ax(2).Position([1 2 4]) = [0.5 0.25 0.7];

	set(ax(1),'XLim',[0.5 n_2_plot+0.5],'YLim',ax(2).YLim,'XTickLabel',[],'YTickLabel',[]);
	% ax(1).Position([2 3 4]) = ax(2).Position([2 3 4]);
	xlabel(ax(1),{'Cells in superficial ACC'});
	ylabel(ax(1),'Delay Activity');
	text(ax(1),7-2.3,3,'Left Choice','Color',colors(1,:),'FontSize',13);
	text(ax(1),6.7-2.3,2.3,'Right Choice','Color',colors(2,:),'FontSize',13);

	f = gcf;
	f.Position([3 4]) = [200 267];

% figure C2
	% calculate B myself
	load fig3
	B_ = my_RRR(X,Y);

	% cue activity direction
	ind = {Epoch.correct, Epoch.incorrect};
	choice_ax = (mean(X(ind{1},:),1)-mean(X(ind{2},:),1))';
	theta_dist = arrayfun(@(ideep) vec_theta(B_(:,ideep),choice_ax), 1:size(B_,2));
	% vec_theta(choice_ax,comm_ax);

	% ax = np(2,1); 
	ax = np;
	clear h;
	h(1) = histogram(ax(2),arrayfun(@(i) vec_theta(mean(Y(ind{1},:),1)-mean(Y(ind{2},:),1),sum(B_(i,:),1)), 1:size(B_,1)),'BinEdges',0:5:90,'FaceColor',[0.7 0.7 0.7]);
	h(2) = histogram(ax(4),arrayfun(@(i) vec_theta(mean(X(ind{1},:),1)-mean(X(ind{2},:),1),sum(B_(:,i),2)'), 1:size(B_,2)),'BinEdges',0:5:90,'FaceColor',[0.7 0.7 0.7]);


	set(ax(2),'XTick',[],'XLim',[0 90]);
	% ylabel(ax(2),'# superficial cells');

	set(ax(4),'XTick',0:30:90,'XLim',[0 90]);
	
	xlabel(ax(4),'\theta_{cue}-\theta_{comm}');
	% ylabel(ax(4),'# deep cells');

	f = gcf;
	f.Position([3 4]) = [400 267];

	ax(2).Position([2 3 4]) = ax(1).Position([2 3 4]);
	ax(4).Position([2 3 4]) = ax(3).Position([2 3 4]);
	% ax(3).Position([2 4]) 

%% example of random
	Dims = [2 3 10 50 110];
	ax = np(numel(Dims));
	for id = 1:numel(Dims)
		nd = Dims(id);
		v = zeros(nd,1); 
		v(1) = 1;

		randv = randn(1000,nd);

		theta = arrayfun(@(i) vec_theta(randv(i,:),v'), 1:size(randv,1));

		histogram(ax(id),theta);
		title(ax(id),[num2str(nd) ' dim']);
	end

	arrayfun(@(h) set(h,'XLim',[0 90]), ax);

function theta = vec_theta(u,v)

	if prod(size(u)==size(v)) == 0
		error('size mismatch');
	end
	if size(u,1) ~= 1
		u = u';
		v = v';
	end
	theta = acosd(u*v' / (norm(u)*norm(v)));
	if theta>90 theta = 180-theta; end
end