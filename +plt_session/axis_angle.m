function m = axis_angle(Epoch,Deep)
	% concatenate activity
	Xdeep = catcell(Epoch.Spike_Epoch_sm(Deep==1));
	Xsup  = catcell(Epoch.Spike_Epoch_sm(Deep==0));

	% RRR across multiple period
	ops = struct;
	ops.twin = getOr(ops,'twin',[1 5;6 25;26 40;41 47;48 77]);
	ops.win_name = getOr(ops,'win_name',{'start','cue','delay','arm','outcome'});

	
	% extract activity during delay period
	iwin = 3;
	X = squeeze(mean(Xsup(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));
	Y = squeeze(mean(Xdeep(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));

	X = X - mean(X,1);
	Y = Y - mean(Y,1);

	% calculate communication axis
	B_ = my_RRR(X,Y);

	% choice vs. outcome
	ind = {{Epoch.left, Epoch.right},{Epoch.correct, Epoch.incorrect}};

	ax = np(1,2);
	for iPeriod = 1:2
		% calculate choice / outcome axis
		X_svm = X([ind{iPeriod}{1} ind{iPeriod}{2}],:);
		choice = [ones(size(ind{iPeriod}{1})) 2*ones(size(ind{iPeriod}{2}))];
		Mdl = fitcsvm(X_svm,choice);
		choice_ax = Mdl.Beta;

		% calculate angle between communication axes and choice axis 
		theta_dist = arrayfun(@(ideep) vec_theta(B_(:,ideep),choice_ax), 1:size(B_,2));

		% calculate angle between communication plane and choice axis 
		[U,S,V] = svd(B_); % A = U*S*V', U(:,1) U(:,2) column vectors as eigenbasis
		alpha = 0:.01:1;
		theta_base_choice = arrayfun(@(a) vec_theta(U(:,1)*a + U(:,2)*(1-a),choice_ax), alpha);
		% take min angle
		[m(iPeriod),I] = min(theta_base_choice);
		fprintf('min angle as %.1f when alpha = %.2f\n',m(iPeriod),alpha(I));


		% plot distribution and min angle
		histogram(ax(iPeriod),theta_dist,'BinEdges',colon_right(min(theta_dist),5,90),'FaceColor',[0.7695 0.4688 0.0273]);
		vline(m(iPeriod),ax(iPeriod),'linespec','k--','linewidth',0.7);

	end

	% figure setting
	xlabel(ax(1),'\theta_{comm - choice}'); xlabel(ax(2),'\theta_{comm - outcome}');
	ax(1).XLim(1) = min([ax(1).XLim(1) ax(2).XLim(1)]); ax(2).XLim(1) = ax(1).XLim(1);


end