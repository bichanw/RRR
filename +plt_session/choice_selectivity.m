function choice_selectivity(Xsup,Xdeep,ind_left,ind_right)

ops = struct;
ops.twin = getOr(ops,'twin',[6 25;26 40]);
ops.win_name = getOr(ops,'win_name',{'cue','delay'});

nwins = numel(ops.win_name);
ax = np(2,2); 
for iwin = 1:nwins
	X = squeeze(mean(Xsup(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));
	Y = squeeze(mean(Xdeep(:,ops.twin(iwin,1):ops.twin(iwin,2),:),2));


	sig = {X,Y};
	for isig = 1:numel(sig)
		% examine the choice selectivity of each session
		L = mean(sig{isig}(ind_left,:),1);
		R = mean(sig{isig}(ind_right,:),1);

		choice = (L-R) ./ (L+R);
		dp = my_dp(sig{isig}(ind_left,:),sig{isig}(ind_right,:)); 

		histogram(ax(sub2ind([2 2],iwin,isig)),choice);
		histogram(ax(sub2ind([2 2],iwin,isig)),dp);
	end
end

% figure setting
title(ax(1),'choice selectivity'); title(ax(2),"d'");
ylabel(ax(1),'X'); ylabel(ax(3),'Y');
legend(ax(end),ops.win_name);

% plots to examine value
	% % examine scatter plots
	% 	ax = np; 
	% 	scatter(dp,choice); 
	% 	xlabel('dp'); ylabel('choice');
	% 	export_fig tmp.pdf

	
	% % examine dprime, choice selectivity and its distribtuion
	% 	ax = np(size(X,2));
	% 	for icell = 1:size(X,2)
	% 		histogram(ax(icell),X(ind_left,icell),'FaceColor',colors(1,:),'FaceAlpha',0.3);
	% 		histogram(ax(icell),X(ind_right,icell),'FaceColor',colors(2,:),'FaceAlpha',0.3);
	% 		title(ax(icell),sprintf('choice=%.2f dp=%.2f',choice(icell),dp(icell)));
	% 	end
	% 	export_fig tmp.pdf
end
