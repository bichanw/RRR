% set up
clear;clc;addpath(genpath('communication-subspace'));addpath(genpath('helpfun'));

Folders = {'/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_7/07302021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_7/08022021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_7/08042021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_7/08052021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_9/07272021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_9/07302021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_45/11022021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_45/11042021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_47/11042021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_47/11092021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_47/11112021_T10',...
		   '/jukebox/Bezos/Ryan/Bay2/ACC/CaMKIIa_tetO_51/11022021_T10'};

cvloss = [];theta_plane_choice = [];
for folder = Folders

	% initiation
	clear Epoch; load([folder{1} '/Results.mat']);
	session_name = get_ses_name(folder{1});

	% concatenate activity
	% matrix size - trials * time points * cells
	Xdeep = catcell(Epoch.Spike_Epoch_sm(Deep==1));
	Xsup  = catcell(Epoch.Spike_Epoch_sm(Deep==0));

	% --- Plotting Functions ---
	% parameter variable, change if default is not desired
	% this is the ops used for 7 bins of 11 bin width
		% ops = struct('dim',0:10,'if_plot',false,'twin',[0 10]+[1:11:77]'); ops.win_name = arrayfun(@(i) num2str(i), 1:7,'UniformOutput',false);% average bin
	% ops used for ridge regression
		% ops = struct('method','ridge','if_plot',true,'lamdba',[0:100:1000 1e4]); % ridge

	% contains analysis you do session-by-session
	try
		% ind_left = Epoch.left; ind_right = Epoch.right;
		% return
		[B,choice,dp] = plt_session.choice_selectivity(Xsup,Xdeep,Epoch.left,Epoch.right);
		% plt_session.PC_RRR(Xsup,Xdeep);
		% plt_session.averaged_RRR(Xdeep,Xsup);
		% [ops,cvloss(end+1,:,:,:)] = plt_session.averaged_RRR(Xsup,Xdeep);
		% theta_plane_choice(end+1,:) = plt_session.axis_angle(Epoch,Deep);
	catch ME
		save('error.mat','ME','folder');
		continue;
	end

	% save(sprintf('%s.mat',session_name),'B','choice','dp');
	sgtitle(session_name,'Interpreter','none');
	export_fig(sprintf('choice_vs_B_%s.pdf',session_name));
	fprintf('%s\n',session_name);

end
return

%% ------ extra functions ------ 
function session_name = get_ses_name(folder)
	
	ind = find(folder=='/',1,'last');
	session_name = folder(ind-2:end);
	session_name(session_name=='/') = '_';
end

