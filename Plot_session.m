% set up
clear;clc;addpath(genpath('communication-subspace'));

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
	Xdeep = catcell(Epoch.Spike_Epoch_sm(Deep==1));
	Xsup  = catcell(Epoch.Spike_Epoch_sm(Deep==0));

	% --- Plotting Functions ---
	% ops = struct('dim',0:10,'if_plot',false,'twin',[0 10]+[1:11:77]'); ops.win_name = arrayfun(@(i) num2str(i), 1:7,'UniformOutput',false);% average bin
	% ops = struct('method','ridge','if_plot',true,'lamdba',[0:100:1000 1e4]); % ridge
	try
		% plt_session.choice_selectivity(Xsup,Xdeep,Epoch.left,Epoch.right);
		% plt_session.PC_RRR(Xsup,Xdeep);
		% [ops,cvloss(end+1,:,:,:)] = plt_session.averaged_RRR(Xsup,Xdeep);
		% plt_session.averaged_RRR(Xdeep,Xsup);
		theta_plane_choice(end+1,:) = plt_session.axis_angle(Epoch,Deep);
	catch ME
		save('error.mat','ME','folder');
		continue;
	end


	title(session_name,'Interpreter','none');
	export_fig(sprintf('angle_%s.pdf',session_name));
	fprintf('%s\n',session_name);

end
return

%% ------ extra functions ------ 
function session_name = get_ses_name(folder)
	
	ind = find(folder=='/',1,'last');
	session_name = folder(ind-2:end);
	session_name(session_name=='/') = '_';
end

