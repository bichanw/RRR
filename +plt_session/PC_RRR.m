function PC_RRR(Xsup,Xdeep)

% PCA on original activity
NPC = 10;
Xsup  = permute(catcell(arrayfun(@(iepo) RR(squeeze(Xsup(:,iepo,:)),20), 1:77,'UniformOutput',false)), [1 3 2]);
Xdeep = permute(catcell(arrayfun(@(iepo) RR(squeeze(Xdeep(:,iepo,:)),10), 1:77,'UniformOutput',false)), [1 3 2]);

% run RRR
ops = struct('method','RRR','dim',0:10);
cvloss = catcell(arrayfun(@(iepo) RRR(squeeze(Xsup(:,iepo,:)),squeeze(Xdeep(:,iepo,:)),ops), 1:77,'UniformOutput',false));

% plot cvloss lines
colors = cbrewer2('spectral',size(cvloss,2));
ax = np;
% cvloss lines
arrayfun(@(i) errorbar(1:77, squeeze(1-cvloss(1,i,:)), squeeze(cvloss(2,i,:)), 'o--', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6),1:size(cvloss,2));
% vertical dashed lines for events
vline([5.5 25.5 40.5 47.5],ax,'linewidth',0.7,'linespec','k--'); % mark events
hline(0,ax,'linewidth',0.7,'linespec','k--');
colormap(colors); c = colorbar; c.Label.String = '# rank'; set(c,'Ticks',[0 1],'TickLabels',{'0','full'});
set(gcf,'Position',[0 0 350 200]);

