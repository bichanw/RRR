function [ax,r,c] = np(n,nc)

		
	close all;

	if nargin == 0;
		ax = subplot(1,1,1,'NextPlot','add');
		r = 1; c = 1;
	elseif nargin == 1
		[r,c] = cal_subplot_numbers(n);
		ax = arrayfun(@(i) subplot(r,c,i,'NextPlot','add'), 1:n);
	elseif nargin == 2
		r = n; c = nc;
		ax = arrayfun(@(i) subplot(n,nc,i,'NextPlot','add'), 1:(n*nc));
	end

	set(gcf,'Position',[0 0 c*150 r*100],'Color','w');

end