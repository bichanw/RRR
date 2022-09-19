colormap(flip(cbrewer2('RdBu')));
arrayfun(@(h) set(h,'CLim',[-1 1]*max(abs(h.CLim)),'YDir','reverse'), ax);
c = arrayfun(@(h) colorbar(h), ax);
	