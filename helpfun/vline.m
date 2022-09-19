function h=vline(x,ax,varargin)

if nargin < 2 || isempty(ax)
    ax = gca; ax.NextPlot = 'add';
end

% parser
p = inputParser;
addParameter(p,'y',[]);
addParameter(p,'linespec','k-');
addParameter(p,'linewidth',2);
parse(p,varargin{:});

% prepare input
if isempty(p.Results.y)
	y = repmat(ax.YLim',[1 numel(x)]);
else
	y = repmat(p.Results.y,[1 numel(x)]);
end
if size(x,1)>1 x=x'; end
x = repmat(x, [2 1]);
h = plot(ax,x,y,p.Results.linespec,'LineWidth',p.Results.linewidth);

