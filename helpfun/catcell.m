function x = catcell(x,idim)
% empty cell
if isempty(x) 
    x = [];
    return
end

% finding the dimension
if nargin < 2

    % matrix in cell is full rank
    if isempty(find(size(x{1})==1))
        x = cat(ndims(x{1})+1,x{:});
    else
    % matrix has a dim is 1
        x = cat(find(size(x{1})==1,1),x{:});
    end

% given the dimension
else
    x = cat(idim,x{:});

end

end
