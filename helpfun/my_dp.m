function dp = my_dp(Y1,Y2)
% check if the # of neurons is the same
if size(Y1,2) ~= size(Y2,2)
	error('Feature number mismatch');
end

% calculate d-prime
M = [mean(Y1,1); mean(Y2,1)];
V = std([Y1-M(1,:); Y2-M(2,:)],[],1);
dp = diff(M,1,1) ./ V;