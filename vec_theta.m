function theta = vec_theta(u,v)

if prod(size(u)==size(v)) == 0
	error('size mismatch');
end
if size(u,1) ~= 1
	u = u';
	v = v';
end
theta = acosd(u*v' / (norm(u)*norm(v)));

if theta>90
	theta = 180-theta;
end
