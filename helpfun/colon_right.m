function a = colon_right(j,i,k)
	a = j:i:k;
	if a(end)<k
		a(end+1) = a(end) + i;
	end
end