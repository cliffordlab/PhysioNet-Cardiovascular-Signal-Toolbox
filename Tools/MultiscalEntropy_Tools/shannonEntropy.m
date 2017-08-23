function se = shannonEntropy(data)
	denom = sum(data);
	temp = 0;
	for i = 1:length(data)
		temp = temp + (data(i) / denom) * log10(data(i));
	end
	se = -temp;
end