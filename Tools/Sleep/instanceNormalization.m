function dlY = instanceNormalization(dlX,fmt)

reductionDims = find(fmt == 'S');
mu = mean(dlX,reductionDims);
sigmaSq = var(dlX,1,reductionDims);

epsilon = 1e-5;
dlY = (dlX-mu) ./ sqrt(sigmaSq+epsilon);

end