function W = normw(W)

n = ndims(W);
nterm = sum(W, n);
nterm = repmat(nterm,[ones(1,n-1) size(W,n)]);
nterm = nterm + (nterm==0); % protect against zeros before dividing
W = W ./ nterm; % use the same W as output here to save RAM

end
