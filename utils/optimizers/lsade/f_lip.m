function Y = f_lip(Yb,Xb,k,X)
% Lipschitz-based surrogate value at point x
[~,N] = size(X);
Y = zeros(N,1);
for n = 1:N
    vals = zeros(size(Yb));
    for i=1:length(Yb)
        vals(i) = min(Yb(i) + k*norm(Xb(:,i)-X(:,n)));
    end
    Y(n) = min(vals);
end