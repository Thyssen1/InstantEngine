function X = lm_mat(x)

[m,n] = size(x);

X = zeros(m,1+2*n+nchoosek(n,2));

X(:,1) = 1;

k = 1;
for i = 1:n
    k = k+1;
    X(:,k) = x(:,i);
end

for i = 1:(n-1)
    for j = (i+1):n
        k = k + 1;
        X(:,k) = prod(x(:,[i j]),2);
    end
end

for i = 1:n
    k = k+1;
    X(:,k) = prod(x(:,[i i]),2);
end

end