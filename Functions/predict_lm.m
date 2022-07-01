function y = predict_lm(x,p)
% x: Input vector
% p: Parameter vector

p = p(:);

% Create design matrix
X = lm_mat(x);

% Compute response
y = X * p;

% 
% x = x(:)';
% 
% n = length(x);
% k = 1;
% for i = 1:n
%     k = k+1;
%     P{k} = i;
% end
% 
% for i = 1:(n-1)
%     for j = (i+1):n
%         k = k + 1;
%         P{k} = [i j];
%     end
% end
% 
% for i = 1:n
%     k = k+1;
%     P{k} = [i i];
% end
% 
% y = p(1);
% for i = 2:length(P)
%     y = y + p(i) * prod(x(P{i}));
% end

end