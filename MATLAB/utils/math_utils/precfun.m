function y = precfun(x)

global P
global PT
global D
global L
global LT

% T = transpose(P) * transpose(L_inv) * D * L_inv * P;
% y = T * x; (but never form T explicitly)
y = P * (LT \ (D * (L \ (PT * x))));
end

