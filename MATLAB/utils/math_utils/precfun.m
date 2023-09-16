function y = precfun(x)

global P
global PT
global D
global L
global LT

% A = scale * S;
% A = A(p, :);
% A = L \ A;
% A = D_orig \ A;
% A = L' \ A;
% A = scale * A(:, p);
% A = P * A;

has_scale = 0;
% T = transpose(P) * transpose(L_inv) * D * L_inv * P;
% y = T * x; (but never form T explicitly)
if has_scale
    % x1 = scale*x;
    % x2 = (LT \ (D* (L \ x1(p, :))));
    % y = scale * (x2(p, :));
else
    % x1 = (LT \ (D * (L \ x(p, :))));
    % y = x1(p, :);
    y = P * (LT \ (D * (L \ (PT * x))));
end
end

