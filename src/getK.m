function [ks_d, ks_p] = getK(Vmax, kcat, n, r, A, D, u)
%% Solve for ks at the diffusive limit
ks_d = (n.*kcat)./ (r .*D); % moles m-3

%% Solve for ks at the porter limit
% mass transfer coefficient
mtc = D./r + u/2; % m s-1
% capture probability
alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless

ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3



end