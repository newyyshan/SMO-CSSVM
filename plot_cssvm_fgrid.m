%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fhat = plot_cssvm_fgrid(~, y, Kgrid, ~, alpha, b)
%function fhat = plot_cssvm_fgrid(x, y, Kgrid, lambda, alpha, alpha0, blinear, gnum)
%	calculate a decision boundary f(x)
%	only 2 dim inputs
%
 gnum = 100;

% if (blinear)
% 	beta = (alpha.*y)'*x;
% 	fhat = -beta(1)/beta(2)*minmax(x(:,1)') - alpha0/beta(2);
% 	fhat(2,:) = fhat(1,:) + lambda/beta(2);
% 	fhat(3,:) = fhat(1,:) - lambda/beta(2);
% else
% 	fhat = (alpha.*y)'*Kgrid + b;
% 	fhat = reshape(fhat, gnum, gnum);
% end
	fhat = (alpha.*y)'*Kgrid + b;
	fhat = reshape(fhat, gnum, gnum);
