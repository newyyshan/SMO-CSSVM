%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kgrid, X1, X2] = util_knlgrid(x)
%function [Kgrid X1 X2] = util_knlgrid(x, hkernel, kernelparam, gnum)
%	generate kernel grid matrix
%	only 2 dim inputs
%
% if nargin < 2, hkernel = @knl_normal; end
% if nargin < 3, kernelparam = 1; end
 gnum = 100;

x1 = linspace(min(x(:,1)), max(x(:,1)), gnum)';
x2 = linspace(min(x(:,2)), max(x(:,2)), gnum)';
[X1, X2] = ndgrid(x1,x2);
XG = [X1(:), X2(:)];
Kgrid =CSSVM.Kernel(x(:,1:2), XG);
