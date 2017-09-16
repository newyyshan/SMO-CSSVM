function [] = plotcssvm( x,y,~ )
%PLOTCSSVM Summary of this function goes here
%   Detailed explanation goes here
plot(x(y==-1,1), x(y==-1,2), 'ob', 'MarkerSize', 7, 'Linewidth', 2);%¸ºÀý
hold on
plot(x(y== 1,1), x(y== 1,2), '+r', 'MarkerSize', 8, 'Linewidth', 2);%ÕýÀý
%     [C h] = contour(X1, X2, fhat, [0 0], 'k', 'LineWidth', 2);
hold off
xr = minmax(x(:,1:2)')';
axis(xr(:) + 0.1*[-1 1 -1 1]')
set(gca, 'xtick', []), set(gca, 'ytick', []);
% if (~bfile)
	drawnow
end

