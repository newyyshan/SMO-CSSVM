%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plot_cssvm_step(x, y, X1, X2, fhat)
%function h = plotstep(x, y, fhat, blinear, bfile, step)
%	plot data points and the decision boundary
%	only 2 dim inputs
%


plot(x(y==-1,1), x(y==-1,2), 'ob', 'MarkerSize', 7, 'Linewidth', 2);
hold on
plot(x(y== 1,1), x(y== 1,2), '+r', 'MarkerSize', 8, 'Linewidth', 2);
% if (blinear)
% 	plot(minmax(x(:,1)'), fhat(1,:), 'k', 'LineWidth', 2)
% 	plot(minmax(x(:,1)'), fhat(2,:), ':k', 'LineWidth', 1)
% 	plot(minmax(x(:,1)'), fhat(3,:), ':k', 'LineWidth', 1)
% else
    [~, h] = contour(X1, X2, fhat, [0 0], 'k', 'LineWidth', 2);
%     countour(X1, X2, fhat, [1 -1], ':k', 'LineWidth', 1);
% end
hold off
xr = minmax(x(:,1:2)')';
axis(xr(:) + 0.1*[-1 1 -1 1]')
set(gca, 'xtick', []), set(gca, 'ytick', []);
% if (~bfile)
	drawnow
% else
% %	saveas(gcf, [filename '.fig'], 'fig');
% 	print('-dpng', '-r96', [filename '.png'] );
% end
