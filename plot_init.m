%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_init()
%function plot_init()
%	plot initialization
%
clf;
figsz = get(gcf, 'Position');
set(gcf, 'Position', [figsz(1) 250 540 540]);
set(gca, 'Position', [0 0 1 1]);
set(gca, 'xtick', []), set(gca, 'ytick', []);
box on;
set(gcf, 'PaperPositionMode', 'auto')
