clear
clc
home
global KTYPE KSCALE fake_zero size_training x y epsilon
load banana_1;
size_training=400;
data=train.data;
labels=train.labels;
x=data(1:size_training,:);
y=labels(1:size_training);
KTYPE = 6;
KSCALE =0.5;
fake_zero=10^-9;
epsilon=10^-3;
tic;
objcs=CSSVM(2,2,0.8);
toc;
[Kgrid, X1, X2]= util_knlgrid(x);
plot_init();
alpha=objcs.alpha;
b=objcs.b;
fhat = plot_cssvm_fgrid(x, y, Kgrid, 1, alpha, b);
plot_cssvm_step(x, y, X1, X2, fhat);