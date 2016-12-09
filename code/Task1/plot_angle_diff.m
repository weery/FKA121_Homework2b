clear all; clc; close all

angle_diff = importdata('angle_diff.dat');

x = cos(angle_diff);

[h,bin]=histcounts(x,'Normalization','pdf');

plot(bin(1:end-1),h)
axis([-1 1 0 1])

xlabel('$x = \cos \theta$','interpreter','latex', 'fontsize', 19)
ylabel('$P(x)$','interpreter','latex', 'fontsize', 19)
title('Distribution of $\Delta \theta$ for the electrons',...
    'interpreter','latex', 'fontsize', 19)