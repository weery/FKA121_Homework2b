clear all, clc, close all

angle_diff = importdata('angle_diff.dat');


[h,bin]=histcounts(angle_diff,'Normalization','pdf');

plot(bin(1:end-1),h)

xlabel('Och jag med','interpreter','latex')
ylabel('Jag \"ar fel','interpreter','latex')
title('En plott som \"ar fel och ska r\"attas till','interpreter','latex')