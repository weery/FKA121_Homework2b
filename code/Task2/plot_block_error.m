clear all, close all, clc

block_error=importdata('block_length.dat');

plot(block_error)
hold on
plot([0, length(block_error)],[11, 11],'--')

ylabel('Statistical inefficiency, $s$','interpreter','latex')
xlabel('Block Length, $B$ [$\#$]','interpreter','latex')
title('Estimate statistical inefficiency $s$ using blocks','interpreter','latex')

