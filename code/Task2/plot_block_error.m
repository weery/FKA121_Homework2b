clear all, close all, clc

block_error=importdata('block_length.dat');

plot(block_error)

ylabel('Statistical inefficiency','interpreter','latex')
xlabel('Block Length [$\#$]','interpreter','latex')
title('Estimate statistical inefficiency $s$ using blocks','interpreter','latex')

