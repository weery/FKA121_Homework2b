clear all, close all, clc

alpha_energy=importdata('alpha_avg_energy.dat');

plot(alpha_energy(:,1),alpha_energy(:,2))

ylabel('Statistical inefficiency','interpreter','latex')
xlabel('Block Length [$\#$]','interpreter','latex')
title('Estimate statistical inefficiency $s$ using blocks','interpreter','latex')

