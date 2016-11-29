clear all, close all, clc

alpha_energy=importdata('alpha_avg_energy.dat');

plot(alpha_energy(:,1),mean(alpha_energy(:,2:end),2))

ylabel('Energy [$\#$]','interpreter','latex')
xlabel('$\alpha$ [$\#$]','interpreter','latex')
title('Estimation of lowest energy, vilket nu ser väldigt dåligt ut','interpreter','latex')

