clear all, close all, clc

alpha_energy=importdata('alpha_avg_energy.dat');

semilogx(alpha_energy(:,1),mean(alpha_energy(:,2:end),2))

ylabel('Energy [$a.u.$]','interpreter','latex')
xlabel('$\alpha$ [$\#$]','interpreter','latex')
title('Estimation of lowest energy based on $\alpha$','interpreter','latex')

