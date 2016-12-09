clear all, close all, clc

energies=importdata('energy.dat');

x=(1:length(energies))';
mean = cumsum(energies)./x;

plot(mean(1:10000))

ylabel('Energy [$a.u$]','interpreter','latex','fontsize',19)
xlabel('Number of trials [$\#$]','interpreter','latex','fontsize',19)
title('Estimated ground state energy $E_0$','interpreter','latex','fontsize',19)

