clear all, close all, clc

energies=importdata('energy.dat');

x=(1:length(energies))';
mean = cumsum(energies)./x;

plot(mean)

ylabel('Energy [$a.u$]','interpreter','latex')
xlabel('Number of trials [$\#$]','interpreter','latex')
title('Estimated ground state energy $E_0$','interpreter','latex')

