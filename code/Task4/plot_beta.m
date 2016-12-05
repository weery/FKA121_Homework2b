close all, clear all, clc

delimiterIn = '\t';
headerlinesIn = 1;
data = importdata('beta.dat',delimiterIn,headerlinesIn);

plot(data.data(:,1), data.data(:,3))

xlabel('Beta, $\beta$','interpreter','latex')
ylabel('Energy E [$a.u.$]','interpreter','latex')
title('Estimating energy, $E$, using different $\beta$ as damping parameters.','interpreter','latex')
