close all, clear all, clc

delimiterIn = '\t';
headerlinesIn = 1;
data = importdata('beta.dat',delimiterIn,headerlinesIn);

plot(data.data(:,1), data.data(:,3))

xlabel('Beta','interpreter','latex')
ylabel('Energy E [$a.u.$]','interpreter','latex')
title('Estimation of $\alpha$ which yeild lowest energy','interpreter','latex')
