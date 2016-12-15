close all; clear all; clc

data = importdata('auto_correlation.dat');

plot(data(:,1), data(:,2))
hold on
plot(data(12,1),data(12,2),'r*')

xlabel('$k$','interpreter','latex','fontsize',19)
ylabel('Auto-correlation function $\Phi_k$','interpreter','latex','fontsize',19)
title('Auto-correlation decay','interpreter','latex','fontsize',19)
L = legend('$\Phi_k$','Correlation decay threshold');
set(L, 'interpreter', 'latex', 'fontsize', 14)
