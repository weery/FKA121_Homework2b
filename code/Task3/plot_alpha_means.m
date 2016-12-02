clear all, clc, clf, close all

alpha_mean_var=importdata('alpha_mean_var.dat');

plot(alpha_mean_var(:,1),alpha_mean_var(:,2))
hold on
errorbar(alpha_mean_var(:,1),alpha_mean_var(:,2),alpha_mean_var(:,3)*sqrt(11/100000))

ylabel('Energy [$a.u.$]','interpreter','latex')
xlabel('$\alpha$ [$\#$]','interpreter','latex')
title('Estimation of lowest energy based on $\alpha$','interpreter','latex')
