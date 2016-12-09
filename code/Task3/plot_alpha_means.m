clear all, clc, clf, close all

alpha_mean_var=importdata('alpha_mean_var.dat');

%plot(alpha_mean_var(:,1),alpha_mean_var(:,2))
hold on
errorbar(alpha_mean_var(:,1),alpha_mean_var(:,2),alpha_mean_var(:,3)*sqrt(11/1000000))

xlim([alpha_mean_var(1,1),alpha_mean_var(end,1)])

ylabel('Energy [$a.u.$]','interpreter','latex','fontsize',19)
xlabel('$\alpha$ [$\#$]','interpreter','latex','fontsize',19)
title('Estimation of lowest energy based on one $\alpha$ run','interpreter','latex','fontsize',19)




[min_alpha, min_alpha_idx]=min(alpha_mean_var(:,2));
plot(alpha_mean_var(min_alpha_idx,1),alpha_mean_var(min_alpha_idx,2),'rx')

min_energy_plot=-2.879;

plot([alpha_mean_var(min_alpha_idx,1) alpha_mean_var(min_alpha_idx,1)],[min_energy_plot,alpha_mean_var(min_alpha_idx,2)],'k--')
plot([alpha_mean_var(1,1) alpha_mean_var(min_alpha_idx,1)],[alpha_mean_var(min_alpha_idx,2),alpha_mean_var(min_alpha_idx,2)],'k--')

alpha_txt = sprintf('%f', alpha_mean_var(min_alpha_idx,1));
alpha_txt = strcat('$\leftarrow \alpha = ',alpha_txt, '$');

text(alpha_mean_var(min_alpha_idx,1)+0.001,alpha_mean_var(min_alpha_idx,2)+(min_energy_plot-alpha_mean_var(min_alpha_idx,2))/2,alpha_txt,'interpreter','latex')

energy_txt = sprintf('%f', alpha_mean_var(min_alpha_idx,2));
energy_txt = strcat('$\uparrow E = ',energy_txt, '$');

text((alpha_mean_var(min_alpha_idx,1)+alpha_mean_var(1,1))/2-0.03,alpha_mean_var(min_alpha_idx,2)-0.00025,energy_txt,'interpreter','latex')

xlim([alpha_mean_var(1,1),alpha_mean_var(end,1)])

