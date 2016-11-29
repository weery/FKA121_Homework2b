clear all, close all, clc

alpha_energy=importdata('alpha.dat');

size_alpha=size(alpha_energy);


mean_min=0;
mean_min_alpha=0;
tot_min =0;
tot_min_alpha=0;

for i=1:(size_alpha(2)/2)
    j=2*i-1;
    [min_val,min_idx]=min(alpha_energy(:,j+1));
    plot(alpha_energy(:,j),alpha_energy(:,j+1))
    mean_min = mean_min+ min_val;
    mean_min_alpha = mean_min_alpha +alpha_energy(min_idx,j);
    if (tot_min > min_val)
        tot_min = min_val;
        tot_min_alpha=alpha_energy(min_idx,j);
    end
end

mean_min=mean_min/(size_alpha(2)/2);
mean_min_alpha=mean_min_alpha/(size_alpha(2)/2);

ylabel('Energy [$\#$]','interpreter','latex')
xlabel('$\alpha$ [$\#$]','interpreter','latex')
title('Estimation of lowest energy','interpreter','latex')

