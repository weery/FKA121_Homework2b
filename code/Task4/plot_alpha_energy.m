clear all, close all, clc

alpha_energy=importdata('alpha.dat');

size_alpha=size(alpha_energy);



tot_min =0;
tot_min_alpha=0;
idx=0;
tot_min_idx=0;
% find the run and index of the alpha that gives the lowest energy
for i=1:(size_alpha(2)/2)
    j=2*i-1;
    [min_val,min_idx]=min(alpha_energy(:,j+1));
    if (tot_min > min_val)
        tot_min = min_val;
        tot_min_alpha=alpha_energy(min_idx,j);
        idx=j;
        tot_min_idx=min_idx;
    end
end


plot(alpha_energy(:,idx))
hold on
plot(tot_min_idx,alpha_energy(tot_min_idx,idx),'*')

xlabel('Run $p$ [$\#$]','interpreter','latex')
ylabel('$\alpha$ [$\#$]','interpreter','latex')
title('Estimation of $\alpha$ which yeild lowest energy','interpreter','latex')

figure
plot(alpha_energy(:,idx+1))
hold on
plot(tot_min_idx,alpha_energy(tot_min_idx,idx+1),'*')
