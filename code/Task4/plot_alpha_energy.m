clear all, close all, clc

alpha_energy=importdata('alpha.dat');

size_alpha=size(alpha_energy);

sim_length = 500;

% make new loop that starts at k=1:4
% and start in alpha_energy ((k-1)*200)

tot_min =zeros(1,4);
tot_min_alpha=zeros(1,4);
idx=zeros(1,4);
tot_min_idx=zeros(1,4);
% find the run and index of the alpha that gives the lowest energy
for k=1:4
    lastIdx = sim_length*k;
    for i=1:(size_alpha(2)/2)
        j=2*i-1;
        min_val = alpha_energy(lastIdx,j+1);
        if (tot_min(k) > min_val)
            tot_min(k) = min_val;
            tot_min_alpha(k)=alpha_energy(lastIdx,j);
            idx(k)=j;
        end
    end
end


figure
hold on


for k = 1:4
    startIdx=(k-1)*sim_length+1;
    endIdx = k*sim_length;
    startShift=0;
    plot((startShift+1):sim_length,alpha_energy((startIdx+startShift):endIdx,idx(k)))
end

legend({'$\beta = 0.5$','$\beta = 0.625$','$\beta = 0.75$','$\beta = 0.875$'},'interpreter','latex')

xlabel('Run $p$ [$\#$]','interpreter','latex','fontsize',19)
ylabel('$\alpha$ [$\#$]','interpreter','latex','fontsize',19)
title('Estimation of $\alpha$ which yeild lowest energy','interpreter','latex','fontsize',19)
xlim([startShift, sim_length])
