clear all, close all, clc

alpha_energy=importdata('alpha.dat');

size_alpha=size(alpha_energy);

% make new loop that starts at k=1:4
% and start in alpha_energy ((k-1)*200)

tot_min =zeros(1,4);
tot_min_alpha=zeros(1,4);
idx=zeros(1,4);
tot_min_idx=zeros(1,4);
% find the run and index of the alpha that gives the lowest energy
for k=1:4
    for i=1:(size_alpha(2)/2)
        j=2*i-1;
        startIdx=(k-1)*200+1;
        endIdx = k*200;
        [min_val,min_idx]=min(alpha_energy(startIdx:endIdx,j+1));
        if (tot_min(k) > min_val)
            tot_min(k) = min_val;
            tot_min_alpha(k)=alpha_energy(min_idx+startIdx,j);
            idx(k)=j;
            tot_min_idx(k)=min_idx;
        end
    end
end

hold on

for k = 1:4
    startIdx=(k-1)*200+1;
    endIdx = k*200;
    startShift=0;
    plot(startShift+1:200,alpha_energy(startIdx+startShift:endIdx,idx(k)))
end

xlabel('Run $p$ [$\#$]','interpreter','latex','fontsize',19)
ylabel('$\alpha$ [$\#$]','interpreter','latex','fontsize',19)
title('Estimation of $\alpha$ which yeild lowest energy','interpreter','latex','fontsize',19)
xlim([startShift, 200])

figure
plot(alpha_energy(10:end,idx+1))
hold on
plot(tot_min_idx,alpha_energy(tot_min_idx,idx+1),'*')
