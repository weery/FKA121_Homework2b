clf
radius = importdata('rads.dat');

r=0:0.01:max(radius);

p=@(r,Z)Z.^3*4*r.^2.*exp(-2*Z*r);

[h,b]=histcounts(radius,100,'Normalization','pdf');

hold on
ax = gca;
ax.ColorOrderIndex = 2;
plot(r,p(r,2),'-','LineWidth',2)
plot(r,p(r,27/16),'-','LineWidth',2)
ax.ColorOrderIndex = 1;
plot((b(1:end-1)+b(2:end))/2,h,'*')

xlim([0,max(radius)])
legend('Simulated distribution','Unscreened nucleus','Variationally optimized')

xlabel('$Radius / [{\AA}]$','interpreter','latex', 'fontsize', 15)
ylabel('$RDF / [{\AA}^{-1}]$','interpreter','latex', 'fontsize', 15)
title('Radial distribution function of helium electrons',...
    'interpreter','latex', 'fontsize', 15)