clf
radius = importdata('rads.dat');

r=0:0.01:max(radius);

p=@(r,Z)Z.^3*4*r.^2.*exp(-2*Z*r);

[h,b]=histcounts(radius,100,'Normalization','pdf');
plot((b(1:end-1)+b(2:end))/2,h)
hold on
plot(r,p(r,2))
plot(r,p(r,27/16))

xlim([0,max(radius)])
legend('Simulated distribution','Unscreened nucleus','Variationally optimized')

xlabel('$r / [\AA]$','interpreter','latex')
ylabel('$RDF / [\AA^{-1}]$','interpreter','latex')
title('Radial distribution function of helium','interpreter','latex')