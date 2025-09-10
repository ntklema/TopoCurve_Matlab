function SpecFiltPlot(C)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
kfilt=1./C.Filter;
km=C.DEM.km;
sigma=abs(kfilt(2)-kfilt(1))/3;
F=exp(-(km-kfilt(1)).^2/(2*sigma^2));
F=F/max(F(:));
F(km<kfilt(1))=1;


figure
scatter(km(:),C.DEM.ZP(:),1,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',0.1); hold on
scatter(km(:),C.DEM.ZP(:).*F(:),1,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',0.1); 
set(gca,'xscale','log','yscale','log')
 ylim([1e-5 1e10])
yyaxis("right")
plot(km(:),F(:),'LineWidth',2)
set(gca,'yscale','log')
ylim([1e-5 1])
end