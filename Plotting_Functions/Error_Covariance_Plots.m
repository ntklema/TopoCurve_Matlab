function Error_Covariance_Plots(DM,DS,km,sl)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% DS(isnan(DM))=nan;

cmap=(crameri('lajolla'));
Alph = ones(size(DM));
Alph(isnan(DM))=0;

DSD=DS;
DSD(isnan(DM))=nan;

figure
subplot(2,2,2)
imagesc(sl,km,DM,'AlphaData',Alph)
colorbar
set(gca,'ydir','normal','colorscale','linear')
colormap(cmap)
clim([0 300])

subplot(2,2,1)
imagesc(sl,km,DS)
colorbar
set(gca,'ydir','normal','colorscale','linear')
% clim([0 11e-3])

sm=50;
subplot(2,2,3); hold on
patch([km,fliplr(km)],[mean(DM,2,"omitmissing")'-mad(DM,0,2)', fliplr(mean(DM,2,"omitmissing")'+mad(DM,0,2)')],[0,0,0],'FaceAlpha',0.1)
plot(km,mean(DM,2,"omitmissing"),'k','LineWidth',1)
% plot(km,mean(DM,2,"omitmissing")+mad(DM,0,2),'k')
% plot(km,mean(DM,2,"omitmissing")-mad(DM,0,2),'k')
% plot(km,min(DM,[],2))
% plot(km,max(DM,[],2))
% plot(km,mean(DM,2,"omitmissing")+2.*std(DM,0,2,"omitmissing"))
% plot(km,mean(DM,2,"omitmissing")-2.*std(DM,0,2,"omitmissing"))
plot(km,mean(DS,2,"omitmissing"),'color',[0,0,0],'LineWidth',1,'LineStyle','--')
plot(km,mean(DSD,2,"omitmissing"),'color',[0.8,0,0],'LineWidth',1,'LineStyle','-')
set(gca,'yscale','linear')

subplot(2,2,4); hold on
patch([sl,fliplr(sl)],[mean(DM,1,"omitmissing")-mad(DM,0,1), fliplr(mean(DM,1,"omitmissing")+mad(DM,0,1))],[0,0,0],'FaceAlpha',0.1)
plot(sl,mean(DM,1,"omitmissing"),'k','LineWidth',1)
plot(sl,-(1-(1+sl.^2).^(3/2)).*100,'m')
% scatter(sl,mean(DM,1,"omitmissing"),'k','filled','MarkerEdgeColor','w','LineWidth',1)
% plot(sl,mean(DM,1,"omitmissing")+mad(DM,0,1),'k')
% plot(sl,mean(DM,1,"omitmissing")-mad(DM,0,1),'k')
% plot(sl,min(DM,[],1))
% plot(sl,max(DM,[],1))
plot(sl,mean(DM,1,"omitmissing")+1.96.*std(DM,0,1,"omitmissing")./sqrt(numel(sl)))
plot(sl,mean(DM,1,"omitmissing")-1.96.*std(DM,0,1,"omitmissing")./sqrt(numel(sl)))

plot(sl,mean(DS,1,"omitmissing"),'color',[0,0,0],'LineWidth',1,'LineStyle','--')
plot(sl,mean(DSD,1,"omitmissing"),'color',[0.8,0,0],'LineWidth',1,'LineStyle','-')
set(gca,'yscale','linear')





end