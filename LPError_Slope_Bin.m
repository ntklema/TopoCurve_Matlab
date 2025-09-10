function PE = LPError_Slope_Bin(C,S,BS,nb)
in=find(BS.Z==1);

PE=(C.CMAP.LP.*2-C.CMAP.KM)./(C.CMAP.KM).*100.*BS.Z;
b=BIN(C.CMAP.Sl(in),PE(in),nb);
x=b(:,1);
y=b(:,3);
se=b(:,7)./2;

US=CurveObj(S);
US=CurveCalc(US,0);

PES=(US.CMAP.LP.*2-US.CMAP.KM)./(US.CMAP.KM).*100;

in=find(~isnan(PES));

b=BIN(US.CMAP.Sl(in),PES(in),2000);
xs=b(:,1);
ys=b(:,2);

figure; hold on
fig = figure('units','inch','position',[5,5,8.5,7.8]); hold on
errorbar(x,y,se,'sq','Color','k','MarkerFaceColor',[0.4 0 0.8],'MarkerSize',20,'LineWidth',1)
plot(xs,ys,'r','LineStyle','--','LineWidth',3)

xlim([0 1.5])
ylim([0,400])

% 
% ae=-(1-(1+sl.^2).^(3/2)).*100;
% y=mean(DM,1,'omitmissing');
% sd=std(DM,0,1,'omitmissing');
% ci=1.96.*sd./sqrt(numel(sd));
% 
% figure
% hold on
% % patch([sl,fliplr(sl)],[y-ci, fliplr(y+ci)],[0,0,0],'FaceAlpha',0.1)
% % scatter(sl,y);
% errorbar(sl,y,ci,'sq','Color','k','MarkerFaceColor',[0.4 0 0.8],'MarkerSize',20,'LineWidth',1)
% 
% plot(sl,ae,'r','LineStyle','--','LineWidth',3)


end