function PED=Slope_Error(C,BS,R,nb)

cmap=(crameri('vik'));

S=Sphere_Curvature(1001,R);
US=CurveObj(S);
US=CurveCalc(US,0);
D=GRIDobj(S.X,S.Y,S.Z);
G=gradient8(D);
[Fx,Fy]=gradient(S.Z,S.dx,S.dx);
F=sqrt(Fx.^2+Fy.^2);
[Y,X]=meshgrid(S.X,S.Y);
sl=sqrt(X.^2+Y.^2)./S.Z;

PE=(G.Z-sl)./abs(sl).*100;
PET=abs(US.CMAP.Sl-sl)./abs(sl).*100;


bs=BIN(US.CMAP.az(~isnan(US.CMAP.az)),PE(~isnan(US.CMAP.az)),1000);
bsc=BIN(US.CMAP.az(~isnan(US.CMAP.az)),PET(~isnan(US.CMAP.az)),1000);

Sl=C.CMAP.Sl;
DEM=C.DEMGRIDobj;
DEM.Z=C.DEM.ZFilt;
G=gradient8(DEM);


PED=(G.Z-Sl)./abs(Sl).*100;
PED=PED.*BS.Z;
in=find(PED~=0);
bd=BIN(C.CMAP.az(in),PED(in),nb);
sed=bd(:,7)./2;

% bin by area
A=(C.Stream.A.Z);
A2=(C.Stream.A2.Z);

LP=C.CMAP.LP.*2;
PEA=(A2-A)./(A).*100;
PEL=abs(LP-C.CMAP.KM)./abs(C.CMAP.KM).*100;
in=find(BS.Z==1);

baa=BIN(log10(A(in)),PEA(in),nb);
aa=10.^baa(:,1);
sea=baa(:,7)./2;

bal=BIN(log10(A(in)),PEL(in),nb);
al=10.^bal(:,1);
sel=bal(:,7)./2;

bas=BIN(log10(A(in)),PED(in),nb);
as=10.^bas(:,1);
ses=bas(:,7)./2;
% figure
% imAlpha = ones(size(PE));
% imAlpha(isnan(PE)) = 0;

% imagesc(PE, 'AlphaData', imAlpha);
% colorbar
% colormap(cmap)
% clim([-7 7])

% figure
% plot(bs(:,1),bs(:,3),'LineWidth',2,'Color','k')
% hold on
% plot(bsc(:,1),bsc(:,3),'b','LineWidth',2,'LineStyle','-.')
% errorbar(bd(:,1),bd(:,3),sed,'^','MarkerFaceColor','c', ...
%     'Color','k','MarkerSize',10,'LineWidth',0.2);
% xlim([0 pi/2])
% ylim([-7.5 1])
% 
figure
fig = figure('units','inch','position',[5,5,8,7.8]);
% errorbar(aa,baa(:,8),sea)
patch([aa' fliplr(aa')],[(baa(:,3)-sea)' fliplr((baa(:,3)+sea)')],'k','FaceAlpha',0.4)
hold on


errorbar(as,bas(:,3),ses,'^','MarkerFaceColor','c','Color','k','MarkerSize',12,'LineWidth',0.5);
% patch([as' fliplr(as')],[(bas(:,8)-ses)' fliplr((bas(:,8)+ses)')],'c','FaceAlpha',0.4)

errorbar(al,bal(:,3),sel,'sq','Color','k','MarkerFaceColor',[0.4 0 0.8],'MarkerSize',12,'LineWidth',0.5)
% patch([al' fliplr(al')],[(bal(:,8)-sel)' fliplr((bal(:,8)+sel)')],[0.4 0 0.8],'FaceAlpha',0.4)

set(gca,'xscale','log')
xlim([10.^1.5 10.^7.3])
ylim([-50 70])