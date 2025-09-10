
filePath = matlab.desktop.editor.getActiveFilename;
addpath(genpath(fileparts(filePath)));

% Load DEM
demname="Umpqua_10m_2";
DEM=GRIDobj(demname+".tif");
% DEM=resample(DEM,5);
%DEM=crop(DEM,[4.19e5, 4.35e5],[4.834e6, 4.848e6]);

%DEM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/GIS_Repositories/Coast Range/DEMs/Lidar/Umpqua_North_Lidar_1.tif')

% Load ROI Shapefile
bs=shaperead(demname+"_Outline.shp");

%% Generate CurveObj and compute curvature
filter=[200 150];
C=CurveObj(DEM);
C=SpecFilt(C,filter,'lowpass');
% C.DEM.ZFilt=imgaussfilt(C.DEM.Z,5);
C=CurveCalc(C,0);


% Route flow using Dinf algorithm
C=RouteFlow(C,6000);

LogA=C.Stream.A;
LogA.Z=log10(LogA.Z);
clear filter demname

% Bin Data by drainage area for area defined by bs polygon
[BS,BD] = bindata(C,bs,100);

%% Area binned distribution plots
Area_Dist(BD,'save','APlots_FullArea') % Full Area
% Area_Dist(BD);

%% Cropped Curvature Map
athresh=[0 log10(375)];
cmap=[0.4,0,0;
    1,0,0;
    1,1,1;
    0,1,1;
    0,0,1];
D=DEM;
D.Z=C.DEM.ZFilt;
SM=C.SMAP;
SM(or(LogA.Z<=athresh(1),LogA.Z>athresh(2)))=nan;
SM(BS.Z==0)=nan;
% hillshade(DEM)

figure
imageschs(DEM,SM,'colormap',cmap)

% colorbar
xlim([4.20100 4.26600].*1e5)
ylim([4.8398 4.8456].*1e6)

per=numel(find(~isnan(SM)))/numel(find(BS.Z==1))

%%  Make shapefile of zoomed in region
S=struct;
S.Geometry='Polygon';
S.X=[4.20100 4.20100 4.26600 4.26600].*1e5;
S.Y=[4.8398 4.8456 4.8456 4.8398 ].*1e6;
shapewrite(S,'/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Shapefiles/Map_Focus_Region.shp')

%% Make shapefile of map overview footprint
S=struct;
S.Geometry='Polygon';
S.X=[419000.000 419000.000 435000.000 435000.000];
S.Y=[4834000.000 4849000.000 4849000.000 4834000.000];
shapewrite(S,'/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Shapefiles/Map_Overview_Region.shp')
%% Make Pie Plot

X=[numel(find(SM==-3)) numel(find(SM==-2)) numel(find(SM==2)) numel(find(SM==3))]./numel(find(~isnan(SM)));

cmap=[0.4,0,0;
    1,0,0;
    0,1,1;
    0,0,1];
figure
ax=axes;
pie(X,'%.3f%%');
colormap(ax,cmap)

%% Laplacian Error of Synthetic Surface Compared to Topographic Surface
[DM,DS,km,sl] = Error_Covariance_Mat(C,BS,50);

%
LPError_Slope_Bin(DM,sl)
%%

Error_Covariance_Plots(DM,DS,km,sl)


%%
S=Sphere_Curvature(1001,1);
US=CurveObj(S);
US=CurveCalc(US,0);
D=GRIDobj(S.X,S.Y,S.Z);
G=gradient8(D);
[Fx,Fy]=gradient(S.Z,S.dx,S.dx);
F=sqrt(Fx.^2+Fy.^2);
[Y,X]=meshgrid(S.X,S.Y);
Sl=sqrt(X.^2+Y.^2)./S.Z;



PE=abs(G.Z-US.CMAP.Sl)./abs(US.CMAP.Sl).*100;
b=bin(US.CMAP.az(~isnan(PE)),PE(~isnan(PE)),2000);
scatter(b(:,1),b(:,2))
hold on

%%

D=DEM;
D.Z=C.DEM.ZFilt;
G=gradient8(D);
PE=abs(C.CMAP.Sl-G.Z)./abs(C.CMAP.Sl).*100;
b2=bin(C.CMAP.az(~isnan(PE)),PE(~isnan(PE)),200);
scatter(b2(:,1),b2(:,8))
b3=bin(C.CMAP.az(~isnan(PE)),C.CMAP.KG(~isnan(PE)),200);
% scatter(b3(:,2),b2(:,2)-b(:,2))

%%
figure
plot(S.X,S.Z(500,:))
yyaxis('right')
yline(1/S.R)
hold on
plot(S.X,US.CMAP.LP(500,:).*2)

PE=abs(US.CMAP.LP(500,:).*2-US.CMAP.KM(500,:))./abs(US.CMAP.KM(500,:)).*100;
mean(PE(:))

plot(S.X,US.CMAP.KM(500,:))
set(gca,'yscale','log')
%%
% dz=linspace(0,1.2,100);
% ae = -(1-(1+dz.^2).^(3/2)).*100;
figure; hold on
% plot(dz,ae); hold on
scatter(bs(:,1),bs(:,2))
plot(bss(:,1),bss(:,2))
% scatter(bs(:,1),bs(:,2)+bs(:,3))
% scatter(bs(:,1),bs(:,2)-bs(:,3))
ylim([0 200])
xlim([0 1.2])

%%
bm=bin(C.CMAP.KM(in),abs(C.CMAP.LP(in)./2-C.CMAP.KM(in)),1000);
figure; hold on
% plot(dz,ae); hold on
scatter(bm(:,1),bm(:,2))
% ylim([50 80])
%% Shape PDF
figure
smwin=10;
bp=BD.P_b; %./BD.P_a; bp(bp>1)=1;
bd=BD.P_d;% ./BD.P_a; bd(bd>1)=1;
bss=BD.P_ss;% ./BD.P_a; bss(bss>1)=1;
bas=BD.P_as;% ./BD.P_a; bas(bas>1)=1;
b = plot(10.^BD.area,smoothdata(bp,'gaussian',smwin),'color',cmap(4,:),'LineWidth',2); hold on
d = plot(10.^BD.area,smoothdata(bd,'gaussian',smwin),'color',cmap(1,:),'LineWidth',2); hold on
ss = plot(10.^BD.area,smoothdata(bss,'gaussian',smwin),'color',cmap(3,:),'LineWidth',2); hold on
as = plot(10.^BD.area,smoothdata(bas,'gaussian',smwin),'color',cmap(2,:),'LineWidth',2); hold on

%%
B=bin(C.CMAP.Sl,abs(C.CMAP.KM),20);
B2=bin(C.CMAP.Sl,abs(del2(C.DEM.ZFilt,C.DEM.dx,C.DEM.dy)),20);

%% Calculate curvature along given channels
S=C.Stream.S;
A=C.Stream.A;

%%\
figure

FD=FLOWobj(DEM,'preprocess','fill');
A=flowacc(FD).*FD.cellsize^2;
S=STREAMobj(FD,A>5700);
St=modify(S,'interactive','polyselect');
St=trunk(klargestconncomps(St));

%%
SA=struct;
[~,~,~,~,SA.a]=STREAMobj2array(St,DEM,A);
K=DEM;

K.Z=C.CMAP.K2;
[SA.x,SA.y,SA.d,SA.z,SA.k2]=STREAMobj2array(St,DEM,K);

K.Z=C.CMAP.K1;
[~,~,~,~,SA.k1]=STREAMobj2array(St,DEM,K);

K.Z=C.CMAP.KG;
[~,~,~,~,SA.kg]=STREAMobj2array(St,DEM,K);

K.Z=C.CMAP.KM;
[~,~,~,~,SA.km]=STREAMobj2array(St,DEM,K);

K.Z=C.CMAP.Sl;
[~,~,~,~,SA.sl]=STREAMobj2array(St,DEM,K);

%% Generate pdfs and map for given area threshold
Meanpdfs(DEM,C,BS,log10(740))


%% Generate filter scale comparison figure
Stack=Scale_Space(DEM,[50 0],50,30);
Scale_Comparison(Stack,LogA.Z,BS,100);

%%
f=50;
for i=1:30
    KM=Stack.KG(:,:,i);
    km(i)=sum(KM(BS.Z==1));
    ft(i)=f;
    f=f+50;
end
figure

plot(ft,1./sqrt(abs(km)).*km./abs(km))
%% Dependencies 

% Flow routing using precarved drainage network
function [A,LogA,S] = a_multi(DEM)
%Computes area grid for DEM using multiflow algorithm
    FD1  = FLOWobj(DEM,'preprocess','carve','mex',true);
    DEM = imposemin(FD1,DEM,0.0001);
    S   = STREAMobj(FD1,'minarea',1e6/DEM.cellsize^2);
    DEMc = DEM;
    DEMc.Z(S.IXgrid) = DEMc.Z(S.IXgrid)-100; 
    FD  = FLOWobj(DEMc,'multi');
    A  = flowacc(FD).*DEM.cellsize^2;    
    LogA=A;
    LogA.Z=log10(A.Z);
end

% Bin Data by drainage area for area defined by bs polygon

% function [BS,BD] = bindata(C,LogA,DEM,bs,nb)
% 
%     % generate a data structure holding binned topographic geometry
%     % metrics.
%     A=LogA.Z;
%     BS = polygon2GRIDobj(DEM,bs);
%     in = find(BS.Z==1);
% 
%     % define vector of area values
%     av = linspace(min(LogA.Z(in)),max(LogA.Z(in)),nb+1);
% 
%     for i = 1:nb
%         ain = find(BS.Z==1 & A>=av(i) & A< av(i+1));
%         a(i) = mean([av(i),av(i+1)]);
%         km(i) = mean(C.CMAP.KM(ain));
%         lp(i) = mean(C.CMAP.LP(ain));
%         kg(i) = mean(C.CMAP.KG(ain));
%         k1(i) = mean(C.CMAP.K1(ain));
%         k2(i) = mean(C.CMAP.K2(ain));
%         alpha(i) = mean(C.CMAP.A(ain));
%         sl(i) = mean(C.CMAP.Sl(ain));
%         az(i) = mean(C.CMAP.az(ain));
%         p_a(i) = numel(ain)/numel(find(BS.Z==1)); % Drainage area pdf
%         S = C.SMAP(ain);
%         p_b(i) = numel(find(S==3))/numel(ain);
%         p_d(i) = numel(find(S==-3))/numel(ain);
%         p_as(i) = numel(find(S==-2))/numel(ain);
%         p_ss(i) = numel(find(S==2))/numel(ain);
%         chi(i)=trapz(a(1:i).^-0.45);
% 
% 
%     end
% 
%     BD = struct('area',a,'KM',km,'KG',kg,'K1',k1,'K2',k2,'A',alpha,'Sl',sl,'Az',az, ...
%         'P_a',p_a,'P_b',p_b,'P_d',p_d,'P_as',p_as,'P_ss',p_ss,'Chi',chi,'LP',lp);
% end

