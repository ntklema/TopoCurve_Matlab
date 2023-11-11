addpath('/Users/nathanielklema/OneDrive - University Of Oregon/Cascades_Curvature_Final/Data/DEMS');
addpath('/Users/nathanielklema/Library/Mobile Documents/com~apple~CloudDocs/Cascades_Curvature/Coast_Range/Curvature_Shapes');
addpath('/Users/nathanielklema/OneDrive - Fort Lewis College/Research_Projects/Curvature/Data/DEMS')
%% Error in Laplacian figure
% DEM=GRIDobj('LIdar_Clip_1.tif');
% DEM.Z=DEM.Z.*0.3048;

[rad,c]=size(DEM.Z);
in=round(c/2);
x=linspace(DEM.georef.SpatialRef.XWorldLimits(1),DEM.georef.SpatialRef.XWorldLimits(2),c);
y=linspace(DEM.georef.SpatialRef.YWorldLimits(1),DEM.georef.SpatialRef.YWorldLimits(2),rad);
cross=DEM.Z(in,:);
 cross=smoothdata(cross,'movmean',400);
dx=x(2)-x(1);
LP=del2(cross,dx);
%  LP=smoothdata(LP,'movmean',200);
curve=(del2(cross,dx))./((1+gradient(cross,dx).^2).^(3/2));
%  curve=smoothdata(curve,'movmean',200);


% 
% subplot(2,2,1)
% imagesc(DEM2); hold on
% yline(y(in),'LineWidth',2)

figure
plot(x,cross); xlim([x(1) max(x)])

yyaxis('right')
plot(x,LP); hold on
plot(x,curve,'k'); 




%%
bn=100;
SB=bin(log10(U.Stream.A.Z(:)),(U.Stream.G.Z(:)),bn);
K1B=bin(log10(U.Stream.A.Z(:)),(U.CMAP.K1(:)),bn);
K2B=bin(log10(U.Stream.A.Z(:)),(U.CMAP.K2(:)),bn);
KMB=bin(log10(U.Stream.A.Z(:)),(U.CMAP.KM(:)),bn);
KGB=bin(log10(U.Stream.A.Z(:)),((U.CMAP.KG(:))),bn);
figure
subplot(1,4,1); hold on
scatter(10.^SB(:,1),SB(:,2));
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Area'); ylabel('Slope'); 

subplot(1,4,2); hold on

scatter(10.^KGB(:,1),KGB(:,2)); 
set(gca,'xscale','log')
yline(0)
xlabel('Area'); ylabel('Gaussian Curvature');

%
% subplot(1,4,2)
% scatter(10.^KMB(:,1),KMB(:,2))
% set(gca,'xscale','log')
% yline(0)
% xlabel('Area')

subplot(1,4,3); hold on
scatter(10.^K1B(:,1),K1B(:,2));
set(gca,'xscale','log')
yline(0)
xlabel('Area'); ylabel('K1')

subplot(1,4,4); hold on
scatter(10.^K2B(:,1),K2B(:,2));
set(gca,'xscale','log')
yline(0)
xlabel('Area'); ylabel('K2')

%% Evolution of curvature metrics along profile
% U=StreamExt(U);
 nb=200;
 b=bin(log10(U.Stream.a),U.Stream.k1,nb);

figure

scatter(10.^b(:,1),b(:,2));
set(gca,'xscale','log')
set(gca,'yscale','log')

 %% Slope-Area plot
Filter=[90,100];
filtertype='lowpass';
DEM1=GRIDobj("Umpqua_10m_1.tif");
C1=CurveObj2(DEM1);
C1=SpecFilt(C1,Filter,filtertype);
C1=CurveCalc(C1,0);
C1=RouteFlow(C1,1e5);


%
DEM2=GRIDobj("Coast_10m_5.tif");
C2=CurveObj2(DEM2);
C2=SpecFilt(C2,Filter,filtertype);
C2=CurveCalc(C2,0);
C2=RouteFlow(C2,1e6);


%%
bn=20;

MG1=mean(C1.CMAP.KG,'all');
SC1=bin(log10(C1.Stream.A.Z(:)),(C1.Stream.G.Z(:)),bn);
GB1=bin(log10(C1.Stream.A.Z(:)),(C1.CMAP.KG(:)),bn);

MG2=mean(C2.CMAP.KG,'all');
SC2=bin(log10(C2.Stream.A.Z(:)),(C2.Stream.G.Z(:)),bn);
GB2=bin(log10(C2.Stream.A.Z(:)),(C2.CMAP.KG(:)),bn);

%
figure; hold on
scatter(10.^SC1(:,1),SC1(:,2))
set(gca,'xscale','log')
set(gca,'yscale','log')
scatter(10.^SC2(:,1),SC2(:,2))

figure; hold on
scatter(10.^GB1(:,1),GB1(:,2))
set(gca,'xscale','log')
scatter(10.^GB2(:,1),GB2(:,2))
yline(0)

%%  Angle diff between curvature and slope vectors

FVA=atan2d(U.Stream.FV.Z,U.Stream.FU.Z);
FVC=atan2d(U.CMAP.K1V,U.CMAP.K1U);
[rad,c]=size(FVA);

for i=1:rad
    for j=1:c
        u=[U.Stream.FU.Z(i,j), U.Stream.FV.Z(i,j)];
        v=[U.CMAP.K1V(i,j), U.CMAP.K1U(i,j)];
        CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
        An(i,j) = real(acosd(CosTheta));

        if 180-An(i,j)<90
            An(i,j)=180-An(i,j);
        end
    end
end



%
B=bin(log10(U.Stream.A.Z(:)),An,100);
figure
scatter(10.^B(:,1),B(:,2))
set(gca,'xscale','log')
%% Scale Dependence of curvature metrics (comparitive)
DEM1=GRIDobj('Umpqua_10m_1.tif');
DEM2=GRIDobj('King_10m_1.tif');

C1=CurveObj2(DEM1);
C2=CurveObj2(DEM2);

Filter=[90,100];
filtertype='lowpass';

for i=1:20
    C1=SpecFilt(C1,Filter,filtertype);
    C2=SpecFilt(C2,Filter,filtertype);

end

%% Bin by stream order
so=U.Stream.so;
IXgrid=U.Stream.S.IXgrid;
sov=min(so):1:max(so);


mn=[];
mna=[];
for i=1:numel(sov)
    in=IXgrid(so==sov(i));
    mn(i)=mean(U.CMAP.KG(in));
    mna(i)=mean(U.Stream.A.Z(in));
%     mna(i)=min(U.Stream.A.Z(in));
end

figure
% scatter(sov,mna,'filled')
% set(gca,'yscale','log')
% yyaxis('right')
scatter(sov,mn,'filled')
%  set(gca,'yscale','log')
yyaxis('right')
scatter(sov,mna,'filled')
set(gca,'yscale','log')
%% Area CDF

nx=100;
loga=log10(U.Stream.A.Z);
na=numel(loga);
maxa=max(max(loga));
mina=min(min(loga));
avec=linspace(mina,maxa,nx);

for i=1:nx
    cdf(i)=numel(find(loga<=avec(i)))./na.*100;
end

% figure
scatter(10.^avec,cdf)
set(gca,'xscale','log')

%% Shape class distribution versus dranage area

nx=30;
loga=log10(U.Stream.A.Z);
na=numel(loga);
maxa=max(max(loga));
mina=min(min(loga));
avec=linspace(mina,maxa,nx);
Abin=[]; D=[]; ASS=[]; B=[]; SSS=[]; Aperc=[]; cdf=[];
for i=1:nx-1
    in=find(loga>=avec(i) & loga<avec(i+1));
    Smat=U.SMAP(in);
    Abin(i)=(avec(i)+avec(i+1))/2;
    Aperc(i)=numel(in)/na;
    D(i)=numel(find(Smat==-3))/numel(Smat)*Aperc(i)*100;
    ASS(i)=numel(find(Smat==-2))/numel(Smat)*Aperc(i)*100;
    B(i)=numel(find(Smat==3))/numel(Smat)*Aperc(i)*100;
    SSS(i)=numel(find(Smat==2))/numel(Smat)*Aperc(i)*100;
    K2(i)=mean(U.CMAP.K2(in));
    KG(i)=mean(U.CMAP.KG(in));
    KGstd(i)=std(U.CMAP.KG(in));
    KM(i)=mean(U.CMAP.KM(in));
    K1(i)=mean(U.CMAP.K1(in));
    S(i)=mean(U.Stream.G.Z(in));
%    Tt(i)=mean(abs(A(in)));
%     P(i)=numel(find(Smat==0))/numel(Smat)*Aperc(i);

%     D(i)=numel(find(Smat==-3))*Aperc(i);
%     AS(i)=numel(find(Smat==-2))*Aperc(i);
%     B(i)=numel(find(Smat==3))*Aperc(i);
%     SS(i)=numel(find(Smat==2))*Aperc(i);
%     P(i)=numel(find(Smat==0))/numel(Smat);

%     D(i)=numel(find(Smat==-3))/numel(Smat)*100;
%     AS(i)=numel(find(Smat==-2))/numel(Smat)*100;
%     B(i)=numel(find(Smat==3))/numel(Smat)*100;
%     SS(i)=numel(find(Smat==2))/numel(Smat)*100;

    Gmat=U.CMAP.KG(in);
    G(i)=mean(Gmat)*Aperc(i);
    cdf(i)=numel(find(loga<=Abin(i)))./na;
    


end

figure
plot(10.^Abin,D,'r','LineWidth',2); hold on
plot(10.^Abin,ASS,'k','LineWidth',2);
plot(10.^Abin,SSS,'c','LineWidth',2);
plot(10.^Abin,B,'b','LineWidth',2);
xlabel('Area (m^2)'); ylabel('Weighted % Composition');

% scatter(10.^Abin,P,'b','LineWidth',2);
% scatter(10.^Abin,Aperc,'m','LineWidth',2);
set(gca,'xscale','log')
% scatter(10.^Abin,Aperc,'m','LineWidth',2);
set(gca,'yscale','log')
yyaxis('right')
plot(10.^Abin,cdf,'LineWidth',2)
ylabel('% of landscape')
legend('Domes','Antiformal Saddles','Synformal Saddles','Basins','Drainage Area CDF')
%% Curvature Spectra

filt=[40 50];
nx=10;
df=50;

for i=1:nx
    U=SpecFilt(U,filt,'lowpass');
    U=CurveCalc(U,0);
    GM(i)=mean(abs(U.CMAP.KG(:)));
    MM(i)=mean(abs(U.CMAP.KM(:)));
    F(i)=filt(2);
    filt=filt+df;
end

figure
scatter(F,GM)

%% Fourier Method Issue
[rad,c]=size(U.DEM.Z);
x=linspace(min(U.DEM.X(:)),max(U.DEM.X(:)),c);
y=ones(size(x)).*350;
for i=2:100
    y(i)=y(i-1)-2;
end
for i=101:400
    y(i)=y(i-1)+1.25;
end
for i=401:c-100
    y(i)=y(i-1)-1.25;
end

 y(c-99:end)=0;
for i=1:rad
    U.DEM.Z(i,:)=y;

end
Filter=[200 300];
U=SpecFilt(U,Filter,filtertype);
U=CurveCalc(U,0);


% Now use gaussian convolution filter
U2=U;
U2.DEM.ZFilt=imgaussfilt(U2.DEM.Z,2);
% U2.DEM.ZFilt = butterworthbpf(U2.DEM.Z,30,120,4);
U2=CurveCalc(U2,0);
%%
figure

subplot(3,2,1);

imagesc(U.DEM.ZFilt)
title('Spectral Filter')
colorbar

subplot(3,2,3)
imagesc(U.CMAP.K1)
title('K1')
colorbar

subplot(3,2,5)
plot(U.DEM.Z(1,:))
hold on
plot(U.DEM.ZFilt(1,:))
ylim([-100 600])
yyaxis('right')
plot(U.CMAP.K1(100,:))
xlim([0 870])


subplot(3,2,2)

imagesc(U2.DEM.ZFilt)
title('Gaussian Convolution')
colorbar

subplot(3,2,4)
imagesc(U2.CMAP.K1)
title('K1')
colorbar

subplot(3,2,6)
plot(U2.DEM.Z(1,:))
hold on
plot(U2.DEM.ZFilt(1,:))
ylim([-100 600])
yyaxis('right')
plot(U2.CMAP.K1(100,:))
xlim([0 870])

%%
C=U.CMAP.K1;

mb=mean(C(U.SMAP==-3))
md=mean(C(U.SMAP==3))
mss=mean(C(U.SMAP==-2))
mas=mean(C(U.SMAP==2))

%% Gaussian Pyramid calculations
DEM=GRIDobj("Trans_30m_1.tif");
dx=DEM.cellsize;
X=[DEM.georef.SpatialRef.XWorldLimits(1) DEM.georef.SpatialRef.XWorldLimits(2)];
Y=[DEM.georef.SpatialRef.YWorldLimits(1) DEM.georef.SpatialRef.YWorldLimits(2)];
xo=X(1)+dx/2:dx:X(2)-dx/2;
yo=fliplr(Y(1)+dx/2:dx:Y(2)-dx/2);
Zo=(DEM.Z);
nx=512;
x=xo(1:nx);
y=fliplr(yo(1:nx));
Z=flipud(Zo(1:nx,1:nx));

D=struct; D.X=x; D.Y=y; D.Z=Z;
U(1)=CurveObj2(D);
U(1)=CurveCalc(U(1),0);

U(1)=RouteFlow(U(1),1e6);
A=U(1).Stream.A;
A.Z=log10(A.Z);
A(1)=A;

%%
np=8;
for i=1:np
    i
    x=(x(1)+x(2))/2:2^i*dx:(x(end-1)+x(end))/2;
    y=(y(1)+y(2))/2:2^i*dx:(y(end-1)+y(end))/2;
    Z=impyramid(Z,'reduce');
    D=struct; D.X=x; D.Y=y; D.Z=Z;
    U(i+1)=CurveObj2(D);
    U(i+1)=CurveCalc(U(i+1),0);
    Az=A(1);
    Az.Z=impyramid(A(i).Z,'reduce');
    A(i+1)=Az;

end

%%
figure
for i=1:10
subplot(2,5,i)
C=U(i).CMAP.K1;
% C(C>0)=U(i).SMAP(C>0);
imagesc((C))
% title("\sigma="+round(sv(i),2)+"")
% colorbar
% scatter(i,numel(find((U(i).SMAP==3))));hold on
end

%%
% figure
rad=1;
C=-U(rad).CMAP.KG;
C=im2bw(C,0);
S1=double(bwskel(C));
C2=U(rad).CMAP.K1;
C2=im2bw(C2,0);
S2=double(bwskel(C2));

M=zeros(size(C));
M(S1==1)=1;
M(S2==1)=-1;
% M(U(r).Stream.A.Z<1e3)=0;
figure
imagesc(M)
% D=DEM;
% D.Z=U(r).DEM.ZFilt;
% imageschs(D,M,'colormap',[1,0,0;0,0,1])
% colorbar('off')
%% Gaussian Scale Space
clear U
DEM=GRIDobj("Umpqua_10m_1.tif");
s=1;
time=tic;
for i=1:40   
    U(i)=CurveObj3(DEM);
    U(i).DEM.ZFilt=imgaussfilt(DEM.Z,s);
    U(i)=CurveCalc(U(i),0);
    sv(i)=s*DEM.cellsize;
    s=s+0.5;
    if i==i
        U(1)=RouteFlow(U(1),1e6);
    end
    disp("Completed "+i+" iterations at time t = "+round(toc(time)/60,2)+" minutes")
end

%% Curvature Spectra of each pixel
clear U f
DEM=GRIDobj('Umpqua_30m_1.tif');
ns=20;
filtertype='lowpass';
Filter=[100 90];
time=tic;
for i=1:ns
    U(i)=CurveObj3(DEM);
    
    
    U(i)=SpecFilt(U(i),Filter,filtertype);
    U(i)=CurveCalc(U(i),0);

    if i==1
        U(i)=RouteFlow(U(i),1e6);
    end
    
    f(i)=Filter(1);
    Filter=Filter+20;
    disp("Completed "+i+" iterations at time t = "+round(toc(time)/60,2)+" minutes")
end


%%
for i=1:ns
    k1m(i)=mean(U(i).CMAP.K1(:));
    k2m(i)=mean(U(i).CMAP.K2(:));
    kgm(i)=mean((U(i).CMAP.KG(:)));

end

figure

subplot(1,3,1)
plot(f,k1m)
subplot(1,3,2)
plot(f,k2m)
subplot(1,3,3)
plot(f,kgm)

%%
s=[];
for i=1:numel(U)
    s(:,:,i)=U(i).CMAP.KG;
end
figure; hold on
[rad,c]=size(U(1).CMAP.K1);

for i=1:rad
    for j=1:c
    v=[s(i,j,1), s(i,j,2),s(i,j,3),s(i,j,4),s(i,j,5),s(i,j,6),s(i,j,7),s(i,j,8),s(i,j,9),s(i,j,10)];
    plot(f,v,'k','LineWidth',0.02)
    end
end

%%
DEM1=GRIDobj("Umpqua_1m_1.tif");
DEM2=resample(DEM1,30);
U2=CurveObj3(DEM2);
U2.DEM.ZFilt=imgaussfilt(U2.DEM.Z,3,'Padding','symmetric');
for i=1:10
    U(i)=CurveObj2(DEM1);
    U(i).DEM.ZFilt=imgaussfilt(U(i).DEM.Z,2*i,'Padding','symmetric');
    s(i)=i;
    DEM=DEM1;
    DEM.Z=U(i).DEM.ZFilt;
    DEM=resample(DEM,30);
    RMS(i)=sum((DEM.Z(~isnan(U2.DEM.ZFilt))-U2.DEM.ZFilt(~isnan(U2.DEM.ZFilt))).^2);

end

%%  Gorge Area Excess Figure - This is garbage code
DF=GRIDobj('/Users/nathanielklema/OneDrive - University Of Oregon/Gorge_Project_Final/Data/DEMs/Projected_Flexure.tif');
DF=resample(DF,DEM);

AK=U.CMAP.A;
AK(DEM.Z>DF.Z)=nan;
% DF.Z=
AK(DEM.Z<50)=nan;
% M=createmask(DEM);
%
DF=resample(DF,DEM);
AD=DEM;
AD.Z=AK;
AD=crop(AD,M,NaN);
DF=crop(DF,M,NaN);
DF.Z(isnan(AD.Z))=nan;
b=bin(DF.Z(~isnan(DF.Z)),AD.Z(~isnan(DF.Z)),30);
figure
scatter(b(:,1),b(:,2));

%%
[rad,c]=size(AD.Z);
DF.Z(isnan(DF.Z))=0;
for i=1:c
    u=DF.Z(:,i)';
    a=AD.Z(:,i)';
    uv(i)=mean(u(u>0));
    av(i)=mean(a(~isnan(a)));
end


%% NUmber of basins vs domes
clear SSN ASN BN DN NO
for i=1:numel(U)
    S=U(i).SMAP;
    SB=S; SD=S; ASS=S; SSS=S;
    SB(SB>-3)=0; SD(SD<3)=0; SSS(SSS~=-2)=0; ASS(ASS~=2)=0;
    SB=im2bw(abs(SB)); SD=im2bw(abs(SD)); SSS=im2bw(abs(SSS)); ASS=im2bw(abs(ASS));
    BN(i)=bwconncomp(SB).NumObjects;
    DN(i)=bwconncomp(SD).NumObjects;
    SSN(i)=bwconncomp(SSS).NumObjects;
    ASN(i)=bwconncomp(ASS).NumObjects;
    
    ba=regionprops(bwconncomp(ASS),'Area');
    ba = struct2cell(ba);

    BA(i)= sum([ba{:}])./numel(S);
    
end
NO=BN+DN+SSN+ASN;



%% Drainage Parameter Comparison
DEM=GRIDobj("Umpqua_10m_1.tif"); %DEM.Z=DEM.Z.*0.3048;
DEM.Z(isnan(DEM.Z))=0;
% DEM=resample(DEM,20);
% ClipRange=[7.18e5,7.34e5,3.956e6,3.968e6]; % ROI Gabilan
% ClipRange=[5.826e5,5.826e5+7e3,4.912e6,4.912e6+7e3]; % ROI Santiam
% ClipRange=[4.23e5,4.256e5,4.842e6,4.844e6]; % ROI Umqua Lidar
% ClipRange=[4.075e5,4.256e5,4.842e6,4.844e6]; % ROI King 1
% ClipRange=[5.74e5,5.78e5,5.047e6,5.05e6]; % ROI Ainsworth
U=CurveObj3(DEM);
% %
% U=CurveObj2(DEM,'clip',ClipRange);
Filter=[200,190];
filtertype='lowpass';
U=SpecFilt(U,Filter,filtertype);

% U.DEM.ZFilt=imgaussfilt(U.DEM.Z,4);
tic
U=CurveCalc(U,0);
toc

U=RouteFlow(U,1e7);

%% Spherical Maps
% Find elevation angle to calculate slope from normal
[az,el,~] = cart2sph(U.CMAP.NX,U.CMAP.NY,U.CMAP.NZ);
Sl=(((pi/2)-el));
D=DEM;
D.Z=U.DEM.ZFilt;

figure
subplot(2,2,1)
imageschs(D,S)

subplot(2,2,2)
x=[];
y=[];
nb=50;
[N,Xedges,Yedges] = histcounts2(U.CMAP.NY(:),U.CMAP.NX(:),nb,'Normalization','pdf');
N(N==0)=nan;
for i=1:nb
    x(i)=(Xedges(i)+Xedges(i+1))/2;
    y(i)=(Yedges(i)+Yedges(i+1))/2;
end
h=imagesc(x,y,N);
set(gca,'ydir','normal')
axis image off
set(h,'AlphaData', ~isnan(N))
colorbar
hold on
th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);
plot(xunit,yunit,'k','LineWidth',2)



subplot(2,2,3)
[phi,rad]=cart2pol(U.CMAP.NX,U.CMAP.NY);

[Np,Redges,Pedges] = histcounts2(S(:),phi(:),nb,'Normalization','pdf');
Np(Np==0)=nan;%
p=[];
r=[];

for i=1:nb
    p(i)=(Pedges(i)+Pedges(i+1))/2;
    r(i)=(Redges(i)+Redges(i+1))/2;
end
hp=imagesc(p,r,Np);
set(gca,'ydir','normal')
% ylim([0 1])
% axis image off
set(hp,'AlphaData', ~isnan(Np))




%% Slope/Curvature Alignment

DEM=GRIDobj("Umpqua_10m_1.tif"); %DEM.Z=DEM.Z.*0.3048;
DEM.Z(isnan(DEM.Z))=0;

U2=CurveObj2(DEM);
U=CurveObj3(DEM);
% %
% U=CurveObj2(DEM,'clip',ClipRange);
Filter=[100,90];
filtertype='lowpass';
U=SpecFilt(U,Filter,filtertype);

% U.DEM.ZFilt=imgaussfilt(U.DEM.Z,4);
% U2.DEM.ZFilt=imgaussfilt(U.DEM.Z,4);

tic
U2=CurveCalc(U2,0);
U=CurveCalc(U,0);
toc

U=RouteFlow(U,1e7);

%% Make unit normal slope vectors 
SVNX=U.CMAP.NX./(sqrt(U.CMAP.NX.^2+U.CMAP.NY.^2));
SVNY=U.CMAP.NY./(sqrt(U.CMAP.NX.^2+U.CMAP.NY.^2));
A=U.Stream.A.Z;
[thn,~]=cart2pol(SVNX,SVNY);
[thk1,~]=cart2pol(U2.CMAP.K1U,U2.CMAP.K1V);
[thk2,~]=cart2pol(U2.CMAP.K2U,U2.CMAP.K2V);

dt1=abs(cos(thn-thk1));
dt2=abs(cos(thn-thk2));

 b1=bin(log10(A(:)),dt1(:),1000); 
 b2=bin(log10(A(:)),dt2(:),1000);
 figure
 plot(10.^b1(:,1),b1(:,2)); hold on
 plot(10.^b2(:,1),b2(:,2));
set(gca,'xscale','log')

%% Gradient Comparison
nb=100;
[az,el,~] = cart2sph(U.CMAP.NX,U.CMAP.NY,U.CMAP.NZ);
S=(((pi/2)-el));
[Gx, Gy]=gradient(U.DEM.ZFilt,DEM.cellsize,DEM.cellsize);
G=abs(Gx)+abs(Gy);
figure
b=bin(log10(U.Stream.A.Z(:)),G(:),nb);
scatter(10.^b(:,1),b(:,2)); hold on
set(gca','xscale','log','yscale','log')

b=bin(log10(U.Stream.A.Z(:)),U.Stream.G.Z(:),nb);
scatter(10.^b(:,1),b(:,2))

b=bin(log10(U.Stream.A.Z(:)),S(:),nb);
scatter(10.^b(:,1),b(:,2))
xlabel('Area (m^2)')
ylabel('Slope')
legend('\nabla z','8-point gradient','Tangential Slope')