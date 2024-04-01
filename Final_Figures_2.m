addpath(genpath('/Users/nathanielklema/OneDrive - Fort Lewis College/Research_Projects/Curvature'))
% Build a stack of Curveobj for a range of low pass filter cutoffs
% clear 
DEM=GRIDobj('Umpqua_10m_1.tif');%DEM.Z=DEM.Z.*0.3048;

filtertype='lowpass';
Filter=[200 190];


U=CurveObj3(DEM);
U2=CurveObj2(DEM);

U=SpecFilt(U,Filter,filtertype);

U2=SpecFilt(U2,Filter,filtertype);
%U.DEM.ZFilt=imgaussfilt(U.DEM.Z,3,'Padding','symmetric');
U=CurveCalc(U,0);
U2=CurveCalc(U2,0);


U=RouteFlow(U,1e6);
clear Filter filtertype
    
% Define shape class color map
cmap=[0.4,0,0;
    1,0,0;
    1,1,1;
    1,1,1;
    0,1,1;
    0,0,1];

%% Bin Everything - Just one profile
clear a nb nd svec in ns ss as S Sl Slv ast sst nbt ndt kg kmn k1 k2 dt1p dt2p st s8 sg
na=200;

[az,el,~] = cart2sph(U.CMAP.NX,U.CMAP.NY,U.CMAP.NZ);
Sl=(tan((pi/2)-el));

A=log10(U.Stream.A.Z); %log10(U(1).Stream.A.Z);
av=linspace(min(A(:)),max(A(:)),na);
da=av(2)-av(1);

[Gx, Gy]=gradient(U.DEM.ZFilt,DEM.cellsize,DEM.cellsize);
Gs=abs(Gx)+abs(Gy);

S=U.SMAP;
    
% All for alignment 
SVNX=U.CMAP.NX./(sqrt(U.CMAP.NX.^2+U.CMAP.NY.^2));
SVNY=U.CMAP.NY./(sqrt(U.CMAP.NX.^2+U.CMAP.NY.^2));

[thn,~]=cart2pol(SVNX,SVNY);
[thk1,~]=cart2pol(U2.CMAP.K1U,U2.CMAP.K1V);
[thk2,~]=cart2pol(U2.CMAP.K2U,U2.CMAP.K2V);

dt1=abs(cos(thn-thk1));
dt2=abs(cos(thn-thk2));

     
    for j=1:na
        in=find(A>=(av(j)-da) & A<(av(j)+da));
        svec=S(in);
        nb(j)=numel(find(svec==3))./numel(in);
        nd(j)=numel(find(svec==-3))./numel(in);
        ss(j)=numel(find(svec==2))./numel(in);
        as(j)=numel(find(svec==-2))./numel(in);
        a(j)=numel(in)./numel(S);
        mnb(j)=mean(U.CMAP.KG(svec==3));
        Slv(j)=mean(Sl(in));
        El(j)=mean(U.DEM.ZFilt(in));
%         E(i,j)=mean(A(in).*s(in));

%         dt1p(i,j)=numel(find(dt1(in)>0.5))./numel(in);
%         dt2p(i,j)=numel(find(dt2(in)>0.5))./numel(in);
        st(j)=mean(Sl(in));
        sg(j)=mean(Gs(in));
        s8(j)=mean(U.Stream.G.Z(in));

        % Alignement
        dt1p(j)=mean(dt1(in));
        dt2p(j)=mean(dt2(in));
        
        nbt(j)=numel(find(svec==3));
        ndt(j)=numel(find(svec==-3));
        sst(j)=numel(find(svec==2));
        ast(j)=numel(find(svec==-2));

        kmn(j)=mean(U.CMAP.KM(in));
        kg(j)=mean(U.CMAP.KG(in)); kms(j)=min(U.CMAP.KM(in));
        k1(j)=(mean(U.CMAP.K1(in)));
        k2(j)=(mean(U.CMAP.K2(in)));
        lp(j)=mean(U.CMAP.LP(in)./2);

    end

    winsz=15;
    window='gaussian';

    nb=smoothdata(nb,window,winsz);
    ss=smoothdata(ss,window,winsz);
    as=smoothdata(as,window,winsz);
    nd=smoothdata(nd,window,winsz);
    a=smoothdata(a,window,winsz);

    nbt=smoothdata(nbt,window,winsz);
    sst=smoothdata(sst,window,winsz);
    ast=smoothdata(ast,window,winsz);
    ndt=smoothdata(ndt,window,winsz);

    kmn=smoothdata(kmn,window,winsz);
    kg=smoothdata(kg,window,winsz);
    k1=smoothdata(k1,window,winsz);
    k2=smoothdata(k2,window,winsz);
    lp=smoothdata(lp,window,winsz);

    dt2p=smoothdata(dt2p,window,winsz);
    dt1p=smoothdata(dt1p,window,winsz);

    %st=smoothdata(st,window,winsz);
    %sg=smoothdata(sg,window,winsz);
    %s8=smoothdata(s8,window,winsz);


% DEfine thresholds

% clear th1 th2 th3 th4 th5 th6
% kgv=[];
% for i=2:numel(kg)-1
%     if abs(kg(i-1))>abs(kg(i)) && abs(kg(i+1))>abs(kg(i)) 
%         kgv=[kgv i];
%     end
% 
%     if abs(kmn(i-1))>abs(kmn(i))  && abs(kmn(i+1))>abs(kmn(i))
%         th2=10^av(i);
%     end
%     
% end
% kgv=10.^av(kgv(abs(kg(kgv))<0.4e-5));
th1 = 334;
th2=750;
th3=1317;
th4=4488;
th5=83000;

th6=980000;

% th1 = 108;
% th2=365;
% th3=614;
% th4=2920;
% 
% th6=980000;


%
figure

subplot(5,1,1)
area(10.^av,a.*100)
set(gca,'xscale','log','xticklabels',[])
xlim([10.^av(1) 10^7])
% ylim([0 18]);
xline(th1,'--','LineWidth',0.5);
xline(th2,'-','LineWidth',2);
xline(th3,'--','LineWidth',0.5);
xline(th4,'--','LineWidth',0.5);
xline(th6,'--','LineWidth',0.5);


subplot(5,1,2)
St=plot(10.^av,st,'Color',[0,0,0,0.5]);
% St.MarkerFaceAlpha=0.3;
% St.MarkerFaceColor='k';
St.LineWidth=4;
set(gca,'xscale','log','xticklabels',[])
set(gca,'yscale','log')
hold on

S8=plot(10.^av,s8,'--','Color',[1,0.6,0]);
% S8.MarkerFaceColor='r';
% S8.MarkerFaceAlpha=0.3;
S8.LineWidth=2;

Sg=plot(10.^av,sg,'-.','Color',[0.4,0,0]);
% Sg.MarkerFaceColor='y';
% Sg.MarkerFaceAlpha=0.3;
Sg.LineWidth=2;

% D=DEM;
% % D.Z=U.DEM.ZFilt;
% D=fillsinks(D);
% FD=FLOWobj(D,'preprocess','carve');
% AST=flowacc(FD).*D.cellsize^2;
% S=STREAMobj(FD,AST>1e3);
% SA=slopearea(S,DEM,AST./D.cellsize^2,'areabins',500,'streamgradient','robust','plot',0);
% sa=plot(SA.a,SA.g,'LineWidth',1,'Color','b');
xlim([10.^av(1) 10^7])

ylabel('Slope');

xline(th1,'--','LineWidth',0.5);
xline(th2,'-','LineWidth',2);
xline(th3,'--','LineWidth',0.5);
xline(th4,'--','LineWidth',0.5);
xline(th6,'--','LineWidth',0.5);

legend([St S8 Sg],{'Slope of tangent plane','8-point neighborhood gradient','\nabla z'})
xlim([min(10.^av) 10^7])
   
% ylim([0.01 1])


% Principal Curvatures
subplot(5,1,4)
plot(10.^av,k1,'LIneWidth',2)
set(gca,'xscale','log','xticklabels',[])
hold on
plot(10.^av,k2,'LIneWidth',2)
yline(0)
xlim([min(10.^av) 10^7])
ylim([-0.015 0.015])
% xline(th1); xline(th2); xline(th3); xline(th4); xline(th5); xline(th6);
yyaxis('right')
plot(10.^av,dt1p)
set(gca,'xscale','log')
hold on
plot(10.^av,dt2p)
ylim([0 1])

xline(th1,'--','LineWidth',0.5);
xline(th2,'-','LineWidth',2);
xline(th3,'--','LineWidth',0.5);
xline(th4,'--','LineWidth',0.5);
xline(th6,'--','LineWidth',0.5);


subplot(5,1,3)
plot(10.^av,kg)
set(gca,'xscale','log','xticklabels',[])
hold on
ylim([-1.5 1.5].*1e-5)
yline(0)

xline(th1,'--','LineWidth',0.5);
xline(th2,'-','LineWidth',2);
xline(th3,'--','LineWidth',0.5);
xline(th4,'--','LineWidth',0.5);
xline(th6,'--','LineWidth',0.5);
yyaxis('right')
plot(10.^av,kmn)

xlim([10.^av(1) 10^7])
ylim([-7 7].*1e-3)
%%open 
% xline(th1); xline(th2); xline(th3); xline(th4); xline(th5); xline(th6);


subplot(5,1,5)
plot(10.^av,as,'r','LineWidth',3);
set(gca,'xscale','log')
hold on


% Plot Basins as function of Area
plot(10.^av,nb,'b','LineWidth',3);


% Plot Synformal Saddles as function of Area
plot(10.^av,ss,'c','LineWidth',3);

% Plot Domes as function of Area
plot(10.^av,nd,'Color',[0.4 0 0],'LineWidth',3);

xlabel('Upstream Drainage Area (m^2)')
ylabel('P(C|A)')

ylim([0 0.7])
xlim([10.^av(1) 10^7])

xline(th1,'--','LineWidth',0.5);
xline(th2,'-','LineWidth',2);
xline(th3,'--','LineWidth',0.5);
xline(th4,'--','LineWidth',0.5);
xline(th6,'--','LineWidth',0.5);

%%
figure
subplot(5,1,1)

plot(10.^av,a,'k','LineWidth',2)
set(gca,'xscale','log','xticklabels',[]);
xlim([10.^av(1) 10^7])
ylabel('P(A)')

subplot(5,1,[2,3])
 % Plot Antiformal Saddles as function of Area
asf=plot(10.^av,ast./numel(U.SMAP),'r','LineWidth',3);
set(gca,'xscale','log')
hold on

% Plot Basins as function of Area
bf=plot(10.^av,nbt./numel(U.SMAP),'b','LineWidth',3);

% Plot Synformal Saddles as function of Area
ssf=plot(10.^av,sst./numel(U.SMAP),'c','LineWidth',3);

% Plot Domes as function of Area
df=plot(10.^av,ndt./numel(U.SMAP),'Color',[0.4 0 0],'LineWidth',3);

xlim([10.^av(1) 10^7])
ylabel('P(C)')
legend([bf ssf asf df],{'Antiformal Saddles','Basins','Synformal Saddles','Domes'},...
    'location','northeast')

subplot(5,1,[4 5])
plot(10.^av,as,'r','LineWidth',3);
set(gca,'xscale','log')
hold on


% Plot Basins as function of Area
plot(10.^av,nb,'b','LineWidth',3);


% Plot Synformal Saddles as function of Area
plot(10.^av,ss,'c','LineWidth',3);

% Plot Domes as function of Area
plot(10.^av,nd,'Color',[0.4 0 0],'LineWidth',3);

xlabel('Upstream Drainage Area (m^2)')
ylabel('P(C|A)')
xlim([10.^av(1) 10^7])
ylim([0 0.8])

% xline(th1); xline(th2); xline(th3); xline(th4); xline(th5); xline(th6);



%%
 D=DEM;
D.Z=U.DEM.ZFilt;
D=fillsinks(D);
FD=FLOWobj(D,'preprocess','carve');
Ad=flowacc(FD).*D.cellsize^2;
S=STREAMobj(FD,Ad>1e3);
ST=trunk(klargestconncomps(modify(S,'interactive','polyselect')));
order = ST.orderednanlist;
            dist  = ST.distance;
           
            zz = getnal(ST,DEM); % gives the GRIDobj Z values of grid sqaures crossed by the streamobj (referenced in the IXgrid)
            I     = ~isnan(order);
            d     = nan(size(order));
            d(I)  = dist(order(I)); % Gives distance along streamlength of each point
            z     = nan(size(order)); %  This makes all of the arrays the same size as the orderednanlist
            z(I)  = zz(order(I));
            x     = nan(size(order));
            x(I)  = ST.x(order(I));
            y     = nan(size(order));
            y(I)  = ST.y(order(I));
            aa    = getnal(ST,Ad);
            a     = nan(size(order));
            a(I)  = aa(order(I));
            
%%
   K1=DEM;
   K1.Z=U.CMAP.KG;
   k1=getnal(ST,K1);
   k1p=nan(size(order));
   k1p(I)=k1(order(I));
% Kmp=nan(size(order));
% Kmp(I)=U.CMAP.K2M(order(I));
%%

x=linspace(min(d),max(d),1000);
zi=interp1(d(~isnan(d)),z(~isnan(d)),x);
ai=interp1(d(~isnan(d)),a(~isnan(d)),x);
zs=smoothdata(zi,'gaussian',20);
as=smoothdata(ai,'gaussian',20);
dx=x(2)-x(1);
g=gradient(zs,dx);
b=bin(log10(ai),g,20);

%%
plot(10.^av,min(st))
asf=fill(x2, [max(st), fliplr(min(st))], 'b');
asf.FaceAlpha=0.5;
plot(10.^av,st(ex,:),'b','LineWidth',2);

plot(10.^av,max(s8))
set(gca,'xscale','log')
set(gca,'yscale','log')
hold on
plot(10.^av,min(s8))
asf=fill(x2, [max(s8), fliplr(min(s8))], 'b');
asf.FaceAlpha=0.5;
plot(10.^av,s8(ex,:),'b','LineWidth',2);

%% Alignment figure
figure
plot(10.^av,dt1p)
set(gca,'xscale','log')
hold on
plot(10.^av,dt2p)
xlim([min(10.^av) 10^7])
ylim([0 1])
% xline(th1); xline(th2); xline(th3); xline(th4); xline(th5); xline(th6);


%% Map colluvial hollows
S=U.SMAP;
A=U.Stream.A.Z.*U.CMAP.A;
SB=S; 
SB(SB<3)=0;
SB=im2bw(abs(SB));
BN=bwconncomp(SB);
BM=nan(size(S));
for i=1:BN.NumObjects
   
    if max(A(cell2mat(BN.PixelIdxList(i))))<=1e4 %&& min(A(cell2mat(BN.PixelIdxList(i))))>th3
        BM(cell2mat(BN.PixelIdxList(i)))=1;
    end

end
CH=DEM;
CH.Z=nan(size(DEM.Z));
BM(isnan(BM))=0;
BM(BM==-1)=0;
BM=imbinarize(BM);
[B,L]=bwboundaries(BM,'holes');
for k = 1:length(B)
   b = B{k};
   [r,~]=size(b);
   for i=1:r
       CH.Z(b(i,1),b(i,2))=1;
   end
end
[~,x,y] = GRIDobj2polygon(CH,'multipart',1);

%% Region Selection
D.Z=U.DEM.ZFilt;
D=fillsinks(D);
FD=FLOWobj(D,'multi');
AST=flowacc(FD).*D.cellsize^2;

%%
S=STREAMobj(FD,AST>th6);
SA =slopearea(S,DEM,AST./D.cellsize^2,'areabins',10000,'streamgradient','robust','plot',0);


D=DEM; D.Z=U.DEM.ZFilt;
A=U.Stream.A.Z;
St1=U.SMAP;
in=find(A>0); % & AST<th6);
St=nan(size(St1));
St(in)=St1(in);

% Map
figure
imageschs(D,St,'colormap',cmap)
hold on
%  plot(x,y,'w','LineWidth',2)
% xlim([4.2 4.27].*10^5)

%% Slope
bn=50;
sh=[];
[sh,et] = histcounts(Sl(in),bn, 'Normalization', 'probability');
[gsh,eg] = histcounts(Gs(in),bn, 'Normalization', 'probability');
[s8,e8] = histcounts(U.Stream.G.Z(in),bn, 'Normalization', 'probability');
%[scn,ecn] = histcounts(SA.g,bn, 'Normalization', 'probability');
bc1=[]; bcg=[]; bc8=[]; bcn=[];
for i=1:numel(sh)
    bc1(i)=mean([et(i), et(i+1)]);
    bcg(i)=mean([eg(i), eg(i+1)]);
    bc8(i)=mean([e8(i), e8(i+1)]);
    %bcn(i)=mean([ecn(i), ecn(i+1)]);

end
% 
figure
smwin=10;
as8=area(bc8,smoothdata(s8,'gaussian',smwin),'FaceColor',[1 0.6 0],'FaceAlpha',0.1); hold on
plot(bc8,smoothdata(s8,'gaussian',smwin),'LIneWidth',2,'Color',[1 0.6 0]);
asg=area(bcg,smoothdata(gsh,'gaussian',smwin),'FaceColor',[0.4 0 0],'FaceAlpha',0.1);
plot(bcg,smoothdata(gsh,'gaussian',smwin),'LIneWidth',2,'Color',[0.4 0 0]);
ast=area(bc1,smoothdata(sh,'gaussian',smwin),'FaceColor','k','FaceAlpha',0.1); hold on
% area(bcn,smoothdata(scn,'gaussian',smwin),'FaceColor','b','FaceAlpha',0.1);
% plot(bcn,smoothdata(scn,'gaussian',smwin),'LIneWidth',2,'Color',[0 0 1]);
plot(bc1,smoothdata(sh,'gaussian',smwin),'k','LineWidth',2); hold on
legend([ast,asg,as8],{'Slope of tangent plane','\nabla z','8-point neighorhood gradient'})
xlim([0 1.6]); ylim([0 0.08])
% 
% Pie
X=[numel(find(St1(in)==-3)), numel(find(St1(in)==-2)),0,0, numel(find(St1(in)==2)), numel(find(St1(in)==3))]./numel(in);
figure
pie(X);
colormap(cmap)

%% Mean Curvature Inflection point
th2=934;
D=DEM; D.Z=U.DEM.ZFilt;
A=U.Stream.A.Z.*U.CMAP.A;
in1=find(A<th2);
in2=find(A>=th2);
St=nan(size(A));
St(in1)=1;
St(in2)=-1;

cmap2=[0,0,1;1,0,0];
figure
hillshade(D)
% xlim([4.2 4.27].*10^5)

DS=D;
DS.Z=St;
figure
imagesc(DS);
colormap(cmap2)
% xlim([4.2 4.27].*10^5)
hold on
figure

plot(x,y,'w','LineWidth',2)
% xlim([4.2 4.27].*10^5)
%%
bn=100;
sh1=[];
[sh1,et1] = histcounts(Sl(in1),bn, 'Normalization', 'probability');
sh2=[];
[sh2,et2] = histcounts(Sl(in2),bn, 'Normalization', 'probability');


kmh1=[];
[kmh1,em1] = histcounts(U.CMAP.KG(in1),bn, 'Normalization', 'probability');
kmh2=[];
[kmh2,em2] = histcounts(U.CMAP.KG(in2),bn, 'Normalization', 'probability');

bc1=[]; bc2=[]; emc1=[]; emc2=[];
for i=1:numel(sh1)
    bc1(i)=mean([et1(i), et1(i+1)]);
    bc2(i)=mean([et2(i), et2(i+1)]);

    emc1(i)=mean([em1(i), em1(i+1)]);
    emc2(i)=mean([em2(i), em2(i+1)]);


end

figure; hold on
smwin=30;
plot(bc1,smoothdata(sh1,'gaussian',smwin))
plot(bc2,smoothdata(sh2,'gaussian',smwin))

figure; hold on
plot(emc1,smoothdata(kmh1,'gaussian',smwin))
plot(emc2,smoothdata(kmh2,'gaussian',smwin))

%%
UB=imbinarize(U.CMAP.KM);
