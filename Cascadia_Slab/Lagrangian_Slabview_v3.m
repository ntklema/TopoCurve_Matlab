DEM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Bathy_15as_UTM.tif');

CAS=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/CASIE_TOC_UTM.tif');

FS=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/CASIE_FS_UTM.tif');

BG=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Gravity/Bouguer_UTM10m2.tif');

domainx=CAS.georef.XWorldLimits;
domainy=CAS.georef.YWorldLimits;
% DEM=resample(DEM,BG,'bilinear');
% FS = resample(FS, S2, 'bilinear');
% CAS = resample(CAS, S2, 'bilinear');

% 
% FS=resample(FS,DEM,'bilinear');
% CAS=resample(CAS,DEM,'bilinear');
% S2 = GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_UTM10N.tif');
% CAS = resample(CAS, S2, 'bilinear');
% FS = resample(FS, S2, 'bilinear');
% S2=resample(S2,DEM,'bilinear');
% 
% D=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Dip_Slab2_UTM10N.tiff');
% D=resample(D,DEM,'bilinear');
% 
% H=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/HeatFlow.tif');
% H=resample(H,DEM,'bilinear');
% 
% KM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/KM_UTM10N.tiff');
% KM=resample(KM,DEM,'bilinear');
% 
% K1=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_K1.tif');
% K1=resample(K1,DEM,'bilinear'); K1.Z=K1.Z./1e-6;
% 
% K2=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_K2.tif');
% K2=resample(K2,DEM,'bilinear');K2.Z=K2.Z./1e-6;
% 
% KG=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/KG_UTM10N.tif');
% KG=resample(KG,DEM,'bilinear');
% 
% SMAP=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_SMAP.tif');
% SMAP=resample(SMAP,DEM,'nearest');

% Import Shapefile Datasets
TS_Poly=shaperead('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Shapefiles/Casie_Wedgefront_Swath.shp');

Out=Slab_Project(CAS,TS_Poly,90-52,40e-3);

%% Make depth swath plot
FS.Z(FS.Z<CAS.Z)=CAS.Z(FS.Z<CAS.Z);
DEMr=resample(DEM,CAS);
dist=Out.dist.Z;
n=10;
cmap=(colormap(copper(n)));

time=Out.time.Z;
d=[0,0,1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4];
lab={'0 km','0 - 5 km (0 - 250 ky)','10 - 20 km (250 - 500 ky)','20 - 30 km (500 - 750 ky)','30 - 40 km (750 - 1 My)','40 - 50 km (1 My - 1.25 My)',...
    '50 - 60 km (1.25 My - 1.5 My)','60 - 70 km (1.5 My - 1.75 My)','70 - 80 km (1.75 My - 2 My)','80 - 90 km (2 My - 2.25 My)'};

nb=50;
f1=figure; hold on
set(f1,'units','inches','Position',[3,3,18,10],'color','w')


Th=[];
subplot(2,3,[1,2,4,5])
for i=1:n

    if i==1
        [X,Y]=getcoordinates(DEMr,'matrix');
        in=find(Out.S.Z==1);
        b=bin(Y(in),DEMr.Z(in),nb);
        bc=bin(Y(in),X(in),nb);

        [lat, lon] = projinv(projcrs(26910), bc(:,2)', bc(:,1)');
        plot(lat,smoothdata(b(:,2)./1e3,'gaussian',3),'k','LineWidth',4,'DisplayName','Bathymetric Depth'); hold on

        in=find(and(and(dist>2500,dist<=11000),~isnan(CAS.Z)));
        bc=bin(Out.trench_y(in),Out.trench_x(in),nb);
        b=bin(Out.trench_y(in),CAS.Z(in),nb);
        b2=bin(Out.trench_y(in),FS.Z(in),nb);
        [lat, ~] = projinv(projcrs(26910), bc(:,2)', bc(:,1)');


        toc=b(:,8)';
        fs=b2(:,8)';
        nin=find(toc~=0);
        toc=smoothdata(toc(nin),'gaussian',5)./1e3;
        fs=smoothdata(fs(nin),'gaussian',5)./1e3;

        fill([lat(nin) fliplr(lat(nin))],[toc fliplr(fs)],cmap(i,:),'FaceAlpha',0.8,'DisplayName',lab{i}); 
        plot(lat(nin),toc,'k','LineWidth',1,'HandleVisibility', 'off'); 
        plot(lat(nin),fs,'k','LineWidth',1, 'HandleVisibility', 'off')
        
        
    else

        in=find(and(and(dist>d(i),dist<=d(i+1)),~isnan(CAS.Z)));
        
        bc=bin(Out.trench_y(in),Out.trench_x(in),nb);
        b=bin(Out.trench_y(in),CAS.Z(in),nb);
        b2=bin(Out.trench_y(in),FS.Z(in),nb);
        [lat, ~] = projinv(projcrs(26910), bc(:,2)', bc(:,1)');
    
        
        toc=b(:,8)';
        fs=b2(:,8)';
        nin=find(toc~=0);
        toc=smoothdata(toc(nin),'gaussian',5)./1e3;
        fs=smoothdata(fs(nin),'gaussian',5)./1e3;

        fill([lat(nin) fliplr(lat(nin))],[toc fliplr(fs)],cmap(i,:),'FaceAlpha',0.8,'DisplayName',lab{i}); 
        plot(lat(nin),toc,'k','LineWidth',1,'HandleVisibility', 'off'); 
        plot(lat(nin),fs,'k','LineWidth',1,'HandleVisibility', 'off'); 
        xlabel('Latitude (^\circN)')
        ylabel('Depth (km)')
    end
end

% Annotation for CFTB
an_height=-4;
plot([45,48.3],[an_height an_height],'k','LineWidth',1,'HandleVisibility', 'off');
plot([45 45],[an_height-0.2 an_height+0.2],'k','LineWidth',1,'HandleVisibility', 'off');
plot([48.3 48.3],[an_height-0.2 an_height+0.2],'k','LineWidth',1,'HandleVisibility', 'off');
text(45.6,an_height+0.4,'Cascadia Fold-Thrust Belt','FontSize',14)

annotation('textarrow',[0.3 0.36], [0.9 0.86], ...
    'String', 'Astoria Fan  ', 'FontSize', 14);
annotation('textarrow',[0.53 0.51], [0.90 0.865], ...
    'String', ' Nitinat Fan', 'FontSize', 14);
set(gca,'fontsize', 14)
ylim([-17 -1])
xlim([42,49.8])

% Time Binning Plot
% in=find(and(~isnan(Out.time.Z),~isnan(CAS.Z)));
% bt=bin(Out.time.Z(in),FS.Z(in)-CAS.Z(in),20);
% 
% subplot(2,3,6)
% scatter(bt(:,1),bt(:,2))
legend()
%% make Time history plot v1
nb=21;
tv=linspace(0,2.5e6,nb);
dt=tv(2)-tv(1);
t=[];
diff=[];
d_25=[];
[~, n1, ~] = deg2utm(42, 127);
[~, n2, ~] = deg2utm(45.5, 127);
[~, n3, ~] = deg2utm(49, 127);

for i=1:nb-1
    t(i)=mean([tv(i), tv(i+1)]);
    in=find(and(and(time>=tv(i),time<tv(i+1)),~isnan(FS.Z-CAS.Z)));
    in2=find(and(and(and(time>=tv(i),time<tv(i+1)),~isnan(FS.Z-CAS.Z)),and(Out.trench_y>n1,Out.trench_y<=n2)));
    in3=find(and(and(and(time>=tv(i),time<tv(i+1)),~isnan(FS.Z-CAS.Z)),and(Out.trench_y>n2,Out.trench_y<=n3)));
    % diff(i)=mean(FS.Z(in)-CAS.Z(in),'omitnan')*40e-3*1.5*1000;
    % diff_s(i)=mean(FS.Z(in2)-CAS.Z(in2),'omitnan')*40e-3*1.5*1000;
    % diff_n(i)=mean(FS.Z(in3)-CAS.Z(in3),'omitnan')*40e-3*1.5*1000;

    diff(i)=mean(FS.Z(in)-CAS.Z(in),'omitnan')*1e-3*40e-6;
    diff_s(i)=mean(FS.Z(in2)-CAS.Z(in2),'omitnan')*1e-3*40e-6;
    diff_n(i)=mean(FS.Z(in3)-CAS.Z(in3),'omitnan')*1e-3*40e-6;

end

figure; hold on
% plot(t,d_25)
% plot(t,d_75)
t=t/1e6;
plot(t,diff,'LineWidth',2)
plot(t,diff_s,'LineWidth',2)
plot(t,diff_n,'LineWidth',2)

xlabel('Time before present (Ma)')
ylabel('Average subducted volume (km^3/km/yr)')
legend('Full trench: 41^\circ - 49^\circ N',"Southern segment: 41^\circ - 45.5^\circ N","Northern segment: 45.5^\circ - 49^\circ N",'location','northwest')

%% make Time history plot 21
nb=21;
tv=linspace(0,2.5e6,nb);
dt=tv(2)-tv(1);
t=[];
diff=[];
d_25=[];
[~, n1, ~] = deg2utm(42, 127);
[~, n2, ~] = deg2utm(45.5, 127);
[~, n3, ~] = deg2utm(49, 127);

for i=1:nb-1
    t(i)=mean([tv(i), tv(i+1)]);
    in=find(and(and(time>=tv(i),time<tv(i+1)),~isnan(FS.Z-CAS.Z)));
    in2=find(and(and(and(time>=tv(i),time<tv(i+1)),~isnan(FS.Z-CAS.Z)),and(Out.trench_y>n1,Out.trench_y<=n2)));
    in3=find(and(and(and(time>=tv(i),time<tv(i+1)),~isnan(FS.Z-CAS.Z)),and(Out.trench_y>n2,Out.trench_y<=n3)));
    diff(i)=mean(FS.Z(in)-CAS.Z(in),'omitnan')*40e-3*1.2;
    diff_s(i)=mean(FS.Z(in2)-CAS.Z(in2),'omitnan')*40e-3*1.2;
    diff_n(i)=mean(FS.Z(in3)-CAS.Z(in3),'omitnan')*40e-3*1.2;

    diff2(i)=mean(FS.Z(in)-CAS.Z(in),'omitnan')*40e-3*2.2;
    diff_s2(i)=mean(FS.Z(in2)-CAS.Z(in2),'omitnan')*40e-3*2.2;
    diff_n2(i)=mean(FS.Z(in3)-CAS.Z(in3),'omitnan')*40e-3*2.2;

    % diff(i)=mean(FS.Z(in)-CAS.Z(in),'omitnan')*1e-3*40e-6;
    % diff_s(i)=mean(FS.Z(in2)-CAS.Z(in2),'omitnan')*1e-3*40e-6;
    % diff_n(i)=mean(FS.Z(in3)-CAS.Z(in3),'omitnan')*1e-3*40e-6;

end

figure; hold on
% plot(t,d_25)
% plot(t,d_75)
t=t/1e6;
fill([t fliplr(t)],[diff fliplr(diff2)],'b',FaceAlpha=0.5)
fill([t fliplr(t)],[diff_s fliplr(diff_s2)],[1,0.8,0.8],FaceAlpha=0.8)
fill([t fliplr(t)],[diff_n fliplr(diff_n2)],[0.8,0.8,0.8],FaceAlpha=0.5)

xlabel('Time before present (Ma)')
ylabel('Average subducted mass (kT/km/yr)')
legend('Full trench: 41^\circ - 49^\circ N',"Southern segment: 41^\circ - 45.5^\circ N","Northern segment: 45.5^\circ - 49^\circ N",'location','northwest')
%% Gravity 
BG=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Gravity/Bouguer_UTM10m2.tif');
BG=crop(BG,domainx,domainy);
Out_G=Slab_Project(BG,TS_Poly,90-52,40e-3);

%%

dist=Out_G.dist.Z;
n=9;
cmap=(colormap(copper(n)));

time=Out_G.time.Z;
d=[0,0,1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4];

nb=100;
f1=figure; hold on
set(f1,'units','inches','Position',[3,3,18,10],'color','w')


Th=[];
subplot(2,3,[1,2,4,5]); hold on
for i=1:n

    if i==1
        in=find(and(and(dist>2500,dist<=5000),~isnan(BG.Z)));
        bc=bin(Out_G.trench_y(in),Out_G.trench_x(in),nb);
        b=bin(Out_G.trench_y(in),BG.Z(in),nb);
        [lat, ~] = projinv(projcrs(26910), bc(:,2)', bc(:,1)');


        % g=b(:,2)';
        % 
        % nin=find(g~=0);
        % g=smoothdata(g(nin),'gaussian',5);

        % fill([lat(nin) fliplr(lat(nin))],[g_min(nin) fliplr(g_max(nin))],cmap(i,:),'FaceAlpha',0.8); hold on
        % plot(lat(nin),g,'k','LineWidth',1); 


    else

        in=find(and(and(dist>d(i),dist<=d(i+1)),~isnan(BG.Z)));

        bc=bin(Out_G.trench_y(in),Out_G.trench_x(in),nb);
        b=bin(Out_G.trench_y(in),BG.Z(in),nb);

        [lat, ~] = projinv(projcrs(26910), bc(:,2)', bc(:,1)');


        g_mn=smoothdata(b(:,7)','gaussian',5);
        g_mx=smoothdata(b(:,6)','gaussian',5);

        nin=find(g~=0);
        g=smoothdata(g(nin),'gaussian',5);
   
        fill([lat(nin), fliplr(lat(nin))],[g_mx, fliplr(g_mn)],cmap(i,:),'FaceAlpha',0.8)
       
        plot(lat(nin),g_mn,'k','LineWidth',1); 
        plot(lat(nin),g_mx,'k','LineWidth',1); 

        xlabel('Latitude (^\circN)')
        ylabel('Depth (km)')
    end
end

%%

CASr=resample(CAS,BG);
FSr=resample(FS,BG);
in=find(and(~isnan(CASr.Z),abs(CASr.Z)<=1e4));
p=polyfit(abs(CASr.Z(in)),log(BG.Z(in)),1);
d=linspace(0,1e4,10);
% figure; hold on
% scatter(abs(CASr.Z(in)),log(BG.Z(in)))
% plot(d,p(2)+p(1).*d)

D=CASr;
D.Z=(abs(CASr.Z).*p(1)+p(2));

in=find(and(and(dist>40000,dist<=80000),~isnan(CASr.Z)));
bc=bin(Out_G.trench_y(in),Out_G.trench_x(in),nb);

bss=bin(Out_G.trench_y(in),FSr.Z(in)-CASr.Z(in),nb);
bg=bin(Out_G.trench_y(in),BG.Z(in),nb);
bss_m=bin(Out_G.trench_y(in),D.Z(in),nb);

[lat, ~] = projinv(projcrs(26910), bc(:,2)', bc(:,1)');

figure
plot(lat,bss(:,2)); hold on
yyaxis('right')
plot(lat,(bg(:,2))); hold on
plot(lat,exp(bss_m(:,2))); hold on

% plot(lat,-bss_m(:,2));
%% Slab 2.0
S2 = GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_UTM10N.tif');
S2=crop(S2,domainx,domainy);
Out_S=Slab_Project(S2,TS_Poly,90-52,40e-3);


%%

K2=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_K2.tif');
K1=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_K1.tif');
KG=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_KG.tif');
K2=crop(K2,domainx,domainy);
K1=crop(K1,domainx,domainy);
KG=crop(KG,domainx,domainy);

BGD=BGr;
BGD.Z=(log(BGr.Z)-p(2))./p(1);

dist=Out_S.dist;
nb=100;

in=find(and(and(dist>2e4,dist<=4e4),~isnan(K2.Z)));

b1=bin(Out_S.trench_y(in),K1.Z(in),nb);
b2=bin(Out_S.trench_y(in),K2.Z(in),nb);
bs2=bin(Out_S.trench_y(in),S2.Z(in),nb);
bgs2=bin(Out_S.trench_y(in),BGD.Z(in),nb);
bg=bin(Out_S.trench_y(in),KG.Z(in),nb);
bc=bin(Out_S.trench_y(in),Out_S.trench_x(in),nb);

[lat, ~] = projinv(projcrs(26910), bc(:,2)', bc(:,1)');

k1=b1(:,2);
k2=b2(:,2);
s2=bs2(:,2);
bgs2=bgs2(:,2);
nin=find(k2~=0);

figure; hold on
% plot(lat(nin),k1(nin))
% plot(lat(nin),k2(nin))
plot(lat(nin),s2(nin))
plot(lat(nin),-bgs2(nin))

%%
BGr=resample(BG,S2);
in=find(and(and(~isnan(BGr.Z),~isnan(S2.Z)),abs(S2.Z)<1.5e4));
p=polyfit(abs(S2.Z(in)),log(BGr.Z(in)),1);

d=linspace(0,1.5e4,10);

figure; hold on
scatter(abs(S2.Z(in)),log(BGr.Z(in)))
plot(d,d.*p(1)+p(2))