DEM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Bathy_15as_UTM.tif');

CAS=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/CASIE_TOC_UTM.tif');

FS=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/CASIE_FS_UTM.tif');

BG=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Gravity/Bouguer_UTM10m2.tif');

% DEM=resample(DEM,BG,'bilinear');
% FS = resample(FS, DEM, 'bilinear');
% CAS = resample(CAS, DEM, 'bilinear');

% 
% FS=resample(FS,DEM,'bilinear');
% CAS=resample(CAS,DEM,'bilinear');
% S2 = GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_UTM10N.tif');
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
dist=Out.dist.Z;
cmap=(colormap(copper(5)));

time=Out.time.Z;
d=[0 2e4];

in=find(and(and(dist>d(1),dist<d(2)),~isnan(CAS.Z)));

b=bin(Out.trench_y(in),CAS.Z(in),200);
b2=bin(Out.trench_y(in),FS.Z(in),200);

figure
plot(b(:,1),b(:,8)); hold on
plot(b2(:,1),b2(:,8));

%% Time binning plot
Diff=CAS;
Diff.Z=FS.Z-CAS.Z;
Diff.Z(Diff.Z<0)=0;

T=Out.time.Z;
N=Out.trench_y;

nb=13;

tv=linspace(0,1.5e6,nb);
nv=linspace(4.6e6,5.5e6,nb);
t=[];
v=[];
n=[];
for i=1:numel(tv)-1
    t(i)=mean([tv(i),tv(i+1)]);
    for j=1:numel(nv)-1
        in=find(and(and(and(T>tv(i),T<=tv(i+1)),and(N>nv(j),N<=nv(j+1))),~isnan(Diff.Z)));
        n(j)=mean([nv(j),nv(j+1)]);

    
    v(i,j)=sum(Diff.Z(in))/numel(in);
    end
end
[lat,~]=projinv(projcrs(26910), ones(size(n)).*4e5, n);

figure
subplot(2,2,3)
imagesc(lat,t,v)

t_v=[];
r_nan=[];
l_v=[];
c_nan=[];

for i=1:nb-1
    data1=v(i,:);
    t_v(i)=sum(data1,'omitnan')/sum(~isnan(data1));
    r_nan(i)=sum(~isnan(data1));
    
    data2=v(:,i);
    l_v(i)=sum(data2,'omitnan')/sum(~isnan(data2));
    c_nan(i)=sum(~isnan(data2));
end

subplot(2,2,1)
plot(lat,l_v)
yyaxis('right')
plot(lat,c_nan)

subplot(2,2,4)
plot(t,t_v)
yyaxis('right')
plot(t,r_nan)
