DEM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Bathy_15as_UTM.tif');
S2 = GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_UTM10N.tif');

CR_Poly=shaperead('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Shapefiles/Coast_Range_2.shp');
TS_Poly=shaperead('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Shapefiles/Casie_Wedgefront_Swath.shp');
AP_Poly=shaperead('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Shapefiles/Offshore_Acc_Prism.shp');
A_Poly=shaperead('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Shapefiles/Arc_Broad_Outline_2.shp');
SH_Poly=shaperead('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Shapefiles/Shelf_Swath_20.shp');
KM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/KM_UTM10N.tiff');


K1=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_K1.tif');
K2=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab2_K2.tif');
Kc=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab_K_c.tif');
Kp=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Slab_K_p.tif');

K1_s=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Surface_K1.tif');
K2_s=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Surface_K2.tif');
KM_s=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Surface_KM.tif');
Kc_s=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Surface_k_c.tif');
Kp_s=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Surface_k_p.tif');

%%

pole=[-32.0 68.3 1.081];
DEMc=resample(DEM,2000);
OutK = Slab_Rotate(DEMc, AP_Poly, pole, Mask = ~isnan(DEMc.Z), TMax = 5e6);


%% Curvature of continental slope
figure

R(1) = raster_swath_stats(K1_s, AP_Poly, Out = OutK, Scale = 1);
m1=smoothdata(R(1).median./1e-3,'gaussian',smwin);


R(2) = raster_swath_stats(K2_s, AP_Poly, Out = OutK, Scale = 1);
m2=smoothdata(R(2).median./1e-3,'gaussian',smwin);

R(3) = raster_swath_stats(KM_s, AP_Poly, Out = OutK, Scale = 1);
m3=smoothdata(R(3).median./1e-3,'gaussian',smwin);
%
subplot(4,1,1); hold on
in=find(~isnan(m1));
x=R(1).lat(in);

fill([x, fliplr(x)], [m1(in) fliplr(m2(in))],[0.2,0.2,1],'FaceAlpha',0.2)
plot(R(1).lat,m1,'LineWidth',2,'Color',[0.2,0.2,1]);
plot(R(1).lat,m2,'LineWidth',2,'Color',[0.2,0.2,1]);
plot(R(1).lat,m3,'LineWidth',2,'Color','k','LineStyle','--');
yline(0)

xlabel('Trench latitude')
ylabel('Curvature (km^{-1})')
title('Continental shelf (K1 and K2)')

% Profile and contour of bathymetric surface
R(1) = raster_swath_stats(Kp_s, AP_Poly, Out = OutK, Scale = 1);
m1=smoothdata(R(1).median./1e-3,'gaussian',smwin);


R(2) = raster_swath_stats(Kc_s, AP_Poly, Out = OutK, Scale = 1);
m2=smoothdata(R(2).median./1e-3,'gaussian',smwin);

R(3) = raster_swath_stats(KM_s, AP_Poly, Out = OutK, Scale = 1);
m3=smoothdata(R(3).median./1e-3,'gaussian',smwin);
%
subplot(4,1,2); hold on
in=find(~isnan(m1));
x=R(1).lat(in);

fill([x, fliplr(x)], [m1(in) fliplr(m2(in))],[0.2,0.2,0.5],'FaceAlpha',0.2)
plot(R(1).lat,m1,'LineWidth',2,'Color',[0.2,0.2,0.5]);
plot(R(1).lat,m2,'LineWidth',2,'Color',[0.2,0.2,0.5]);
plot(R(1).lat,m3,'LineWidth',2,'Color','k','LineStyle','--');
yline(0)

xlabel('Trench latitude')
ylabel('Curvature (km^{-1})')
title('Continental shelf (profile and contour curvatures)')

%
% K1 and K2 of slap
R(1) = raster_swath_stats(K1, AP_Poly, Out = OutK, Scale = 1);
m1=smoothdata(R(1).median./1e-3,'gaussian',smwin);


R(2) = raster_swath_stats(K2, AP_Poly, Out = OutK, Scale = 1);
m2=smoothdata(R(2).median./1e-3,'gaussian',smwin);

R(3) = raster_swath_stats(KM, AP_Poly, Out = OutK, Scale = 1);
m3=smoothdata(R(3).median./1e3,'gaussian',smwin);
%
subplot(4,1,3); hold on
in=find(~isnan(m1));
x=R(1).lat(in);

fill([x, fliplr(x)], [m1(in) fliplr(m2(in))],[0.5,0.2,0.5],'FaceAlpha',0.2)
plot(R(1).lat,m1,'LineWidth',2,'Color',[0.5,0.2,0.5]);
plot(R(1).lat,m2,'LineWidth',2,'Color',[0.5,0.2,0.5]);
plot(R(1).lat,m3,'LineWidth',2,'Color','k','LineStyle','--');
yline(0)

xlabel('Trench latitude')
ylabel('Curvature (km^{-1})')
title('Slab (K1 and K2)')

% Slab profile and contour curvature
R(1) = raster_swath_stats(Kp, AP_Poly, Out = OutK, Scale = 1);
m1=smoothdata(R(1).median./1e-3,'gaussian',smwin);


R(2) = raster_swath_stats(Kc, AP_Poly, Out = OutK, Scale = 1);
m2=smoothdata(R(2).median./1e-3,'gaussian',smwin);

R(3) = raster_swath_stats(KM, AP_Poly, Out = OutK, Scale = 1);
m3=smoothdata(R(3).median./1e3,'gaussian',smwin);
%
subplot(4,1,4); hold on
in=find(~isnan(m1));
x=R(1).lat(in);

fill([x, fliplr(x)], [m1(in) fliplr(m2(in))],[0.5,0,0.2],'FaceAlpha',0.2)
plot(R(1).lat,m1,'LineWidth',2,'Color',[0.5,0,0.2]);
plot(R(1).lat,m2,'LineWidth',2,'Color',[0.5,0,0.2]);
plot(R(1).lat,m3,'LineWidth',2,'Color','k','LineStyle','--');
yline(0)

xlabel('Trench latitude')
ylabel('Curvature (km^{-1})')
title('Slab (profile and countour curvatures)')
