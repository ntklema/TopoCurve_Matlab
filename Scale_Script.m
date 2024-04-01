dem_path='/Users/nathanielklema/OneDrive - Fort Lewis College/GIS_Repositories/Coast Range/DEMs/Lidar/Umpqua_North_Lidar_3.tif';
addpath('/Users/nathanielklema/OneDrive - Fort Lewis College/Research_Projects/Curvature/Matlab_Code/TopoCurve_Matlab');
DEM=GRIDobj(dem_path);


FD=FLOWobj(DEM,'preprocess','carve');
A=flowacc(FD).*DEM.cellsize^2;

nb=200;

for i=1:nb
    
end