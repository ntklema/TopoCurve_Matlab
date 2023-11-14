dem_path='/Users/nathanielklema/OneDrive - Fort Lewis College/Research_Projects/Curvature/Data/DEMS/Umpqua_Lidar_1.tif';

DEM=GRIDobj(dem_path);
M=createmask(DEM);


%%
DEM=crop(DEM,M,nan);

%%

FD=FLOWobj(DEM,'preprocess','carve');
A=flowacc(FD).*DEM.cellsize.^2;
S=STREAMobj(FD,'minarea',100);
ST=trunk(klargestconncomps(S));

DEM10=resample(DEM,10);
FD10=FLOWobj(DEM10,'preprocess','carve');
A10=flowacc(FD10).*DEM10.cellsize.^2;
S10=STREAMobj(FD10,'minarea',100);
ST10=trunk(klargestconncomps(S10));

%%
DEM20=resample(DEM,20);
FD20=FLOWobj(DEM20,'preprocess','carve');
A20=flowacc(FD20).*DEM20.cellsize.^2;
S20=STREAMobj(FD20,'minarea',100);
ST20=trunk(klargestconncomps(S20));