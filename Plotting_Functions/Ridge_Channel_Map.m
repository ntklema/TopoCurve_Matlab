function Ridge_Channel_Map(U,BS,Ridge,Channel,Junctions)
cmap=[0.4,0,0;
    1,0,0;
    1,1,1;
    0,1,1;
    0,0,1];
J=load(Junctions);
R=shaperead(Ridge);
C=shaperead(Channel);
DEM=U.DEMGRIDobj;
x = [425500 432500]; 
y = [4834400 4838500];
S=DEM;
S.Z=U.SMAP.*BS.Z;
S.Z(S.Z==0)=nan;

DEM=crop(DEM,x,y);
S=crop(S,x,y);

figure
imageschs(DEM,S,'colormap',cmap); 
hold on

RP=struct;
CP=struct;

for i=1:numel(J.Junctions.RidgeX)
    RP(i).Geometry='Point';
    RP(i).X=J.Junctions.RidgeX(i);
    RP(i).Y=J.Junctions.RidgeY(i);

    CP(i).Geometry='Point';
    CP(i).X=J.Junctions.ChanX(i);
    CP(i).Y=J.Junctions.ChanY(i);
end

shapewrite(RP,'/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Shapefiles/Franklin_Ridge_Points.shp')
shapewrite(CP,'/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Shapefiles/Franklin_Channel_Points.shp')


end