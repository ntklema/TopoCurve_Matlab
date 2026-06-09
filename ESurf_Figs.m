% Code written by Nathaniel Klema to accompany 

% Klema, N., Karlstrom, L., and Roering, J.: Discrete differential geometry 
% of fluvial landscapes, EGUsphere [preprint], 
% https://doi.org/10.5194/egusphere-2025-4431, 2025

% All uses of this code should cite this manuscript.  Questions can be
% directed to ntklema@fortlewis.edu


% Load filepath of script.  DEM and plotting functions should be stored
% in a sub directory in the same folder.
filePath = matlab.desktop.editor.getActiveFilename;
addpath(genpath(fileparts(fileparts(filePath))));

% Load DEM as a GRIDobj object class of the TopoToolbox library
demname="Umpqua_10m_2";
DEM=GRIDobj(demname+".tif");

% % Uncomment the following lines to downsample or crop the DEM.
% DEM=resample(DEM,5); Downsample DEM
% DEM=crop(DEM,[4.19e5, 4.35e5],[4.834e6, 4.848e6]); Crop DEM

% Pass DEM into 'CurveObj' object class.  The DEM will be stored under C.DEM.
C=CurveObj(DEM);

% Load ROI Shapefile - this will be used later to select cells within study
% area
bs=shaperead(demname+"_Outline.shp");

%% Filter DEM and compute curvatures

% Assign filter cutoff wavelengths.  Order does not matter.  Wavelengths
% greater than the large number will be included in the filtered DEM.
% Wavelengths smaller than smaller number will not be included.  For 
% transitional wavelengths amplitudes will be supressed.
filter=[150 200];

% Low-pass filter DEM using filter defined above.  Default is to mirror the 
% DEM and apply a 2-D tukey window. 
C=SpecFilt(C,filter,'lowpass');

% Calculate Curvatures  
% - Second input is curvature threshold (k_t of Mynatt et al., 2007).  
% Principal curvatures below this value will be assigned a value of zero.  
% Setting this to zero will maintain raw curvature values.
C=CurveCalc(C,0);

% The shape class distribution is stored in the C.SMAP matrix where numbers
% correspond to shape classes via the following assignments:
%
% Perfect Saddles: -4
% Domes: -3
% Antiformal Saddles: -2
% Antiform: -1
% Planes: 0
% Synforms: 1
% Synformal Saddles: 2
% Basins: 3
%
% NOTE: If you do not assign a curvature threshold you will only extract
% Domes, Antiformal Saddles, Synformal Saddles, and Basins.

% Surface Geometry Metrics are stored in C.CMAP with the following labels:
%
% NX,NY,NZ: x, y, and z components of unit normal vector.
% A,E,G: Metric coefficients
% e,f,g: Curvature coefficients
% V1,V2: Horizontal vector orientations (dv/du) of principal curvature vectors.  
% KG: Gaussian curvature
% KM: Mean curvature
% K1: First principal curvature
% K2: Second principal curvature
% az: Azimuth angle of tangent plane
% el: Elevation angle of tangent plane
% Sl: Slope of tangent plane
% LP: Projected Laplacian curvature

clear demname filePath DEM filter
%% Flow Routing and binning by area
% Route flow over DEM using Dinf algorithm from TopoToolbox.  Second input is drainage area
% threshold for assigning a channel in the STREAMobj object class.  This
% value does not influence the area binning procedure.  For fastest
% calculation compile .mex function in accordance with TopoToolbox
% documentation.
C=RouteFlow(C,6000);

%% Bin Data by drainage area for area defined by bs polygon.  
% Generates struct (Area_Binned_Metrics) within the CurveObj object
% with the following fields:

% area - Mean of area bin
% KM - Mean curvature
% KG - Gaussian curvature
% K1 - 1st principal curvature
% K2 - 2nd principal curvature
% A - Area expansion factor (alpha in manuscript)
% Sl - Tangent slope
% Az - Azimuth angle
% P_a - Probability aof area value
% P_b - Probability of basin occurance
% P_d - Probability of dome occurance
% P_as - Probability of antiformal saddle occurance
% P_ss - Probability of synformal saddle occurance
% LP - Projected Laplacian curvature

nb = 100; % Define number of area bins

% Perform data binning.  C is the modified CurvObj, while BS is a GRIDobj
% with logical values indicating the area of interest.
[C,BS] = bindata(C,bs,nb);

% Generate plots of geometry metrics binned by area
Area_Dist(C) 

clear bs nb

%% Generate map of shape classes for a defined range of drainage areas

% Define range of drainage areas to plot.  We use log10 to calculate log
% drainage areas.
a_range=[0 log10(350)];

% Create raster with only pixels within defined area range included.  Set
% others to nan.
LogA=C.Stream.A;
LogA.Z=log10(LogA.Z);
D=C.DEMGRIDobj;
D.Z=C.DEM.ZFilt;
SM=C.SMAP;
SM(or(LogA.Z<=a_range(1),LogA.Z>a_range(2)))=nan;

% Select only grid cells within ROI
SM(BS.Z==0)=nan;

% Define colormap
cmap=[0.4,0,0;
    1,0,0;
    1,1,1;
    0,1,1;
    0,0,1];

% Plot map
figure
imageschs(C.DEMGRIDobj,SM,'colormap',cmap)

% colorbar
xlim([4.20100 4.26600].*1e5)
ylim([4.8398 4.8456].*1e6)

clear a_range D cmap LogA

%% Make Pie Plot of shape composition

X=[numel(find(SM==-3)) numel(find(SM==-2)) numel(find(SM==2)) numel(find(SM==3))]./numel(find(~isnan(SM)));

% Define colormap
cmap=[0.4,0,0;
    1,0,0;
    0,1,1;
    0,0,1];
figure
ax=axes;
pie(X,'%.3f%%');
colormap(ax,cmap)

clear SM X ax cmap


%% Calculate curvature along given channels

% Interactively select stream network to analyze
figure
DEM=C.DEMGRIDobj;
FD=FLOWobj(DEM,'preprocess','fill');
A=flowacc(FD).*FD.cellsize^2;
S=STREAMobj(FD,A>5700);
St=modify(S,'interactive','polyselect');

% Extract trunk stream
St=trunk(klargestconncomps(St));

% Write values sampled along trunk stream to struct SA with the following
% attritubes:

% SA.a - Upstream area (m^2)
% SA.x - Easting (m)
% SA.y - Northing (m)
% SA.d - Upstream distance (m)
% SA.z - Elevation (m)
% SA.k1 - 1st principal curvature
% SA.k2 - 2nd principal curvature
% SA.kg - Gaussian curvature
% SA.km - Mean curvature
% SA.sl - Tangent slope

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

% plot first an second principal curvatures vs upstream distanc
figure
plot(SA.d,SA.k1); hold on
plot(SA.d,SA.k2)
xlabel ('Upstream Distance (m)')
ylabel('Curvature (1/m)')
legend('k1','k2')

clear FD St S A DEM
%% Generate pdfs of geometry metrics above and below a given area threshold
thresh=740; % Area threshold
Meanpdfs(C.DEMGRIDobj,C,BS,log10(thresh))

clear thresh
%% Generate filter scale comparison figure
LogA=C.Stream.A;
LogA.Z=log10(LogA.Z);
% Scale_Space generates grids of geometry metrics for a range of filter
% cutoofs, which are stored in the Stack struct.
Stack=Scale_Space(C.DEMGRIDobj,[50 0],50,4);

% Scale_Comparison bins slope, KM, and KG for each of the filter cutoffs
% stored in Struct.  Depends on bin.m written by Taylor Perron (% Copyright
% (C) 2004-2010 Taylor Perron <perron@mit.edu>) which is stored within the
% repository.
Scale_Comparison(Stack,LogA.Z,BS,100);

clear LogA