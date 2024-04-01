
addpath('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/TopoCurve_Matlab')
% Can use this to check dependencies of script
[functionList,packageList] = matlab.codetools.requiredFilesAndProducts('Curvature_SandBox.m');
 % compilemexfiles % will need to run this command if you have not run topotoolbox mex functions before
%%

CoastPath="Umpqua_10m_1.tif"; 
HighPath="Santiam_Pass_10m_1.tif";

CoastDEM=GRIDobj(CoastPath);
HighDEM=GRIDobj(HighPath);
clear CoastPath HighPath
%% Build out CurveObj 
% Start by putting the DEM into the CurveObj class.  This automatically
% puts it into a data format in which it is easy to perform spectral
% filtering and calculate curvatures.

%% 1) Define object using CurveObj command
    % Inputs
        % DEM - Elevation raster in GRIDobj format from topotoolbox
        % 
        % 'location'(optional) - Name to attribute to raster to be referenced in
        %       mapping functions.
        %
        % 'clip'(optional) - Region of interest defined in a vector of UTM 
        %       coordinates of the form [x1 x2 y1 y2].
DEM=CoastDEM;   
Location='Oregon Coast Range'; % location name
ClipRange=[4.22e5,4.22e5+5e3,4.84e6,4.84e6+5e3]; % ROI
Coast=CurveObj(DEM,'location',Location,'clip',ClipRange);  %Create obj

% Plot DEM to see clipped region in context of full raster
TopoPlot(Coast,'clipped',true)

clear DEM Location ClipRange 

%% 2) Low pass filter topography using SpecFilt command

        % Inputs:
            % obj - CurveObj to be operated on
            %
            % Filter - two componant vector that defines transition 
            %       wavelengths in the filter.  In 'lowpass'mode the second 
            %       number is the wavelengt at which the filter increases 
            %       significantly above 0 while the first number is where 
            %       it reaches 1.
            %
            % filtertype - either 'lowpass' or 'highpass'.

Filter=[190,200];
filtertype='lowpass';
Coast=SpecFilt(Coast,Filter,filtertype); %  Apply filter

% Now if we use the TopoPlot command it will add panels that show the 
    %lowpass filtered and "differenced" (highpass filtered) topography.
TopoPlot(Coast,'clipped',true)

clear Filter filtertype 

% 3) Calculate invariant curvature metrics using CurveCalc command
%       Inputs:
            % obj - CurveObj to be operated on
            % kt - CUrvature threshold value below which all principal
            %   curatures will be set to 0.
kt=0.0000;
Coast=CurveCalc(Coast,kt); % Calculate curvatures
clear kt

% 4) Visualize curvature distributions using CurvePlot command

% Use 'topocompare' to look at the distribution of shape classes alongside
% topography
CurvePlot(Coast,'topocompare',true)

%% Use 'invariants' to see maps of each calculated curvature metric
CurvePlot(Coast,'invariants',true)

%% 5) use CurvePlot 'shapecompare' to do side by side comparisons
    % Input: CurvePlot(obj1,'shapecompare',[obj2,obj3,...])
        % obj1 - Object associated with a reference map
        % [obj2,obj3,...] - array containing any number of CurveObjs with 
            % filtering and calculation of curvature metrics already done.

%% Example 1 : Comparison of two curvature between two locations
% If I want to look at how curvature varies between two landscapes, or
% parts of the same landscape I can define to objects.  For this example I
% will compare the Oregon Coast Range outside Reedsport to the High
% Cascades around Santiam Pass using the same filter and curvature
% threshold.

% Create object holding Santiam Pass data
Location='Santiam Pass';
ClipRange=[5.826e5,5.826e5+5e3,4.912e6,4.912e6+5e3];
Filter=[290,300];
kt=0.00001;
High=CurveObj(HighDEM,'location',Location,'clip',ClipRange);
High=SpecFilt(High,Filter,'lowpass');
High=CurveCalc(High,kt);
TopoPlot(High,'clipped',true)% show region of high cascades selected
CurvePlot(High,'topocompare',true) % Show comparison between topography and shapes
clear Filter kt ClipRange Location
%%
% Create object holding Oregon Coast Range data
Location='Oregon Coast Range';
ClipRange=[4.22e5,4.22e5+5e3,4.84e6,4.84e6+5e3];
Filter=[190,200];
kt=0.00001;
Coast=CurveObj(CoastDEM,'location',Location,'clip',ClipRange); 
Coast=SpecFilt(Coast,Filter,'lowpass');
Coast=CurveCalc(Coast,kt);

% Use Curveplot to compare the two locations
CurvePlot(Coast,'shapecompare',High)

clear Filter kt ClipRange Location
%% Example 2: Comparison of same location with varied inputs
% If I want to see how varying an input (like the curvature threshold) I
% can again use the 'shapecompare' input.  This time I will try to see how
% the distribution of shape classes varies in the coast range for a range
% of kt values by creating a few objects, each with a different kt value
% then inputting them into the 'shapecompare' argument as an array.

% I will create 4 instances of the coast obj with shapes calculated for
% different values of kt

Coast1=CurveCalc(Coast,0);
Coast2=CurveCalc(Coast,0.000025);
Coast3=CurveCalc(Coast,0.00005);
Coast4=CurveCalc(Coast,0.0001);

CurvePlot(Coast1,'shapecompare',[Coast2,Coast3,Coast4])

%% Or could do the same mixing different landscapes

High1=CurveCalc(High,0);
High2=CurveCalc(High,0.00005);

CurvePlot(Coast1,'shapecompare',[High1,Coast3,High2])

%%  Synthetic surface generated by Leif Karlstrom that superimposes a quadratic basin 
%   with sinusoidal channels

Lx = 1e4; %length in X dir of domain in meters
Ly = 1e4;
H = 2e3; %desired scale height of basin (m)
Hch = 100; %desired amplitude of channels (m)

nx = 1000; %number of grid points
ny = 1000;

%spatial vectors 
X = linspace(-Lx/2,Lx/2,nx);
Y = linspace(-Ly/2,Ly/2,ny);

%generate basin, using function in bergbauer 2003
Zbasin = zeros(length(X),length(Y));

for i=1:length(X)
    Zbasin(i,:) = 3*X(i).^2 + Y.^2;
end

%scale height 
Zbasin = Zbasin./max(Zbasin(:)) * H;

%generate sinusoidal channels
Zchannel = zeros(length(X),length(Y)); 

for i=1:length(X)
    Zchannel(i,:) = Hch/2*(cos(X.*40/Lx)+1);
end

%composite surface is sum of both
Ztot = Zbasin + Zchannel;

figure(1)
subplot(1,3,1)
imagesc(X,Y,Zbasin)
title('basin')
colorbar
subplot(1,3,2)
imagesc(X,Y,Zchannel)
title('channels')
colorbar
subplot(1,3,3)
imagesc(X,Y,Ztot)
title('composite surface: basin+channels')
colorbar
%% Generate CurveObjs of each synthetic surface.  First we will plot all of 
%   them with kt=0 so see the full distribution of shape classes on the 
%   synthetic landscape.

DEM=struct;
DEM.X=X;
DEM.Y=Y;

DEM.Z=Zbasin;
SynthBasin=CurveObj(DEM,'location','Synthetic Basin');

DEM.Z=Zchannel;
SynthChannel=CurveObj(DEM,'location','Synthetic Channels');

DEM.Z=Ztot;
SynthTot=CurveObj(DEM,'location','Total Synthetic');

kt=0;
SynthBasin=CurveCalc(SynthBasin,kt);
SynthChannel=CurveCalc(SynthChannel,kt);
SynthTot=CurveCalc(SynthTot,kt);

CurvePlot(SynthBasin,'shapecompare',[SynthChannel, SynthTot])

%% WE can use the 'invariants' plot command for the composite landscape to
% see the magnitudes of principal curvatures associated with each
% "landform",  which gives a clue as to what curvature threshold will
% eliminate the background basin shape.  

CurvePlot(SynthTot,'invariants',1)

%%  Going of this we can set a curvature threshold of kt=4e-6 and see how 
%   it changes the shape distribution.  We see that with kt=4e-6 we
%   completely eliminate the effects of the basin including saddle shapes 
%   associated with cylindrical structures superimposed on the longer 
%   wavelength struncture fo the basin, so the composite map looks exactly 
%   like the channel map.

kt=4e-6;
STkt=CurveCalc(SynthTot,kt);
SCkt=CurveCalc(SynthChannel,kt);
SBkt=CurveCalc(SynthBasin,kt);

CurvePlot(SBkt,'shapecompare',[SCkt,STkt]);

%% While with a curvature threshold of 1e-6 the basin, which is slightly 
%   asymmetric appears synformal leading to a more complex shape
%   distribution in the composite landscape

kt=1e-6;
STkt=CurveCalc(SynthTot,kt);
SCkt=CurveCalc(SynthChannel,kt);
SBkt=CurveCalc(SynthBasin,kt);

CurvePlot(SBkt,'shapecompare',[SCkt,STkt]);

%% Try FLow routing to look for signal of channel heads
DEM=GRIDobj(Coast.DEM.X,Coast.DEM.Y,Coast.DEM.ZFilt);
FD=FLOWobj(DEM,'preprocess','carve');
a=flowacc(FD).*FD.cellsize^2;
SMAP=Coast.SMAP;
nx=20;
av=10.^(linspace(3,7.2,nx));
fv=linspace(100,600,nx);

bd=nan(nx,nx);
ssd=nan(nx,nx);
asd=nan(nx,nx);
dd=nan(nx,nx);
kg=zeros(nx,nx);
for j=1:nx

    Filter=[fv(j)-30,fv(j)];
    filtertype='lowpass';
    Coast=SpecFilt(Coast,Filter,filtertype); %  Apply filter
    Coast=CurveCalc(Coast,0); % Calculate curvatures
    SMAP=Coast.SMAP;
    KG=Coast.CMAP.KG; %./mean(Coast.CMAP.KG,'all');


    for i=2:numel(av)-1
        ain=find((a.Z>av(i-1)) & (a.Z<av(i+1)));
        shps=SMAP(ain);
        kg(i,j)=median(KG(ain));
        bd(i,j)=numel(find(shps==3))./numel(shps); 
        ssd(i,j)=numel(find(shps==2))./numel(shps); 
        asd(i,j)=numel(find(shps==-2))./numel(shps); 
        dd(i,j)=numel(find(shps==-3))./numel(shps);
    end
    j
end 
bd(isnan(bd))=0;
ssd(isnan(ssd))=0;
asd(isnan(asd))=0;
 dd(isnan(dd))=0;
%%
figure
imagesc(fv,av,dd)
colorbar
set(gca,'yscale','log'); hold on

%%
SMAP=Coast.SMAP;
in=find((a.Z>1e6) & (a.Z<5e6) & (SMAP==3));
kg=Coast.CMAP.KG(in);
figure
histogram(kg)