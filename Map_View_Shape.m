function Map_View_Shape(C,DEM,BS,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
    p.FunctionName = 'Map_View_Shape';
    addRequired(p,'C',@(x) isa(x,'CurveObj'));
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'BS',@(x) isa(x,'GRIDobj'));
    addParameter(p,'save',false)
    addParameter(p,'loc',[])
    addParameter(p,'fig_name',{})
    parse(p,C,DEM,BS,varargin{:});


cmap=[0.4,0,0;
    1,0,0;
    1,1,1;
    1,1,1;
    0,1,1;
    0,0,1];

D=DEM;
D.Z=C.SMAP;
D.Z(BS.Z==0)=nan;

figure
imageschs(DEM,D,'colormap',cmap,'colorbar',0)

if ~isempty(p.Results.loc)
    xlim([p.Results.loc(1) p.Results.loc(1)+1e4 ])
    ylim([p.Results.loc(2) p.Results.loc(2)+1e4 ])
end
end