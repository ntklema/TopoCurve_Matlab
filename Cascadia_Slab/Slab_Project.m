function [Out] = Slab_Project(grid,shape,angle,velocity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    grid GRIDobj
    shape struct
    angle 
    velocity
end

S=polygon2GRIDobj(grid,shape);
[X,Y]=getcoordinates(grid,'matrix');
Xo=X;
X=X+abs(min(X(:)));
[r,c]=size(X);
D=nan(size(grid.Z));
TY=nan(size(grid.Z));
TX=nan(size(grid.Z));
time=nan(size(grid.Z));
lat=nan(size(grid.Z));
lon=nan(size(grid.Z));


for i=1:r

    for j=1:c
        R=sqrt((X-X(i,j)).^2+(Y-Y(i,j)).^2);
        TH=atand((Y-Y(i,j))./(X-X(i,j)));
        in=find(and(S.Z==1,abs(TH-angle)==min(abs(TH-angle))));
        
        if ~isempty(in)
            in2=find(R(in)==min(R(in)));
            
            D(i,j) = R(in(in2));
            TY(i,j)=Y(in(in2));
            TX(i,j)=Xo(in(in2));

        else
            D(i,j) = nan; 
    end
    
    s=S.Z(i,:);
    ti=find(s==1,1);
    D(i,1:ti)=nan;
end

D(S.Z==1)=0;
Dist=grid;
Dist.Z=D;
Out.dist=Dist;
Out.time=Dist./velocity;
Out.S=S;
Out.trench_y=TY;
Out.trench_x=TX;
end