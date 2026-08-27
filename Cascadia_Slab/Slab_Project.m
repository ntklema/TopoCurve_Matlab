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
X=X+abs(min(X(:)));
Y=Y+abs(min(Y(:)));
[r,c]=size(X);

for i=1:r
    for j=1:c
        R=sqrt((X-X(i,j)).^2+(Y-Y(i,j)).^2);
        TH=atand((Y-Y(i,j))./(X-X(i,j)));
        in=find(and(S.Z==1,abs(TH-angle)==min(abs(TH-angle))));

        if ~isempty(in)
            D(i,j) = min(R(in));
        else
            D(i,j) = nan; 
    

    end
end

Out=D;
end