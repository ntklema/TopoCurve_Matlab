function Slope_Compare(C,BS,bn)
DEM=C.DEMGRIDobj;
DEM.Z=C.DEM.ZFilt;
A=log10(C.Stream.A.Z);
in=find(BS.Z==1);

[Gx,Gy]=gradient(DEM.Z,C.DEM.dx,C.DEM.dy);
G=sqrt(Gx.^2+Gy.^2);
G8=gradient8(DEM);
Sl=C.CMAP.Sl;

bg=bin(A(in),G(in),bn);
bg8=bin(A(in),G8.Z(in),bn);
bsl=bin(A(in),Sl(in),bn);

figure
bgs=scatter(10.^bg(:,1),bg(:,2),'s'); hold on
bg8s=scatter(10.^bg8(:,1),bg8(:,2),'o'); 
bsls=scatter(10.^bsl(:,1),bsl(:,2),'d'); 

set(gca,'xscale','log')



end