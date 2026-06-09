function Error_Covariance_Plot(C,BS,nx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

km=linspace(0,0.008,nx);
sl=linspace(0,1.5,nx);

dk=km(2)-km(1);
ds=sl(2)-sl(1);


E=abs(C.CMAP.KM-C.CMAP.LP./2).*BS.Z;
PE=abs(C.CMAP.KM-C.CMAP.LP./2)./abs(C.CMAP.KM).*100.*BS.Z;

KM=abs(C.CMAP.KM.*BS.Z);
% KM(C.CMAP.KG<0)=0;
Sl=C.CMAP.Sl.*BS.Z;

for i=1:nx
    S=Sphere_Curvature(1000,1/km(i));
    US=CurveObj(S);
    US=CurveCalc(US,0);
    ES=abs(US.CMAP.KM-US.CMAP.LP./2);

    for j=1:nx
        DM(i,j)=mean(E(and(and(KM>=km(i)-dk,KM<km(i)+dk),and(Sl>=sl(j)-ds,Sl<sl(j)+ds))));
        DS(i,j)=mean(ES(and(US.CMAP.Sl>=sl(j)-ds,US.CMAP.Sl<sl(j)+ds)));
        
    end
end

DM(DM<1e-4)=nan;
figure
subplot(2,2,2)
imagesc(sl,km,DM)
colorbar
set(gca,'ydir','normal','colorscale','linear')

subplot(2,2,3)
imagesc(sl,km,DS)
colorbar
set(gca,'ydir','normal','colorscale','linear')

subplot(2,2,1); hold on
plot(km,mean(DM,2,"omitmissing"))
plot(km,min(DM,[],2))
plot(km,max(DM,[],2))
% plot(km,mean(DM,2,"omitmissing")+2.*std(DM,0,2,"omitmissing"))
% plot(km,mean(DM,2,"omitmissing")-2.*std(DM,0,2,"omitmissing"))
plot(km,mean(DS,2,"omitmissing"))
set(gca,'yscale','log')

subplot(2,2,4); hold on
plot(sl,mean(DM,1,"omitmissing"))
plot(sl,min(DM,[],1))
plot(sl,max(DM,[],1))
% plot(sl,mean(DM,1,"omitmissing")+2.*std(DM,0,1,"omitmissing"))
% plot(sl,mean(DM,1,"omitmissing")-2.*std(DM,0,1,"omitmissing"))
plot(sl,mean(DS,1,"omitmissing"))
set(gca,'yscale','log')

