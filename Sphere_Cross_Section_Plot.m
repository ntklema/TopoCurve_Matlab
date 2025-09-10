function Sphere_Cross_Section_Plot(R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% SPhere cross section figure
fig = figure('units','inch','position',[5,5,11,8]);

for i=1
    S=Sphere_Curvature(1001,R);
    DEM=GRIDobj(S.X,S.Y,S.Z);
    G=gradient8(DEM);
    US=CurveObj(S);
    US=CurveCalc(US,0);
    ES=abs(1-US.CMAP.LP.*2);

    [Y,X]=meshgrid(S.X,S.Y);
    sl=sqrt(X.^2+Y.^2)./S.Z;
    r=sqrt(X.^2+Y.^2);
    sl(isnan(sl))=100;
    in=find(sl<90);
    PEG=(G.Z-sl)./abs(sl).*100;
    imagesc(PEG)
    bg=bin(r(in),abs(sl(in)),100);
    


    % hold on
    % plot(S.X./S.R,US.CMAP.LP(500,:).*2,'k','LineWidth',4,'LineStyle','-')
    % plot(S.X./S.R,US.CMAP.KM(500,:),'b','LineWidth',4,'LineStyle','-.')
    % plot(S.X./S.R,ES(500,:),'r','LineWidth',1,'LineStyle','--')
    % % plot(S.X,(1+(S.X./S.Z(500,:)).^2).^(3/2)-1)
    % 
    % set(gca,'yscale','linear')
    % ylim([0 2])
    % 
    % yyaxis('right')
    % plot(S.X./S.R,US.CMAP.Sl(500,:),'LineWidth',4,'LineStyle',':')
    % hold on
    % plot(bg(:,1),bg(:,2))
    % % plot(S.X,-G(500,:))
    % xlim([0 R])

end