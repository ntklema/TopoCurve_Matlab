function LPError(CurveObj)
    C=CurveObj;
    LP=del2(C.DEM.ZFilt,C.DEM.dx,C.DEM.dy);
    LPE=abs(LP./2-C.CMAP.KM)./abs(C.CMAP.KM);

    % bin by slope and Curvature
    bc=bin((C.CMAP.KM),LPE,20);
    bg=bin((C.CMAP.KG),LPE,20);
    bs=bin(C.CMAP.Sl,LPE,20);

    figure
    subplot(1,3,1)
    scatter(bc(:,1),bc(:,2))
    subplot(1,3,2)
    scatter(bg(:,1),bg(:,2))

    subplot(1,3,3)
    scatter(bs(:,1),bs(:,2)); hold on
    plot(linspace(0,1,100),abs(1-(1+linspace(0,1,100).^2).^(3/2)))
    % 
    % figure
    % histogram(LPE)

    
end