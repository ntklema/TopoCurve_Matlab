function LPError(CurveObj)
    C=CurveObj;
    LP=del2(C.DEM.ZFilt,C.DEM.dx,C.DEM.dy);
    LPE=abs(LP./2-C.CMAP.KM);

    % bin by slope and Curvature
    bc=bin(C.CMAP.KM,LPE,20);
    bs=bin(C.CMAP.Sl,LPE,20);

    figure
    subplot(1,2,1)
    scatter(bc(:,1),bc(:,2))
    scatter(bs(:,1),bs(:,2))

    
end