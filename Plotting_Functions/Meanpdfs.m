function Meanpdfs(DEM,C,BS,thresh)
    DEM.Z=C.DEM.ZFilt;
    KM=C.CMAP.KM;
    KM=KM.*BS.Z;
    A=log10(C.Stream.A.Z);
    M=nan(size(KM));

    M(and(A<thresh,BS.Z==1))=0;
    M(and(A>thresh,BS.Z==1))=1;
    
    % Hillshade of ROI
    % figure
    % hillshade(DEM)
    % xlim([4.20100 4.26600].*1e5)
    % ylim([4.8398 4.8456].*1e6)
    
    % Mean curvature bimodal colormap
    % figure
    % DM=DEM;
    % DM.Z=M;
    % imagesc(DM);
    % colormap([1,0,0;0,0,1]);
    % xlim([4.20100 4.26600].*1e5)
    % ylim([4.8398 4.8456].*1e6)
    f1=figure; hold on
    f1.Units = 'inches';
    f1.Position=[3 3 5.5 9];
    %% PDF of slopes in each domain
    nb=1000;
    smwin=100;
    sedges = linspace(0,1.4,nb);
    
    in=find(and(A>thresh,BS.Z==1));
    [SNpdf,~]=histcounts(C.CMAP.Sl(in),sedges,'Normalization','probability');
    in=find(and(A<thresh,BS.Z==1));
    [SPpdf,~]=histcounts(C.CMAP.Sl(in),sedges,'Normalization','probability');
    
    Svec=nan(size(SNpdf));
    for i=1:numel(sedges)-1
        Svec(i)=mean([sedges(i),sedges(i+1)]);
    end
    t1 =  tiledlayout(3,1,'TileSpacing','Compact');
    nexttile(t1) % Plot Area pdf
    plot(Svec,smoothdata(SNpdf,'movmean',smwin),'b');
    hold on
    
    plot(Svec,smoothdata(SPpdf,'movmean',smwin),'r');
    ylabel('Probability')
    xlabel('Slope (m/m)')
    legend("A > "+10^thresh+" m^2","A < "+10^thresh+" m^2")
  

    %% PDF of Mean Curvatures in each domain
    nb = 1000;
    smwin=100;
    cedges=linspace(-0.015,0.015,nb);

    in=find(and(A>thresh,BS.Z==1));
    [KMPpdf,~]=histcounts(C.CMAP.KM(in),cedges,'Normalization','probability');

    in=find(and(A<thresh,BS.Z==1));

    [KMNpdf,sedges]=histcounts(C.CMAP.KM(in),cedges,'Normalization','probability');
    Svec=nan(size(KMNpdf));
    for i=1:numel(sedges)-1
        Svec(i)=mean([sedges(i),sedges(i+1)]);
    end

    nexttile(t1)
    plot(Svec,smoothdata(KMPpdf,'movmean',smwin),'b');
    hold on
    plot(Svec,smoothdata(KMNpdf,'movmean',smwin),'r');
    ylabel('Probability')
    xlabel('Mean curvature (1/m)')
    legend("A > "+10^thresh+" m^2","A < "+10^thresh+" m^2")

    %% PDF of Gaussian Curvatures in each domain
    nb = 1000;
    smwin=100;

    gedges=linspace(-1e-4,1e-4,nb);
    in=find(and(A>thresh,BS.Z==1));
    [KGNpdf,~]=histcounts(C.CMAP.KG(in),gedges,'Normalization','probability');
    
    in=find(and(A<thresh,BS.Z==1));
    [KGPpdf,~]=histcounts(C.CMAP.KG(in),gedges,'Normalization','probability');
    
    KGPvec=nan(size(KGNpdf));
    for i=1:numel(gedges)-1
        KGvec(i)=mean([gedges(i),gedges(i+1)]);
    end

    nexttile(t1)
    plot(KGvec,smoothdata(KGNpdf,'movmean',smwin),'b');
    hold on    
    plot(KGvec,smoothdata(KGPpdf,'movmean',smwin),'r');
    ylabel('Probability')
    xlabel('Gaussian curvature (1/m^2)')
    legend("A > "+10^thresh+" m^2","A < "+10^thresh+" m^2")
end