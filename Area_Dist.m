function Area_Dist(BD,varargin)
    % Creates stacked plot of curvature metrics binned by area with side by
    % side shape class map.
    

    % Define shape class color map
    cmap=[0.4,0,0;
    1,0,0;
    1,1,1;
    1,1,1;
    0,1,1;
    0,0,1];

    p = inputParser;
    p.FunctionName = 'Area_Dist';
    addRequired(p,'BD',@(x) isa(x,'struct'));
    addParameter(p,'save',false)
    addParameter(p,'Area_Range',[])
    addParameter(p,'fig_name',{})
    parse(p,BD,varargin{:});
    
    fontl = 14; % fontsize of labels
    fontt = 12; % fontsize of labels

    smwin=10;
    f1=figure; hold on
    f1.Units = 'inches';
    f1.Position=[3 3 5.5 9];
    
    t1 =  tiledlayout(3,1,'TileSpacing','Compact');
    nexttile(t1) % Plot Area pdf
    area(10.^BD.area,smoothdata(BD.P_a.*100,'gaussian',smwin)); hold on
    ylabel('Area PDF (%)','FontName','helvetica','fontsize',fontl,'Color','k')
    set(gca,'xscale','log','xticklabels',[],'xlim',[10.^2 10.^7.3],'ylim',[0 7],'FontSize',fontt)
    yyaxis('right')
    plot(10.^BD.area,smoothdata(BD.Sl,'gaussian',smwin),'Color','k','LineWidth',2,'LineStyle','--')
    set(gca,'xscale','log','xticklabels',[],'xlim',[10.^2 10.^7.3],'ylim',[0 0.8],'FontSize',fontt,'YColor','k')
    
    ylabel('Slope (m/m)','FontName','helvetica','fontsize',fontl,'Color','k')

    if ~isempty(p.Results.Area_Range)
        fill([1e2 1e2 p.Results.Area_Range(1) p.Results.Area_Range(1)],[0 8 8 0],'w','FaceAlpha',0.9,'Linestyle','none')
        fill([p.Results.Area_Range(2) p.Results.Area_Range(2) 10^7.3 10^7.3],[0 8 8 0],'w','FaceAlpha',0.9,'Linestyle','none')
        fill([p.Results.Area_Range(1) p.Results.Area_Range(1) p.Results.Area_Range(2) p.Results.Area_Range(2)],[0 8 8 0],...
            'w','FaceAlpha',0,'Linestyle','-','EdgeColor','k','LineWidth',2)
    else
    end

    
    nexttile(t1) % Mean and Gaussian Curvature plots
    yyaxis('left')
    kg=plot(10.^BD.area,smoothdata(BD.KG,'gaussian',smwin),'color',[0.2 0.4 0],'LInewidth',2); hold on
    xline(700)
    yline(0)
    set(gca,'ylim',[-1.5e-5 1.5e-5],'YColor',[0.2 0.4 0],'FontSize',fontt)
    ylabel('Gaussian Curv.','FontName','helvetica','fontsize',fontl,'Color','k')
    
    
    yyaxis('right')
    km = plot(10.^BD.area,smoothdata(BD.KM,'gaussian',smwin),'color',[0.4 0 0.6],'LInewidth',2,'LineStyle','--'); hold on
    lp = plot(10.^BD.area,smoothdata(BD.LP./2,'gaussian',smwin),'color','k','Linewidth',1,'LineStyle','-'); hold on
    set(gca,'xscale','log','xticklabels',[],'xlim',[10.^2 10.^7.3],'ylim',[-8e-3 8e-3],'YColor',[0.4 0 0.6],'FontSize',fontt)
    ylabel('Mean Curv.','FontName','helvetica','fontsize',fontl,'Color','k')
    legend([kg,km,lp],{'K_G','K_M','LP'},'location','northeast')

    if ~isempty(p.Results.Area_Range)
        fill([1e2 1e2 p.Results.Area_Range(1) p.Results.Area_Range(1)],[-6.5e-3 6.5e-3 6.5e-3 -6.5e-3],'w','FaceAlpha',0.9,'Linestyle','none')
        fill([p.Results.Area_Range(2) p.Results.Area_Range(2) 10^7.3 10^7.3],[-6.5e-3 6.5e-3 6.5e-3 -6.5e-3],'w','FaceAlpha',0.9,'Linestyle','none')
        fill([p.Results.Area_Range(1) p.Results.Area_Range(1) p.Results.Area_Range(2) p.Results.Area_Range(2)],[-6.5e-3 6.5e-3 6.5e-3 -6.5e-3],...
            'w','FaceAlpha',0,'Linestyle','-','EdgeColor','k','LineWidth',2)

    else
    end
    

    nexttile(t1) % Shape Class distribution plot
    fill([0 1e3],[1 1],'r'); hold on
    bp=BD.P_b.*BD.P_a./BD.P_a; bp(bp>1)=1;
    bd=BD.P_d.*BD.P_a./BD.P_a; bd(bd>1)=1;
    bss=BD.P_ss.*BD.P_a./BD.P_a; bss(bss>1)=1;
    bas=BD.P_as.*BD.P_a./BD.P_a; bas(bas>1)=1;
    b = plot(10.^BD.area,smoothdata(bp,'gaussian',smwin),'color',cmap(6,:),'LineWidth',2); hold on
    d = plot(10.^BD.area,smoothdata(bd,'gaussian',smwin),'color',cmap(1,:),'LineWidth',2); hold on
    ss = plot(10.^BD.area,smoothdata(bss,'gaussian',smwin),'color',cmap(5,:),'LineWidth',2); hold on
    as = plot(10.^BD.area,smoothdata(bas,'gaussian',smwin),'color',cmap(2,:),'LineWidth',2); hold on
    
    if ~isempty(p.Results.Area_Range)
        fill([1e2 1e2 p.Results.Area_Range(1) p.Results.Area_Range(1)],[0 1 1 0],'w','FaceAlpha',0.9,'Linestyle','none')
        fill([p.Results.Area_Range(2) p.Results.Area_Range(2) 10^7.3 10^7.3],[0 1 1 0],'w','FaceAlpha',0.9,'Linestyle','none')
        fill([p.Results.Area_Range(1) p.Results.Area_Range(1) p.Results.Area_Range(2) p.Results.Area_Range(2)],[0 1 1 0],...
            'w','FaceAlpha',0,'Linestyle','-','EdgeColor','k','LineWidth',2)
    else
    end


    set(gca,'xscale','log','xlim',[10.^2 10.^7.3],'FontSize',fontt,'ylim',[0 1])
    xlabel('Upstream Drainage Area (m^2)','FontName','helvetica','fontsize',fontl)
    ylabel('P(C|A)','FontName','helvetica','fontsize',fontl)
    legend([b,ss,as,d],{'B','SS','AS','D'},'location','southeast')

    if p.Results.save == 1
        if isempty(p.Results.fig_name)
            saveas(gcf,'Area_Dist.pdf')
        else
            saveas(gcf,strcat(p.Results.fig_name,'.pdf'))
        end
    end

    




end