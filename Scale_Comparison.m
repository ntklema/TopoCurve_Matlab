function Scale_Comparison(Stack,LogArea,BS,nb)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
[~,~,n]=size(Stack.K1);
smwin=10;
figure; hold on
cmap=crameri('romaO',n);
in=find(BS.Z==1);
for i=1:n
    subplot(1,3,1); hold on
    Sl=Stack.Sl(:,:,i);
    b=bin(LogArea(in),Sl(in),nb);
    plot(10.^b(:,1),smoothdata(b(:,2),'movmean',smwin),'Color',cmap(i,:),'LineWidth',1);
    set(gca,'xscale','log')
    title('Sl')
    xlim([min(10.^b(:,1)) max(10.^b(:,1))])

    subplot(1,3,2); hold on
    KM=Stack.KM(:,:,i);
    b=bin(LogArea(in),KM(in),nb);
    plot(b(:,1),smoothdata(b(:,2),'movmean',smwin),'Color',cmap(i,:),'LineWidth',1);
    yline(0)
    title('KM')
    xlim([min(b(:,1)) max(b(:,1))])
    
    subplot(1,3,3); hold on
    KG=Stack.KG(:,:,i);
    b=bin(LogArea(in),KG(in),nb);
    plot(b(:,1),smoothdata(b(:,2),'movmean',smwin),'Color',cmap(i,:),'LineWidth',1);
    yline(0)
    title('KG')
    xlim([min(b(:,1)) max(b(:,1))])
end

end