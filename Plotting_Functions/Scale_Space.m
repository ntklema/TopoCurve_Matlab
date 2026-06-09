function Stack=Scale_Space(DEM,InitialFilt,StepSize,num)
% InitialFilt - Lowest filter cutoff
% StepSize - Increment to increase lowpass filter cutoffs by


C=CurveObj(DEM);
Filt=InitialFilt;
Stack=struct;
time=tic;
for n=1:num
    C=SpecFilt(C,Filt,'lowpass');
    C=CurveCalc(C,0);
    
    Stack.K1(:,:,n)=C.CMAP.K1;
    Stack.K2(:,:,n)=C.CMAP.K2;
    Stack.KM(:,:,n)=C.CMAP.KM;
    Stack.KG(:,:,n)=C.CMAP.KG;
    Stack.Sl(:,:,n)=C.CMAP.Sl;
    Stack.A(:,:,n)=C.CMAP.A;
    
    Filt=Filt+StepSize;
    disp("Number of Iterations: "+n+"")
    disp("Elapsed Time: "+round(toc(time)./60,2)+" minutes")


end
    
end