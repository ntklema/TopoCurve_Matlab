function SS = Synth_Surf(r,c,dx)
    x=linspace(0,c*dx,c);
    y=linspace(0,r*dx,r);

    xc=200.*sin(8*2*pi/max(x).*x);
    xc2=900.*sin(pi/max(x).*x);
    yc=700.*sin(2*pi/max(y).*y);
    yc2=100.*sin(16*pi/max(y).*y);

    for i=1:r
        for j=1:c
            XD(i,:)=xc;
            XD2(i,:)=xc2;
            YD(:,j)=yc;
            YD2(:,j)=yc2;
        end
    end
    SS=XD+XD2;
    % SS=XD+YD+XD2+YD2;



end