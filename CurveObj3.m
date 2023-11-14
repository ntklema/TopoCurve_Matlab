classdef CurveObj3
    % Inputs:
    %   DEM - GRIDobj containing elevation raster
    %   
    %   'location' - String with specific regional name.  This will be
    %           stored and referenced in comparitive plotting functions.
    %
    %   'clip' - Clips original raster by extent.  Input is a vector fo the
    %           form [x1,x2,y1,y3]



    properties
        Location
        DEMGRIDobj
        DEM
        SMAP
        SDist
        CMAP
        kt
        Filter
        Stream
        
    end
    
    methods
        function obj = CurveObj3(DEM,varargin)
            set(0,'defaultfigurecolor',[1 1 1],'defaultfigureposition',[400 250 900 750])
           if isa(DEM,'GRIDobj')
               obj.DEM=struct('FullExtent',[],'X',[],'Y',[],'Z',[],'georef',[],'dx',[],'dy',[],'DTPlane',[],'Zmn',[],'ZFilt',[],'ZDiff',[]);
               
               X=DEM.georef.SpatialRef.XWorldLimits;
               Y=DEM.georef.SpatialRef.YWorldLimits;
               georef=[X,Y];
               

               p = inputParser;
               p.FunctionName = 'CurveObj3';
               addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));               
               addParameter(p,'clip',georef)
               addParameter(p,'location',[])
               parse(p,DEM,varargin{:});
               
               if ~isempty(p.Results.location)
                    obj.Location=p.Results.location;
               end
               obj.DEMGRIDobj=DEM;
               [obj.DEM.FullExtent.X, obj.DEM.FullExtent.Y]=getcoordinates(DEM);
               obj.DEM.FullExtent.Z=DEM.Z;
               obj.DEM.FullExtent.georef=georef;
               
               if isequal(p.Results.clip,georef)
                    obj.DEM.X=obj.DEM.FullExtent.X;
                    obj.DEM.Y=obj.DEM.FullExtent.Y;
                    obj.DEM.Z=obj.DEM.FullExtent.Z;
                    obj.DEM.georef=georef;

               else 
                   if p.Results.clip(1)<georef(1) || p.Results.clip(2)>georef(2) || p.Results.clip(3)<georef(3) || p.Results.clip(4)>georef(4)
                       error('User defined clipping range is outside spatial extent of raster')

                   elseif p.Results.clip(1)>p.Results.clip(2) || p.Results.clip(3)>p.Results.clip(4)
                        error('Clipping Extents must be of form [x1,x2,y1,y2]')
                   end

                    [~,xin1]=min(abs(p.Results.clip(1)-obj.DEM.FullExtent.X));
                    [~,xin2]=min(abs(p.Results.clip(2)-obj.DEM.FullExtent.X));
                    [~,yin1]=min(abs(p.Results.clip(3)-obj.DEM.FullExtent.Y));
                    [~,yin2]=min(abs(p.Results.clip(4)-obj.DEM.FullExtent.Y));
                    
                    obj.DEM.X=obj.DEM.FullExtent.X(xin1:xin2);
                    obj.DEM.Y=obj.DEM.FullExtent.Y(yin2:yin1);
                    obj.DEM.Z=obj.DEM.FullExtent.Z(yin2:yin1,xin1:xin2);
                    obj.DEM.georef=p.Results.clip;
               end
           obj.DEM.dx=DEM.cellsize;
           obj.DEM.dy=DEM.cellsize;
           
           elseif isa(DEM,'struct')
               obj.DEM=struct('FullExtent',[],'X',[],'Y',[],'Z',[],'georef',[],'dx',[],'dy',[],'DTPlane',[],'Zmn',[],'ZFilt',[],'ZDiff',[]);
               p = inputParser;
               p.FunctionName = 'CurveObj3';
               addRequired(p,'DEM',@(x) isa(x,'struct'));               
%                addParameter(p,'clip',georef)
               addParameter(p,'location',[])
               parse(p,DEM,varargin{:});
               
               obj.DEM.X=DEM.X;
               obj.DEM.Y=DEM.Y;
               obj.DEM.Z=DEM.Z;
               obj.DEM.dx=DEM.X(2)-DEM.X(1);
               obj.DEM.dy=DEM.Y(2)-DEM.Y(1);
               obj.Location=p.Results.location;



           else
              error('Error: Invalid Data Type')
           end

        end

        function obj=SpecFilt(obj,Filter,filtertype)
                    % Use spectral filtering to smooth data
                    % Inputs:
                        % obj - CurveObj to be operated on
                        %
                        % Filter - two componant vector that defines transition
                        %       wavelengths in the filter.  In 'lowpass'
                        %       mode the second number is the wavelength
                        %       at which the filter increases significantly 
                        %       above 0 while the first number is where it
                        %       reaches 1.
                        %
                        % filtertype - either 'lowpass' or 'highpass'.

  
                        
                    Zo=obj.DEM.Z; dx=obj.DEM.dx; dy=obj.DEM.dy;
                    
                    % detrend original dem
                    [ny,nx]=size(Zo);
                    [X,Y]=meshgrid(1:nx,1:ny);
                    [centroid, cosines] = lsplane([X(:) Y(:) Zo(:)]);

                    % for a plane with equation z = ax + by + c
                    a = -cosines(1)/cosines(3);
                    b = -cosines(2)/cosines(3);
                    c = centroid(3) + ((cosines(1)*centroid(1) + cosines(2)*centroid(2))/cosines(3));
                    Z = Zo - (a*X + b*Y + c);
                    plane=Zo-Z;
                    obj.DEM.DTPlane=plane;
                    clear nx ny
                    
                    [r,c]=size(Z);

                    % mirror DEM
                    Zm=[rot90(Z,2),flipud(Z),rot90(Z,2);
                        fliplr(Z),Z,fliplr(Z);
                        rot90(Z,2),flipud(Z),rot90(Z,2)];
                    [rm,cm]=size(Zm);
                    obj.DEM.Zm=Zm;
                    zw=window2(rm,cm,@tukeywin); %Window using flat top cosine taper
                    obj.DEM.Zmw=Zm.*zw; 
                    
                    zmw=obj.DEM.Zmw;
                    [ny, nx]=size(zmw); % number of rows and columns in mirrored dem
%                     [Pmat, fmat] = fft2D(Zmw, dx, dy, pad, window);
                    Lx=2^nextpow2(nx);  Ly=2^nextpow2(ny); % pad data to nearest power of 2
                    
                    ZMW = fftshift(fft2(zmw,Ly,Lx));
                    
                    dkx = 1/(dx*Lx); dky = 1/(dy*Ly); %Define wavenumber increments
                    ZMW(Ly/2 + 1, Lx/2 + 1)=0; %Zero out DC component to build radial frequencies
                    
                    % Make radial wavenumber matrix
                    xc = Lx/2+1; yc = Ly/2+1; % matrix indices of zero wavenumber
                    [cols, rows] = meshgrid(1:Lx,1:Ly); % matrices of column and row indices
                    km = sqrt((dky*(rows-yc)).^2 + (dkx*(cols-xc)).^2); % matrix of radial wavenumbers
                   
%                     ZMW = ZMW .* conj(ZMW) / (Lx * Ly * zw);
                    obj.DEM.ZP=ZMW;
                    obj.DEM.km=km;
                
                    switch filtertype

                        case 'lowpass'
                            kfilt=1./Filter;
                        
                            sigma=abs(kfilt(2)-kfilt(1))/3;
                            F=exp(-(km-kfilt(1)).^2/(2*sigma^2));
                            F=F/max(F(:));
                            F(km<kfilt(1))=1;

                        case 'highpass'
                            kfilt=1./Filter;
                            sigma=abs(kfilt(2)-kfilt(1))/3;
                            F=exp(-(km-kfilt(2)).^2/(2*sigma^2));
                            F(km>=kfilt(2))=1;
                    end

                
                    ZMWF = real(ifft2(ifftshift(ZMW.*F)));
                    obj.DEM.ZFilt=ZMWF(r+1:ny-r,c+1:nx-c)+plane;
                    obj.DEM.ZDiff=obj.DEM.Z-obj.DEM.ZFilt;
                    obj.Filter=Filter;
                end
        
        function obj = CurveCalc(obj,kt)
            %  Inputs:
                    % obj - CurveObj that has been filtered using the SpecFilt
                    % command
                    %
                    % kt - curvture thrshold below which curvature is set to zero.
            obj.kt=kt;
            if isempty(obj.DEM.ZFilt)
                SZ=obj.DEM.Z;
                obj.DEM.ZFilt=SZ;
                obj.Filter=[nan,nan];
            else
                SZ=obj.DEM.ZFilt;
            end
            
            [SX,SY]=meshgrid(obj.DEM.X,obj.DEM.Y);
            du=SX(1,2)-SX(1,1); dv=SY(2,1)-SY(1,1);
            [m,n]=size(SZ); % Number of rows and cols for grid on parameter plane.
            
            SXU=ones(size(SX)); SXV=zeros(size(SX)); % partial derivatives of SX
            SYU=zeros(size(SY)); SYV=ones(size(SY)); % partial derivatives of SY
            [SZU,SZV]=gradient(SZ,du,dv); % partial derivatives of SZ

            % Construct two 3D arrays, each of size (m,n,3).
            % Each page in 3rd dimension contains a component of the tangent vectors
            % SU and SV are ordered as SXU, SYU, SZU and SXV, SYV, SZV respectively.
            SU(:,:,1)=SXU; SU(:,:,2)=SYU; SU(:,:,3)=SZU;
            SV(:,:,1)=SXV; SV(:,:,2)=SYV; SV(:,:,3)=SZV;
            % Calculate the coefficients of the first fundamental form.
            E=dot(SU,SU,3); F=dot(SU,SV,3); G=dot(SV,SV,3); % Dividing by du and dv is an alteration of the mybnatt code

            % Calculate components of the unit normal vector, N, to the surface s.
            % Numerator is vector (cross) product of tangent vectors.
            % Denominator is absolute value of this cross product.
            CUV = cross(SU,SV,3);
            AC = sqrt(CUV(:,:,1).^2 + CUV(:,:,2).^2 + CUV(:,:,3).^2);
            NX = CUV(:,:,1)./AC; NY = CUV(:,:,2)./AC; NZ = CUV(:,:,3)./AC; 
            obj.CMAP.NX=NX; obj.CMAP.NY=NY; obj.CMAP.NZ=NZ;
            % Calculate the first partial derivatives of N w/re to u and v.
            % These are used to calculate the second fundamental form coefficients.
            [NXU,NXV]=gradient(NX,du,dv); [NYU,NYV]=gradient(NY,du,dv); 
            [NZU,NZV]=gradient(NZ,du,dv); 

            % Construct two 3D arays each of size (m,n,3).
            % Each page in 3rd dimension contains a component of the normal vectors.
            NU(:,:,1)=NXU; NU(:,:,2)=NYU; NU(:,:,3)=NZU;
            NV(:,:,1)=NXV; NV(:,:,2)=NYV; NV(:,:,3)=NZV;
            % Calculate the coefficients of the second fundamental form.
            e=-dot(NU,SU,3); f=-0.5*(dot(NU,SV,3)+dot(NV,SU,3)); g=-dot(NV,SV,3);

            % Preallocate arrays that are built inside the following loops.
            K1=zeros(size(SZ)); K1U=zeros(size(SZ)); K1V=zeros(size(SZ));
            K2=zeros(size(SZ)); K2U=zeros(size(SZ)); K2V=zeros(size(SZ));

            % The first and second fundamental form coefficients are used to
            % construct two 2 by 2 matrices, I and II, which then are used to
            % calculate the 2 by 2 shape operator matrix, SO, which then is used to
            % calculate the principal curvature magnitudes and apparent directions.
%             F=zeros(size(K1)); M=zeros(size(K1));
%             for i=1:m % Loop over rows of grid points on parameter plane.
%                 for j=1:n % Loop over columns of grid points on parameter plane.
%                     I=[E(i,j) F(i,j);F(i,j) G(i,j)]; % Matrix of first fund. form.
%                     II=[e(i,j) f(i,j);f(i,j) g(i,j)]; % Matrix of second fund. form.
%                     SO=I\II; % 2 by 2 shape operator matrix: same as inv(I)*II.
%                     A=max(max(isnan(SO))); % Test for nan elements in SO.
%                     if A==1 % if any are nan set all elements of KD and KM to nan
%                     KD=[NaN NaN;NaN NaN]; KM=[NaN NaN;NaN NaN];
%                     elseif A==0 % if not nans compute eigenvalues and vectors of SO
%                         [KD,KM]=eig(SO); % KD and KM are 2 by 2 matrices.
%                     % Cols of KD are eigenvectors; trace of KM are eigenvalues.
%                         if abs(KM(1,1))<abs(KM(2,2)) % Order principal curvatures so K1 > K2.
%                             KM=rot90(KM,2); KD=fliplr(KD);
%                         end
%                     % KD(:,1) and KM(1,1) are eigenvector and value for K1.
%                     % KD(:,2) and KM(2,2) are eigenvector and value for K2.
%                     
%                     end
% 
%                 % Extract principal curvature magnitudes, K1 and K2.
% %                   K1(i,j)=KM(1,1); K2(i,j)=KM(2,2); %As in Mynatt code
%                   K1(i,j)=-KM(1,1); K2(i,j)=-KM(2,2); % I think more consistent with convention?
%                 
% %                % High magnitude curvature as K1
% %                 PC=[-KM(1,1) -KM(2,2)];  
% %                 if PC(1)==0 && PC(2)==0
% %                     K1(i,j)=PC(1);
% %                     K2(i,j)=PC(2);
% %                 end
% %                 if PC(1)
% %                         K1(i,j)=PC(abs(PC)==max(abs(PC)));
% %                         K2(i,j)=PC(abs(PC)==min(abs(PC)));
% %                 end
% 
%                 % Extract components of apparent principal curvature directions.
%                 % Note: these are 2D vectors in the parameter plane which must be
%                 % projected onto planes tangent to s for actual directions.
%                     K1U(i,j)=KD(1,1); K1V(i,j)=KD(2,1);
%                     K2U(i,j)=KD(1,2); K2V(i,j)=KD(2,2);
% 
                    obj.CMAP.A=sqrt(E.*G-F.^2);
                    obj.CMAP.E=E; obj.CMAP.G=G; obj.CMAP.F=F;
                    obj.CMAP.e=e; obj.CMAP.f=f; obj.CMAP.g=g;
% 
% 
%                 end
%             end


%           Calculate the Principal Curvatures

%             if abs(F)<Fthresh
%                 
% 
%             else    
            a=E.*G-F.^2;
            b=-(g.*E-2.*f.*F+e.*G);
            c=e.*g-f.^2;
            K1=-(-b+sqrt(abs(b.^2-4.*a.*c)))./(2.*a);
            K2=-(-b-sqrt(abs(b.^2-4.*a.*c)))./(2.*a);

%           Calculate the principal directions
            al=F.*g-G.*f;
            bl=E.*g-G.*e;
            cl=E.*f-F.*e;
            obj.CMAP.V1=(-bl+sqrt(bl.^2-4.*al.*cl))./(2.*al);
            obj.CMAP.V2=(-bl-sqrt(bl.^2-4.*al.*cl))./(2.*al);
            obj.CMAP.al=al;
            obj.CMAP.bl=bl;
            obj.CMAP.cl=cl;
%         end
%             obj.CMAP.KG2=obj.CMAP.K12.*obj.CMAP.K22;
%              obj.CMAP.KM2=0.5.*(obj.CMAP.K12+obj.CMAP.K22);
            %K1X=K1U; K1Y=K1V; K1Z=K1U.*SZU+K1V.*SZV;
            %K2X=K2U; K2Y=K2V; K2Z=K2U.*SZU+K2V.*SZV;
            %K1M = sqrt(K1X.^2+K1Y.^2+K1Z.^2);
            %K1X = K1X./K1M; K1Y = K1Y./K1M; K1Z = K1Z./K1M;
           % K2M = sqrt(K2X.^2+K2Y.^2+K2Z.^2);
            %K2X = K2X./K2M; K2Y = K2Y./K2M; K2Z = K2Z./K2M;

            % Set K1 to 0 if |K1| < kt and set K2 to 0 if |K2| < kt. 
            K1((abs(K1)<=kt))=0; K2((abs(K2)<=kt))=0;
            % Calculate Gaussian curvature, KG, and mean curvature, KM.
            KG=K1.*K2; KM=0.5*(K1+K2);
            obj.SMAP=nan(size(KG));

            obj.SDist=cell(2,9);


            % Identify perfect saddle points.
            in=find(KG<0 & abs(KM)<=0); 
            obj.SMAP(in)=-4; 
            obj.SDist{1,1}={'Perfect Saddles'};
            obj.SDist{2,1}=numel(in)/numel(obj.SMAP);

            % Identify domal points.
            in=find(KG>0 & KM<0); 
            obj.SMAP(in)=3;
            obj.SDist{1,2}={'Peaks'};
            obj.SDist{2,2}=numel(in)/numel(obj.SMAP);


            % Identify antiformal saddle points.
            in=find(KG<0 & KM<0); 
            obj.SMAP(in)=2; 
            obj.SDist{1,3}={'Antiformal Saddles'};
            obj.SDist{2,3}=numel(in)/numel(obj.SMAP);
            
            % Identify cylindrical antiformal points.
            in=find(KG==0 & KM<0); 
            obj.SMAP(in)=1;
            obj.SDist{1,4}={'Antiforms'};
            obj.SDist{2,4}=numel(in)/numel(obj.SMAP);

            % Identify planar points.
            in=find(KG==0 & abs(KM)<=0); 
            obj.SMAP(in)=0;
            obj.SDist{1,5}={'Planes'};
            obj.SDist{2,5}=numel(in)/numel(obj.SMAP);

            % Identify cylindrical synformal points.
            in=find(KG==0 & KM>0); 
            obj.SMAP(in)=-1;
            obj.SDist{1,6}={'Synforms'};
            obj.SDist{2,6}=numel(in)/numel(obj.SMAP);

            % Identify synformal saddle points.
            in=find(KG<0 & KM>0); 
            obj.SMAP(in)=-2;
            obj.SDist{1,7}={'Synformal Saddles'};
            obj.SDist{2,7}=numel(in)/numel(obj.SMAP);

            
            % Identify Basin points.
            in=find(KG>0 & KM>0); 
            obj.SMAP(in)=-3;
            obj.SDist{1,8}={'Basins'};
            obj.SDist{2,8}=numel(in)/numel(obj.SMAP);

            % Store Gaussian, Mean, and principal curvatures in addition to
            % orientation vectors
            obj.CMAP.KG=KG; obj.CMAP.KM=KM;
            obj.CMAP.K1=K1; obj.CMAP.K2=K2;
%             obj.CMAP.K1U=K1U; obj.CMAP.K1V=K1V;
%             obj.CMAP.K2U=K2U; obj.CMAP.K2V=K2V;

            % Store slope vectors form filtered topographic map
            [obj.DEM.GU, obj.DEM.GV]=gradient(obj.DEM.ZFilt,obj.DEM.dx, obj.DEM.dy);

            % Store principal Curvatures sorted by magnitude
            obj.CMAP.K1M=K1; obj.CMAP.K2M=K2;
            obj.CMAP.K1M(abs(K1)<abs(K2))=K2(abs(K1)<abs(K2));


            
            
        end
                
        function TopoPlot(obj,varargin)


            p = inputParser;
            p.FunctionName = 'PlotTopo';
            addRequired(p,'obj',@(x) isa(x,'CurveObj2'));
            addParameter(p,'clipped',false)
            parse(p,obj,varargin{:});

            if true(p.Results.clipped)
                
                if ~isempty(obj.DEM.ZFilt)
                figure
                set(gcf,'units','inches','Position',[5,5,15,10]) 
                subplot(2,2,1);
                imagesc(obj.DEM.FullExtent.X,obj.DEM.FullExtent.Y,obj.DEM.FullExtent.Z)
                xlabel('Easting (km)')
                ylabel('Northing (km)')
                title('Full DEM')
                set(gca,'ydir','normal')
                a=colorbar;
                a.Label.String=('Elevation (m)');
                rectangle('position',[obj.DEM.georef(1),obj.DEM.georef(3),obj.DEM.georef(2)-obj.DEM.georef(1),obj.DEM.georef(4)-obj.DEM.georef(3)],'LineWidth',2)


                subplot(2,2,2)
                imagesc(obj.DEM.X,obj.DEM.Y,obj.DEM.ZFilt)
                xlabel('Easting (km)')
                ylabel('Northing (km)')
                title('Filtered Topography')
                set(gca,'ydir','normal')
                a=colorbar;
                a.Label.String=('Elevation (m)');

                subplot(2,2,3)
                imagesc(obj.DEM.X,obj.DEM.Y,obj.DEM.Z)
                xlabel('Easting (km)')
                ylabel('Northing (km)')
                title('Clipped DEM')
                set(gca,'ydir','normal')
                a=colorbar;
                a.Label.String=('Elevation (m)');

                subplot(2,2,4)
                imagesc(obj.DEM.X,obj.DEM.Y,obj.DEM.ZDiff)
                xlabel('Easting (km)')
                ylabel('Northing (km)')
                set(gca,'ydir','normal')
                title('Differenced Topography')
                a=colorbar;
                a.Label.String=('Elevation (m)');

                else 
                    figure
                    set(gcf,'units','inches','Position',[5,5,6,10])
                    subplot(2,1,1);
                    imagesc(obj.DEM.FullExtent.X,obj.DEM.FullExtent.Y,obj.DEM.FullExtent.Z)
                    xlabel('Easting (km)')
                    ylabel('Northing (km)')
                    title('Full DEM')
                    set(gca,'ydir','normal')
                    a=colorbar;
                    a.Label.String=('Elevation (m)');
                    rectangle('position',[obj.DEM.georef(1),obj.DEM.georef(3),obj.DEM.georef(2)-obj.DEM.georef(1),obj.DEM.georef(4)-obj.DEM.georef(3)],'LineWidth',2)


                    subplot(2,1,2)
                    imagesc(obj.DEM.X,obj.DEM.Y,obj.DEM.Z)
                    xlabel('Easting (km)')
                    ylabel('Northing (km)')
                    title('Filtered Topography')
                    set(gca,'ydir','normal')
                    a=colorbar;
                    a.Label.String=('Elevation (m)');
                end
            else
%                 X=obj.DEM.FullExtent.X;
%                 Y=obj.DEM.FullExtent.Y;
%                 Z=obj.DEM.FullExtent.Z;
                figure
                imagesc(obj.DEMGRIDobj)
                xlabel('Easting (km)'),ylabel('Northing (km)')
                set(gca,'ydir','normal')
                a=colorbar;
                a.Label.String=('Elevation (m)');
            end

         
        end
        
        function CurvePlot(obj,varargin)

            p = inputParser;
            p.FunctionName = 'PlotShapes';
            addRequired(p,'obj',@(x) isa(x,'CurveObj2'));
            addParameter(p,'topocompare',false)
            addParameter(p,'shapecompare',[])
            addParameter(p,'invariants',false)
            addParameter(p,'pc1',false)
            addParameter(p,'pc2',false)
            addParameter(p,'mean',false)
            addParameter(p,'gaussian',false)

            parse(p,obj,varargin{:});

            %Define shape colormap
             map=[0.6,0.1,0.2; 
                1.0,0.0,0.0;
                1.0,0.4,0.2;
                0,0,0;
                0.1,1.0,1.0;
                0.3,0.7,0.9;
                0.0,0.0,1.0];

%% topocompare plot function
            if true(p.Results.topocompare)
                figure
                r=subplot(2,2,1); 
                imagesc(obj.DEM.X,obj.DEM.Y,obj.DEM.Z);
                ylabel('Northing (m)'); title('Unfiltered DEM');
                set(gca,'ydir','normal')
                colormap(r,'gray')
                
                t=subplot(2,2,2);
                pie([obj.SDist{2,2},obj.SDist{2,3},obj.SDist{2,4},...
                    obj.SDist{2,5},obj.SDist{2,6},obj.SDist{2,7},obj.SDist{2,8}])
                colormap(t,map)
                colorbar
                colorbar('Ticks',[1+6/14,1+3*6/14, 1+5*6/14,1+7*6/14,1+9*6/14,1+11*6/14,1+13*6/14],...
                    'TickLabels',[obj.SDist{1,2:end}])%
                title("kt = "+obj.kt+", Filter Cutoff = "+obj.Filter(2)+"");

                ct=subplot(2,2,3);
                imagesc(obj.DEM.X,obj.DEM.Y,obj.DEM.ZFilt);
                xlabel('Easting (m)'); title('Filtered DEM');
                set(gca,'ydir','normal')
                colormap(ct,'gray')
                
                cm=subplot(2,2,4);
                imagesc(obj.DEM.X,obj.DEM.Y,obj.SMAP);
                xlabel('Easting (m)'); title('Shape Classifications');
                set(gca,'ydir','normal')
                colormap(cm,map)
                clim([-3 3])
            end

%% shapecompare plot function  
            if ~isempty(p.Results.shapecompare)
                if isa(p.Results.shapecompare,'CurveObj2')

                count=numel(p.Results.shapecompare);
            
                figure
                set(gcf,'units','inches','Position',[10,10,10,5])
                tiledlayout(2,count+1,'TileSpacing','loose');
                
                no=nexttile(1);
                imagesc(obj.DEM.X,obj.DEM.Y,obj.SMAP);
                xlabel('Easting (m)'); ylabel('Northing (m)')
                set(gca,'ydir','normal')
                colormap(no,map)
                clim([-3 3])
                if ~isempty(obj.Location)
                    title(""+obj.Location+"");
                else
                    title("DEM 1");
                end

                pn=nexttile(2+count);
                pie([obj.SDist{2,2},obj.SDist{2,3},obj.SDist{2,4},...
                    obj.SDist{2,5},obj.SDist{2,6},obj.SDist{2,7},obj.SDist{2,8}])
                colormap(pn,map)
%                 colorbar('Ticks',[1+6/14,1+3*6/14, 1+5*6/14,1+7*6/14,1+9*6/14,1+11*6/14,1+13*6/14],...
%                     'TickLabels',[obj.SDist{1,2:end}])
                title("kt = "+obj.kt+", Filter \lambda = ["+obj.Filter(1)+","+obj.Filter(2)+"] m")

                for i=1:count
                    o=p.Results.shapecompare(i);
                    co=nexttile(1+i);
                    imagesc(o.DEM.X,o.DEM.Y,o.SMAP);
                    xlabel('Easting (m)'); 
                    set(gca,'ydir','normal')
                    colormap(co,map)
                    clim([-3 3])
                    if ~isempty(o.Location)
                        title(""+o.Location+"");
                    else
                        title("DEM 2");
                    end
                    
    
                    p2=nexttile(count+2+i);
                    pie([o.SDist{2,2},o.SDist{2,3},o.SDist{2,4},...
                        o.SDist{2,5},o.SDist{2,6},o.SDist{2,7},o.SDist{2,8}])
                    colormap(p2,map)

                    title("kt = "+o.kt+", Filter \lambda = ["+obj.Filter(1)+","+obj.Filter(2)+"] m")
                end
                lgd=colorbar('Ticks',[1+6/14,1+3*6/14, 1+5*6/14,1+7*6/14,1+9*6/14,1+11*6/14,1+13*6/14],...
                'TickLabels',[o.SDist{1,2:end}]);
                lgd.Layout.Tile='east';
                else 
                    error('shapecompare function requires input of one or more CurveObjs directly following the argument')
            

                end
            end

             %% map invariants function
             if true(p.Results.invariants)
                 figure
                tiledlayout(2,2,'TileSpacing','loose');
                nexttile(1);
                imagesc(obj.DEM.X,obj.DEM.Y,obj.CMAP.KG)
                title('Gaussian Curvature')
                set(gca,'ydir','normal')
                colorbar
            
                nexttile(2);
                imagesc(obj.DEM.X,obj.DEM.Y,obj.CMAP.KM)
                title('Mean Curvature')
                set(gca,'ydir','normal')
                colorbar
            
                nexttile(3);
                imagesc(obj.DEM.X,obj.DEM.Y,obj.CMAP.K1)
                title('First Principal Curvature')
                set(gca,'ydir','normal')
                colorbar
            
                nexttile(4);
                imagesc(obj.DEM.X,obj.DEM.Y,obj.CMAP.K2)
                title('Second Principal Curvature')
                set(gca,'ydir','normal')
                colorbar
                     
             end
        end       
    
        function obj=RouteFlow(obj,thresh)

            D=obj.DEMGRIDobj;
             D.Z=obj.DEM.ZFilt;
%             DS=D;

            D=fillsinks(D);
            FD=FLOWobj(D,'Dinf');
%             FD=FLOWobj(D,'preprocess','carve');

            A=flowacc(FD).*D.cellsize^2;
            [obj.Stream.FU, obj.Stream.FV]=flowvec(FD);
            obj.Stream.S=STREAMobj(FD,A>thresh);
            obj.Stream.A=A;
            obj.DEM.Clipped=D;
            obj.Stream.SGRID=STREAMobj2GRIDobj(obj.Stream.S);
            obj.Stream.G=gradient8(D);
            obj.Stream.so=streamorder(obj.Stream.S,'strahler');
            CL=-del2(obj.DEM.ZFilt,obj.DEM.dx,obj.DEM.dx);
            obj.CMAP.LP=CL;
            C=curvature(D);
            obj.CMAP.TT=C.Z;
        end
    
        function obj=StreamExt(obj)
            ST=modify(obj.Stream.S,'interactive','polyselect');
            ST=trunk(klargestconncomps(ST));
% ST=obj.Stream.S;
            obj.Stream.ST=ST;
            order = ST.orderednanlist;
            dist  = ST.distance;
            zz = getnal(ST,obj.DEM.Clipped); % gives the GRIDobj Z values of grid sqaures crossed by the streamobj (referenced in the IXgrid)
            I     = ~isnan(order);
            d     = nan(size(order));
            d(I)  = dist(order(I)); % Gives distance along streamlength of each point
            z     = nan(size(order)); %  This makes all of the arrays the same size as the orderednanlist
            z(I)  = zz(order(I));
            x     = nan(size(order));
            x(I)  = ST.x(order(I));
            y     = nan(size(order));
            y(I)  = ST.y(order(I));
            aa    = getnal(ST,obj.Stream.A);
            a     = nan(size(order));
            a(I)  = aa(order(I));

            D=obj.DEM.Clipped;
            D.Z=obj.CMAP.KG;
            gg=getnal(ST,D);
            obj.Stream.g = nan(size(order));
            obj.Stream.g(I)=gg(order(I));
            obj.Stream.g=obj.Stream.g(~isnan(d));

            D.Z=obj.CMAP.KM;
            mm=getnal(ST,D);
            obj.Stream.m = nan(size(order));
            obj.Stream.m(I)=mm(order(I));
            obj.Stream.m=obj.Stream.m((~isnan(d)));

            D.Z=obj.CMAP.K1;
            k1=getnal(ST,D);
            obj.Stream.k1 = nan(size(order));
            obj.Stream.k1(I)=k1(order(I));
             obj.Stream.k1=obj.Stream.k1((~isnan(d)));
            obj.Stream.k1=obj.Stream.k1;

            D.Z=obj.CMAP.K2;
            k2=getnal(ST,D);
            obj.Stream.k2 = nan(size(order));
            obj.Stream.k2(I)=k2(order(I));
            obj.Stream.k2=obj.Stream.k2(~isnan(d));
            
            obj.Stream.x=x(~isnan(d)); obj.Stream.y=y(~isnan(d)); obj.Stream.z=z(~isnan(d));
            obj.Stream.a=a(~isnan(d));
            obj.Stream.d=d(~isnan(d)); 
           
        end
    
        function obj=AmpSpec(obj)
            
        end
    end
end

%% 
