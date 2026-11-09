%% Slab_Swath_MORVEL.m
% Depth swaths of the subducted Juan de Fuca plate reconstructed by rotation
% about the MORVEL JF-NA Euler pole (DeMets et al., 2010, GJI 181).
%
% Requires in workspace: CAS (slab-top GRIDobj), FS (fault-surface GRIDobj),
% DEM (GRIDobj), TS_Poly (trench polygon mapstruct, UTM 10N).
% Requires on path: Slab_Rotate.m, bin.m, TopoToolbox, Mapping Toolbox.

%% ===================== Inputs =====================
% clearvars -except CAS FS DEM TS_Poly     % optional: clears stale/shadowing variables
DEM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/Bathy_15as_UTM.tif');
KM=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/KM_UTM10N.tiff');
CAS=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/CASIE_TOC_UTM.tif');

FS=GRIDobj('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Rasters/CASIE_FS_UTM.tif');

TS_Poly=shaperead('/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Cascadia_Slab/GIS/Shapefiles/Casie_Wedgefront_Swath.shp');

pole      = [-32.0 68.3 1.081];            % MORVEL JF relative to NA [lat lon deg/Myr]
proj      = projcrs(26910);                % UTM 10N / NAD83
binBy     = 'time';                        % 'time' or 'dist'
edgesT    = 0:2.5e5:2.25e6;                % time-bin edges (yr)
edgesD    = 0:1e4:9e4;                     % distance-bin edges (m)
nb        = 50;                            % latitude bins for swath statistics

faultFile = 'Surface_Faults.shp';            % wedge fault traces
nameField = 'Name';                        % fault label attribute ('' for none)
fcol      = [0.1 0.3 0.8];                 % fault color
Ls        = 3e3;                           % fault-trace smoothing length (m)
smoothAge = false;                         % smooth bin variable beneath faults
demMinCount = 5;                           % min cells per latitude bin for DEM medians

swathOut  = 'slab_swaths_MORVEL.shp';
flowOut   = 'flowlines_MORVEL.shp';
tickOut   = 'flowline_ticks_MORVEL.shp';
latList   = 42:50;                         % flowline start latitudes
TMaxFlow  = 3e6;  dtFlow = 1e4;            % flowline duration and vertex spacing (yr)

switch binBy                               % reconstruction must cover the last bin
    case 'time', TMax = edgesT(end) + 2.5e5;
    case 'dist', TMax = 1.1 * edgesD(end) / 0.028;   % slowest margin rate ~28-30 mm/yr
end

%% ===================== Reconstruction =====================
FS.Z(FS.Z < CAS.Z) = CAS.Z(FS.Z < CAS.Z);
valid = ~isnan(CAS.Z) & ~isnan(FS.Z);
Out   = Slab_Rotate(CAS, TS_Poly, pole, Mask = valid, TMax = TMax);
DEMr  = resample(DEM, CAS);

%% ===================== Bins and labels =====================
switch binBy
    case 'time', Vb = Out.time.Z;  edges = edgesT;
    case 'dist', Vb = Out.dist.Z;  edges = edgesD;
end
n = numel(edges) - 1;
masks = cell(1, n);  lab = cell(1, n);

for i = 1:n
    m = Vb > edges(i) & Vb <= edges(i+1) & valid;
    masks{i} = m;
    switch binBy
        case 'time'
            dr = prctile(Out.dist.Z(m), [5 95]) / 1e3;              % km
            tr = timeRange(edges(i), edges(i+1));
            if any(isnan(dr))
                lab{i} = tr;
            else
                lab{i} = sprintf('%.0f - %.0f km (%s)', dr(1), dr(2), tr);
            end
        case 'dist'
            tq = round(prctile(Out.time.Z(m), [5 95]) / 1e4) * 1e4;  % nearest 10 ky
            dl = sprintf('%g - %g km', edges(i)/1e3, edges(i+1)/1e3);
            if any(isnan(tq))
                lab{i} = dl;
            else
                lab{i} = sprintf('%s (%s)', dl, timeRange(tq(1), tq(2)));
            end
    end
end
cmap = copper(n);

%% ===================== Swath figure =====================
f1 = figure('units','inches','Position',[3 3 18 10],'color','w');
ax = subplot(2,3,[1 2 4 5]);  hold(ax,'on')

% Bathymetry along trench polygon
[X, Y] = getcoordinates(DEMr, 'matrix');
in  = Out.S.Z == 1;
b   = bin(Y(in), DEMr.Z(in), nb);
bc  = bin(Y(in), X(in), nb);
lat = projinv(proj, bc(:,2)', bc(:,1)');
hB  = plot(ax, lat, smoothdata(b(:,2)'/1e3,'gaussian',5), 'k', ...
           'LineWidth', 4, 'DisplayName', 'Bathymetric Depth');

% Slab swaths, plotted at trench-entry latitude
hF = gobjects(0);
bandLa = cell(1, n);  bandTop = cell(1, n);  bandBot = cell(1, n);
for i = 1:n
    in = masks{i};
    if ~any(in(:)), continue, end
    bc    = bin(Out.trench_y(in), Out.trench_x(in), nb);
    b     = bin(Out.trench_y(in), CAS.Z(in), nb);
    b2    = bin(Out.trench_y(in), FS.Z(in),  nb);
    lat   = projinv(proj, bc(:,2)', bc(:,1)');
    zSlab = b(:,8)';  zFS = b2(:,8)';
    k     = isfinite(zSlab) & isfinite(zFS) & isfinite(lat) & zSlab ~= 0;
    if ~any(k), continue, end
    zSlab = smoothdata(zSlab(k), 'gaussian', 5) / 1e3;
    zFS   = smoothdata(zFS(k),   'gaussian', 5) / 1e3;
    la    = lat(k);
    bandLa{i} = la;  bandTop{i} = zFS;  bandBot{i} = zSlab;
    hF(end+1) = fill(ax, [la fliplr(la)], [zSlab fliplr(zFS)], cmap(i,:), ...
                     'FaceAlpha', 0.8, 'DisplayName', lab{i}); %#ok<SAGROW>
    plot(ax, la, zSlab, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(ax, la, zFS,   'k', 'LineWidth', 1, 'HandleVisibility', 'off');
end
xlabel(ax, 'Latitude (^\circN)');  ylabel(ax, 'Depth (km)')

%% ===================== Median surface elevation per swath =====================
% Median DEM elevation of the cells in each swath, binned by trench-entry
% latitude (same projection as the slab swaths). Dashed, in the swath colour.
[Xg, Yg] = getcoordinates(CAS, 'matrix');          % grid vectors (also used by exports)
xv = Xg(1,:);  cs = CAS.cellsize;
Zdem = double(DEMr.Z) / 1e3;                       % km
zEmax = -Inf;
for i = 1:n
    in = masks{i} & ~isnan(Zdem) & ~isnan(Out.trench_x);
    if ~any(in(:)), continue, end
    la  = projinv(proj, Out.trench_x(in), Out.trench_y(in));
    ze  = Zdem(in);
    eL  = linspace(min(la), max(la), nb + 1);
    g   = discretize(la, eL);  ok = ~isnan(g);
    cnt = accumarray(g(ok), 1, [nb 1]);
    med = accumarray(g(ok), ze(ok), [nb 1], @median, NaN);
    med(cnt < demMinCount) = NaN;
    med = smoothdata(med, 'gaussian', 5, 'omitnan');
    med(cnt < demMinCount) = NaN;
    lc  = (eL(1:end-1) + eL(2:end)).' / 2;
    plot(ax, lc, med, '--', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off');
    zEmax = max(zEmax, max(med));
end
yline(ax, 0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');   % sea level
hE = plot(ax, NaN, NaN, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, ...
          'DisplayName', 'Median surface elevation');                  % legend proxy

% ---- Wedge fault section disabled; remove the %{ and %} lines to restore ----
%{
%% ===================== Wedge faults on the swaths =====================
% Faults are placed at the trench-entry latitude of the slab beneath them,
% at a vertical position interpolated smoothly between band centre-lines
% according to the age (or distance) of that slab. Vertical position is a
% band coordinate, not a fault depth.
Lf = shaperead(faultFile);

% Grid vectors for interp2 (must be increasing)
[Xg, Yg] = getcoordinates(CAS, 'matrix');
xv = Xg(1,:);  yv = Yg(:,1);  cs = CAS.cellsize;
if yv(2) < yv(1)
    yv = flipud(yv);  flipG = @flipud;
else
    flipG = @(G) G;
end
Gtx = flipG(Out.trench_x);
Gty = flipG(Out.trench_y);
Gv  = flipG(Vb);                                   % binning variable (time or dist)
F   = @(G, qx, qy) interp2(xv, yv, G, qx, qy, 'linear');   % NaN outside grid

hL = gobjects(0);
for j = 1:numel(Lf)
    lx = Lf(j).X;  ly = Lf(j).Y;
    keep = ~isnan(lx) & ~isnan(ly);
    lx = lx(keep);  ly = ly(keep);
    lx = lx(:).';   ly = ly(:).';
    if numel(lx) < 2, continue, end
    if all(abs(lx) <= 180) && all(abs(ly) <= 90)          % geographic -> UTM 10N
        [lx, ly] = projfwd(proj, ly, lx);
    end

    % Densify to half-cell spacing
    s = [0 cumsum(hypot(diff(lx), diff(ly)))];
    [s, iu] = unique(s);  lx = lx(iu);  ly = ly(iu);
    if numel(s) < 2, continue, end
    si = 0:cs/2:s(end);
    lx = interp1(s, lx, si);  ly = interp1(s, ly, si);

    % Slab beneath each point: trench-entry latitude and bin variable
    la = projinv(proj, F(Gtx, lx, ly), F(Gty, lx, ly));
    vv = F(Gv, lx, ly);
    if all(isnan(la) | isnan(vv)), continue, end

    % Smooth along the trace, per continuous segment
    good = ~isnan(la) & ~isnan(vv);
    edgeIdx = diff([false good false]);
    st = find(edgeIdx == 1);  en = find(edgeIdx == -1) - 1;
    for q = 1:numel(st)
        idx = st(q):en(q);
        if numel(idx) < 3, continue, end
        la(idx) = smoothdata(la(idx), 'gaussian', Ls, 'SamplePoints', si(idx));
        if smoothAge
            vv(idx) = smoothdata(vv(idx), 'gaussian', Ls, 'SamplePoints', si(idx));
        end
    end

    % Band centre-line depths at each trace latitude
    Zm = nan(n, numel(la));
    for kb = 1:n
        if isempty(bandLa{kb}), continue, end
        [lu, iu2] = unique(bandLa{kb});
        zmid = (bandTop{kb}(iu2) + bandBot{kb}(iu2)) / 2;
        Zm(kb,:) = interp1(lu, zmid, la, 'linear', NaN);
    end

    % Smooth curve through band centre-lines in fractional-bin coordinate
    u  = interp1(edges, 0:n, vv);                  % 0 at trench, k at k-th bin edge
    cK = (1:n) - 0.5;                              % band centres
    zb = nan(size(la));
    for p = find(~isnan(u))
        zc  = Zm(:, p);
        okc = ~isnan(zc);
        if nnz(okc) == 0, continue, end
        ck = cK(okc);
        if u(p) < min(ck) - 1 || u(p) > max(ck) + 1
            continue                               % don't extrapolate more than one bin
        end
        if nnz(okc) == 1
            zb(p) = zc(okc);
        else
            zb(p) = interp1(ck, zc(okc), u(p), 'pchip', 'extrap');
        end
    end
    if all(isnan(zb)), continue, end

    % Light smoothing of vertical position along the trace, per segment
    goodz = ~isnan(zb) & ~isnan(la);
    edgeIdx = diff([false goodz false]);
    st = find(edgeIdx == 1);  en = find(edgeIdx == -1) - 1;
    for q = 1:numel(st)
        idx = st(q):en(q);
        if numel(idx) < 3, continue, end
        zb(idx) = smoothdata(zb(idx), 'gaussian', Ls, 'SamplePoints', si(idx));
    end

    hL(end+1) = plot(ax, la, zb, '-', 'Color', fcol, 'LineWidth', 2, ...
                     'HandleVisibility', 'off'); %#ok<SAGROW>

    if ~isempty(nameField) && isfield(Lf, nameField)
        kEnd = find(~isnan(la) & ~isnan(zb), 1, 'last');
        text(ax, la(kEnd), zb(kEnd) + 0.3, string(Lf(j).(nameField)), ...
             'Color', fcol, 'FontSize', 11, 'HorizontalAlignment', 'center');
    end
end
%}

%% ===================== Annotations, limits, legend =====================
an_height = -4;
plot(ax, [45 48.3],   [an_height an_height],         'k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(ax, [45 45],     [an_height-0.2 an_height+0.2], 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(ax, [48.3 48.3], [an_height-0.2 an_height+0.2], 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
text(ax, 45.6, an_height+0.4, 'Cascadia Fold-Thrust Belt', 'FontSize', 14)
annotation('textarrow', [0.3 0.36],  [0.9 0.86],   'String', 'Astoria Fan  ', 'FontSize', 14);
annotation('textarrow', [0.53 0.51], [0.90 0.865], 'String', ' Nitinat Fan',  'FontSize', 14);
set(ax, 'FontSize', 14)
yTop = max(-1, ceil(2*zEmax)/2 + 0.5);            % room for onshore elevations
ylim(ax, [-17 yTop])
xlim(ax, [42 49.8])

legend(ax, [hB hF hE], 'Location', 'southwest');

%% ===================== Export map-view swath polygons =====================
xs = xv;  ys = Yg(:,1);                            % original (unflipped) row order
warning('off', 'MATLAB:polyshape:repairedBySimplify')
S = struct('Geometry', {}, 'BoundingBox', {}, 'X', {}, 'Y', {}, ...
           'Bin', {}, 'T_lo_ky', {}, 'T_hi_ky', {}, ...
           'D_lo_km', {}, 'D_hi_km', {}, 'Label', {});
for i = 1:n
    M = masks{i};
    if ~any(M(:)), continue, end
    P = mask2poly(M, xs, ys, cs);
    [bx, by] = boundary(P);
    switch binBy
        case 'time'
            tl = edges(i)/1e3;  th = edges(i+1)/1e3;
            dq = prctile(Out.dist.Z(M), [5 95]) / 1e3;  dl = dq(1);  dh = dq(2);
        case 'dist'
            dl = edges(i)/1e3;  dh = edges(i+1)/1e3;
            tq = prctile(Out.time.Z(M), [5 95]) / 1e3;  tl = tq(1);  th = tq(2);
    end
    S(end+1) = struct('Geometry', 'Polygon', ...
        'BoundingBox', [min(bx) min(by); max(bx) max(by)], ...
        'X', [bx(:).' NaN], 'Y', [by(:).' NaN], ...
        'Bin', i, 'T_lo_ky', tl, 'T_hi_ky', th, ...
        'D_lo_km', dl, 'D_hi_km', dh, 'Label', lab{i}); %#ok<SAGROW>
end
shapewrite(S, swathOut);
writePrj(swathOut, proj);

%% ===================== Flowlines from the trench =====================
Re = 6371e3;
e  = [cosd(pole(1))*cosd(pole(2)); cosd(pole(1))*sind(pole(2)); sind(pole(1))];
w  = deg2rad(pole(3)) / 1e6;
px = [TS_Poly.X];  py = [TS_Poly.Y];
tPath  = (0:dtFlow:TMaxFlow).';
tTicks = edgesT;

FLs = struct('Geometry', {}, 'BoundingBox', {}, 'X', {}, 'Y', {}, ...
             'Lat0', {}, 'Lon0', {}, 'Rate_mm', {}, 'Azim_deg', {}, 'TMax_ky', {});
Pt  = struct('Geometry', {}, 'X', {}, 'Y', {}, 'Lat0', {}, 'T_ky', {}, 'D_km', {});

for phi = latList
    lonS = -132:0.005:-120;
    [qx, qy] = projfwd(proj, phi*ones(size(lonS)), lonS);
    [xi, yi] = polyxpoly(qx, qy, px, py);
    if isempty(xi)
        warning('No trench crossing at %g N; skipped.', phi);  continue
    end
    [~, k] = max(xi);                              % landward (eastern) edge
    [la0, lo0] = projinv(proj, xi(k), yi(k));
    p0 = [cosd(la0)*cosd(lo0), cosd(la0)*sind(lo0), sind(la0)];

    P = eulerRot(p0, e, w*tPath);
    [fx, fy] = projfwd(proj, asind(P(:,3)), atan2d(P(:,2), P(:,1)));

    v  = Re * w * cross(e.', p0);
    nh = [-sind(la0)*cosd(lo0), -sind(la0)*sind(lo0), cosd(la0)];
    eh = [-sind(lo0), cosd(lo0), 0];
    vn = v*nh.';  ve = v*eh.';

    FLs(end+1) = struct('Geometry', 'Line', ...
        'BoundingBox', [min(fx) min(fy); max(fx) max(fy)], ...
        'X', [fx(:).' NaN], 'Y', [fy(:).' NaN], ...
        'Lat0', phi, 'Lon0', lo0, 'Rate_mm', 1e3*hypot(vn, ve), ...
        'Azim_deg', mod(atan2d(ve, vn), 360), 'TMax_ky', TMaxFlow/1e3); %#ok<SAGROW>

    sinD = norm(cross(e.', p0));
    Pk = eulerRot(p0, e, w*tTicks(:));
    [kxp, kyp] = projfwd(proj, asind(Pk(:,3)), atan2d(Pk(:,2), Pk(:,1)));
    for q = 1:numel(tTicks)
        Pt(end+1) = struct('Geometry', 'Point', 'X', kxp(q), 'Y', kyp(q), ...
            'Lat0', phi, 'T_ky', tTicks(q)/1e3, ...
            'D_km', Re*sinD*w*tTicks(q)/1e3); %#ok<SAGROW>
    end
end
shapewrite(FLs, flowOut);  writePrj(flowOut, proj);
shapewrite(Pt,  tickOut);  writePrj(tickOut, proj);

%% ===================== Local functions =====================
function s = timeRange(a, b)
% Format a time interval in years as in the original legend
if b <= 1e6
    s = sprintf('%g - %g ky', a/1e3, b/1e3);
elseif a >= 1e6
    s = sprintf('%g My - %g My', a/1e6, b/1e6);
else
    s = sprintf('%g ky - %g My', a/1e3, b/1e6);
end
end

function P = mask2poly(M, xs, ys, cs)
% Exact cell-edge outline of a logical mask, built in integer index space
% so shared edges match exactly, then mapped to coordinates.
sgn   = sign(ys(2) - ys(1));
rects = polyshape.empty;
for r = find(any(M, 2)).'
    edgeIdx = diff([false M(r,:) false]);
    st = find(edgeIdx == 1);
    en = find(edgeIdx == -1) - 1;
    v1 = sgn*(r-1) - 0.5;  v2 = sgn*(r-1) + 0.5;
    for k = 1:numel(st)
        u1 = (st(k)-1) - 0.5;  u2 = (en(k)-1) + 0.5;
        rects(end+1) = polyshape([u1 u2 u2 u1], [v1 v1 v2 v2]); %#ok<AGROW>
    end
end
P = union(rects);
P = translate(scale(P, cs, [0 0]), xs(1), ys(1));
end

function P = eulerRot(p0, e, th)
% Rotate unit vector p0 (1x3) about unit axis e (3x1) by angles th (rad).
th = th(:);
P  = p0.*cos(th) + cross(repmat(e.', numel(th), 1), repmat(p0, numel(th), 1), 2).*sin(th) ...
     + (p0*e) * e.' .* (1 - cos(th));
end

function writePrj(shpFile, proj)
% Write a single-line WKT1 .prj next to a shapefile
[pth, nm] = fileparts(shpFile);
fid = fopen(fullfile(pth, [nm '.prj']), 'w');
fprintf(fid, '%s', wktstring(proj, 'Version', 'wkt1', 'Format', 'compact'));
fclose(fid);
end
