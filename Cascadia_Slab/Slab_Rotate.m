function Out = Slab_Rotate(grid, shape, pole, opts)
%SLAB_ROTATE Back-rotate grid cells about an Euler pole to their trench-entry point.
%
%   Out = Slab_Rotate(grid, shape, pole) back-rotates every cell landward of
%   the trench polygon about the Euler pole until it first enters the
%   polygon. Returns elapsed time, flowline arc length, local plate rate,
%   and the trench-entry coordinates for each cell.
%
%   Inputs
%     grid   GRIDobj defining the output grid (projected, e.g. UTM 10N)
%     shape  polygon mapstruct of the trench (same CRS as grid)
%     pole   [lat lon omega] of the subducting plate relative to the
%            overriding plate; degrees and deg/Myr, positive = CCW
%
%   Name-value options
%     Proj          projcrs of the grid (default EPSG:26910)
%     Mask          logical mask of cells to solve (default all)
%     TMax          maximum back-rotation time, yr ([] = automatic, default)
%     TMaxCap       upper limit for automatic TMax, yr (default 10e6)
%     TMaxMargin    factor on the straight-line time estimate (default 1.3)
%     Dt            time step, yr (default: half a cell of travel at max rate)
%     LandwardOnly  exclude cells seaward (west) of the trench (default true)
%
%   Outputs (struct)
%     dist      GRIDobj, flowline arc length from trench (m)
%     time      GRIDobj, time since trench entry (yr)
%     rate      GRIDobj, local plate rate (m/yr)
%     S         GRIDobj, rasterized trench polygon
%     trench_x  trench-entry easting  (matrix, m)
%     trench_y  trench-entry northing (matrix, m)
%     TMax      back-rotation time limit used (yr)

arguments
    grid GRIDobj
    shape struct
    pole (1,3) double
    opts.Proj = projcrs(26910)
    opts.Mask logical = true(grid.size)
    opts.TMax double = []
    opts.TMaxCap (1,1) double = 10e6
    opts.TMaxMargin (1,1) double = 1.3
    opts.Dt double = []
    opts.LandwardOnly logical = true
end

Re  = 6371e3;
S   = polygon2GRIDobj(grid, shape);
inS = S.Z == 1;

[X, Y] = getcoordinates(grid, 'matrix');
xv = X(1,:);  yv = Y(:,1);  cs = grid.cellsize;
x0 = xv(1);   y0 = yv(1);   dy = yv(2) - yv(1);

% Euler axis (unit vector) and rate (rad/yr)
e = [cosd(pole(1))*cosd(pole(2)); cosd(pole(1))*sind(pole(2)); sind(pole(1))];
w = deg2rad(pole(3)) / 1e6;

% Landward side of the trench: per row, east of the easternmost trench cell.
% Rows without trench cells use the nearest row's value.
if opts.LandwardOnly
    xT = nan(grid.size(1), 1);
    for r = find(any(inS, 2)).'
        xT(r) = max(X(r, inS(r,:)));
    end
    xT = fillmissing(xT, 'nearest');
    landward = X > xT;
else
    landward = true(grid.size);
end

todo = find(opts.Mask & ~inS & landward);
N = numel(todo);
[lat, lon] = projinv(opts.Proj, X(todo), Y(todo));
P = [cosd(lat).*cosd(lon), cosd(lat).*sind(lon), sind(lat)];
E = repmat(e.', N, 1);
sinD = vecnorm(cross(E, P, 2), 2, 2);              % sin(angular distance to pole)

if ~any(inS(:))
    error('Trench polygon does not overlap the grid (check that it is in the grid CRS).');
end
if N == 0
    error('No cells to reconstruct: check Mask and that cells lie landward of the trench.');
end

% Automatic TMax: farthest cell from the trench at its local plate rate.
% Straight-line distance is a lower bound on flowline length, hence the margin.
if isempty(opts.TMax)
    dTr   = bwdist(inS) * cs;                      % m to nearest trench cell
    tNeed = dTr(todo) ./ (Re * w * sinD);          % yr
    TMax  = opts.TMaxMargin * max(tNeed);
    if TMax > opts.TMaxCap
        warning('Automatic TMax %.1f My capped at %.1f My; farthest cells may stay NaN.', ...
                TMax/1e6, opts.TMaxCap/1e6);
        TMax = opts.TMaxCap;
    end
else
    TMax = opts.TMax;
end

if isempty(opts.Dt)
    dt = 0.5*cs / (w*Re*max(sinD));
else
    dt = opts.Dt;
end
nt = ceil(TMax / dt);
fprintf('Slab_Rotate: %d cells, TMax = %.2f My, %d steps\n', N, TMax/1e6, nt);

T = nan(N,1);  TX = T;  TY = T;  open = true(N,1);
for k = 1:nt
    th  = -w*k*dt;                                 % backward in time
    idx = find(open);
    p   = P(idx,:);  ek = E(idx,:);
    pr  = p*cos(th) + cross(ek, p, 2)*sin(th) + (p*e).*ek*(1 - cos(th));
    [xr, yr] = projfwd(opts.Proj, asind(pr(:,3)), atan2d(pr(:,2), pr(:,1)));
    c  = round((xr - x0)/cs) + 1;
    r  = round((yr - y0)/dy) + 1;
    ok = r >= 1 & r <= grid.size(1) & c >= 1 & c <= grid.size(2);
    hit = false(numel(idx), 1);
    hit(ok) = inS(sub2ind(grid.size, r(ok), c(ok)));
    T(idx(hit))  = k*dt;
    TX(idx(hit)) = xr(hit);
    TY(idx(hit)) = yr(hit);
    open(idx(hit)) = false;
    if ~any(open), break, end
end

D = nan(grid.size);  Tm = D;  TXg = D;  TYg = D;  V = D;
D(todo)   = Re * w * sinD .* T;                    % flowline arc length (m)
Tm(todo)  = T;                                     % yr
TXg(todo) = TX;  TYg(todo) = TY;
V(todo)   = Re * w * sinD;                         % m/yr
D(inS) = 0;  Tm(inS) = 0;  TXg(inS) = X(inS);  TYg(inS) = Y(inS);

Out.dist = grid;  Out.dist.Z = D;
Out.time = grid;  Out.time.Z = Tm;
Out.rate = grid;  Out.rate.Z = V;
Out.S = S;
Out.trench_x = TXg;
Out.trench_y = TYg;
Out.TMax = TMax;
end
