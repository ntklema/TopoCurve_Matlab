function [hL, hB] = plot_raster_swaths(ax, G, swaths, opts)
%PLOT_RASTER_SWATHS Swath statistics of any raster within swath polygons.
%
%   [hL, hB] = plot_raster_swaths(ax, G, swaths)
%   For each polygon in a swath shapefile (e.g. slab_swaths_MORVEL.shp),
%   samples raster G inside the polygon, bins the cells by latitude, and
%   plots the median (line) and interquartile range (band).
%
%   Latitude is either the cell's own map latitude (default) or, if a
%   Slab_Rotate result is passed as Out, the trench-entry latitude, which
%   reproduces the projection used in the depth-swath figure.
%
%   Inputs
%     ax      target axes
%     G       GRIDobj to sample (any variable: depth, curvature, heat flow...)
%     swaths  shapefile path or polygon mapstruct (same CRS as G)
%
%   Name-value options
%     Out        Slab_Rotate output; if given, plot at trench-entry latitude
%     Proj       projcrs of the grids               (default EPSG:26910)
%     BinField   attribute used for order/colour    (default 'Bin')
%     LabelField attribute used for legend text     (default 'Label')
%     Cmap       colormap name, handle, RGB triplet or m x 3 matrix
%                (resampled to the number of swaths; default copper)
%     LatEdges   latitude bin edges, common to all  (default 42:0.15:50)
%     MinCount   minimum cells per latitude bin     (default 10)
%     Smooth     Gaussian window in bins, 0 = off   (default 3)
%     Scale      multiplier applied to G values     (default 1)
%     YLabel     y-axis label                       (default '')
%     ShowIQR    draw interquartile band            (default true)
%     FaceAlpha  band transparency                  (default 0.3)
%     Legend     draw legend                        (default true)
%
%   Outputs
%     hL   median line handles, one per plotted swath
%     hB   IQR band handles

arguments
    ax (1,1) matlab.graphics.axis.Axes
    G GRIDobj
    swaths
    opts.Out = []
    opts.Proj = projcrs(26910)
    opts.BinField char = 'Bin'
    opts.LabelField char = 'Label'
    opts.Cmap = []
    opts.LatEdges double = 42:0.15:50
    opts.MinCount (1,1) double = 10
    opts.Smooth (1,1) double = 3
    opts.Scale (1,1) double = 1
    opts.YLabel char = ''
    opts.ShowIQR logical = true
    opts.FaceAlpha (1,1) double = 0.3
    opts.Legend logical = true
end

% ---- Swath polygons, ordered by bin attribute ----
if ischar(swaths) || isstring(swaths)
    S = shaperead(char(swaths));
else
    S = swaths;
end
if isfield(S, opts.BinField)
    [~, ord] = sort([S.(opts.BinField)]);
    S = S(ord);
end
n = numel(S);
cmap = makeCmap(opts.Cmap, n);

% ---- Sampling grid and latitude source ----
useEntry = ~isempty(opts.Out);
if useEntry
    ref = opts.Out.dist;                         % reconstruction grid
    G = alignTo(G, ref);
    LX = opts.Out.trench_x;  LY = opts.Out.trench_y;
else
    ref = G;
    [LX, LY] = getcoordinates(G, 'matrix');
end
V = double(G.Z) * opts.Scale;

eL = opts.LatEdges(:).';
nL = numel(eL) - 1;
lc = (eL(1:end-1) + eL(2:end)) / 2;
q25 = @(v) prctile(v, 25);
q75 = @(v) prctile(v, 75);

hold(ax, 'on')
hL = gobjects(0);  hB = gobjects(0);
for i = 1:n
    M  = polygon2GRIDobj(ref, S(i));
    in = M.Z == 1 & ~isnan(V) & ~isnan(LX) & ~isnan(LY);
    if ~any(in(:)), continue, end

    la = projinv(opts.Proj, LX(in), LY(in));
    v  = V(in);
    g  = discretize(la, eL);
    ok = ~isnan(g);
    if ~any(ok), continue, end

    cnt = accumarray(g(ok), 1, [nL 1]).';
    med = accumarray(g(ok), v(ok), [nL 1], @median, NaN).';
    lo  = accumarray(g(ok), v(ok), [nL 1], q25, NaN).';
    hi  = accumarray(g(ok), v(ok), [nL 1], q75, NaN).';
    bad = cnt < opts.MinCount;
    med(bad) = NaN;  lo(bad) = NaN;  hi(bad) = NaN;

    if opts.Smooth > 0
        med = smoothdata(med, 'gaussian', opts.Smooth, 'omitnan');
        lo  = smoothdata(lo,  'gaussian', opts.Smooth, 'omitnan');
        hi  = smoothdata(hi,  'gaussian', opts.Smooth, 'omitnan');
        med(bad) = NaN;  lo(bad) = NaN;  hi(bad) = NaN;
    end

    if isfield(S, opts.LabelField)
        lab = char(string(S(i).(opts.LabelField)));
    else
        lab = sprintf('Swath %d', i);
    end

    % IQR band, drawn per continuous run so gaps stay open
    if opts.ShowIQR
        good = ~isnan(lo) & ~isnan(hi);
        e  = diff([false good false]);
        st = find(e == 1);  en = find(e == -1) - 1;
        for k = 1:numel(st)
            idx = st(k):en(k);
            hB(end+1) = fill(ax, [lc(idx) fliplr(lc(idx))], [lo(idx) fliplr(hi(idx))], ...
                cmap(i,:), 'FaceAlpha', opts.FaceAlpha, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off'); %#ok<AGROW>
        end
    end
    hL(end+1) = plot(ax, lc, med, '-', 'Color', cmap(i,:), 'LineWidth', 2.5, ...
                     'DisplayName', lab); %#ok<AGROW>
end

if useEntry
    xlabel(ax, 'Trench-entry latitude (^\circN)')
else
    xlabel(ax, 'Latitude (^\circN)')
end
if ~isempty(opts.YLabel), ylabel(ax, opts.YLabel); end
if opts.Legend && ~isempty(hL)
    legend(ax, hL, 'Location', 'best');
end
end

function cmap = makeCmap(c, n)
% Build an n x 3 colormap from: [] (copper), a colormap name ('parula'),
% a function handle (@turbo), a single RGB triplet, or an m x 3 matrix
% (interpolated to n rows if m ~= n).
if isempty(c)
    cmap = copper(n);
elseif ischar(c) || isstring(c)
    cmap = feval(char(c), n);
elseif isa(c, 'function_handle')
    cmap = c(n);
else
    cmap = double(c);
end
if size(cmap, 2) ~= 3 || any(cmap(:) < 0 | cmap(:) > 1) || any(isnan(cmap(:)))
    error(['Cmap must be a colormap name, function handle, RGB triplet, ' ...
           'or m x 3 matrix with values in [0, 1].']);
end
m = size(cmap, 1);
if m == 1
    cmap = repmat(cmap, n, 1);
elseif m ~= n
    cmap = interp1(linspace(0, 1, m), cmap, linspace(0, 1, n));
end
end

function G = alignTo(G, ref)
% Resample G onto ref's grid if the two are not aligned
try
    validatealignment(G, ref);
catch
    G = resample(G, ref);
end
end
