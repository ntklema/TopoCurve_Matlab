function R = raster_swath_stats(G, swaths, opts)
%RASTER_SWATH_STATS Swath statistics of any raster within swath polygons.
%
%   R = raster_swath_stats(G, swaths)
%   For each polygon in a swath shapefile, samples raster G inside the
%   polygon, bins the cells by latitude, and returns the median and
%   interquartile range per latitude bin. Nothing is plotted.
%
%   Latitude is the cell's own map latitude (default) or, if a Slab_Rotate
%   result is passed as Out, the trench-entry latitude.
%
%   Inputs
%     G       GRIDobj to sample
%     swaths  shapefile path or polygon mapstruct (same CRS as G)
%
%   Name-value options
%     Out        Slab_Rotate output; if given, bin by trench-entry latitude
%     Proj       projcrs of the grids               (default EPSG:26910)
%     BinField   attribute used for ordering        (default 'Bin')
%     LabelField attribute used for labels          (default 'Label')
%     LatEdges   latitude bin edges, common to all  (default 42:0.15:50)
%     MinCount   minimum cells per latitude bin     (default 10)
%     Smooth     Gaussian window in bins, 0 = off   (default 0)
%     Scale      multiplier applied to G values     (default 1)
%
%   Output: 1 x n struct array, one element per swath polygon, with fields
%     label    swath label (from LabelField, or 'Swath i')
%     bin      value of BinField (or i)
%     lat      1 x nL latitude bin centres (trench-entry or map latitude)
%     median   1 x nL median of G
%     q25      1 x nL 25th percentile
%     q75      1 x nL 75th percentile
%     iqr      1 x nL interquartile range (q75 - q25)
%     count    1 x nL number of cells per latitude bin
%     latType  'trench-entry' or 'map'
%   Bins with count < MinCount are NaN in median, q25, q75 and iqr.
%
%   Example
%     R = raster_swath_stats(K2, 'slab_swaths_MORVEL.shp', Out = Out, Scale = 1e3);
%     plot(R(1).lat, R(1).median)

arguments
    G GRIDobj
    swaths
    opts.Out = []
    opts.Proj = projcrs(26910)
    opts.BinField char = 'Bin'
    opts.LabelField char = 'Label'
    opts.LatEdges double = 42:0.15:50
    opts.MinCount (1,1) double = 10
    opts.Smooth (1,1) double = 0
    opts.Scale (1,1) double = 1
end

% ---- Swath polygons, ordered by bin attribute ----
if ischar(swaths) || isstring(swaths)
    S = shaperead(char(swaths));
else
    S = swaths;
end
hasBin = isfield(S, opts.BinField);
if hasBin
    [~, ord] = sort([S.(opts.BinField)]);
    S = S(ord);
end
n = numel(S);

% ---- Sampling grid and latitude source ----
useEntry = ~isempty(opts.Out);
if useEntry
    ref = opts.Out.dist;
    G = alignTo(G, ref);
    LX = opts.Out.trench_x;  LY = opts.Out.trench_y;
    latType = 'trench-entry';
else
    ref = G;
    [LX, LY] = getcoordinates(G, 'matrix');
    latType = 'map';
end
V = double(G.Z) * opts.Scale;

eL = opts.LatEdges(:).';
nL = numel(eL) - 1;
lc = (eL(1:end-1) + eL(2:end)) / 2;
q25f = @(v) prctile(v, 25);
q75f = @(v) prctile(v, 75);
blank = nan(1, nL);

R = repmat(struct('label', '', 'bin', NaN, 'lat', lc, 'median', blank, ...
                  'q25', blank, 'q75', blank, 'iqr', blank, ...
                  'count', zeros(1, nL), 'latType', latType), 1, n);

for i = 1:n
    if isfield(S, opts.LabelField)
        R(i).label = char(string(S(i).(opts.LabelField)));
    else
        R(i).label = sprintf('Swath %d', i);
    end
    if hasBin
        R(i).bin = S(i).(opts.BinField);
    else
        R(i).bin = i;
    end

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
    lo  = accumarray(g(ok), v(ok), [nL 1], q25f, NaN).';
    hi  = accumarray(g(ok), v(ok), [nL 1], q75f, NaN).';
    bad = cnt < opts.MinCount;
    med(bad) = NaN;  lo(bad) = NaN;  hi(bad) = NaN;

    if opts.Smooth > 0
        med = smoothdata(med, 'gaussian', opts.Smooth, 'omitnan');
        lo  = smoothdata(lo,  'gaussian', opts.Smooth, 'omitnan');
        hi  = smoothdata(hi,  'gaussian', opts.Smooth, 'omitnan');
        med(bad) = NaN;  lo(bad) = NaN;  hi(bad) = NaN;
    end

    R(i).median = med;
    R(i).q25    = lo;
    R(i).q75    = hi;
    R(i).iqr    = hi - lo;
    R(i).count  = cnt;
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
