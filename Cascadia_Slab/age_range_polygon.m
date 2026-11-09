function [PS, P] = age_range_polygon(Out, edges, opts)
%AGE_RANGE_POLYGON Polygon mapstruct of cells within age range(s) from Slab_Rotate.
%
%   PS = age_range_polygon(Out, [tMin tMax])
%   Returns an in-memory polygon mapstruct (same structure shaperead returns
%   and shapewrite accepts) outlining every cell of the Slab_Rotate output
%   whose age lies in (tMin, tMax]. Polygons follow exact cell edges, with
%   holes and multiple parts preserved.
%
%   PS = age_range_polygon(Out, edges) with a longer edges vector returns
%   one feature per consecutive range, e.g. edges = 0:2.5e5:2.25e6.
%
%   Inputs
%     Out    Slab_Rotate output (uses Out.time, Out.dist)
%     edges  range limits in years (or metres if Var = 'dist'), ascending
%
%   Name-value options
%     Var         'time' (default) or 'dist'
%     Mask        extra logical mask, e.g. ~isnan(CAS.Z)   (default none)
%     MinCells    drop connected parts smaller than this   (default 0 = keep all)
%     Geographic  return lon/lat vertices instead of projected (default false)
%     Proj        projcrs of the grid, used if Geographic  (default EPSG:26910)
%
%   Outputs
%     PS   1 x k polygon mapstruct with fields Geometry, BoundingBox, X, Y,
%          Bin, Lo, Hi (range in input units), T_lo_ky, T_hi_ky,
%          D_lo_km, D_hi_km (5th-95th percentile of the other variable), NCells
%     P    1 x k polyshape array (projected coordinates)
%
%   Example
%     PS = age_range_polygon(Out, [0 1.25e6], Mask = ~isnan(CAS.Z), MinCells = 20);
%     S  = polygon2GRIDobj(CAS, PS);      % use like a shapefile read with shaperead

arguments
    Out struct
    edges (1,:) double {mustBeNonempty}
    opts.Var char {mustBeMember(opts.Var, {'time', 'dist'})} = 'time'
    opts.Mask = []
    opts.MinCells (1,1) double = 0
    opts.Geographic logical = false
    opts.Proj = projcrs(26910)
end

if numel(edges) < 2 || any(diff(edges) <= 0)
    error('edges must contain at least two ascending values.');
end

G = Out.(opts.Var);                         % GRIDobj of the binning variable
V = double(G.Z);
[X, Y] = getcoordinates(G, 'matrix');
xs = X(1,:);  ys = Y(:,1);  cs = G.cellsize;

base = ~isnan(V);
if ~isempty(opts.Mask)
    base = base & logical(opts.Mask);
end

warning('off', 'MATLAB:polyshape:repairedBySimplify')
nR = numel(edges) - 1;
PS = struct('Geometry', {}, 'BoundingBox', {}, 'X', {}, 'Y', {}, ...
            'Bin', {}, 'Lo', {}, 'Hi', {}, 'T_lo_ky', {}, 'T_hi_ky', {}, ...
            'D_lo_km', {}, 'D_hi_km', {}, 'NCells', {});
P = polyshape.empty;

for i = 1:nR
    M = base & V > edges(i) & V <= edges(i+1);
    if opts.MinCells > 0
        M = bwareaopen(M, opts.MinCells, 8);
    end
    if ~any(M(:)), continue, end

    Pi = mask2poly(M, xs, ys, cs);
    [bx, by] = boundary(Pi);
    if opts.Geographic
        okv = ~isnan(bx);
        [la, lo] = deal(nan(size(bx)));
        [la(okv), lo(okv)] = projinv(opts.Proj, bx(okv), by(okv));
        bx = lo;  by = la;
    end

    tq = prctile(Out.time.Z(M), [5 95]) / 1e3;
    dq = prctile(Out.dist.Z(M), [5 95]) / 1e3;

    PS(end+1) = struct('Geometry', 'Polygon', ...
        'BoundingBox', [min(bx) min(by); max(bx) max(by)], ...
        'X', [bx(:).' NaN], 'Y', [by(:).' NaN], ...
        'Bin', i, 'Lo', edges(i), 'Hi', edges(i+1), ...
        'T_lo_ky', tq(1), 'T_hi_ky', tq(2), ...
        'D_lo_km', dq(1), 'D_hi_km', dq(2), 'NCells', nnz(M)); %#ok<AGROW>
    P(end+1) = Pi; %#ok<AGROW>
end

if isempty(PS)
    warning('No cells fall within the requested range(s).');
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
