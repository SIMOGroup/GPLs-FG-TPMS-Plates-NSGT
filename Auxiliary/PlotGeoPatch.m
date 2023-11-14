function PlotGeoPatch(NURBS, varargin)
% PlotGeo(NURBS)
% Plot geometry
% ------------------------------------------------------------

if NURBS.Dim == 1
    xi = linspace(0, 1, 101); % parametric points
    Cw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, {xi});
    w = Cw(4, :);
    C = bsxfun(@rdivide, Cw, w);
    if (nargin == 2)
        color = varargin{1};
        plot3 (C(1, :), C(2, :), C(3, :), '--', 'Linewidth', 1, 'color', color);
    else
        plot3 (C(1, :), C(2, :), C(3, :), 'Linewidth', 1, 'color', 'k');
    end
elseif NURBS.Dim == 2
    xi = linspace(0, 1, 101); % parametric points
    eta = linspace(0, 1, 101); % parametric points
    % elvaluate the parametric points to plot the surface
    Sw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, {xi, eta});
    % ----------------------------------------------------------
    
    % Plot the NURBS surface
    bcol=[0.76,0.83,0.85];
    bcol = repmat(bcol, 3, 1);
    colormap(bcol);
%     colormap('white');
%     light = 'on';
    [~, m, n] = size(Sw);
    w = Sw(4, :, :);
    % parametric points projected into Cartesian 3D space.
    S = bsxfun(@rdivide, Sw, w);
    x = reshape(S(1, :, :), m, n);
    y = reshape(S(2, :, :), m, n);
    z = reshape(S(3, :, :), m, n);
    surfl(x, y, z);
    shading interp
elseif NURBS.Dim == 3
    for dir = 1 : 3
        for side = [0 1]
            Surf = NURBSExtract(NURBS, dir, side);
            PlotGeoPatch(Surf);
        end
    end    
end
end