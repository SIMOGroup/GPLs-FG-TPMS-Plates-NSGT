function NURBS = CreateNURBS(KntVect, CtrlPts, varargin)
% NURBS = CreateNURBS(KntVects, CtrlPts)
% ------------------------------------
% written by Khanh Chau-Nguyen
% last modified: 22-Oct-2018

assert(iscell(KntVect), 'Input knots must be stored in cell format');
assert(numel(KntVect) >= 1);
assert(numel(KntVect) <= 3);
Dim = numel(KntVect); % The dimensionality of the model
for i = 1 : Dim
    assert(size(KntVect{i}, 1) == 1)
    assert(numel(KntVect{i}) >= 4, 'Number of knot values must be equal or greater than 4')
end
% ------------------------------------
W = CtrlPts(4, :, :, :);
CtrlPts3D = bsxfun(@rdivide, CtrlPts(1 : 3, :, :, :), W);

NumPts = size(CtrlPts);
NumCPsDir = zeros(1, Dim);
p = zeros(1, Dim);
uqKntVect = cell(1, Dim);
uqKntsIdcs = cell(1, Dim);
KntMult = cell(1, Dim);
for i = 1 : Dim
    NumCPsDir(i) = NumPts(i + 1);
    p(i) = numel(KntVect{i}) - NumCPsDir(i) - 1; % (m + 1) - (n + 1) - 1 = m - n - 1
    assert(p(i) > 0, 'Degree of the spline is negative, please check the input parameters!');
    uqKntsIdcs{i} = [true, diff(KntVect{i}) > 0];
    uqKntVect{i} = KntVect{i}(uqKntsIdcs{i});
    KntMult{i} = diff([find(uqKntsIdcs{i} == 1), numel(KntVect{i}) + 1]);
    NURBS.NumElemDir(i) = numel(uqKntVect{i}) - 1;
    NURBS.NumElemShpsDir(i) = p(i) + 1;
    ispan = p(i) + 1 : NumCPsDir(i);
    NURBS.SpanDir{i} = ispan(diff(KntVect{i}(p(i) + 1 : NumCPsDir(i) + 1)) > 1e-12);
end
% check the number of space dimensions
if all(abs(CtrlPts(2 : 3, :, :, :)) <= eps)
    NURBS.NumSpaceDim = 1;
elseif all(abs(CtrlPts(3, :, :, :)) <= eps)
    NURBS.NumSpaceDim = 2;
else
    NURBS.NumSpaceDim = 3;
end
NURBS.KntVect = KntVect;
NURBS.uqKntVect = uqKntVect; % unique knot values
NURBS.KntMult = KntMult; % multiplicities of knot values
NURBS.CtrlPts4D = CtrlPts; % control points in 4D space
NURBS.CtrlPts3D = CtrlPts3D; % control points projected into 3D space
NURBS.Weights = W;
NURBS.Dim = Dim; % number of parametric dimension
% number of control points in each direction
NURBS.NumCPsDir = NumCPsDir;
NURBS.Order = p;
% number of total control points
NURBS.NumCPs = prod(NumCPsDir);
% NURBS.NumElemDir = zeros(Dim, 1); % number of NURBS elements per direction
% NURBS.NumElemShpsDir = zeros(Dim, 1); % number of element shape (basis) functions per direction
NURBS.NumElem = prod(NURBS.NumElemDir);
NURBS.NumElemShps = prod(NURBS.NumElemShpsDir); % the total number of element basis functions
% NURBS.SpanDir = cell(Dim, 1);
NURBS = getCtrlPtsConnectivity(NURBS); % control points connectivity of the elements
NURBS = getNURBSBdryData(NURBS);
NURBS.Orient = evalGeoOrientation(NURBS); % orientation of the NURBS patch (= sign of mapping jacobian determinant)
% NURBS.ElemSize = estimateElemSize(NURBS);
end