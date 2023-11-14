function ONURBS = NURBSExtract(NURBS, Dir, val)
% function ONURBS = NURBSExtract(NURBS, Dir, val)
% ------------------------------------------------------------------
% Extract lower dimensional NURBS object.
%-------------------------------------------------------------------
% Input:
%       NURBS: input NURBS structure
%       Dir: direction along which to extract
%       (xi: dir = 1, eta: dir = 2, zeta: dir = 3)
%       val: parametric value from which to exact
%-------------------------------------------------------------------
% Output:
%       ONURBS: output NURBS structure
%-------------------------------------------------------------------

% if 'Dir'th dimension of CtrlPts4D matrix is not the last
% dimension, then permute it until it lies at last dimension
if Dir ~= NURBS.Dim
    Dirs = 1 : NURBS.Dim + 1;
    Dirs(Dir + 1) = [];
    Dirs = [Dirs, Dir + 1];
    temp = permute(NURBS.CtrlPts4D, Dirs);
else
    temp = NURBS.CtrlPts4D;
end
dim = size(temp);
temp = reshape(temp, [], dim(end));

CtrlPts = CurvPntByCornerCut(NURBS.NumCPsDir(Dir), NURBS.Order(Dir), NURBS.KntVect{Dir}, temp, val);

CtrlPts = reshape(CtrlPts, dim(1 : end - 1));

KntVect = NURBS.KntVect;
KntVect{Dir} = [];
KntVect = KntVect(~cellfun(@isempty, KntVect));
ONURBS = CreateNURBS(KntVect, CtrlPts);
end