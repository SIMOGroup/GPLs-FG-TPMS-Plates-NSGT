function N0 = BasisFuns(Idx, Pts, p, KntVect)
% N0 = BasisFuns(Idx, Pts, p, KntVect)
% -------------------------------------------------------------------
% Compute the nonvanishing B-splines basis functions.
%--------------------------------------------------------------------
% Input:
%      idx: knot span index
%      pts: parametric points
%      p: order of basis
%      KntVect: knot vector
%--------------------------------------------------------------------
% Output:
%      N0: B-spline basis functions
%--------------------------------------------------------------------
% Based on Algorithm A2.2 [The NURBS BOOK, p.70]
%--------------------------------------------------------------------

assert(numel(Idx) == numel(Pts));
N0 = zeros(numel(Pts), p + 1);
for i = 1 : numel(Pts)
    index = Idx(i);
    XiI = Pts(i);
    left = zeros(1, p + 1);
    right = zeros(1, p + 1);
    Ni = zeros(1, p + 1);
    Ni(1) = 1;
    for j = 1 : p
        left(j + 1) = XiI - KntVect(index + 1 - j);
        right(j + 1) = KntVect(index + j) - XiI;
        saved = 0;
        for r = 0 : j - 1
            temp = Ni(r + 1)/(right(r + 2) + left(j - r + 1));
            Ni(r + 1) = saved + right(r + 2) * temp;
            saved = left(j - r + 1) * temp;
        end
        Ni(j + 1) = saved;
    end
    N0(i, :) = Ni;
end
end