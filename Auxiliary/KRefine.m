function oNURBS = KRefine(iNURBS, NumElems, Order, Continuity)
% function ONURBS = KRefine(INURBS, NEl, Order, Continuity)
% -----------------------------------------------------------------
% Refine a NURBS object (curve, surface, volume) by knot refinement
% and order (degree) elevation
% -----------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       NEl: number of elements per direction
%       Order: order of basis functions per direction
%       Continuity: continuity of basis functions per direction
% -------------------------------------------------------------
% Output:
%       ONURBS: output NURBS structure
% -------------------------------------------------------------

oNURBS = iNURBS;
% order elevation
assert(all(Order - oNURBS.Order >= 0));
if all(Order)
    for d = 1 : oNURBS.Dim
        p = Order(d);
        t = p - oNURBS.Order(d);
        if t > 0
            oNURBS = PRefine(oNURBS, d, t);
        end
    end
end
% knot refinement
assert(all(Continuity >= 0));
assert(all(Continuity < Order));
assert(all(NumElems > 0));
for d = 1 : oNURBS.Dim
    knts = oNURBS.uqKntVect{d};
    mlts = oNURBS.KntMult{d};
    N = NumElems(d);
    % compute the number of knots to be inserted
    dxi = diff(knts) / N;
    knts = knts(1 : end - 1);
    xi = repmat(knts, N, 1); 
    step = 1 : N - 1;    
    xi(2 : end, :) = xi(2 : end, :) + dxi .* step';
    knts = xi(:);
    knts = knts(2 : end);
        
    mlts = mlts(1 : end - 1);
    repsMat = zeros(N, numel(mlts));
    repsMat(1, :) = oNURBS.Order(d) - Continuity(d) - mlts;
    repsMat(2 : end, :) = oNURBS.Order(d) - Continuity(d);  
    repsMat = repsMat(:);
    repsVect = repsMat(2 : end);
    
    insKnts = knts(repsVect > 0); % inserted knots
    repsVect = repsVect(repsVect > 0);
    
    idcs = zeros(sum(repsVect), 1);
    idcs(cumsum([1; repsVect(1 : end - 1)])) = 1;
    idcs = cumsum(idcs);
    
    if ~isempty(insKnts)
        KntsMult = insKnts(idcs);
        % insert multiple knots
        oNURBS = HRefine(oNURBS, d, KntsMult);
    end
end
end