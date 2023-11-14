function [InterCPs, SharedElems] = getInteractingCPs(p, KntVect, KntSpan)
% Infer interacting control points for each control point
%==========================================================================
% Input:
%       p: order
%       KntVect: knot vector
%       KntSpan: knot span index
%--------------------------------------------------------------------------
% Output:
%       InterCPs [NumCPs x 2]: matrix of indices of interacting control 
%       points, each row contains the beginning and ending indices of 
%       influent basis functions
%       SharedElems [NumCPs x 2]: matrix of indices of support knot spans, 
%       each row contains the beginning and ending element indices of 
%       the support domain
%==========================================================================

ext = mexext;
switch ext
    case 'mexw64'
        if exist(['mexGetInteractingCPs', '.', mexext], 'file')
            % Call the mex routine mexGetInteractingCPs.
            [InterCPs, SharedElems] = mexGetInteractingCPs(p, KntVect, KntSpan);
        else
            error('.MEX-file not found on path.');
        end
    otherwise
        error('Unsupported Architecture')
end
end