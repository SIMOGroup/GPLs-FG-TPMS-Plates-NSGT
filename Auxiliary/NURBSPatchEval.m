function varargout = NURBSPatchEval(NURBS, ParaPts, d)
% function varargout = NURBSPatchEval(NURBS, ParaPts, d)
% ------------------------------------------------------------------
% Interpolate parameter points and field
% ------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       ParaPts: parameter points
%       varargin: field variable (optional)
% ------------------------------------------------------------------
% Output:
%       varargout: interpolated points and field (optional)
% ------------------------------------------------------------------

assert(iscell(ParaPts), 'ParaPts must be stored in cell format');
if (nargin >= 2)
    if NURBS.Dim == 1
        USpan = uint32(FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), ParaPts{1}, NURBS.KntVect{1}));
        NU = mDersBasisFuns(USpan, ParaPts{1}, NURBS.Order(1), 0, NURBS.KntVect{1}); % [p + 1, d + 1, NGPs]
        C = mNURBSEval1D(NU, NURBS.CtrlPts4D, USpan, NURBS.NumSpaceDim, 0);
    elseif NURBS.Dim == 2
        assert(numel(ParaPts) == 2, 'For 2D NURBS patch, two parameter point arrays should be passed as the inputs!')
        USpan = uint32(FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), ParaPts{1}, NURBS.KntVect{1}));
        VSpan = uint32(FindSpan(NURBS.NumCPsDir(2), NURBS.Order(2), ParaPts{2}, NURBS.KntVect{2}));
        
        NU = mDersBasisFuns(USpan, ParaPts{1}, NURBS.Order(1), 0, NURBS.KntVect{1}); % [p + 1, d + 1, NGPs]
        NV = mDersBasisFuns(VSpan, ParaPts{2}, NURBS.Order(2), 0, NURBS.KntVect{2}); % [q + 1, d + 1, NGPs]
        
        C = mNURBSEval2D(NU, NV, NURBS.CtrlPts4D, USpan, VSpan, NURBS.NumSpaceDim, 0);
    elseif NURBS.Dim == 3
        assert(numel(ParaPts) == 3, 'For 3D NURBS patch, two parameter point arrays should be passed as the inputs!')
        USpan = uint32(FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), ParaPts{1}, NURBS.KntVect{1}));
        VSpan = uint32(FindSpan(NURBS.NumCPsDir(2), NURBS.Order(2), ParaPts{2}, NURBS.KntVect{2}));
        WSpan = uint32(FindSpan(NURBS.NumCPsDir(3), NURBS.Order(3), ParaPts{3}, NURBS.KntVect{3}));
        
        NU = mDersBasisFuns(USpan, ParaPts{1}, NURBS.Order(1), 0, NURBS.KntVect{1}); % [p + 1, d + 1, NGPs]
        NV = mDersBasisFuns(VSpan, ParaPts{2}, NURBS.Order(2), 0, NURBS.KntVect{2}); % [q + 1, d + 1, NGPs]
        NW = mDersBasisFuns(WSpan, ParaPts{3}, NURBS.Order(3), 0, NURBS.KntVect{3}); % [r + 1, d + 1, NGPs]
        
        C = mNURBSEval3D(NU, NV, NW, NURBS.CtrlPts4D, USpan, VSpan, WSpan, 0);
    end
    if (exist('d', 'var'))
        Field = reshape(d, NURBS.NumCPs, [])';
        NumComp = size(Field, 1);
        Field = reshape(Field, [NumComp, NURBS.NumCPsDir]);
        Fw = zeros(size(NURBS.CtrlPts4D));
        Fw(4, :, :, :) = NURBS.Weights;
        aux = bsxfun(@times, Field, NURBS.Weights);
        
        NumPnts = cellfun(@numel, ParaPts);
        F = zeros([NumComp, NumPnts]);
        for c = 1 : NumComp
            Fw(1, :, :, :) = aux(c, :, :, :);
            if NURBS.Dim == 1
                tmp = mNURBSEval1D(NU, Fw, USpan, NumComp, 0);
            elseif NURBS.Dim == 2
                tmp = mNURBSEval2D(NU, NV, Fw, USpan, VSpan, NumComp, 0);
            elseif NURBS.Dim == 3
                tmp = mNURBSEval3D(NU, NV, NW, Fw, USpan, VSpan, WSpan, 0);
            end
            F(c, :, :, :) = tmp(1, :, :, :);
        end
    end
end

if (nargout == 1)
    varargout{1} = C;
elseif (nargout == 2)
    if (exist('d', 'var'))
        varargout{1} = C;
        varargout{2} = F;
    else
        error('You must input vector of temperature or displacement field');
    end
end
end