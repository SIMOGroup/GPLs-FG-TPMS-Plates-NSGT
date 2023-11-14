function orient = evalGeoOrientation(NURBS)

% evaluate orientation of the geometry
if NURBS.Dim == 2
    knt1 = 0.5;
    knt2 = 0.5;
    USpan = FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), knt1, NURBS.KntVect{1});
    NU = mDersBasisFuns(uint32(USpan), knt1, NURBS.Order(1), 1, NURBS.KntVect{1});
    
    VSpan = FindSpan(NURBS.NumCPsDir(2), NURBS.Order(2), knt2, NURBS.KntVect{2});
    NV = mDersBasisFuns(uint32(VSpan), knt2, NURBS.Order(2), 1, NURBS.KntVect{2});
    
    [~, jab] = mNURBSEval2D(NU, NV, NURBS.CtrlPts4D, uint32(USpan), uint32(VSpan), NURBS.NumSpaceDim, 1);
    
    if NURBS.NumSpaceDim == 2
        jacDet = det(jab);
        orient = (0 < jacDet) - (jacDet < 0);
    else
        orient = 1;
    end
    
elseif NURBS.Dim == 3
    knt1 = 0.5;
    knt2 = 0.5;
    knt3 = 0.5;
    USpan = FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), knt1, NURBS.KntVect{1});
    NU = mDersBasisFuns(uint32(USpan), knt1, NURBS.Order(1), 1, NURBS.KntVect{1});
    
    VSpan = FindSpan(NURBS.NumCPsDir(2), NURBS.Order(2), knt2, NURBS.KntVect{2});
    NV = mDersBasisFuns(uint32(VSpan), knt2, NURBS.Order(2), 1, NURBS.KntVect{2});
    
    WSpan = FindSpan(NURBS.NumCPsDir(3), NURBS.Order(3), knt3, NURBS.KntVect{3});
    NW = mDersBasisFuns(uint32(WSpan), knt3, NURBS.Order(3), 1, NURBS.KntVect{3});
    
    [~, jab] = mNURBSEval3D(NU, NV, NW, NURBS.CtrlPts4D, uint32(USpan), uint32(VSpan), uint32(WSpan), 1);
    
    jacDet = det(jab);
    
    orient = (0 < jacDet) - (jacDet < 0);
else
    orient = 1;
end
end