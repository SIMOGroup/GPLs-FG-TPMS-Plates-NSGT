function Elem = getElemFromParaPt(NURBS, u, v, w)

if NURBS.Dim == 1 && exist('u', 'var')
    USpan = FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), u, NURBS.KntVect{1});
    % List of non-zero basis functions in a given knot-span
    XiConn = (USpan - NURBS.Order(1) + (0 : NURBS.Order(1)));
    Elem.ElemCPs = XiConn;
    Elem.ElemIND = find(NURBS.SpanDir{1} == USpan);
elseif NURBS.Dim == 2 && exist('u', 'var') && exist('v', 'var')
    USpan = FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), u, NURBS.KntVect{1});
    VSpan = FindSpan(NURBS.NumCPsDir(2), NURBS.Order(2), v, NURBS.KntVect{2});
    
    % List of non-zero basis functions in a given knot-span
    XiConn = (USpan - NURBS.Order(1) + (0 : NURBS.Order(1)));
    Elem.ElemCPsDir{1} = XiConn;
    XiConn = repmat(XiConn', NURBS.Order(2) + 1, 1);
    
    EtaConn = (VSpan - NURBS.Order(2) + (0 : NURBS.Order(2)));
    Elem.ElemCPsDir{2} = EtaConn;
    EtaConn = repmat(EtaConn', 1, NURBS.Order(1) + 1)';
    EtaConn = EtaConn(:);
    
    LConn = sub2ind(NURBS.NumCPsDir, XiConn, EtaConn);
    Elem.ElemCPs = LConn;
    
    ElXiInd = find(NURBS.SpanDir{1} == USpan);
    ElEtInd = find(NURBS.SpanDir{2} == VSpan);
    
    Elem.ElemIND = sub2ind(NURBS.NumElemDir, ElXiInd, ElEtInd); % element index
elseif NURBS.Dim == 3 && exist('u', 'var') && exist('v', 'var') && exist('w', 'var')
    USpan = FindSpan(NURBS.NumCPsDir(1), NURBS.Order(1), u, NURBS.KntVect{1});
    VSpan = FindSpan(NURBS.NumCPsDir(2), NURBS.Order(2), v, NURBS.KntVect{2});
    WSpan = FindSpan(NURBS.NumCPsDir(3), NURBS.Order(3), w, NURBS.KntVect{3});
    
    % List of non-zero basis functions in a given knot-span
    XiConn = (USpan - NURBS.Order(1) + (0 : NURBS.Order(1)));
    Elem.ElemCPsDir{1} = XiConn;
    
    XiConn = repmat(XiConn', prod(NURBS.Order(2 : 3) + 1), 1);
    
    EtConn = (VSpan - NURBS.Order(2) + (0 : NURBS.Order(2)));
    Elem.ElemCPsDir{2} = EtConn;
    EtConn = repmat(EtConn', 1, NURBS.Order(1) + 1)';
    EtConn = EtConn(:);
    EtConn = repmat(EtConn, 1, NURBS.Order(3) + 1);
    EtConn = EtConn(:);
    
    ZeConn = (WSpan - NURBS.Order(3) + (0 : NURBS.Order(3)));
    Elem.ElemCPsDir{3} = ZeConn;
    ZeConn = repmat(ZeConn', 1, prod(NURBS.Order(1 : 2) + 1))';
    ZeConn = ZeConn(:);
    
    LConn = sub2ind(NURBS.NumCPsDir, XiConn, EtConn, ZeConn);
    Elem.ElemCPs = LConn;
    
    ElXiInd = find(NURBS.SpanDir{1} == USpan);
    ElEtInd = find(NURBS.SpanDir{2} == VSpan);
    ElZeInd = find(NURBS.SpanDir{3} == WSpan);
    
    Elem.ElemIND = sub2ind(NURBS.NumElemDir, ElXiInd, ElEtInd, ElZeInd); % element index
end
end