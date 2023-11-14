function NURBS = getCtrlPtsConnectivity(NURBS)
% ElemCPs = getCtrlPtsConnectivity(NURBS)
% -------------------------------------------------------------------------
% return control point connectivities for all NURBS elements
% size(Element) = NEN * NEl
% -------------------------------------------------------------------------

% written by Khanh Chau-Nguyen
% last modified: 14-June-2018

NURBS.ElemCPs = zeros(NURBS.NumElemShps, NURBS.NumElem, 'uint32');

for i = 1 : NURBS.Dim
    NURBS.ElemCPsDir{i} = zeros(NURBS.NumElemShpsDir(i), NURBS.NumElemDir(i), 'uint32');
end
% disp('getCtrlPtsConnectivity: This function is under testing!')
for el = 1 : NURBS.NumElem
    subidcs = cell(NURBS.Dim, 1);
    [subidcs{:}] = ind2sub([NURBS.NumElemDir, 1], el);
    aux = cell(NURBS.Dim, 1);
    for d = 1 : NURBS.Dim
        aux{d} = (NURBS.SpanDir{d}(subidcs{d}) - NURBS.Order(d) + (0 : NURBS.Order(d)));
        NURBS.ElemCPsDir{d}(:, subidcs{d}) = aux{d};
    end
    conn = cell(NURBS.Dim, 1);
    [conn{:}] = ndgrid(aux{:});
    idcs = sub2ind ([NURBS.NumCPsDir, 1], conn{:});
    NURBS.ElemCPs(:, el) = idcs(:);
end

% if NURBS.Dim == 1
%     USpan = NURBS.SpanDir{1};
%     Elem = USpan - NURBS.Order; % element number
%     for el = 1 : NURBS.NumElem
%         % List of non-zero basis functions in a given knot-span
%         XiConn = (Elem(el) + (0 : NURBS.Order))';
%         NURBS.ElemCPsDir{1}(:, el) = XiConn;
%         NURBS.ElemCPs(:, el) = XiConn;
%     end
% elseif NURBS.Dim == 2
%     USpan = NURBS.SpanDir{1};
%     VSpan = NURBS.SpanDir{2};
%     for el = 1 : NURBS.NumElem
%         [ue, ve] = ind2sub(NURBS.NumElemDir, el);
%         % List of non-zero basis functions in a given knot-span
%         XiConn = (USpan(ue) - NURBS.Order(1) + (0 : NURBS.Order(1)))';
%         NURBS.ElemCPsDir{1}(:, ue) = XiConn;
%         XiConn = repmat(XiConn, NURBS.Order(2) + 1, 1);
%         
%         EtConn = (VSpan(ve) - NURBS.Order(2) + (0 : NURBS.Order(2)))';
%         NURBS.ElemCPsDir{2}(:, ve) = EtConn;
%         EtConn = repmat(EtConn, 1, NURBS.Order(1) + 1)';
%         EtConn = EtConn(:);
%         
%         NURBS.ElemCPs(:, el) = sub2ind(NURBS.NumCPsDir, XiConn, EtConn);
%     end
% elseif NURBS.Dim == 3
%     USpan = NURBS.SpanDir{1};
%     VSpan = NURBS.SpanDir{2};
%     WSpan = NURBS.SpanDir{3};
%     for el = 1 : NURBS.NumElem
%         [ue, ve, we] = ind2sub(NURBS.NumElemDir, el);
%         % List of non-zero basis functions in a given knot-span
%         XiConn = (USpan(ue) - NURBS.Order(1) + (0 : NURBS.Order(1)))';
%         NURBS.ElemCPsDir{1}(:, ue) = XiConn;
%         XiConn = repmat(XiConn, prod(NURBS.Order(2 : 3) + 1), 1);
%         
%         EtConn = (VSpan(ve) - NURBS.Order(2) + (0 : NURBS.Order(2)))';
%         NURBS.ElemCPsDir{2}(:, ve) = EtConn;
%         EtConn = repmat(EtConn, 1, NURBS.Order(1) + 1)';
%         EtConn = EtConn(:);
%         EtConn = repmat(EtConn, 1, NURBS.Order(3) + 1);
%         EtConn = EtConn(:);
%         
%         ZeConn = (WSpan(we) - NURBS.Order(3) + (0 : NURBS.Order(3)))';
%         NURBS.ElemCPsDir{3}(:, we) = ZeConn;
%         ZeConn = repmat(ZeConn, 1, prod(NURBS.Order(1 : 2) + 1))';
%         ZeConn = ZeConn(:);
%         
%         NURBS.ElemCPs(:, el) = sub2ind(NURBS.NumCPsDir, XiConn, EtConn, ZeConn);
%     end
% end

end