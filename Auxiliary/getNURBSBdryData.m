function NURBS = getNURBSBdryData(NURBS)

if NURBS.Dim == 1
    NURBS.Boundary(1).LocalToGlobalCPs = 1;
    NURBS.Boundary(2).LocalToGlobalCPs = NURBS.NumCPsDir(1);
elseif NURBS.Dim == 2
    for iSide = 1 : 4
        Boundary = NURBSBoundary(NURBS, iSide);
        s = int32(iSide);
        if mod(s + idivide(s + 1, 2), 2) == 0
            sgn = 1; % side 1, 4 (clock wise)
        else
            sgn = -1; % side 2, 3 (counter clock wise)
        end
        Boundary.SideOrient = sgn;
        axis = floor((iSide - 1) / 2);
        Boundary.Axis = axis + 1;
        NURBS.Boundary(iSide) = Boundary;
        clear bdry
    end
    NURBS.Boundary(1).LocalToGlobalCPs = uint32(sub2ind(NURBS.NumCPsDir, ones(1, NURBS.NumCPsDir(2)), 1 : NURBS.NumCPsDir(2))');
    NURBS.Boundary(2).LocalToGlobalCPs = uint32(sub2ind(NURBS.NumCPsDir, NURBS.NumCPsDir(1) * ones(1, NURBS.NumCPsDir(2), 'uint32'), 1 : NURBS.NumCPsDir(2))');
    NURBS.Boundary(3).LocalToGlobalCPs = uint32(sub2ind(NURBS.NumCPsDir, 1 : NURBS.NumCPsDir(1), ones(1, NURBS.NumCPsDir(1)))');
    NURBS.Boundary(4).LocalToGlobalCPs = uint32(sub2ind(NURBS.NumCPsDir, 1 : NURBS.NumCPsDir(1), NURBS.NumCPsDir(2) * ones(1, NURBS.NumCPsDir(1), 'uint32'))');
    
    NURBS.Boundary(1).AdjacentCPs = uint32(sub2ind(NURBS.NumCPsDir, 2 * ones(1, NURBS.NumCPsDir(2)), 1 : NURBS.NumCPsDir(2))');
    NURBS.Boundary(2).AdjacentCPs = uint32(sub2ind(NURBS.NumCPsDir, (NURBS.NumCPsDir(1) - 1) * ones(1, NURBS.NumCPsDir(2), 'uint32'), 1 : NURBS.NumCPsDir(2))');
    NURBS.Boundary(3).AdjacentCPs = uint32(sub2ind(NURBS.NumCPsDir, 1 : NURBS.NumCPsDir(1), 2 * ones(1, NURBS.NumCPsDir(1), 'uint32'))');
    NURBS.Boundary(4).AdjacentCPs = uint32(sub2ind(NURBS.NumCPsDir, 1 : NURBS.NumCPsDir(1), (NURBS.NumCPsDir(2) - 1) * ones(1, NURBS.NumCPsDir(1), 'uint32'))');
elseif NURBS.Dim == 3
    for iSide = 1 : 6
        Boundary = NURBSBoundary(NURBS, iSide);
        ind1 = setdiff(1 : 3, ceil(iSide / 2)); %ind = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2]
        ind2 = floor((iSide + 1) / 2); % ind2 = [1 1 2 2 3 3];
        Idcs = ones(3, Boundary.NumCPs);
        [Idcs(ind1(1), :), Idcs(ind1(2), :)] = ind2sub(Boundary.NumCPsDir, 1 : Boundary.NumCPs);
        if (rem(iSide, 2) == 0)
            Idcs(ind2, :) = NURBS.NumCPsDir(ind2);
        end
        Boundary.LocalToGlobalCPs = uint32(sub2ind(NURBS.NumCPsDir, Idcs(1, :), Idcs(2, :), Idcs(3, :))'); % control point numbering
        s = int32(iSide);
        if mod(s + idivide(s + 1, 2), 2) == 0
            sgn = 1; % side 1, 4 (clock wise)
        else
            sgn = -1; % side 2, 3 (counter clock wise)
        end
        Boundary.SideOrient = sgn;
        axis = floor((iSide - 1) / 2);
        Boundary.Axis = axis + 1;
        NURBS.Boundary(iSide) = Boundary;
        clear Boundary
    end
end

end