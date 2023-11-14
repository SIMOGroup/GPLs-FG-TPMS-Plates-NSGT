function Geom = makeGeometry(NURBS, varargin)

% Geom = makeGeometry(Surf, Interfaces, Boundary);

if iscell(NURBS) % mutiple patches geometry
    assert(nargin >= 2, 'For multiple patches geometry, at least 2 input arguments are required (Geom = makeGeometry(NURBS, Boundaries) ) !')
    assert(isfield(varargin{1}(1), 'Patches') && isfield(varargin{1}(1), 'Sides'), 'For multiple patches geometry, the second input argument (Boundary) is mandatory!')
    Boundary = varargin{1};
    Interfaces = [];
    gluedFaces = [];
    if nargin >= 3
        Interfaces = varargin{2};
        if nargin == 4
            gluedFaces = varargin{3};
        end
    end
    Geom = genMultiPatchConn(NURBS, Boundary, Interfaces, gluedFaces);
else
    if nargin == 2 % coupling coincide control points of internal surfaces (curves) by global numbering
        gluedFaces = varargin{1};
        assert(numel(gluedFaces) == 2)
        Geom.gluedFaces = gluedFaces;
        DOFsMapper = zeros(NURBS.NumCPs, 1, 'uint32');
        gluedDofs = union(NURBS.Boundary(gluedFaces(1)).LocalToGlobalCPs, NURBS.Boundary(gluedFaces(2)).LocalToGlobalCPs);
        nonGludedDofs = setdiff(1 : NURBS.NumCPs, gluedDofs);
        
        DOFsMapper(nonGludedDofs) = 1 : numel(nonGludedDofs);
        newDofs = numel(nonGludedDofs) + (1 : numel(NURBS.Boundary(gluedFaces(1)).LocalToGlobalCPs));
        DOFsMapper(NURBS.Boundary(gluedFaces(1)).LocalToGlobalCPs) = newDofs;
        DOFsMapper(NURBS.Boundary(gluedFaces(2)).LocalToGlobalCPs) = newDofs;
        
        %         newDofs = 1 : numel(NURBS.Boundary(gluedFaces(1)).LocalToGlobalCPs);
        %         DOFsMapper(NURBS.Boundary(gluedFaces(1)).LocalToGlobalCPs) = newDofs;
        %         DOFsMapper(NURBS.Boundary(gluedFaces(2)).LocalToGlobalCPs) = newDofs;
        %         DOFsMapper(nonGludedDofs) = numel(newDofs) + (1 : numel(nonGludedDofs));
        
        %         if NURBS.Dim == 2
        %             aux1 = (1 : NURBS.NumCPsDir(NURBS.Boundary(gluedFaces(1)).Axis))';
        %             aux1(end) = 1;
        %             aux2 = (1 : NURBS.Boundary(gluedFaces(1)).NumCPs)';
        %             utmp = repmat(aux1, NURBS.Boundary(gluedFaces(1)).NumCPs, 1);
        %             vtmp = repmat(aux2, 1, NURBS.NumCPsDir(NURBS.Boundary(gluedFaces(1)).Axis))';
        %             vtmp = vtmp(:);
        %
        %             DOFsMapper = uint32(sub2ind([NURBS.NumCPsDir(NURBS.Boundary(gluedFaces(1)).Axis) - 1, NURBS.Boundary(gluedFaces(1)).NumCPs], utmp, vtmp));
        %
        %             NumBaseDOFs = max(DOFsMapper);
        %         end
        
        NumBaseDOFs = numel(nonGludedDofs) + numel(newDofs);
        
        Geom.DOFsMapper{1} = DOFsMapper;
        Geom.Patch{1} = NURBS;
        Geom.NumPatch = 1;
        Geom.ElemDOF = DOFsMapper(NURBS.ElemCPs);
        Geom.Patch{1}.ElemDOF = Geom.ElemDOF;
        Geom.NumBaseDOFs = NumBaseDOFs;
        
        for bdry = 1 : numel(Geom.Patch{1}.Boundary)
            Geom.Patch{1}.Boundary(bdry).LocalToGlobalDOF = DOFsMapper(Geom.Patch{1}.Boundary(bdry).LocalToGlobalCPs);
            if NURBS.Dim == 2
                Geom.Patch{1}.Boundary(bdry).AdjacentDOF = DOFsMapper(Geom.Patch{1}.Boundary(bdry).AdjacentCPs);
            end
        end
    else
        Geom.Patch{1} = NURBS;
        Geom.NumPatch = 1;
        Geom.ElemDOF = NURBS.ElemCPs;
        Geom.Patch{1}.ElemDOF = NURBS.ElemCPs;
        Geom.NumBaseDOFs = NURBS.NumCPs;
        for bdry = 1 : numel(Geom.Patch{1}.Boundary)
            Geom.Patch{1}.Boundary(bdry).LocalToGlobalDOF = Geom.Patch{1}.Boundary(bdry).LocalToGlobalCPs;
            if NURBS.Dim == 2
                Geom.Patch{1}.Boundary(bdry).AdjacentDOF = Geom.Patch{1}.Boundary(bdry).AdjacentCPs;
            end
        end
    end
    % common part
    NumElemsBdry = 0;
    if Geom.Patch{1}.Dim >= 2
        for bdry = 1 : numel(Geom.Patch{1}.Boundary)
            Geom.Patch{1}.Boundary(bdry).ElemDOF = Geom.Patch{1}.Boundary(bdry).LocalToGlobalDOF(NURBS.Boundary(bdry).ElemCPs);
            NumElemsBdry = NumElemsBdry + Geom.Patch{1}.Boundary(bdry).NumElem;
            Geom.Boundary(bdry).NumPatch = 1;
            Geom.Boundary(bdry).Patches(1) = 1;
            Geom.Boundary(bdry).Sides(1) = bdry;
        end
        MaxNumElemShps = max(NURBS.NumElemShpsDir);
        ElemDOFBdry = zeros(MaxNumElemShps, NumElemsBdry, 'uint32'); % dofs connectivity for all 4 boundary sides
        ind = 1;
        for bdry = 1 : numel(Geom.Patch{1}.Boundary)
            NURBSBdry = Geom.Patch{1}.Boundary(bdry);
            ElemDOFBdry(1 : NURBSBdry.NumElemShps, ind : ind + NURBSBdry.NumElem - 1) = NURBSBdry.ElemDOF;
            ind = ind + NURBSBdry.NumElem;
        end
        Geom.MaxNumElemShps = MaxNumElemShps;
        Geom.ElemDOFBdry = ElemDOFBdry;
    end
    if Geom.Patch{1}.Dim >= 1
        [Geom.DOFsGraph{1}.U, Geom.SharedElem{1}.U] = getInteractingCPs(Geom.Patch{1}.Order(1), Geom.Patch{1}.KntVect{1}, uint32(Geom.Patch{1}.SpanDir{1}));
        Geom.DOFsGraph{1}.NumCPsU = diff(Geom.DOFsGraph{1}.U, 1, 2) + 1;
        if Geom.Patch{1}.Dim >= 2
            [Geom.DOFsGraph{1}.V, Geom.SharedElem{1}.V] = getInteractingCPs(Geom.Patch{1}.Order(2), Geom.Patch{1}.KntVect{2}, uint32(Geom.Patch{1}.SpanDir{2}));
            Geom.DOFsGraph{1}.NumCPsV = diff(Geom.DOFsGraph{1}.V, 1, 2) + 1;
            if Geom.Patch{1}.Dim == 3
                [Geom.DOFsGraph{1}.W, Geom.SharedElem{1}.W] = getInteractingCPs(Geom.Patch{1}.Order(3), Geom.Patch{1}.KntVect{3}, uint32(Geom.Patch{1}.SpanDir{3}));
                Geom.DOFsGraph{1}.NumCPsW = diff(Geom.DOFsGraph{1}.W, 1, 2) + 1;
            end
        end
    end
end

% tAllocStart = tic;
[ColPtrBase, RowIndBase] = computeBaseSparsityPattern(Geom);
% disp(['Elapsed allocating time for base sparsity pattern is ', num2str(toc(tAllocStart)), ' seconds.']);
Geom.ColPtrBase = ColPtrBase;
Geom.RowIndBase = RowIndBase;
end

function Geom = genMultiPatchConn(NURBS, Boundary, Interfaces, gluedFaces)
% generate multiple patches connectivity data
% modified from GeoPDEs code

Geom.NumPatch = numel(NURBS);
if ~isempty(Interfaces)
    DOFsMapper = cell(Geom.NumPatch, 1);
    PatchInterfaces = cell(Geom.NumPatch, 1);
    SharedBdryCPs = cell(Geom.NumPatch, numel(Interfaces)); % store the control point numberings of the shared boundaries with respect to patch numberings and interface numberings
    InterfaceCPs = cell(numel(Interfaces), 1);
    
    Geom.MaxNumElemShps = 0; % the maximum number of basis functions per element in the geometry
    Geom.NumElem = 0;
    for iPatch = 1 : Geom.NumPatch
        Geom.Patch{iPatch} = NURBS{iPatch};
        if(Geom.Patch{iPatch}.NumElemShps > Geom.MaxNumElemShps)
            Geom.MaxNumElemShps = Geom.Patch{iPatch}.NumElemShps;
        end
        Geom.NumElem = Geom.NumElem + Geom.Patch{iPatch}.NumElem;
        DOFsMapper{iPatch} = zeros(Geom.Patch{iPatch}.NumCPs, 1, 'uint32');
        PatchInterfaces{iPatch} = union(find([Interfaces.Patch1] == iPatch), find([Interfaces.Patch2] == iPatch)); % find the global indices of interfaces that associated with each patch
    end
    
    for itf = 1 : numel(Interfaces) % loop over the interfaces
        iPatch1  = Interfaces(itf).Patch1;
        iPatch2  = Interfaces(itf).Patch2;
        iSide1 = Interfaces(itf).Side1;
        iSide2 = Interfaces(itf).Side2;
        
        NURBSBdry1 = Geom.Patch{iPatch1}.Boundary(iSide1);
        SharedBdryCPs{iPatch1, itf} = NURBSBdry1.LocalToGlobalCPs;
        
        NURBSBdry2 = Geom.Patch{iPatch2}.Boundary(iSide2);
        
        if Geom.Patch{iPatch1}.Dim == 2
            Ca = NURBSPatchEval(NURBSBdry1, {linspace(0, 1, 5)});
            Cb = NURBSPatchEval(NURBSBdry2, {linspace(0, 1, 5)});
            tfSame = all(all(abs(Ca - Cb) <= 1e-10));
            if (tfSame)
                SharedBdryCPs{iPatch2, itf} = NURBSBdry2.LocalToGlobalCPs;
            else
                SharedBdryCPs{iPatch2, itf} = flipud(NURBSBdry2.LocalToGlobalCPs);
            end
        elseif Geom.Patch{iPatch1}.Dim == 3
            nghbrCPs = reshape(NURBSBdry2.LocalToGlobalCPs, NURBSBdry2.NumCPsDir);
            
            Sa = NURBSPatchEval(NURBSBdry1, {linspace(0, 1, 5), linspace(0, 1, 5)});
            Sb = NURBSPatchEval(NURBSBdry2, {linspace(0, 1, 5), linspace(0, 1, 5)});
            
            %plot3(squeeze( Sa(1, :, :) ), squeeze( Sa(2, :, :) ), squeeze( Sa(3, :, :) ), '.');
            %plot3(squeeze( Sb(1, :, :) ), squeeze( Sb(2, :, :) ), squeeze( Sb(3, :, :) ), '.');
            
            tfSame = all(all(all(abs(Sa - Sb) <= 1e-10)));
            tfSwap = false;
            tfSwap1 = false;
            tfSwap2 = false;
            tfSwap3 = false;
            if (~tfSame)
                tfSame1 = all(all(all(abs(Sa - Sb(:, end : -1 : 1, :)) <= 1e-10))); % if tfSame1 = true, the first parametric direction is reverse
                tfSame2 = all(all(all(abs(Sa - Sb(:, :, end : -1 : 1)) <= 1e-10))); % if tfSame2 = true, the second parametric direction is reverse
                tfSame3 = all(all(all(abs(Sa - Sb(:, end : -1 : 1, end : -1 : 1)) <= 1e-10))); % if tfSame3 = true, the two parametric directions are reverse together
                
                if ~(tfSame1 || tfSame2 || tfSame3)
                    swap = permute(Sb, [1, 3, 2]); % xi -> eta and vice versa
                    tfSwap = all(all(all(abs(Sa - swap) <= 1e-10)));
                    if (~tfSwap)
                        tfSwap1 = all(all(all(abs(Sa - swap(:, end : -1 : 1, :)) <= 1e-10)));
                        tfSwap2 = all(all(all(abs(Sa - swap(:, :, end : -1 : 1)) <= 1e-10)));
                        tfSwap3 = all(all(all(abs(Sa - swap(:, end : -1 : 1, end : -1 : 1)) <= 1e-10)));
                        if ~(tfSwap1 || tfSwap2 || tfSwap3)
                            error('This interface is not conformable, please check the input geometry!')
                        end
                    end
                end
                
                if (tfSwap || tfSwap1 || tfSwap2 || tfSwap3)
                    nghbrCPs = nghbrCPs';
                end
                if (tfSame1 || tfSwap1)
                    nghbrCPs = flipud(nghbrCPs);
                end
                if (tfSame2 || tfSwap2)
                    nghbrCPs = fliplr(nghbrCPs);
                end
                if (tfSame3 || tfSwap3)
                    nghbrCPs = rot90(nghbrCPs, 2);
                end
            end
            nghbrCPs = nghbrCPs(:);
            SharedBdryCPs{iPatch2, itf} = nghbrCPs;
        end
        InterfaceCPs{itf} = zeros(1, Geom.Patch{iPatch1}.Boundary(iSide1).NumCPs);
    end
    
    NumBaseDOFs = 0;
    % We start with the DOFs that do not belong to any interface
    for iPatch = 1 : Geom.NumPatch
        %disp(num2str(ipatch))
        if ~isempty(gluedFaces)
            GluedFace1 = Geom.Patch{iPatch}.Boundary(gluedFaces{iPatch}(1));
            GluedFace2 = Geom.Patch{iPatch}.Boundary(gluedFaces{iPatch}(2));
            
            CPsBdry1 = setdiff(GluedFace1.LocalToGlobalCPs, [SharedBdryCPs{iPatch, :}]);
            CPsBdry2 = setdiff(GluedFace2.LocalToGlobalCPs, [SharedBdryCPs{iPatch, :}]);
            
            newCPsBdry1 = DOFsMapper{iPatch}(CPsBdry1) == 0;
            newCPsBdry2 = DOFsMapper{iPatch}(CPsBdry2) == 0;
            
            newGluedCPs = NumBaseDOFs + (1 : numel(CPsBdry1));
            DOFsMapper{iPatch}(CPsBdry1(newCPsBdry1)) = newGluedCPs;
            DOFsMapper{iPatch}(CPsBdry2(newCPsBdry2)) = newGluedCPs;
            
            NumBaseDOFs = NumBaseDOFs + numel(newCPsBdry1);
            
            GluedCPs = union(GluedFace1.LocalToGlobalCPs, GluedFace2.LocalToGlobalCPs);
            nonInterfaceCPs = setdiff(1 : Geom.Patch{iPatch}.NumCPs, union([SharedBdryCPs{iPatch, :}], GluedCPs));
        else
            nonInterfaceCPs = setdiff(1 : Geom.Patch{iPatch}.NumCPs, [SharedBdryCPs{iPatch, :}]); % [ttform{iPatch, :}]
        end
        NumNonInterfaceCPs = numel(nonInterfaceCPs);
        DOFsMapper{iPatch}(nonInterfaceCPs) = NumBaseDOFs + (1 : NumNonInterfaceCPs);
        NumBaseDOFs = NumBaseDOFs + NumNonInterfaceCPs;
    end
    % Then we set the Interfaces
    for itf = 1 : numel(Interfaces)
        %disp(num2str(itf))
        iPatch = Interfaces(itf).Patch1;
        if ~isempty(gluedFaces)
            BdryIntrfc1 = Geom.Patch{iPatch}.Boundary(gluedFaces{iPatch}(1));
            BdryIntrfc2 = Geom.Patch{iPatch}.Boundary(gluedFaces{iPatch}(2));
            
            GluedCPs = union(BdryIntrfc1.LocalToGlobalCPs, BdryIntrfc2.LocalToGlobalCPs);
            newNonGluedCPsIntrfc = setdiff(SharedBdryCPs{iPatch, itf}, GluedCPs);
            
            SharedCPsInd = DOFsMapper{iPatch}(newNonGluedCPsIntrfc) == 0;
            SharedCPsNumbering = NumBaseDOFs + (1 : numel(newNonGluedCPsIntrfc));
            DOFsMapper{iPatch}(newNonGluedCPsIntrfc(SharedCPsInd)) = SharedCPsNumbering;
            NumBaseDOFs = NumBaseDOFs + numel(newNonGluedCPsIntrfc);
            
            CPsBdry1Intrfc = intersect([SharedBdryCPs{iPatch, itf}], BdryIntrfc1.LocalToGlobalCPs);
            CPsBdry2Intrfc = intersect([SharedBdryCPs{iPatch, itf}], BdryIntrfc2.LocalToGlobalCPs);
            newCPsBdry1Intrfc = find(DOFsMapper{iPatch}(CPsBdry1Intrfc) == 0);
            newCPsBdry2Intrfc = DOFsMapper{iPatch}(CPsBdry2Intrfc) == 0;
            
            newGluedCPsIntrfc = NumBaseDOFs + (1 : numel(newCPsBdry1Intrfc));
            DOFsMapper{iPatch}(CPsBdry1Intrfc(newCPsBdry1Intrfc)) = newGluedCPsIntrfc;
            DOFsMapper{iPatch}(CPsBdry2Intrfc(newCPsBdry2Intrfc)) = newGluedCPsIntrfc;
            
            NumBaseDOFs = NumBaseDOFs + numel(newCPsBdry1Intrfc);
            
            SharedCPsInd = 1 : numel(SharedCPsInd) + 2 * numel(newCPsBdry1Intrfc);
            InterfaceCPs{itf}(SharedCPsInd) = [SharedCPsNumbering, newGluedCPsIntrfc, newGluedCPsIntrfc];
        else
            SharedCPsInd = find(DOFsMapper{iPatch}(SharedBdryCPs{iPatch, itf}) == 0); % find local indices of the (remaining) unnumbered shared control points
            NumSharedCPs = numel(SharedCPsInd); % number of shared control points
            SharedCPsNumbering = NumBaseDOFs + (1 : NumSharedCPs); % numbering these shared control points
            NumBaseDOFs = NumBaseDOFs + NumSharedCPs;
            DOFsMapper{iPatch}(SharedBdryCPs{iPatch, itf}(SharedCPsInd)) = SharedCPsNumbering; % assign the numberings
            InterfaceCPs{itf}(SharedCPsInd) = SharedCPsNumbering;
        end
        [DOFsMapper, InterfaceCPs] = SetSameIntrfc(iPatch, itf, SharedBdryCPs, SharedCPsInd, Interfaces, DOFsMapper, InterfaceCPs, PatchInterfaces);
        [DOFsMapper, InterfaceCPs] = SetSamePatch(iPatch, itf, SharedBdryCPs, SharedCPsInd, Interfaces, DOFsMapper, InterfaceCPs, PatchInterfaces);
    end
else
    NumBaseDOFs = 0;
    DOFsMapper = cell(Geom.NumPatch, 1);
    
    Geom.MaxNumElemShps = 0; % the maximum number of basis functions per element in the geometry
    Geom.NumElem = 0;
    for iPatch = 1 : Geom.NumPatch
        Geom.Patch{iPatch} = NURBS{iPatch};
        if(Geom.Patch{iPatch}.NumElemShps > Geom.MaxNumElemShps)
            Geom.MaxNumElemShps = Geom.Patch{iPatch}.NumElemShps;
        end
        Geom.NumElem = Geom.NumElem + Geom.Patch{iPatch}.NumElem;
        DOFsMapper{iPatch} = uint32(NumBaseDOFs + (1 : Geom.Patch{iPatch}.NumCPs));
        NumBaseDOFs = NumBaseDOFs + Geom.Patch{iPatch}.NumCPs;
    end
end

Geom.DOFsMapper = DOFsMapper;
Geom.NumBaseDOFs = NumBaseDOFs;

Geom.ElemDOF = zeros(Geom.MaxNumElemShps, Geom.NumElem, 'uint32'); % element dofs connectivity of the geometry
ind = 1;
for iPatch = 1 : Geom.NumPatch
    CPsConn = Geom.Patch{iPatch}.ElemCPs; % element control points connectivity of the patch
    Geom.Patch{iPatch}.ElemDOF = Geom.DOFsMapper{iPatch}(CPsConn);
    
    if Geom.Patch{iPatch}.Dim >= 1
        [Geom.DOFsGraph{iPatch}.U, Geom.SharedElem{iPatch}.U] = getInteractingCPs(Geom.Patch{iPatch}.Order(1), Geom.Patch{iPatch}.KntVect{1}, uint32(Geom.Patch{iPatch}.SpanDir{1}));
        Geom.DOFsGraph{iPatch}.NumCPsU = diff(Geom.DOFsGraph{iPatch}.U, 1, 2) + 1;
        if Geom.Patch{iPatch}.Dim >= 2
            [Geom.DOFsGraph{iPatch}.V, Geom.SharedElem{iPatch}.V] = getInteractingCPs(Geom.Patch{iPatch}.Order(2), Geom.Patch{iPatch}.KntVect{2}, uint32(Geom.Patch{iPatch}.SpanDir{2}));
            Geom.DOFsGraph{iPatch}.NumCPsV = diff(Geom.DOFsGraph{iPatch}.V, 1, 2) + 1;
            if Geom.Patch{iPatch}.Dim == 3
               [Geom.DOFsGraph{iPatch}.W, Geom.SharedElem{iPatch}.W] = getInteractingCPs(Geom.Patch{iPatch}.Order(3), Geom.Patch{iPatch}.KntVect{3}, uint32(Geom.Patch{iPatch}.SpanDir{3}));
               Geom.DOFsGraph{iPatch}.NumCPsW = diff(Geom.DOFsGraph{iPatch}.W, 1, 2) + 1;
            end
        end
    end
    
    Geom.ElemDOF(1 : Geom.Patch{iPatch}.NumElemShps, ind : ind + Geom.Patch{iPatch}.NumElem - 1) = Geom.Patch{iPatch}.ElemDOF;
    ind = ind + Geom.Patch{iPatch}.NumElem;
    for bdry = 1 : numel(Geom.Patch{iPatch}.Boundary)
        Geom.Patch{iPatch}.Boundary(bdry).LocalToGlobalDOF = Geom.DOFsMapper{iPatch}(Geom.Patch{iPatch}.Boundary(bdry).LocalToGlobalCPs);
        Geom.Patch{iPatch}.Boundary(bdry).ElemDOF = Geom.Patch{iPatch}.Boundary(bdry).LocalToGlobalDOF(Geom.Patch{iPatch}.Boundary(bdry).ElemCPs);
    end
end
% % process boundaries
Geom.Boundary = Boundary;
NumBdry = numel(Boundary);

MaxNumElemShpsBdry = 0;
NumElemsBdry = 0;
for bdry = 1 : NumBdry
    Geom.Boundary(bdry).NumPatch = numel(Boundary(bdry).Patches);
    for bdrypatch = 1 : Geom.Boundary(bdry).NumPatch
        iPatch = Boundary(bdry).Patches(bdrypatch);
        iSide = Boundary(bdry).Sides(bdrypatch);
        NURBSBdry = Geom.Patch{iPatch}.Boundary(iSide);
        NumElemsBdry = NumElemsBdry + NURBSBdry.NumElem;
        if NURBSBdry.NumElemShps > MaxNumElemShpsBdry
            MaxNumElemShpsBdry = NURBSBdry.NumElemShps;
        end
    end
end

ElemDOFBdry = zeros(MaxNumElemShpsBdry, NumElemsBdry, 'uint32'); % dofs connectivity for all boundary sides
ind = 1;
for bdry = 1 : NumBdry
    for bdrypatch = 1 : Geom.Boundary(bdry).NumPatch
        iPatch = Boundary(bdry).Patches(bdrypatch);
        iSide = Boundary(bdry).Sides(bdrypatch);
        NURBSBdry = Geom.Patch{iPatch}.Boundary(iSide);
        ElemDOFBdry(1 : NURBSBdry.NumElemShps, ind : ind + NURBSBdry.NumElem - 1) = NURBSBdry.ElemDOF;
        ind = ind + NURBSBdry.NumElem;
    end
end
Geom.ElemDOFBdry = ElemDOFBdry;
end

function [DOFsMapper, ppnum] = SetSamePatch(iPatch, intrfc, SharedBdryCPs, SharedCPsInd, Interfaces, DOFsMapper, ppnum, patchIntrfc)
interfaceDofs = SharedBdryCPs{iPatch, intrfc}(SharedCPsInd);
for ii = setdiff(patchIntrfc{iPatch}, intrfc)
    [~, pos1, pos2] = intersect(SharedBdryCPs{iPatch, ii}, interfaceDofs);
    notSet = find(ppnum{ii}(pos1) == 0);
    if (~isempty(notSet))
        ppnum{ii}(pos1(notSet)) = ppnum{intrfc}(pos2(notSet));
        [DOFsMapper, ppnum] = SetSameIntrfc(iPatch, ii, SharedBdryCPs, pos1(notSet), Interfaces, DOFsMapper, ppnum, patchIntrfc);
    end
end
end

function [DOFsMapper, ppnum] = SetSameIntrfc(iPatch, intrfc, SharedBdryCPs, SharedCPsInd, Interfaces, DOFsMapper, ppnum, patchIntrfc)
intrfcDofs = SharedBdryCPs{iPatch, intrfc}(SharedCPsInd);
iPatch2 = setdiff([Interfaces(intrfc).Patch1 Interfaces(intrfc).Patch2],iPatch);
DOFsMapper{iPatch2}(SharedBdryCPs{iPatch2, intrfc}(SharedCPsInd)) = DOFsMapper{iPatch}(intrfcDofs);
[DOFsMapper, ppnum] = SetSamePatch(iPatch2, intrfc, SharedBdryCPs, SharedCPsInd, Interfaces, DOFsMapper, ppnum, patchIntrfc);
end