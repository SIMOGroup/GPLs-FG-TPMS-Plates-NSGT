function [CP,u_knot,v_knot,p,q,nnode,nel,Ien,Inn,B_net,mcp,ncp] = Complex_Domain

[point_coords,~,~] = xlsread('Complicated_shape_FGP_CtrlPts');

CtrlPts = zeros(4, 10, 13);

for j = 1:10
    CtrlPts(1, j, :) = point_coords(13*(j-1)+1:13*j, 1);
    CtrlPts(2, j, :) = point_coords(13*(j-1)+1:13*j, 2);
end
CtrlPts(4, :, :) = 1;

KntVect{1} = [0 0 0 0 0.3006 0.4140 0.5190 0.6149 0.7010 0.7870 1 1 1 1];
KntVect{2} = [0 0 0 0.1285 0.1285 0.2946 0.3765 0.4587 0.5413 0.6235 0.7054 0.8715 0.8715 1 1 1];

SurfAnalysis = CreateNURBS(KntVect, CtrlPts);  % SurfAnalysis now has order [3,2]

% order elevation
p = 3; q = 3;

oNURBS = SurfAnalysis; Order = [p, q]; 
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
SurfAnalysis = oNURBS; % now it has order [3,3]

GeomAnalysis = makeGeometry(SurfAnalysis);
NURBSAnalysis = GeomAnalysis.Patch{1};
u_knot = NURBSAnalysis.KntVect{1};
v_knot = NURBSAnalysis.KntVect{2};
% Ien = NURBSAnalysis.ElemCPs; Ien = Ien';

%% ===== Reshape
Temp = reshape(NURBSAnalysis.CtrlPts3D, 3, []);
CP(:,:,1) = reshape(Temp(1,:),NURBSAnalysis.NumCPsDir(1),NURBSAnalysis.NumCPsDir(2));
CP(:,:,2) = reshape(Temp(2,:),NURBSAnalysis.NumCPsDir(1),NURBSAnalysis.NumCPsDir(2));
CP(:,:,3) = reshape(Temp(3,:),NURBSAnalysis.NumCPsDir(1),NURBSAnalysis.NumCPsDir(2));
CP(:,:,4) = reshape(Temp(3,:)+1,NURBSAnalysis.NumCPsDir(1),NURBSAnalysis.NumCPsDir(2));

nnode = GeomAnalysis.NumBaseDOFs;          % Number of control point
nel = NURBSAnalysis.NumElemDir(1)*NURBSAnalysis.NumElemDir(2);

%============== Input Control Net =====
B_net(:,:,1) = CP(:,:,1);    B_net(:,:,2) = CP(:,:,2);    B_net(:,:,3) = CP(:,:,4);
mcp = length(B_net(:,1,1)); 
ncp = length(B_net(1,:,1));

figure(23)
hold on
set(gcf,'color','white')
Plot_NURBS_Surf(u_knot,v_knot,CP); hold on; axis equal
PlotGeoPatch(SurfAnalysis);
axis off

% plotNURBS_surf_El_CP(p,q,u_knot,v_knot,CP); hold on
% 
daspect([1, 1, 1])
% PlotGeoPatch(VoluDesign);
% PlotKntsPatch(SurfAnalysis);
PlotCtrlPtsPatch(SurfAnalysis, 'r', 1);

% get Ien
for elu = 1:NURBSAnalysis.NumElemDir(1)
    for elv = 1:NURBSAnalysis.NumElemDir(2)
        el = sub2ind(NURBSAnalysis.NumElemDir, elu, elv);
        Ien(el,:) = NURBSAnalysis.ElemCPs(:, el)';
    end    
end

% get Inn
g = 0;
for elu = 1: ncp
    for elv = 1: mcp
        g = g+1;
        Inn(g,1) = elv;
        Inn(g,2) = elu;
    end    
end

