% Static bending and free vibration analyses of GPLs-reinforced FG-TPMS
% using nonlocal strain gradient theory

% Written by: Nam V. Nguyen; Kim Q. Tran, H. Nguyen-Xuan
% Email address: nguyennamkt.2311@gmail.com; nguyennamkt.2311@gmail.com

% Reference: Nonlocal strain gradient-based isogeometric analysis of graphene 
% platelets-reinforced functionally graded triply periodic minimal surface
% nanoplates, Applied Mathematics and Computation, 2023

restoredefaultpath; addpath(genpath('./'));
set(0,'defaulttextinterpreter','latex')

clear all; close all; clc;format long g;

global p q mcp ncp nsd nshl nnode nel;
global u_knot v_knot Ien B_net
global gcoord L sdof ndof h
global E_Dis e_0 t_gpl E_gpl nu_gpl rho_gpl lambda_gpl
global Porous_type Rho_Dis M_mean RD_max RD_0
global GPL_Dis l_gpl w_gpl Mat_Model
%====================================================
FEM.ProbType = 1; % 1: Static analysis; 2: Free vibration;
FEM.Example = 2;
FEM.Load_Type = 1;

%[1]: Free vibration analysis of an FG nanoplate (Verification studies)
% Ref: Composite Structures 251 (2020) 112634 - Table 5

%[2]: FG-TPMS reinforced by GPLs with Rho controling
% Ref: Nonlocal strain gradient-based isogeometric analysis of graphene 
% platelets-reinforced functionally graded triply periodic minimal surface
% nanoplates, Applied Mathematics and Computation, 2023
%====================== Plate geometric properties ========================
if FEM.Example == 1  % Square plate
    L = 10; W = L;             
    ah = 10; h = L/ah;
    ell = 4;
    Mu = 4;
    Lamd = ell^2;
elseif FEM.Example == 2
    L = 10; W = L;
    ah = 10; % h = L/ah;
    FEM.Vol = 1;
    h =  FEM.Vol*W/ah;   % Thickness
    Load = 1;
    ell = 2;
    Mu =  4;
    Lamd = ell^2;
    if FEM.Example == 2
        %====== Porosity Distribution =======
        % 1: A (asymmetric), 2: B (symmetric), 3: C (uniform)
        Rho_Dis = 1;
        M_mean = 0.35; RD_max = 1; RD_0 = 0.8;
        Mat_Model = 2;

    elseif FEM.Example == 3
        % 1: A (symmetric); 2: B (asymmetric); 3: C (uniform)
        E_Dis = 1;
        e_0 = 0.9;
        Mat_Model = 2;
    end
    %====== GPLS Distribution =====
    % 1: A (symmetric); 2: B (asymmetric); 3: C (uniform)
    GPL_Dis =   1; % GPL_ij(ij);
    lambda_gpl =  0.01;  % GPL weight fraction
    l_gpl = 2.5e-6; w_gpl = 1.5e-6; t_gpl = 1.5e-9;
    E_gpl = 1.01e12; nu_gpl = 0.186; rho_gpl = 1062.5;
    %%===== Types of porous core =======
    %  1: Primitive, 2: Gyroid, 3: IWP, 4: Closed-cell, 5: Open-cell
    Porous_type =  5;
end

if FEM.Example == 1
    % ==== Get material properties =====
    FEM.IndexMat = 10;
    FEM.Material = 3;
    FEM.Law = 1;
    FEM.DisFunc = 0; % 0: Using Tuan's model 3 (10.1016/j.ijmecsci.2016.01.012)
    FEM = Cal_FGM_HSDT(FEM);

elseif FEM.Example == 2 || FEM.Example == 3
    FEM.DisFunc = 0; % 0: Using Tuan's model 3 (10.1016/j.ijmecsci.2016.01.012)
    FEM = Cal_TPMS_GPL_Rho_control(FEM);
end
%============== NURBS functions =============================
Deg = 3;    % Choose degree
Ref = 11;   % Refinement
FEM.Ref = Ref;

% ===== Mesh Data ====
[CP,U,V,p,q] = Square_Coarse_Mesh(L,W,Deg); % Coarse mesh with order = p
R1 = Refinement_vec(U,Ref); R2 = Refinement_vec(V,Ref);
[CP,u_knot,v_knot] = Knot_Refine_Surf(p,q,U,V,CP,R1,R2);
%============== Plot Mesh ======
% set(gcf,'color','white')
% Plot_NURBS_Surf(u_knot,v_knot,CP); hold on; axis equal
% axis off
%============== Input Control Net =====
B_net(:,:,1)= CP(:,:,1);    B_net(:,:,2)= CP(:,:,2);    B_net(:,:,3)= CP(:,:,4);
mcp = length(B_net(:,1,1)); ncp = length(B_net(1,:,1));
Numx = mcp-p; Numy = ncp-q;
% ============= Global Coordinate of Control Points ======
gcoord(:,1) = reshape(B_net(:,:,1),mcp*ncp,1);
gcoord(:,2) = reshape(B_net(:,:,2),mcp*ncp,1);
gcoord(:,3) = reshape(B_net(:,:,3),mcp*ncp,1);
% ============= Generate Connectivities ==========
[Ien,Inn] = Gen_IEN_INN_2D(p,q,mcp,ncp);
nodes = Connect_Element_Q4(mcp-1,ncp-1);
% ===== Number of elements, degrees of freedom
nnode = mcp*ncp;          % Number of control point
nshl = (p+1)*(q+1);     % Number of local shape functions
nel = (mcp-p)*(ncp-q);    % Number of element
nsd = 2;                  % Number of spatial dimension
ndof = 7; sdof = nnode*ndof;
ngauss = p + 1;            % Number of gauss point in integration
%============================ K - M matrices ===========================
FEM = KmatNURBS_NSGT_HSDT(FEM,ngauss,nel,Inn,Ien,B_net,Lamd,Mu);
FEM = Bcdof_SSSS_7dof(FEM,ndof,mcp,ncp);     %

if FEM.ProbType == 1 % ===== Static analysis
    FEM = FmatNURBS_Uni_HSDT(FEM,ngauss,nel,Inn,Ien,B_net,Load,Mu);
    %     FEM = FmatNURBS_Sine_HSDT(FEM,ngauss,nel,Inn,Ien,B_net,Load,Mu,L,gcoord);
    [FEM.K,FEM.F] = feaplyc2(FEM.K,FEM.F,FEM.BCDof,FEM.BCVal);
    Disp = FEM.K\FEM.F; % Static
    %============================ Post Processing =============================
    C_ele    = (nel + 1)/2;            % Central element
    C_point  = [L/2,W/2] ;             % Coord of center

    Coord = Calculate_xieta(C_ele,C_point,Ien,gcoord,u_knot,v_knot,B_net);
    sctr = Ien(C_ele,:); sctrW = ndof.*sctr - 4;
    [N,dNdxi,dNdxy,dN2dxy,detj] = Kine_Shape_2nd(C_ele,Coord(1,1),Coord(1,2),u_knot,v_knot,B_net);
    Def_cen = N'*Disp(sctrW);
    % Normalized central deflections
    DD = FEM.Em*h^3/(12*(1 - FEM.Num^2));
    Cen_def = round(100*Def_cen*DD/(Load*L^4),4)

elseif FEM.ProbType == 2 % ===== Free vibration
    FEM = MmatNURBS_NSGT_HSDT(FEM,ngauss,nel,Inn,Ien,B_net,Mu);
    Freedofs = setdiff(1 : sdof, FEM.BCDof);
    [ModeShape,Lamda] = eigs(FEM.K(Freedofs, Freedofs),FEM.M(Freedofs, Freedofs), 20, 'sm');
    [Lamda,id] = sort(diag(Lamda));

    n_mode=1;
    for i=1: length(Lamda)
        if or((Lamda(i)<0),(abs(Lamda(i)-1)<1e-5))
            n_mode=n_mode+1;
        else
            break
        end
    end
    FEM.Freq = Lamda(n_mode:n_mode+1);
    if FEM.Example == 1
        G = FEM.Ec/(2*(1 + FEM.Nuc));
        Freq = round(sqrt(Lamda(1:1))*h*sqrt(FEM.Rhoc/G),4)

    elseif FEM.Example == 2
        rhom = 8960; Em = 130e9; num = 0.34; % Copper
        Freq = round(sqrt(Lamda)*L*sqrt(rhom*(1-num^(2))/Em),4)
    end
end