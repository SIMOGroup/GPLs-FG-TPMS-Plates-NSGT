function FEM = KmatNURBS_NSGT_HSDT(FEM,ngauss,nel,inn,ien,b_net,Lamd,Mu)
global u_knot v_knot ndof sdof L h

FEM.K = zeros(sdof,sdof); 
FEM.KG = zeros(sdof,sdof); 
% Get gaussian points and weights;
[gp,gw] = genGP_GW(ngauss);
nel_nza = 0; % Elements of non-zero area
tol = 1e-8;% To check u_knot(ni)matches u_knot(ni+1)or not loop over elements;

for iel = 1:nel   
    sctr = ien(iel,:) ;         % Element scatter vector
    nn = length(sctr);
    sctrBb(1:ndof:ndof*nn-6) = ndof.*sctr-6 ;
    sctrBb(2:ndof:ndof*nn-5) = ndof.*sctr-5 ;
    sctrBb(3:ndof:ndof*nn-4) = ndof.*sctr-4 ;
    sctrBb(4:ndof:ndof*nn-3) = ndof.*sctr-3 ;
    sctrBb(5:ndof:ndof*nn-2) = ndof.*sctr-2 ;
    sctrBb(6:ndof:ndof*nn-1) = ndof.*sctr-1 ;
    sctrBb(ndof:ndof:ndof*nn-0) = ndof.*sctr-0 ;

    % Check to see if mlv current element has nonzero area;
    ni = inn(ien(iel,1),1); % Get NURBS coordinates
    nj = inn(ien(iel,1),2);
    % Element has positive area in the parametric domain
    if(abs(u_knot(ni)-u_knot(ni+1))>tol) && (abs(v_knot(nj)-v_knot(nj+1))>tol)
        nel_nza = nel_nza + 1;
        % Used in calculating quadrature points. The factor of 4 comes from mapping from the [-1,1]
        % line onto a real segment...(in jacobian det)
        da =(u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
        %--------------------------------------------------------
        % loop over integration points(ngauss in each direction);
        for igauss = 1: ngauss
            for jgauss = 1: ngauss
%                 [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
                [N,dNdxi,dNdxy,dN2dxy,dN3dxy,dN4dxy,detj] = Kine_Shape_4th(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
%                 [N,dNdxi,dNdxy,dN2dxy,dN3dxy,detj] = Kine_Shape_3th(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
                % calculate given element stiffness matrix and force vector;
                gwt = gw(igauss)*gw(jgauss)*da;

                % Membrane matrix
                Bm = zeros(3, ndof * nn);
                Bm(1, 1 : ndof : ndof * nn - 6) = dNdxy(:, 1)';
                Bm(2, 2 : ndof : ndof * nn - 5) = dNdxy(:, 2)';
                Bm(3, 1 : ndof : ndof * nn - 6) = dNdxy(:, 2)';
                Bm(3, 2 : ndof : ndof * nn - 5) = dNdxy(:, 1)';

                % Bending matrix 1
                Bb1 = zeros(3, ndof * nn);
                Bb1(1, 6 : ndof : ndof * nn - 1) = -dNdxy(:, 1)';
                Bb1(2, 7 : ndof : ndof * nn - 0) = -dNdxy(:, 2)';
                Bb1(3, 6 : ndof : ndof * nn - 1) = -dNdxy(:, 2)';
                Bb1(3, 7 : ndof : ndof * nn - 0) = -dNdxy(:, 1)';

                % Bending matrix 2
                Bb2 = zeros(3, ndof * nn);
                Bb2(1, 4 : ndof : ndof * nn - 3) = dNdxy(:, 1)';
                Bb2(2, 5 : ndof : ndof * nn - 2) = dNdxy(:, 2)';
                Bb2(3, 4 : ndof : ndof * nn - 3) = dNdxy(:, 2)';
                Bb2(3, 5 : ndof : ndof * nn - 2) = dNdxy(:, 1)';

                % Shear matrix 1
                Bs1 = zeros(2, ndof * nn);
                Bs1(1, 3 : ndof : ndof * nn - 4) = dNdxy(:, 1)';
                Bs1(1, 6 : ndof : ndof * nn - 1) = -N';
                Bs1(2, 3 : ndof : ndof * nn - 4) = dNdxy(:, 2)';
                Bs1(2, 7 : ndof : ndof * nn - 0) = -N';

                % Shear matrix 2
                Bs2 = zeros(2, ndof * nn);
                Bs2(1, 4 : ndof : ndof * nn - 3) = N';
                Bs2(2, 5 : ndof : ndof * nn - 2) = N';

                %===   
                Bb = [Bm;Bb1;Bb2];
                Bs = [Bs1;Bs2];

                Kbs  = (Bb'*FEM.DDb*Bb + Bs'*FEM.DDs*Bs)*gwt*detj; 
                
                % ==== Strain gradient terms               
                Bm_SGT = zeros(3, ndof * nn);
                Bm_SGT(1, 1 : ndof : ndof * nn - 6) = (dN3dxy(:,1) + dN3dxy(:,4))';
                Bm_SGT(2, 2 : ndof : ndof * nn - 5) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bm_SGT(3, 1 : ndof : ndof * nn - 6) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bm_SGT(3, 2 : ndof : ndof * nn - 5) = (dN3dxy(:,1) + dN3dxy(:,4))';

                % Bending matrix for NSGT
                Bb1_SGT = zeros(3, ndof * nn);
                Bb1_SGT(1, 6 : ndof : ndof * nn - 1) = -(dN3dxy(:,1) + dN3dxy(:,4))';
                Bb1_SGT(2, 7 : ndof : ndof * nn - 0) = -(dN3dxy(:,3) + dN3dxy(:,2))';
                Bb1_SGT(3, 6 : ndof : ndof * nn - 1) = -(dN3dxy(:,3) + dN3dxy(:,2))';
                Bb1_SGT(3, 7 : ndof : ndof * nn - 0) = -(dN3dxy(:,1) + dN3dxy(:,4))';

                % Bending matrix for NSGT
                Bb2_SGT = zeros(3, ndof * nn);
                Bb2_SGT(1, 4 : ndof : ndof * nn - 3) = (dN3dxy(:,1) + dN3dxy(:,4))';
                Bb2_SGT(2, 5 : ndof : ndof * nn - 2) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bb2_SGT(3, 4 : ndof : ndof * nn - 3) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bb2_SGT(3, 5 : ndof : ndof * nn - 2) = (dN3dxy(:,1) + dN3dxy(:,4))';

                % Shear matrix for NSGT
                Bs1_SGT = zeros(2, ndof * nn);
                Bs1_SGT(1, 3 : ndof : ndof * nn - 4) = (dN3dxy(:,1) + dN3dxy(:,4))';
                Bs1_SGT(1, 6 : ndof : ndof * nn - 1) = -(dN2dxy(:,1) + dN2dxy(:,2))'; % N'
                Bs1_SGT(2, 3 : ndof : ndof * nn - 4) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bs1_SGT(2, 7 : ndof : ndof * nn - 0) = -(dN2dxy(:,1) + dN2dxy(:,2))'; % N'

                Bs2_SGT = zeros(2, ndof * nn);
                Bs2_SGT(1, 4 : ndof : ndof * nn - 3) = (dN2dxy(:,1) + dN2dxy(:,2))'; % N'
                Bs2_SGT(2, 5 : ndof : ndof * nn - 2) = (dN2dxy(:,1) + dN2dxy(:,2))'; % N'

                %===  
                Bb_SGT = [Bm_SGT;Bb1_SGT;Bb2_SGT];
                Bs_SGT = [Bs1_SGT;Bs2_SGT];

                Kbs_SGT  = Lamd*(Bb_SGT'*FEM.DDb*Bb + Bs_SGT'*FEM.DDs*Bs)*gwt*detj; 

                FEM.K(sctrBb,sctrBb) = FEM.K(sctrBb,sctrBb) + Kbs - Kbs_SGT;

                % Calculate buckling term
                % Buckling matrix 
                BG = zeros(2, ndof * nn);
                BG(1, 3 : ndof : ndof * nn - 4) = dNdxy(:, 1)';
                BG(2, 3 : ndof : ndof * nn - 4) = dNdxy(:, 2)';

                if FEM.Load_Type == 1
                    N00 = [1 0;0 0]; % Uni-axial load
                elseif FEM.Load_Type == 2
                    N00 = [1 0;0 1]; % Bi-axial load
                end

                BG_SGT = zeros(2, ndof * nn);
                BG_SGT(1, 3 : ndof : ndof * nn - 4) = (dN3dxy(:,1) + dN3dxy(:,4))';
                BG_SGT(2, 3 : ndof : ndof * nn - 4) = (dN3dxy(:,3) + dN3dxy(:,2))';
                FEM.KG(sctrBb,sctrBb) = FEM.KG(sctrBb,sctrBb) + (BG' - Mu*BG_SGT')*N00*BG*gwt*detj;
            end
        end
    end
end
