function [K] = KmatNURBS_NSGT(ngauss,nel,inn,ien,b_net,Db,Ds,Lamd)
global u_knot v_knot ndof sdof

K = zeros(sdof,sdof); 
% Get gaussian points and weights;
[gp,gw] = genGP_GW(ngauss);
nel_nza = 0; % Elements of non-zero area
tol = 1e-8;% To check u_knot(ni)matches u_knot(ni+1)or not loop over elements;

for iel = 1:nel   
    sctr = ien(iel,:) ;         % Element scatter vector
    nn = length(sctr);
    sctrBb(1:ndof:ndof*nn-4) = ndof*sctr-4 ;
    sctrBb(2:ndof:ndof*nn-3) = ndof*sctr-3 ;
    sctrBb(3:ndof:ndof*nn-2) = ndof*sctr-2 ;
    sctrBb(4:ndof:ndof*nn-1) = ndof*sctr-1 ;
    sctrBb(5:ndof:ndof*nn-0) = ndof*sctr ;
    % Check to see if mlv current element has nonzero area;
    ni = inn(ien(iel,1),1); % Get NURBS coordinates
    nj = inn(ien(iel,1),2);
    % Element has positive area in the parametric domain
    if(abs(u_knot(ni)-u_knot(ni+1))>tol)&&(abs(v_knot(nj)-v_knot(nj+1))>tol)
        nel_nza = nel_nza + 1;
        % Used in calculating quadrature points. The factor of 4 comes from mapping from the [-1,1]
        % line onto a real segment...(in jacobian det)
        da =(u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
        %--------------------------------------------------------
        % loop over integration points(ngauss in each direction);
        for igauss = 1: ngauss
            for jgauss = 1: ngauss
% %                 [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
                [N,dNdxi,dNdxy,dN2dxy,dN3dxy,dN4dxy,detj] = Kine_Shape_4th(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
                % calculate given element stiffness matrix and force vector;
                gwt = gw(igauss)*gw(jgauss)*da;
                % Membrane matrix
                Bm = zeros(3, ndof * nn);
                Bm(1, 1 : ndof : ndof * nn - 4) = dNdxy(:, 1)';
                Bm(2, 2 : ndof : ndof * nn - 3) = dNdxy(:, 2)';
                Bm(3, 1 : ndof : ndof * nn - 4) = dNdxy(:, 2)';
                Bm(3, 2 : ndof : ndof * nn - 3) = dNdxy(:, 1)';
                % Bending matrix
                Bb = zeros(3, ndof * nn);
                Bb(1, 4 : ndof : ndof * nn - 1) = dNdxy(:, 1)';
                Bb(2, 5 : ndof : ndof * nn - 0) = dNdxy(:, 2)';
                Bb(3, 4 : ndof : ndof * nn - 1) = dNdxy(:, 2)';
                Bb(3, 5 : ndof : ndof * nn - 0) = dNdxy(:, 1)';

                % Shear matrix
                Bs = zeros(2, ndof * nn);
                Bs(1, 3 : ndof : ndof * nn - 2) = dNdxy(:, 1)';
                Bs(1, 4 : ndof : ndof * nn - 1) = N';
                Bs(2, 3 : ndof : ndof * nn - 2) = dNdxy(:, 2)';
                Bs(2, 5 : ndof : ndof * nn - 0) = N';
                %===   
                Bbd = [Bm;Bb];
                Kb = Bbd'*Db*Bbd*gwt*detj; % Bending part
                Ks = Bs'*Ds*Bs*gwt*detj; % Shear part 
              
                % ==== Strain gradient terms               
                Bm1 = zeros(3, ndof * nn);
                Bm1(1, 1 : ndof : ndof * nn - 4) = (dN3dxy(:,1) + dN3dxy(:,4))';
                Bm1(2, 2 : ndof : ndof * nn - 3) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bm1(3, 1 : ndof : ndof * nn - 4) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bm1(3, 2 : ndof : ndof * nn - 3) = (dN3dxy(:,1) + dN3dxy(:,4))';
                % Bending matrix
                Bb1 = zeros(3, ndof * nn);
                Bb1(1, 4 : ndof : ndof * nn - 1) = (dN3dxy(:,1) + dN3dxy(:,4))';
                Bb1(2, 5 : ndof : ndof * nn - 0) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bb1(3, 4 : ndof : ndof * nn - 1) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bb1(3, 5 : ndof : ndof * nn - 0) = (dN3dxy(:,1) + dN3dxy(:,4))';

                % Shear matrix
                Bs1 = zeros(2, ndof * nn);
                Bs1(1, 3 : ndof : ndof * nn - 2) = (dN3dxy(:,1) + dN3dxy(:,4))';
                Bs1(1, 4 : ndof : ndof * nn - 1) = (dN2dxy(:,1) + dN2dxy(:,2))'; % N'
                Bs1(2, 3 : ndof : ndof * nn - 2) = (dN3dxy(:,3) + dN3dxy(:,2))';
                Bs1(2, 5 : ndof : ndof * nn - 0) = (dN2dxy(:,1) + dN2dxy(:,2))'; % N'
                %===   
                Bbd1 = [Bm1;Bb1];
                Kb1 = Lamd*Bbd1'*Db*Bbd*gwt*detj; % Bending part
                Ks1 = Lamd*Bs1'*Ds*Bs*gwt*detj; % Shear part 
                        
                K(sctrBb,sctrBb) = K(sctrBb,sctrBb) + Kb + Ks - Kb1 - Ks1;
            end
        end
    end
end
