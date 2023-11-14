function FEM = MmatNURBS_NSGT_HSDT(FEM,ngauss,nel,inn,ien,b_net,Mu)                                  
global u_knot v_knot ndof sdof

FEM.M = zeros(sdof,sdof);
% Get gaussian points and weights;
[gp,gw] = genGP_GW(ngauss);
nel_nza = 0; % Elements of non-zero area
tol = 1e-8;% To check u_knot(ni)matches u_knot(ni+1)or not loop over elements;

for iel = 1:nel 
    sctr = ien(iel,:) ;         % Element scatter vector
    nn = length(sctr);
    % Determine the position in the global matrix M
    sctrB = zeros(1, ndof * nn);
    sctrB(1:ndof:ndof*nn-6) = ndof.*sctr-6 ;
    sctrB(2:ndof:ndof*nn-5) = ndof.*sctr-5 ;
    sctrB(3:ndof:ndof*nn-4) = ndof.*sctr-4 ;
    sctrB(4:ndof:ndof*nn-3) = ndof.*sctr-3 ;
    sctrB(5:ndof:ndof*nn-2) = ndof.*sctr-2 ;
    sctrB(6:ndof:ndof*nn-1) = ndof.*sctr-1 ;
    sctrB(ndof:ndof:ndof*nn-0) = ndof.*sctr-0 ;

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
%                 [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
                  [N,dNdxi,dNdxy,dN2dxy,dN3dxy,dN4dxy,detj] = Kine_Shape_4th(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);

                gwt = gw(igauss)*gw(jgauss)*da;
                % Calculate mass matrix corresponding to classical terms
                
                R1 = zeros(3, ndof * nn);
                R1(1, 1 : ndof : ndof * nn-6)  =  N';
                R1(2, 2 : ndof : ndof * nn-5)  =  N';
                R1(3, 3 : ndof : ndof * nn-4)  =  N';
                
                R2 = zeros(3, ndof * nn);
                R2(1, 6 : ndof : ndof * nn-1)  =  -N';
                R2(2, 7 : ndof : ndof * nn-0)  =  -N';


                R3 = zeros(3, ndof * nn);
                R3(1, 4 : ndof : ndof * nn-3)  = N'; % fixed
                R3(2, 5 : ndof : ndof * nn-2)  = N'; % fixed

                % ==== Nonlocal terms
                N1 = zeros(3, ndof * nn);
                N1(1, 1 : ndof : ndof * nn-6)  =  dN2dxy(:, 1)' + dN2dxy(:, 2)';
                N1(2, 2 : ndof : ndof * nn-5)  =  dN2dxy(:, 1)' + dN2dxy(:, 2)';
                N1(3, 3 : ndof : ndof * nn-4)  =  dN2dxy(:, 1)' + dN2dxy(:, 2)';

                N2 = zeros(3, ndof * nn);
                N2(1, 6 : ndof : ndof * nn-1)  =  -(dN2dxy(:, 1)' + dN2dxy(:, 2)');
                N2(2, 7 : ndof : ndof * nn-0)  =  -(dN2dxy(:, 1)' + dN2dxy(:, 2)');
                
                N3=zeros(3, ndof * nn);
                N3(1, 4 : ndof : ndof * nn-3)  =  (dN2dxy(:, 1)' + dN2dxy(:, 2)');
                N3(2, 5 : ndof : ndof * nn-2)  =  (dN2dxy(:, 1)' + dN2dxy(:, 2)');

                R = [R1; R2; R3];
                R_Bar = [N1; N2; N3];
                
%                 m = [I0 zeros(2,4);zeros(2,2) I0 zeros(2,2);zeros(2,4) I0];
%                  MM(sctrB, sctrB) = MM(sctrB, sctrB) + R'*m*R*gwt*detj;               
                FEM.M(sctrB, sctrB) = FEM.M(sctrB, sctrB) + (R' - Mu*R_Bar')*FEM.Im*R*gwt*detj;
            end
        end
    end
end
