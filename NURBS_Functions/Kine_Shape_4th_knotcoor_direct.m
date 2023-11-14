% Subroutine eval_SHAPE.f consumes an element number and the coordinates in the
% parent element of an integration point and returns the vector of all local
% basis functions evaluated at the point and the matrix of gradients for all
% nonzero bais functions with respect to parameters u and v and with respect
% to x and y.
%
%  June 17, 2003
%  J. Austin Cottrell
%  CES Graduate Student
%  Texas Institute for Computational Engineering Science
%  University of Texas at Austin
%
%  Modify to codes Matlab by :
%  Hung Nguyen Xuan
%
%   Faculty of Mathematics & Informatics, University of Natural Sciences
%   Vietnam   National University HCM

function [R,shgradl,shgradg,shgradg2d,shgradg3d,shgradg4d,detj]=Kine_Shape_4th_knotcoor_direct(e,u_hat,v_hat,u_knot,v_knot,b_net)
global nsd nshl mcp ncp p q

shl=zeros(nshl,1);
shgradl=zeros(nshl,nsd);
shgradl2=zeros(nshl,nsd+1);
shgradl3=zeros(nshl,nsd+2);
shgradl4=zeros(nshl,nsd+3);
denom_sum=0;
derv_sum_u=0;
derv_sum_v=0;
derv_sum_uu = 0 ;
derv_sum_vv = 0 ;
derv_sum_uv = 0 ;
derv_sum_uuu = 0 ;
derv_sum_vvv = 0 ;
derv_sum_uuv = 0 ;
derv_sum_uvv = 0 ;
derv_sum_uuuu = 0 ;
derv_sum_vvvv = 0 ;
derv_sum_uuuv = 0 ;
derv_sum_uuvv = 0 ;
derv_sum_uvvv = 0 ;
% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;
[ien,inn]=genIEN_INN_2D(p,q,mcp,ncp);
ni = inn(ien(e,1),1);
nj = inn(ien(e,1),2);
%     get u and v coordinates
u =u_hat;
v =v_hat;

% evaluate 1d size functions and derivatives each direction
% size and M and N = number of derviatives+shape functions, degree of
% poynomial.
% row 1 of M and N => shape functions
% i^{th} row (i > 1) => i^{th} derivative of the shape function
% calculate in u direction
M = dersbasisfuns(ni,p,mcp,u,u_knot) ;

% calculate in v direction
N = dersbasisfuns(nj,q,ncp,v,v_knot) ;

% form basis functions and derivatives dr./du and dr./dv;
icount = 0;

for j = 0:q
    for i = 0:p
        icount = icount+1;
        
        % basis functions
        shl(icount,1) = M(1,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1);
        denom_sum = denom_sum + shl(icount);
        
        % first derivatives
        shgradl(icount,1) = M(2,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1); %u
        derv_sum_u = derv_sum_u + shgradl(icount,1);
        shgradl(icount,2) = M(1,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1); %v
        derv_sum_v = derv_sum_v + shgradl(icount,2);
        
        % second derivatives
        
        % wrt u
        shgradl2(icount,1) = M(3,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uu = derv_sum_uu + shgradl2(icount,1) ;
        
        % wrt v
        shgradl2(icount,2) = M(1,p+1-i)*N(3,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_vv = derv_sum_vv + shgradl2(icount,2) ;
        
        % cross derivative...wrt uv
        shgradl2(icount,3) = M(2,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uv = derv_sum_uv + shgradl2(icount,3) ;
        
        % third derivatives
        
        % wrt uuu
        shgradl3(icount,1) = M(4,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uuu = derv_sum_uuu + shgradl3(icount,1) ;
        
        % wrt vvv
        shgradl3(icount,2) = M(1,p+1-i)*N(4,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_vvv = derv_sum_vvv + shgradl3(icount,2) ;
        
        % wrt uuv
        shgradl3(icount,3) = M(3,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uuv = derv_sum_uuv + shgradl3(icount,3) ;
        
        % wrt uvv
        shgradl3(icount,4) = M(2,p+1-i)*N(3,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uvv = derv_sum_uvv + shgradl3(icount,4) ;
        
        % fourth derivatives
        
        % wrt uuuu
        shgradl4(icount,1) = M(5,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uuuu = derv_sum_uuuu + shgradl4(icount,1) ;
        
        % wrt vvvv
        shgradl4(icount,2) = M(1,p+1-i)*N(5,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_vvvv = derv_sum_vvvv + shgradl4(icount,2) ;
        
        % wrt uuuv
        shgradl4(icount,3) = M(4,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uuuv = derv_sum_uuuv + shgradl4(icount,3) ;
        
        % wrt uuvv
        shgradl4(icount,4) = M(3,p+1-i)*N(3,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uuvv = derv_sum_uuvv + shgradl4(icount,4) ;
        
        % wrt uvvv
        shgradl4(icount,5) = M(2,p+1-i)*N(4,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uvvv = derv_sum_uvvv + shgradl4(icount,5) ;
    end
end

% basis functions
R = shl/denom_sum;

% First derivative.....divide through by denominator
% wrt u
shgradl(:,1) = shgradl(:,1)/denom_sum -(shl(:)*derv_sum_u)/(denom_sum^2);

% wrt v
shgradl(:,2) = shgradl(:,2)/denom_sum -(shl(:)*derv_sum_v)/(denom_sum^2);

% Second derivative....
% wrt uu
shgradl2(:,1) = shgradl2(:,1)/denom_sum - 2*derv_sum_u*shgradl(:,1)/denom_sum ...
    - derv_sum_uu*R(:)/denom_sum ;

% wrt vv
shgradl2(:,2) = shgradl2(:,2)/denom_sum - 2*derv_sum_v*shgradl(:,2)/denom_sum ...
    - derv_sum_vv*R(:)/denom_sum ;

% wrtuv
shgradl2(:,3) = shgradl2(:,3)/denom_sum - shgradl(:,1)*derv_sum_v/denom_sum ...
    - shgradl(:,2)*derv_sum_u/denom_sum ...
    - R(:)*derv_sum_uv/denom_sum ;

% Third derivative....
% wrt uuu
shgradl3(:,1) = shgradl3(:,1)/denom_sum - 3*derv_sum_u*shgradl2(:,1)/denom_sum ...
    - 3*derv_sum_uu*shgradl(:,1)/denom_sum - derv_sum_uuu*R(:)/denom_sum ;

% wrt vvv
shgradl3(:,2) = shgradl3(:,2)/denom_sum - 3*derv_sum_v*shgradl2(:,2)/denom_sum ...
    - 3*derv_sum_vv*shgradl(:,2)/denom_sum - derv_sum_vvv*R(:)/denom_sum ;

% wrt uuv
shgradl3(:,3) = shgradl3(:,3)/denom_sum - derv_sum_v*shgradl2(:,1)/denom_sum ...
    - 2*derv_sum_u*shgradl2(:,3)/denom_sum - 2*derv_sum_uv*shgradl(:,1)/denom_sum ...
    - derv_sum_uu*shgradl(:,2)/denom_sum - derv_sum_uuv*R(:)/denom_sum ;

% wrt uvv
shgradl3(:,4) = shgradl3(:,4)/denom_sum - derv_sum_u*shgradl2(:,2)/denom_sum ...
    - 2*derv_sum_v*shgradl2(:,3)/denom_sum - 2*derv_sum_uv*shgradl(:,2)/denom_sum ...
    - derv_sum_vv*shgradl(:,1)/denom_sum - derv_sum_uvv*R(:)/denom_sum ;

% Fourth derivative....
% wrt uuuu
shgradl4(:,1) = shgradl4(:,1)/denom_sum - 4*derv_sum_u*shgradl3(:,1)/denom_sum ...
    - 6*derv_sum_uu*shgradl2(:,1)/denom_sum - 4*derv_sum_uuu*shgradl(:,1)/denom_sum ...
    - derv_sum_uuuu*R(:)/denom_sum ;

% wrt vvvv
shgradl4(:,2) = shgradl4(:,2)/denom_sum - 4*derv_sum_v*shgradl3(:,2)/denom_sum ...
    - 6*derv_sum_vv*shgradl2(:,2)/denom_sum - 4*derv_sum_uuu*shgradl(:,2)/denom_sum ...
    - derv_sum_vvvv*R(:)/denom_sum ;

% wrt uuuv
shgradl4(:,3) = shgradl4(:,3)/denom_sum - 3*derv_sum_u*shgradl3(:,3)/denom_sum ...
    - derv_sum_v*shgradl3(:,1)/denom_sum - 3*derv_sum_uu*shgradl2(:,3)/denom_sum ...
    - 3*derv_sum_uv*shgradl2(:,1)/denom_sum - 3*derv_sum_uuv*shgradl(:,1)/denom_sum ...
    - derv_sum_uuu*shgradl(:,2)/denom_sum - derv_sum_uuuv*R(:)/denom_sum ;

% wrt uuvv
shgradl4(:,4) = shgradl4(:,4)/denom_sum - 2*derv_sum_u*shgradl3(:,4)/denom_sum ...
    - 2*derv_sum_v*shgradl3(:,3)/denom_sum - 4*derv_sum_uv*shgradl2(:,3)/denom_sum ...
    - derv_sum_uu*shgradl2(:,2)/denom_sum - derv_sum_vv*shgradl2(:,1)/denom_sum ...
    - 2*derv_sum_uuv*shgradl(:,2)/denom_sum - 2*derv_sum_uvv*shgradl(:,1)/denom_sum ...
    - derv_sum_uuvv*R(:)/denom_sum ;

% wrt uvvv
shgradl4(:,5) = shgradl4(:,5)/denom_sum - 3*derv_sum_v*shgradl3(:,4)/denom_sum ...
    - derv_sum_u*shgradl3(:,2)/denom_sum - 3*derv_sum_vv*shgradl2(:,3)/denom_sum ...
    - 3*derv_sum_uv*shgradl2(:,2)/denom_sum - 3*derv_sum_uvv*shgradl(:,2)/denom_sum ...
    - derv_sum_vvv*shgradl(:,1)/denom_sum - derv_sum_uvvv*R(:)/denom_sum ;

% now calculate gradients.;
% calculate dx/dxi;

dxdxi=zeros(nsd,nsd);
icount = 0;
for j = 0: q
    for i = 0: p
        icount = icount + 1;
        dxdxi(1,1) = dxdxi(1,1) + b_net(ni-i,nj-j,1)*shgradl(icount,1);
        dxdxi(1,2) = dxdxi(1,2) + b_net(ni-i,nj-j,1)*shgradl(icount,2);
        dxdxi(2,1) = dxdxi(2,1) + b_net(ni-i,nj-j,2)*shgradl(icount,1);
        dxdxi(2,2) = dxdxi(2,2) + b_net(ni-i,nj-j,2)*shgradl(icount,2);
    end
end

% compute the inverse of deformation gradient and gradient of shapes in physical coordinates;
dxidx=inv(dxdxi);
shgradg = shgradl*dxidx ;

% Note that DetJ resides in common
detj = det(dxdxi);
%  if(detj < 0) % chu y doan code doi dau detj nay khi can thiet.
%  detj = -detj;
%  end
clear icount

% for higher order derivatives
% Second derivative
icount = 0 ;
jac12 = zeros(2,3) ;
jac22 = zeros(3,3) ;
for j = 0:q
    for i = 0:p
        icount = icount + 1 ;
        
        % term 1
        jac12(1,1) = jac12(1,1) + b_net(ni-i,nj-j,1)*shgradl2(icount,1) ;
        jac12(1,2) = jac12(1,2) + b_net(ni-i,nj-j,1)*shgradl2(icount,2) ;
        jac12(1,3) = jac12(1,3) + b_net(ni-i,nj-j,1)*shgradl2(icount,3) ;
        jac12(2,1) = jac12(2,1) + b_net(ni-i,nj-j,2)*shgradl2(icount,1) ;
        jac12(2,2) = jac12(2,2) + b_net(ni-i,nj-j,2)*shgradl2(icount,2) ;
        jac12(2,3) = jac12(2,3) + b_net(ni-i,nj-j,2)*shgradl2(icount,3) ;
        
    end
end

jac22 = [dxdxi(1,1)^2        dxdxi(1,2)^2        dxdxi(1,1)*dxdxi(1,2);...
    dxdxi(2,1)^2        dxdxi(2,2)^2        dxdxi(2,1)*dxdxi(2,2);...
    2*dxdxi(1,1)*dxdxi(2,1)     2*dxdxi(1,2)*dxdxi(2,2)     (dxdxi(2,1)*dxdxi(1,2)+dxdxi(1,1)*dxdxi(2,2))] ;

term12 = shgradg*jac12 ;
term22 = shgradl2 - term12 ;
shgradg2d = term22/jac22 ;

% Third derivative
icount = 0 ;
jac13 = zeros(2,4) ;
jac23 = zeros(3,4) ;
jac33 = zeros(4,4) ;
for j = 0:q
    for i = 0:p
        icount = icount + 1 ;
        
        % term 1
        jac13(1,1) = jac13(1,1) + b_net(ni-i,nj-j,1)*shgradl3(icount,1) ;
        jac13(1,2) = jac13(1,2) + b_net(ni-i,nj-j,1)*shgradl3(icount,2) ;
        jac13(1,3) = jac13(1,3) + b_net(ni-i,nj-j,1)*shgradl3(icount,3) ;
        jac13(1,4) = jac13(1,4) + b_net(ni-i,nj-j,1)*shgradl3(icount,4) ;
        jac13(2,1) = jac13(2,1) + b_net(ni-i,nj-j,2)*shgradl3(icount,1) ;
        jac13(2,2) = jac13(2,2) + b_net(ni-i,nj-j,2)*shgradl3(icount,2) ;
        jac13(2,3) = jac13(2,3) + b_net(ni-i,nj-j,2)*shgradl3(icount,3) ;
        jac13(2,4) = jac13(2,4) + b_net(ni-i,nj-j,2)*shgradl3(icount,4) ;
        
    end
end

% term 2
jac23 = [3*dxdxi(1,1)*jac12(1,1)    3*dxdxi(1,2)*jac12(1,2)...
    dxdxi(1,2)*jac12(1,1)+2*dxdxi(1,1)*jac12(1,3)   dxdxi(1,1)*jac12(1,2)+2*dxdxi(1,2)*jac12(1,3);...
    3*dxdxi(2,1)*jac12(2,1)    3*dxdxi(2,2)*jac12(2,2)...
    dxdxi(2,2)*jac12(2,1)+2*dxdxi(2,1)*jac12(2,3)   dxdxi(2,1)*jac12(2,2)+2*dxdxi(2,2)*jac12(2,3);...
    3*(dxdxi(1,1)*jac12(2,1)+dxdxi(2,1)*jac12(1,1))    3*(dxdxi(1,2)*jac12(2,2)+dxdxi(2,2)*jac12(1,2))     ...
    2*(dxdxi(1,1)*jac12(2,3)+dxdxi(2,1)*jac12(1,3))+dxdxi(2,2)*jac12(1,1)+dxdxi(1,2)*jac12(2,1)    ...
    2*(dxdxi(1,2)*jac12(2,3)+dxdxi(2,2)*jac12(1,3))+dxdxi(2,1)*jac12(1,2)+dxdxi(1,1)*jac12(2,2)] ;

% term 3
jac33 = [(dxdxi(1,1))^3     (dxdxi(1,2))^3      (dxdxi(1,1))^2*dxdxi(1,2)    dxdxi(1,1)*(dxdxi(1,2))^2;...
    (dxdxi(2,1))^3     (dxdxi(2,2))^3      (dxdxi(2,1))^2*dxdxi(2,2)    dxdxi(2,1)*(dxdxi(2,2))^2;...
    3*(dxdxi(1,1))^2*dxdxi(2,1)     3*(dxdxi(1,2))^2*dxdxi(2,2)      ...
    (dxdxi(1,1))^2*dxdxi(2,2)+2*dxdxi(1,1)*dxdxi(1,2)*dxdxi(2,1)     ...
    (dxdxi(1,2))^2*dxdxi(2,1)+2*dxdxi(1,1)*dxdxi(1,2)*dxdxi(2,2);...
    3*(dxdxi(2,1))^2*dxdxi(1,1)     3*(dxdxi(2,2))^2*dxdxi(1,2)      ...
    (dxdxi(2,1))^2*dxdxi(1,2)+2*dxdxi(1,1)*dxdxi(2,1)*dxdxi(2,2)     ...
    (dxdxi(2,2))^2*dxdxi(1,1)+2*dxdxi(1,2)*dxdxi(2,1)*dxdxi(2,2)] ;

term13 = shgradg*jac13+shgradg2d*jac23 ;
term23 = shgradl3 - term13 ;
shgradg3d = term23/jac33 ;

%Fourth derivative
icount = 0 ;
jac14 = zeros(2,5) ;
jac24 = zeros(3,5) ;
jac34 = zeros(4,5) ;
jac44 = zeros(5,5) ;
for j = 0:q
    for i = 0:p
        icount = icount + 1 ;
        
        % term 1
        jac14(1,1) = jac14(1,1) + b_net(ni-i,nj-j,1)*shgradl4(icount,1) ;
        jac14(1,2) = jac14(1,2) + b_net(ni-i,nj-j,1)*shgradl4(icount,2) ;
        jac14(1,3) = jac14(1,3) + b_net(ni-i,nj-j,1)*shgradl4(icount,3) ;
        jac14(1,4) = jac14(1,4) + b_net(ni-i,nj-j,1)*shgradl4(icount,4) ;
        jac14(1,5) = jac14(1,5) + b_net(ni-i,nj-j,1)*shgradl4(icount,5) ;
        jac14(2,1) = jac14(2,1) + b_net(ni-i,nj-j,2)*shgradl4(icount,1) ;
        jac14(2,2) = jac14(2,2) + b_net(ni-i,nj-j,2)*shgradl4(icount,2) ;
        jac14(2,3) = jac14(2,3) + b_net(ni-i,nj-j,2)*shgradl4(icount,3) ;
        jac14(2,4) = jac14(2,4) + b_net(ni-i,nj-j,2)*shgradl4(icount,4) ;
        jac14(2,5) = jac14(2,5) + b_net(ni-i,nj-j,2)*shgradl4(icount,5) ;
        
    end
end
% term 2
jac24 = [3*(jac12(1,1))^2 + 4*dxdxi(1,1)*jac13(1,1)   3*(jac12(1,2))^2 + 4*dxdxi(1,2)*jac13(1,2)  ...
    3*jac12(1,1)*jac12(1,3) + 3*dxdxi(1,1)*jac13(1,3) + dxdxi(1,2)*jac13(1,1)  ...
    2*(jac12(1,3))^2 + 2*dxdxi(1,1)*jac13(1,4) + 2*dxdxi(1,2)*jac13(1,3) + jac12(1,1)*jac12(1,2) ...
    3*jac12(1,3)*jac12(1,2) + 3*dxdxi(1,2)*jac13(1,4) + dxdxi(1,1)*jac13(1,2);...
    3*(jac12(2,1))^2 + 4*dxdxi(2,1)*jac13(2,1)   3*(jac12(2,2))^2 + 4*dxdxi(2,2)*jac13(2,2)  ...
    3*(jac12(2,1)*jac12(2,3)) + 3*dxdxi(2,1)*jac13(2,3) + dxdxi(2,2)*jac13(2,1)  ...
    2*(jac12(2,3))^2 + 2*dxdxi(2,1)*jac13(2,4) + 2*dxdxi(2,2)*jac13(2,3) + jac12(2,1)*jac12(2,2) ...
    3*jac12(2,3)*jac12(2,2) + 3*dxdxi(2,2)*jac13(2,4) + dxdxi(2,1)*jac13(2,2);...
    6*jac12(1,1)*jac12(2,1) + 4*dxdxi(1,1)*jac13(2,1) + 4*dxdxi(2,1)*jac13(1,1)  ...
    6*jac12(1,2)*jac12(2,2) + 4*dxdxi(1,2)*jac13(2,2) + 4*dxdxi(2,2)*jac13(1,2)  ...
    3*jac12(1,3)*jac12(2,1) + 3*jac12(1,1)*jac12(2,3) + 3*dxdxi(1,1)*jac13(2,3) + ...
    3*dxdxi(2,1)*jac13(1,3) + dxdxi(2,2)*jac13(1,1) + dxdxi(1,2)*jac13(2,1) ...
    2*dxdxi(2,1)*jac13(1,4) + 2*dxdxi(2,2)*jac13(1,3) + 2*dxdxi(1,1)*jac13(2,4) + ...
    2*dxdxi(1,2)*jac13(2,3) + 4*jac12(1,3)*jac12(2,3) + jac12(1,1)*jac12(2,2) + jac12(1,2)*jac12(2,1)  ...
    3*jac12(1,3)*jac12(2,2) + 3*jac12(1,2)*jac12(2,3) + 3*dxdxi(1,2)*jac13(2,4) + 3*dxdxi(2,2)*jac13(1,4) + ...
    dxdxi(1,1)*jac13(2,2) + dxdxi(2,1)*jac13(1,2)];
    
% term 3
jac34 = [6*(dxdxi(1,1))^2*jac12(1,1)  6*(dxdxi(1,2))^2*jac12(1,2) ...
    3*((dxdxi(1,1))^2*jac12(1,3) + dxdxi(1,1)*dxdxi(1,2)*jac12(1,1))  ...
    dxdxi(1,1)^2*jac12(1,2) + dxdxi(1,2)^2*jac12(1,1) + 4*dxdxi(1,1)*dxdxi(1,2)*jac12(1,3)  ...
    3*((dxdxi(1,2))^2*jac12(1,3) + dxdxi(1,1)*dxdxi(1,2)*jac12(1,2));  ...
    6*(dxdxi(2,1))^2*jac12(2,1)  6*(dxdxi(2,2))^2*jac12(2,2) ...
    3*((dxdxi(2,1))^2*jac12(2,3) + dxdxi(2,1)*dxdxi(2,2)*jac12(2,1))  ...
    dxdxi(2,1)^2*jac12(2,2) + dxdxi(2,2)^2*jac12(2,1) + 4*dxdxi(2,1)*dxdxi(2,2)*jac12(2,3)  ...
    3*((dxdxi(2,2))^2*jac12(2,3) + dxdxi(2,1)*dxdxi(2,2)*jac12(2,2));  ...
    3*(4*dxdxi(1,1)*dxdxi(2,1)*jac12(1,1) + 2*dxdxi(1,1)^2*jac12(2,1))  ...
    3*(4*dxdxi(1,2)*dxdxi(2,2)*jac12(1,2) + 2*dxdxi(1,2)^2*jac12(2,2))  ...
    3*(2*dxdxi(1,1)*dxdxi(2,1)*jac12(1,3) + dxdxi(1,1)^2*jac12(2,3) + dxdxi(1,1)*dxdxi(2,2)*jac12(1,1) + ...
    dxdxi(1,1)*dxdxi(1,2)*jac12(2,1) + dxdxi(1,2)*dxdxi(2,1)*jac12(1,1)) ...
    4*dxdxi(1,1)*dxdxi(2,2)*jac12(1,3) + 4*dxdxi(1,2)*dxdxi(2,1)*jac12(1,3) + 4*dxdxi(1,1)*dxdxi(1,2)*jac12(2,3) + ...
    2*dxdxi(1,1)*dxdxi(2,1)*jac12(1,2) + 2*dxdxi(1,2)*dxdxi(2,2)*jac12(1,1) + ...
    dxdxi(1,1)^2*jac12(2,2) + dxdxi(1,2)^2*jac12(2,1)   ...
    3*(2*dxdxi(1,2)*dxdxi(2,2)*jac12(1,3) + dxdxi(1,2)^2*jac12(2,3) + dxdxi(1,2)*dxdxi(2,1)*jac12(1,2) + ...
    dxdxi(1,1)*dxdxi(1,2)*jac12(2,2) + dxdxi(1,1)*dxdxi(2,2)*jac12(1,2)); ...
    3*(4*dxdxi(1,1)*dxdxi(2,1)*jac12(2,1) + 2*dxdxi(2,1)^2*jac12(1,1))  ...
    3*(4*dxdxi(1,2)*dxdxi(2,2)*jac12(2,2) + 2*dxdxi(2,2)^2*jac12(1,2))  ...
    3*(2*dxdxi(1,1)*dxdxi(2,1)*jac12(2,3) + dxdxi(2,1)^2*jac12(1,3) + dxdxi(1,2)*dxdxi(2,1)*jac12(2,1) + ...
    dxdxi(2,1)*dxdxi(2,2)*jac12(1,1) + dxdxi(1,1)*dxdxi(2,2)*jac12(2,1)) ...
    4*dxdxi(1,2)*dxdxi(2,1)*jac12(2,3) + 4*dxdxi(1,1)*dxdxi(2,2)*jac12(2,3) + 4*dxdxi(2,1)*dxdxi(2,2)*jac12(1,3) + ...
    2*dxdxi(1,1)*dxdxi(2,1)*jac12(2,2) + 2*dxdxi(1,2)*dxdxi(2,2)*jac12(2,1) + ...
    dxdxi(2,1)^2*jac12(1,2) + dxdxi(1,2)^2*jac12(1,1)   ...
    3*(2*dxdxi(1,2)*dxdxi(2,2)*jac12(2,3) + dxdxi(1,2)^2*jac12(1,3) + dxdxi(1,1)*dxdxi(2,2)*jac12(2,2) + ...
    dxdxi(2,1)*dxdxi(2,2)*jac12(1,2) + dxdxi(1,2)*dxdxi(2,1)*jac12(2,2))];
 
% term 4
jac44 = [dxdxi(1,1)^4  dxdxi(1,2)^4  dxdxi(1,1)^3*dxdxi(1,2)  dxdxi(1,1)^2*dxdxi(1,2)^2  dxdxi(1,1)*dxdxi(1,2)^3;...
    dxdxi(2,1)^4  dxdxi(2,2)^4  dxdxi(2,1)^3*dxdxi(2,2)  dxdxi(2,1)^2*dxdxi(2,2)^2  dxdxi(2,1)*dxdxi(2,2)^3;...
    4*dxdxi(1,1)^3*dxdxi(2,1)  4*dxdxi(1,2)^3*dxdxi(2,2)  3*dxdxi(1,1)^2*dxdxi(1,2)*dxdxi(2,1) + dxdxi(1,1)^3*dxdxi(2,2)  ...
    2*(dxdxi(1,1)^2*dxdxi(1,2)*dxdxi(2,2) + dxdxi(1,1)*dxdxi(1,2)^2*dxdxi(2,1))  ...
    3*dxdxi(1,1)*dxdxi(1,2)^2*dxdxi(2,2) + dxdxi(1,2)^3*dxdxi(2,1);...
    6*dxdxi(1,1)^2*dxdxi(2,1)^2   6*dxdxi(1,2)^2*dxdxi(2,2)^2    3*(dxdxi(1,1)^2*dxdxi(2,1)*dxdxi(2,2) + dxdxi(1,1)*dxdxi(1,2)*dxdxi(2,1)^2)   ...
    dxdxi(1,1)^2*dxdxi(2,2)^2 + dxdxi(1,2)^2*dxdxi(2,1)^2 + 4*dxdxi(1,1)*dxdxi(1,2)*dxdxi(2,1)*dxdxi(2,2)  ...
    3*(dxdxi(1,2)^2*dxdxi(2,1)*dxdxi(2,2) + dxdxi(1,1)*dxdxi(1,2)*dxdxi(2,2)^2);...
    4*dxdxi(2,1)^3*dxdxi(1,1)  4*dxdxi(2,2)^3*dxdxi(1,2)  3*dxdxi(2,1)^2*dxdxi(1,1)*dxdxi(2,2) + dxdxi(2,1)^3*dxdxi(1,2)  ...
    2*(dxdxi(2,1)^2*dxdxi(1,2)*dxdxi(2,2) + dxdxi(1,1)*dxdxi(2,2)^2*dxdxi(2,1))  ...
    3*dxdxi(1,2)*dxdxi(2,2)^2*dxdxi(2,1) + dxdxi(2,2)^3*dxdxi(1,1)];
    
term14 = shgradg*jac14 + shgradg2d*jac24 + shgradg3d*jac34;
term24 = shgradl4 - term14;
shgradg4d = term24/jac44;
end