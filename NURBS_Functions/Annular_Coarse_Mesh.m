function [CP,U,V,p,q] = Annular_Coarse_Mesh(R,r,deg)
p = 2;
q = 2;
U = [0 0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1 1]; %m=9
V = [0 0 0 1 1 1]; %n=2
%Control Point coordinates
CP(:,:,1)=[R (R+r)/2 r;R (R+r)/2 r;0 0 0;-R -(R+r)/2 -r;-R -(R+r)/2 -r;-R -(R+r)/2 -r;0 0 0;R (R+r)/2 r;R (R+r)/2 r];
CP(:,:,2)=[0 0 0;R (R+r)/2 r;R (R+r)/2 r;R (R+r)/2 r;0 0 0;-R -(R+r)/2 -r;-R -(R+r)/2 -r;-R -(R+r)/2 -r;0 0 0];
CP(:,:,3)=[1 1 1;1 1 1;1 1 1;1 1 1 ;1 1 1;1 1 1;1 1 1;1 1 1 ;1 1 1]*5;
%CP(:,:,4)=[ 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2);1 1];
CP(:,:,4)=[ 1 1 1; 1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1];

%=====================================================================
% REFINE
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);


% R1 = refinement_vec_repeated_p2(U,ref);
% R2 = refinement_vec_repeated_p2(V,ref);
% 
% [CP,u_knot,v_knot] = knot_refine_surf(p,q,U,V,CP,R1,R2);
% 
% plotNURBS_surf_El_CP(p,q,u_knot,v_knot,CP); hold on
% plot_ctrlnet(CP);
% view(2)