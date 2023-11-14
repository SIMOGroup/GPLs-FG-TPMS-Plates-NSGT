function [CP,U,V,p,q]=Circle_Coarse_Mesh(R,deg)
% circular shape
p = 2;
q = 2;
U = [0 0 0 1 1 1]; % m=3
V = [0 0 0 1 1 1]; % n=3
theta = 90 * pi / 180;

CP(:, :, 1) = R * [-sqrt(2)/2 -sqrt(2)/1 -sqrt(2)/2; 0 0 0; sqrt(2)/2 sqrt(2)/1 sqrt(2)/2];
CP(:, :, 2) = R * [-sqrt(2)/2 0 sqrt(2)/2; -sqrt(2)/1 0 sqrt(2)/1; -sqrt(2)/2 0 sqrt(2)/2];
CP(:, :, 3) = [0 0 0; 0 0 0; 0 0 0];
CP(:, :, 4) = [1 cos(theta/2) 1; cos(theta/2) 1 cos(theta/2); 1 cos(theta/2) 1];

[CP, U, V, p, q] = degree_elevate_surf(p, q, U, V, CP, deg-p,deg-q);
% [CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);
return
p=1;
q=1;
U=[0 0 1 1];
V=[0 0 1 1];
%Control Point coordinates
CP(:,:,1)=[0 0;a a];
CP(:,:,2)=[0 b;0 b];
CP(:,:,3)=[0 0;0 0];
CP(:,:,4)=[1 1;1 1];

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