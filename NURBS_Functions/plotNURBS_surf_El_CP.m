function plotNURBS_surf_El_CP(p,q,U,V,CP)
% plots the surface, elements and control points

[X,Y,Z] = create_surf(p,q,U,V,CP);
% geometry
%  surf(X/L,Y/L,Z/max(max(abs(Z))),'FaceColor','blue','EdgeColor','none');
 surf(X,Y,Z,'FaceColor','blue','EdgeColor','none');
xlabel('x'); ylabel('y'); zlabel('z');
hold on;
% CP(:,:,1)=[0 0;1 1];
% CP(:,:,2)=[0 1;0 1];
% CP(:,:,3)=[0 0;0 0];
% CP(:,:,4)=[1 1;1 1];

% [CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);
% element edges
create_el_edges(p,q,U,V,CP)

axis equal;

hold off;
