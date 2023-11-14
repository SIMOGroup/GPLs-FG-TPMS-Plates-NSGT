function FEM = Bcdof_SSSS_7dof(FEM,ndof,mcp, ncp)

bcdof=[];
bcval=[];
% u, v, w, betax, betay, thetax, thetay
right1=[mcp:mcp:mcp*ncp];               % on perimeter
right2=[mcp-1:mcp:mcp*ncp-1];           % inside

left1=[1:mcp:mcp*ncp-mcp+1];            % on perimeter
left2=[2:mcp:mcp*ncp-mcp+2];            % inside

lower1=[1:mcp];                         % on perimeter
lower2=[mcp+1:2*mcp];                   % inside

upper1=[mcp*ncp-mcp+1:mcp*ncp];         % on perimeter
upper2=[mcp*ncp-2*mcp+1:mcp*ncp-mcp];   % inside

% right 
bcdof = [bcdof ndof*right1-5];     % v0
bcdof = [bcdof ndof*right1-4];     % W
bcdof = [bcdof ndof*right1-2];       % theta_y
bcdof = [bcdof ndof*right1-0];       % beta_y

% left 
bcdof = [bcdof ndof*left1-5];     % v0
bcdof = [bcdof ndof*left1-4];     % W
bcdof = [bcdof ndof*left1-2];       % theta_y
bcdof = [bcdof ndof*left1-0];       % beta_y

% upper 
bcdof = [bcdof ndof*upper1-6];     % u0
bcdof = [bcdof ndof*upper1-4];     % W
bcdof = [bcdof ndof*upper1-3];     % betax
bcdof = [bcdof ndof*upper1-1];     % betax

% lower 
bcdof = [bcdof ndof*lower1-6];     % u0
bcdof = [bcdof ndof*lower1-4];     % W
bcdof = [bcdof ndof*lower1-3];     % betax
bcdof = [bcdof ndof*lower1-1];     % betax

bcval = [bcval, zeros(1,length(bcdof))];

FEM.BCDof = bcdof;
FEM.BCVal = bcval;
return