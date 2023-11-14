function FEM = Bcdof_CCCC_7dof(FEM,ndof,mcp, ncp)

bcdof=[];
bcval=[];
% u, v, w, betax, betay
right1=[mcp:mcp:mcp*ncp];               % on perimeter
right2=[mcp-1:mcp:mcp*ncp-1];           % inside

left1=[1:mcp:mcp*ncp-mcp+1];            % on perimeter
left2=[2:mcp:mcp*ncp-mcp+2];            % inside

lower1=[1:mcp];                         % on perimeter
lower2=[mcp+1:2*mcp];                   % inside

upper1=[mcp*ncp-mcp+1:mcp*ncp];         % on perimeter
upper2=[mcp*ncp-2*mcp+1:mcp*ncp-mcp];   % inside

% right 
bcdof = [bcdof ndof*right1-6];     % u0
bcdof = [bcdof ndof*right1-5];     % v0
bcdof = [bcdof ndof*right1-4];     % w
% bcdof = [bcdof ndof*right2-4];     % w
bcdof = [bcdof ndof*right1-3];     % 
bcdof = [bcdof ndof*right1-2];     % 
bcdof = [bcdof ndof*right1-1];     % betax
bcdof = [bcdof ndof*right1];       % betay

% left 
bcdof = [bcdof ndof*left1-6];     % u0
bcdof = [bcdof ndof*left1-5];     % v0
bcdof = [bcdof ndof*left1-4];     % w
% bcdof = [bcdof ndof*left2-4];     % w
bcdof = [bcdof ndof*left1-3];     % 
bcdof = [bcdof ndof*left1-2];     % 
bcdof = [bcdof ndof*left1-1];     % betax
bcdof = [bcdof ndof*left1];       % betay

% upper 
bcdof = [bcdof ndof*upper1-6];     % u0
bcdof = [bcdof ndof*upper1-5];     % v0
bcdof = [bcdof ndof*upper1-4];     % w
% bcdof = [bcdof ndof*upper2-4];     % w
bcdof = [bcdof ndof*upper1-3];     % 
bcdof = [bcdof ndof*upper1-2];     % 
bcdof = [bcdof ndof*upper1-1];     % betax
bcdof = [bcdof ndof*upper1];       % betay

% lower 
bcdof = [bcdof ndof*lower1-6];     % u0
bcdof = [bcdof ndof*lower1-5];     % v0
bcdof = [bcdof ndof*lower1-4];     % w
% bcdof = [bcdof ndof*lower2-4];     % w
bcdof = [bcdof ndof*lower1-3];     % 
bcdof = [bcdof ndof*lower1-2];     % 
bcdof = [bcdof ndof*lower1-1];     % betax
bcdof = [bcdof ndof*lower1];       % betay

bcval = [bcval, zeros(1,length(bcdof))];

FEM.BCDof = bcdof;
FEM.BCVal = bcval;
