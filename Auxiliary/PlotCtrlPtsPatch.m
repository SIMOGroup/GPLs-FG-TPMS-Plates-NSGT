function PlotCtrlPtsPatch(NURBS, cflag, varargin)
% PlotCtrlPts(NURBS, varargin)
% Plot control points
% ---------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       varargin (optional): if varargin = 1, then plot control point numbering
% ---------------------------------------------------------------------

CtrlPts3D = reshape(NURBS.CtrlPts3D, 3, []);

plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :), '.', 'MarkerSize', 17, 'color', cflag)
% return
if (nargin >= 3)
    option = varargin{1};    
    if NURBS.Dim == 2
        if(option)
            t = 1;
            for j = 1 : NURBS.NumCPsDir(2)
                for i = 1 : NURBS.NumCPsDir(1)
                    if nargin == 3
                        txt = t;
                    elseif nargin == 4
                        txt = varargin{2}(t);
                    end
                    text(NURBS.CtrlPts3D(1, i, j), NURBS.CtrlPts3D(2, i, j), NURBS.CtrlPts3D(3, i, j), num2str(txt), 'FontWeight', 'bold');
                    t = t + 1;
                end
            end
        end
    elseif NURBS.Dim == 3
        if(option)
            t = 1;
            for k = 1 : NURBS.NumCPsDir(3)
                for j = 1 : NURBS.NumCPsDir(2)
                    for i = 1 : NURBS.NumCPsDir(1)
                        text(NURBS.CtrlPts3D(1, i, j, k), NURBS.CtrlPts3D(2, i, j, k), NURBS.CtrlPts3D(3, i, j, k), num2str(t), 'FontWeight', 'bold');
                        t = t + 1;
                    end
                end
            end
        end
    end
end
end