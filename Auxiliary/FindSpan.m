function ind = FindSpan(n, p, u, U)
% ind = FindSpan(n, p, u, U)
%----------------------------------------------------------
% Determine the knot span index
% ie. find i if u(i) lies in the half-open interval [ u(i), u(i + 1) ), i = 1, 2,..., n.
%----------------------------------------------------------
% Input(s):
%    n: number of control points (or basis functions)
%    p: order (degree) of Bspline basis functions
%    u: parametric points
%    U: clamped knot vector (U = [u_{1} = ... = u_{p + 1} < u_{p + 2} <= ... <= u_{n} < u_{n + 1} = ... = u_{n + p + 1}])
%----------------------------------------------------------
% Output(s):
%    ind: knot span indices
%----------------------------------------------------------
% Based on Algorithm A2.1 [The NURBS BOOK, p.68]
% Written by Khanh Chau-Nguyen
%----------------------------------------------------------

% assert(isinteger(p), 'Integer type is required for degree of basis functions')
% assert(isinteger(n), 'Integer type is required for number of basis functions')

% p = uint32(p);
% n = uint32(n);
ind = zeros(size(u));
for i = 1 : numel(u)
    if (abs(u(i) - U(n + 1)) <= 1e-12)
        ind(i) = n;
        continue
    end
    low = p + 1;
    high = n + 1;
    mid = floor((low + high)/2);
%     mid = idivide((low + high), 2);
    % while ( (U(mid) - u(i)) > 1e-14 || (u(i) - U(mid + 1) + 1e-14) >= 0 )
    while ( u(i) < U(mid)|| u(i) >= U(mid + 1) )
        if (u(i) < U(mid))
            high = mid;
        else
            low = mid;
        end
%         mid = idivide((low + high), 2);
        mid = floor((low + high) / 2);
    end
    ind(i) = mid;
end
end