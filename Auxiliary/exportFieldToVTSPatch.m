function exportFieldToVTSPatch(NURBS, d, vtkPtsPerElement, filename, fieldname)
% exportFieldToVTSPatch(NURBS, d, vtkPtsPerElement, filename, fieldname)
% --------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       d: temperature or displacement field
%       ParaPts: parameter points per direction
%       filename: name of file
%       fieldname: name of field
% ------------------------------------------------------------------
% Output:
%       filename.vts
% ------------------------------------------------------------------

ParaPts = cell(1, NURBS.Dim);
for dir = 1 : NURBS.Dim
    dxi = diff(NURBS.uqKntVect{dir}) / vtkPtsPerElement(dir);
    knts = NURBS.uqKntVect{dir}(1 : end - 1);
    xi = repmat(knts, vtkPtsPerElement(dir), 1);
    step = 1 : vtkPtsPerElement(dir) - 1;
    xi(2 : end, :) = xi(2 : end, :) + dxi .* step';
    xi = xi(:);
    ParaPts{dir} = [xi', 1];
end
[C, F] = NURBSPatchEval(NURBS, ParaPts, d);
exportToVTK(C, F, filename, fieldname);
end