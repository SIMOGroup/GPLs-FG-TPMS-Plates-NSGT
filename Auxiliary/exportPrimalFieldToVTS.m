function exportPrimalFieldToVTS(Geom, d, vtkPtsPerElement, filename, fieldname)
% exportPrimalFieldToVTS(Geom, d, DOFsMapper, ParaPts, filename, fieldname)

if Geom.NumPatch > 1 % multiple patches problem
    str1 = cat (2,'<?xml version="1.0"?> \n', '<VTKFile type="Collection" version="0.1"> \n', '<Collection> \n');
    
    str2 = cat (2, '<DataSet part="%d" file="%s.vts"/> \n');
    
    str3 = cat (2, '</Collection>\n', '</VTKFile> \n');
    
    if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
        pvd_filename = cat (2, filename, '.pvd');
    else
        pvd_filename = filename;
        filename = filename (1:end-4);
    end
    
    fid = fopen (pvd_filename, 'w');
    if (fid < 0)
        error ('exportPrimalFieldToVTS: could not open file %s', pvd_filename);
    end
    
    fprintf (fid, str1);
    for iptc = 1 : Geom.NumPatch
        NURBS = Geom.Patch{iptc};
        filename_patch = cat(2, filename, '_', num2str(iptc));
        fprintf (fid, str2, iptc, filename_patch);
        if Geom.NumBaseDOFs == numel(d) % scalar field
            exportFieldToVTSPatch(NURBS, d(Geom.DOFsMapper{iptc}), vtkPtsPerElement, filename_patch, fieldname);
        else
            if NURBS.NumSpaceDim == 2
                exportFieldToVTSPatch(NURBS, d([Geom.DOFsMapper{iptc}; Geom.DOFsMapper{iptc} + Geom.NumBaseDOFs]), vtkPtsPerElement, filename_patch, fieldname);
            elseif NURBS.NumSpaceDim == 3
                exportFieldToVTSPatch(NURBS, d([Geom.DOFsMapper{iptc}; Geom.DOFsMapper{iptc} + Geom.NumBaseDOFs; Geom.DOFsMapper{iptc} + 2 * Geom.NumBaseDOFs]), vtkPtsPerElement, filename_patch, fieldname);
            else
                error('not implemented yet!')
            end
        end
    end
    fprintf (fid, str3);
    
    fclose (fid);
else
    exportFieldToVTSPatch(Geom.Patch{1}, d, vtkPtsPerElement, filename, fieldname);
end
end