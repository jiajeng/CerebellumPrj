function [matlabbatch,procsNum] = secndLevlresult(filepath,matlabbatch,varargin)
    in = finputcheck(varargin, ...
           {'con'       'string'    []  '';
            'procsNum'  'real'      []  [];
            'threstype' 'string'    []  '';
            'thresh'    'real'      []  [];
           });
    if class(in) ~= "struct", error(in); end
    procsNum = in.procsNum;
    load(fullfile(filepath,'SPM.mat'));
    contrast = {SPM.xCon.name};
    if ~isempty(in.con)
        c = find(cellfun(@(x) string(in.con)==string(x),contrast));
    else
        c = 1:length(contrast);
    end
    for i = 1:length(c)
        if i ~= 1, procsNum = procsNum+1;end
        matlabbatch{procsNum}.spm.stats.results.spmmat = {fullfile(filepath,'SPM.mat')};
        matlabbatch{procsNum}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{procsNum}.spm.stats.results.conspec.contrasts = c(i);
        matlabbatch{procsNum}.spm.stats.results.conspec.threshdesc = in.threstype;
        matlabbatch{procsNum}.spm.stats.results.conspec.thresh = in.thresh;
        matlabbatch{procsNum}.spm.stats.results.conspec.extent = 0;
        matlabbatch{procsNum}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{procsNum}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{procsNum}.spm.stats.results.units = 1;
        matlabbatch{procsNum}.spm.stats.results.export{1}.jpg = true;
    end
end