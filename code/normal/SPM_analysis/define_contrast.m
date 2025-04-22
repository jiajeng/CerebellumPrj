function [matlabbatch,procsNum] = define_contrast(filepath,matlabbatch,varargin)
in = finputcheck(varargin, ...
       {'procsNum'      'real'      []      [];
        'contrastName'  'cell'      []      {};
        'contrast'      'cell'      []      {};
        'weight'        'cell'      []      {};
        'condi'         'real'      []      0;
        'conName'       'cell'      []      {};
       });
procsNum = in.procsNum;

load(fullfile(filepath,'SPM.mat'));

condition = SPM.xX.name;
matlabbatch{procsNum}.spm.stats.con.spmmat = {fullfile(filepath,'SPM.mat')};
conName = cell(size(condition));

for ncond = 1:length(condition)
    weight = zeros(size(condition)); weight(ncond) = 1;
    condName = condition{ncond};
    condName = split(condName,'*'); condName = condName{1};
    condName = split(condName,' '); condName = condName{2};
    conName{ncond} = condName;
    if in.condi
        matlabbatch{procsNum}.spm.stats.con.consess{ncond}.tcon.name = condName;
        matlabbatch{procsNum}.spm.stats.con.consess{ncond}.tcon.weights = weight;
        matlabbatch{procsNum}.spm.stats.con.consess{ncond}.tcon.sessrep = 'none';
        condiN = length(condition);
    else
        condiN = 0;
    end
end


if ~isempty(in.contrast)
    conName = string(conName);
    for ncond = 1:length(in.contrast)
        if length(in.contrast{ncond})~=length(in.weight{ncond})
            error('one contrast map to one weight');
        end 
        matlabbatch{procsNum}.spm.stats.con.consess{ncond+condiN}.tcon.name = char(in.conName{ncond});
        weightPos = arrayfun(@(x) find(conName == x),string(in.contrast{ncond}));
        weight = zeros(size(condition));
        for i = 1:length(weightPos)
            weight(weightPos(i)) = in.weight{ncond}{i};
        end
        
        matlabbatch{procsNum}.spm.stats.con.consess{ncond+condiN}.tcon.weights = weight;
        matlabbatch{procsNum}.spm.stats.con.consess{ncond+condiN}.tcon.sessrep = 'none';
    end
end
matlabbatch{procsNum}.spm.stats.con.delete = 0;

end