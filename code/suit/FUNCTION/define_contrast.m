function define_contrast(filepath,f,varargin)
    % procsNum --> matlabbatch process number, 
    %              [int]
    % contrast --> condition name to set to contrast, only enter "wanted" condition 
    %              this script will get position based the SPM.mat.
    %              ex.  {{condi1,condi2},{condi2,condi3}}
    % weight --> same length with contrast mean this map contrast weight.
    %            ex. {{1,-1},{1,-1}}
    % conName --> same length with contrast, set this contrast name
    %             {char}
    % condi --> set every conditon to a contrast, like condition 1 is 1 other
    %           is 0. if set 1 all condition will set 1 and -1. 
    %           [1,0]
    
    in = finputcheck(varargin, ...
           {'contrast'      'cell'      []      {}; % condition name to set to contrast
            'weight'        'cell'      []      {}; 
            'condi'         'real'      []      0; % every condition is a contrast, 1 yes, 0 no
            'conName'       'cell'      []      {};
            'batchfolder'   'string'    []      [];
           });
    matlabbatch = [];
    load(fullfile(filepath,'SPM.mat'));
    
    SPMCOND = SPM.xX.name;  
    condi = split(SPMCOND',' ');
    condi = condi(:,end)';
    condi = ~(contains(condi,'rp')|contains(condi,'^'));
    tmp = SPMCOND(condi);
    CONNAME = tmp;
    tmp = tmp(contains(tmp,' '));
    tmp = split(tmp',' ');
    tmp = tmp(:,2);
    CONNAME(~contains(tmp,'*')) = tmp(~contains(tmp,'*'));
    tmp = tmp(contains(tmp,'*'));
    tmp = split(tmp,'*');
    tmp = tmp(:,1)';
    CONNAME(contains(CONNAME,'*')) = tmp;
    matlabbatch{end+1}.spm.stats.con.spmmat = {fullfile(filepath,'SPM.mat')};
    condiN = 0;
    
    
    if f.func.fstconP
        for ncond = 1:length(CONNAME)
            wpos = find(condi==1);
            wpos = wpos(contains(SPMCOND(condi),CONNAME{ncond}));
            weight = zeros(size(CONNAME)); weight(wpos) = 1;
            conName = CONNAME{ncond};
            if in.condi
                matlabbatch{end}.spm.stats.con.consess{ncond}.tcon.name = ['pos_',conName];
                matlabbatch{end}.spm.stats.con.consess{ncond}.tcon.weights = weight;
                matlabbatch{end}.spm.stats.con.consess{ncond}.tcon.sessrep = 'none';
                condiN = 1+condiN;
            end
        end
    end

    if f.func.fstconN
        ocondiN = condiN;
        for ncond = 1:length(CONNAME)
            weight = zeros(size(CONNAME)); weight(ncond) = -1;
            condName = CONNAME{ncond+rpL};
            conName{ncond} = condName;
            if in.condi
                matlabbatch{end}.spm.stats.con.consess{ncond+ocondiN}.tcon.name = ['neg_',condName];
                matlabbatch{end}.spm.stats.con.consess{ncond+ocondiN}.tcon.weights = weight;
                matlabbatch{end}.spm.stats.con.consess{ncond+ocondiN}.tcon.sessrep = 'none';
                condiN = 1+condiN;
            end
        end
    end
    
    if ~isempty(in.contrast)
        conName = string(CONNAME);
        for ncond = 1:length(in.contrast)
            if length(in.contrast{ncond})~=length(in.weight{ncond})
                error('one contrast map to one weight');
            end 
            matlabbatch{end}.spm.stats.con.consess{ncond+condiN}.tcon.name = char(in.conName{ncond});
            weight = zeros(size(CONNAME));
            weightPos = arrayfun(@(x) find(conName == x),string(in.contrast{ncond}),'UniformOutput',false);
    
            for i = 1:length(weightPos)
                weight(weightPos{i}) = in.weight{ncond}{i};
            end
            
            matlabbatch{end}.spm.stats.con.consess{ncond+condiN}.tcon.weights = weight;
            matlabbatch{end}.spm.stats.con.consess{ncond+condiN}.tcon.sessrep = 'none';
        end
    end
    matlabbatch{end}.spm.stats.con.delete = 0;
    
    % run job
    % ==================================================
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch);
    % ==================================================

    % save batch
    % ===================================================
    if ~isempty(in.batchfolder)
        if ~exist(in.batchfolder,"dir"), mkdir(in.batchfolder); end
        save(fullfile(in.batchfolder,'func2_def_con.mat'),"matlabbatch");
    end
    % ===================================================
end
