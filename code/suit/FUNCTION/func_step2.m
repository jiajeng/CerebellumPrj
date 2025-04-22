function func_step2(filepath,filename,outputFolder,TR,units,batchfolder,condition,con,conWt,conName,condi,f,varargin)
    % 1st level analysis
    % filepath --> file direction
    % filename --> filename can use * filter
    % outputFolder --> put SPM.mat file direction
    % TR --> TR
    % units --> 'secs' or 'scans'
    % batchfolder --> put batch.mat file, can open in spm batch later
    % condition --> task condition, format 
    %                                   struct('name', {}, 
    %                                          'onset', {}, 
    %                                          'duration', {}, 
    %                                          'tmod', {}, 
    %                                          'pmod', {}, 
    %                                          'orth', {});
    % con --> cell cell array, example {{con1,con2},{con1}}
    % conWt --> cell cell array, example {{1,-1},{1}}
    % conName --> this contrast name, same length as con and conWt 
    % condi --> automatically set all condition contrast 1
    varName = varargin(1:2:end);
    varValue = varargin(2:2:end);
    phyfilepath = varValue(contains(varName,'phyfilepath'));

    matlabbatch = [];

    % get file name and check
    % =================================================
    file = ls(fullfile(filepath,filename));
    Steps.checkfile(file,filepath,filename)

    % get image size
    info = niftiinfo(fullfile(filepath,file));
    ImgSize = info.ImageSize;
    
    % get file frames
    frames = cell(ImgSize(4),1);
    for i = 1:ImgSize(4)
        frames{i} = [fullfile(filepath,file),',',num2str(i)];
    end

    % =================================================
    
    if ~exist(outputFolder,'dir'), mkdir(outputFolder); end

    % 1st Level analysis
    % =================================================
    if f.func.fstL
        matlabbatch{end+1}.spm.stats.fmri_spec.dir = {outputFolder};
        matlabbatch{end}.spm.stats.fmri_spec.timing.units = units;
        matlabbatch{end}.spm.stats.fmri_spec.timing.RT = TR;      
        matlabbatch{end}.spm.stats.fmri_spec.sess.scans = frames;

        % Get realingment regressor
        rpfile = ls(fullfile(filepath,'rp*.txt'));
        rp = readtable(fullfile(filepath,rpfile));
        vn = convertStringsToChars(string([repmat('rp',6,1),num2str([1:6]')]));
        rp.Properties.VariableNames = vn;
        for i = 1:size(rp,2)
            matlabbatch{end}.spm.stats.fmri_spec.sess.regress(i).name = vn{i};
            matlabbatch{end}.spm.stats.fmri_spec.sess.regress(i).val = rp.(vn{i});
        end
        
        matlabbatch{end}.spm.stats.fmri_spec.sess.cond = condition;
        matlabbatch{end+1}.spm.stats.fmri_est.spmmat = {fullfile(outputFolder,'SPM.mat')};
    end
    % ==================================================

    if ~isempty(matlabbatch)
        % run job
        % ==================================================
        spm('defaults', 'FMRI');
        spm_jobman('initcfg')
        spm_jobman('run',matlabbatch);
        % ========================================================
        
        % save batch
        % ===================================================
        if ~exist(batchfolder,"dir"), mkdir(batchfolder); end
        save(fullfile(batchfolder,'func2_1st_level.mat'),"matlabbatch");
        % ===================================================
    end

    if f.func.fstcon
        % define contrast for 1st level analysis
        % ===================================================
        Steps.define_contrast(outputFolder,f, ...
            'contrast',con, ...
            'weight',conWt, ...
            'conName',conName, ...
            'condi',condi, ...
            'batchfolder',batchfolder);
        % ===================================================
    end

end
