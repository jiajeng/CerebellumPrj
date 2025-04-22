function func_step3(T1filepath,T1filename,funcfilepath,funcfilename,batchfolder,f)
    % reslice to structure(cerebullum area) 'c_name_pcereb.nii'
    % smooth functional data
    % T1filepath --> "string", T1 file 'c_name_pcereb.nii' path
    % T1filename --> "string" , 'c_*_pcereb.nii'
    % funcfilepath --> "string", reslice functional data path
    % funcfilename --> "string", function data name can use *, e.g. con*.nii
    %                           if empty then get T_cname.nii
    % batchfolder --> "string", save batch file path
    matlabbatch = [];

    % T1 file
    T1file = ls(fullfile(T1filepath,T1filename));
    Steps.checkfile(T1file,T1filepath,T1filename);
    % func file
    if isempty(funcfilename)
        funcfilename = 'T_*.nii';
    end
    FILE = dir(fullfile(funcfilepath,funcfilename));
    funcfile = cell(1,size(FILE,1));
    for nfile = 1:size(FILE,1)
        funcfile{nfile} = fullfile(FILE(nfile).name);
    end

    % reslice to structure SUIT
    % ========================================================
    if f.func.resliceSUIT
        paramfile = char({dir(fullfile(T1filepath,'mc_*_snc.mat')).name}); 
        for nfile = 1:length(funcfile)
            suit_reslice(fullfile(funcfilepath,funcfile{nfile}),fullfile(T1filepath,paramfile),'mask',fullfile(T1filepath,T1file));
        end

        % smooth normalized functional data
        % ========================================================
        if f.func.smooth
            sfuncfile = ls(fullfile(funcfilepath,'w*'));
            for nfile = 1:size(sfuncfile,1)
                Steps.checkfile(sfuncfile(nfile,:),'w*',funcfilepath)
            end
            for nfile = 1:size(sfuncfile,1)
                file = sfuncfile(nfile,:);
                file(file==' ') = [];
                matlabbatch{end+1}.spm.spatial.smooth.data = {[fullfile(funcfilepath,file),',1']};
                matlabbatch{end}.spm.spatial.smooth.fwhm = [4 4 4];
                matlabbatch{end}.spm.spatial.smooth.dtype = 0;
                matlabbatch{end}.spm.spatial.smooth.im = 0;
                matlabbatch{end}.spm.spatial.smooth.prefix = 's';
            end
        end
    end

    if f.func.resliceMNI
        MNIfuncfile = cell(length(funcfile),1);
        for nfile = 1:length(funcfile)
            MNIfuncfile{nfile} = [fullfile(funcfilepath,funcfile{nfile}),',1'];
        end 
        deformField = {dir(fullfile(T1filepath,'y_*.nii')).name};
        deformField = deformField{~contains(deformField,'c0')};
        
        matlabbatch{end+1}.spm.spatial.normalise.write.subj.def = {fullfile(T1filepath,deformField)};
        matlabbatch{end}.spm.spatial.normalise.write.subj.resample = MNIfuncfile;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                  78 76 85];
        matlabbatch{end}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{end}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix = 'w';

        matlabbatch{end+1}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{end}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{end}.spm.spatial.smooth.dtype = 0;
        matlabbatch{end}.spm.spatial.smooth.im = 0;
        matlabbatch{end}.spm.spatial.smooth.prefix = 's';

    end
    % ========================================================
    
    % ========================================================
    if ~isempty(matlabbatch)
        % run job
        % ==================================================
        spm('defaults', 'FMRI');
        spm_jobman('initcfg')
        spm_jobman('run',matlabbatch);
        % ==================================================

        % save batch
        % ===================================================
        if ~exist(batchfolder,"dir"), mkdir(batchfolder); end
        save(fullfile(batchfolder,'func3_smooth.mat'),"matlabbatch");
        % ===================================================
    end
end
