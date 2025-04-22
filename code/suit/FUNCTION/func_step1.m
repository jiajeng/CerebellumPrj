function func_step1(T1path,T1filename,filepath,filename,TR,sliceorder,batchfolder,f)
    % preprocessing
    %   realignment
    %   slice timing
    %   coregister to anatomical
    % T1 path --> T1 data filepath
    % T1 filename --> T1 data filename can use * filter
    % filepath --> functional data filepath 
    % filename --> functional data filename can use * filter
    % TR --> TR can get in MRinfo.xlsx
    % sliceorder --> sliceorder can get in MRinfo.xlsx
    % batchfolder --> put batch.mat to open in spm later
    
    matlabbatch = [];

    % realignment
    % ========================================================
    % get file name
    file = ls(fullfile(filepath,filename));
    Steps.checkfile(file,filepath,filename)

    mfile = ['mean',file]; % for corigester
    % get image size
    info = niftiinfo(fullfile(filepath,file));
    ImgSize = info.ImageSize;

    % get file frames
    frames = cell(ImgSize(4),1);
    for i = 1:ImgSize(4)
        frames{i} = [fullfile(filepath,file),',',num2str(i)];
    end

    % define job batch
    if f.func.relign
        matlabbatch{end+1}.spm.spatial.realign.estwrite.data = {frames};
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{end}.spm.spatial.realign.estwrite.roptions.prefix = 'u';
    end
    % ========================================================

    % slice Timing 
    % =================================================
    % get file name
    file = ['u',file];

    % get file frames
    frames = cell(ImgSize(4),1);
    for i = 1:ImgSize(4)
        frames{i} = [fullfile(filepath,file),',',num2str(i)];
    end

    % define job parameter
    if f.func.sliTim
        matlabbatch{end+1}.spm.temporal.st.scans = {frames};
        matlabbatch{end}.spm.temporal.st.nslices = ImgSize(3);
        matlabbatch{end}.spm.temporal.st.tr = TR;
        matlabbatch{end}.spm.temporal.st.ta = TR-TR/ImgSize(3);
        matlabbatch{end}.spm.temporal.st.so = sliceorder;
        matlabbatch{end}.spm.temporal.st.refslice = 1;
        matlabbatch{end}.spm.temporal.st.prefix = 'a';
    end
    % ==============================================

    % corigister to anatomical
    % ===============================================
    
    T1file = ls(fullfile(T1path,T1filename));
    Steps.checkfile(T1file,T1path,T1filename);
    cofile = ['a',file];

    % get file frames
    frames = cell(ImgSize(4),1);
    for i = 1:ImgSize(4)
        frames{i} = [fullfile(filepath,cofile),',',num2str(i)];
    end

    if f.func.corig
        matlabbatch{end+1}.spm.spatial.coreg.estimate.ref = {[fullfile(T1path,T1file),',1']};
        matlabbatch{end}.spm.spatial.coreg.estimate.source = {[fullfile(filepath,mfile),',1']};
        matlabbatch{end}.spm.spatial.coreg.estimate.other = frames;
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    end
    % ===============================================

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
        save(fullfile(batchfolder,'func1_au.mat'),"matlabbatch");
        % =================================================== 
    end
end
