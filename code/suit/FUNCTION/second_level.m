function second_level(filepath,round,outputpath,f,varargin)
    in = finputcheck(varargin, ...
           {'procsNum'      'real'      []  1;
            'conName'       'cell'      []  {};
            'method'        'string'    {'onesampleT','flexANOVA'}  'flexANOVA';
            'nAnova'        'real'      []  2;% anova
            'conds'         'real'      []  [];% anova
            'TorF'          'string'    []  'T';% anova
            'contrast'      'cell'      []  {};% anova
            'allcon'        'real'      []  1;% anova
            'weight'        'cell'      []  {};% anova
            'conName_2nd'   'cell'      []  {};% anova
            'analysis'      'string'    []  '';
            'voi'           'string'    []  '';
            'diff'          'string'    []  '';
            'batchfolder'   'string'    []  '';
            'subidntfr'     'string'    []  'SUB';
            'datafolder'      'string'    []  '1st_level';
            'infilepat'     'string'    []  'swcmT*.nii';
            'rmSUB'         'cell'      []  {};
            'SPMorCONN'     'string'    []  'SPM';
            'space'         'string'    []  ''; % CONN
           });
    if isempty(in.infilepat)
        infilepat = [in.TorF,'_*.nii'];
    else
        infilepat = in.infilepat;
    end

    matlabbatch = [];
    subject = string({dir(filepath).name});
    subject = subject(contains(subject,in.subidntfr));
    in.rmSUB = string(in.rmSUB);
    for i = 1:length(in.rmSUB)
        subject(subject==in.rmSUB(i)) = [];
    end

    Nsubject = length(subject);
    
    % find SPM.mat file 
    switch in.SPMorCONN
        case "SPM"
            % check all subject has same SPM.mat contrast Name
            datapath = fullfile(filepath,'**',in.datafolder,'**','SPM.mat');
            datafolder = {dir(datapath).folder};
            load(fullfile(datafolder{1},'SPM.mat'));
            oldcondName = string({SPM.xCon.name});
            
            for nsub = 1:Nsubject
                datapath = datafolder{contains(datafolder,char(subject(nsub)))};
                load(fullfile(datapath,'SPM.mat'));
                conName = string({SPM.xCon.name});
                if ~any(arrayfun(@(x) any(oldcondName==x),conName))
                    error(sprintf("subject%d and subject%d contrast name is not same",nsub,nsub-1));
                end
                oldcondName = conName;
            end
            if ~isempty(in.conName)
                conName = string(in.conName);
            end
            srcName = "";
        case "CONN"
            % check all subject has same scource and same condition
            srcpath = fullfile(filepath,'**',in.datafolder,'**','_list_sources.txt');
            condipath = fullfile(filepath,'**',in.datafolder,'**','_list_conditions.txt');
            datafolder = strcat({dir(condipath).folder},filesep,in.space);

            allcond = struct();
            allsrc = struct();
            for nsub = 1:Nsubject
                % check condition
                condifile = dir(condipath);
                condifile = fullfile(condifile(nsub).folder,condifile(nsub).name);
                conname = readcell(condifile,"Delimiter",'=');
                for i = 1:size(conname,1)
                    allcond(nsub).(conname{i,2}) = conname{i,1};
                end

                % check source
                srcfile = dir(srcpath);
                srcfile = fullfile(srcfile(1).folder,srcfile(1).name);
                srcname = readcell(srcfile,"Delimiter",'=');
                for i = 1:size(srcname,1)
                    srctmp = srcname{i,2};
                    srctmp = split(srctmp,' ');
                    srctmp = srctmp{1};
                    srctmp(srctmp=='.') = '_';
                    allsrc(nsub).(srctmp) = srcname{i,1};
                end
            end

            % select condition that all subject contains
            conName = string(fieldnames(allcond));
            for i = 1:length(conName)
                if any(cellfun(@isempty,{allcond.(conName(i))}))
                    allcond = rmfield(allcond,conName(i));
                end
            end

            % select source that all subject contains
            srcName = string(fieldnames(allsrc));
            for i = 1:length(srcName)
                if any(cellfun(@isempty,{allsrc.(srcName(i))}))
                    allsrc = rmfield(allsrc,srcName(i));
                end
            end
    end
    Ncon = length(conName);


    switch in.method
        case 'onesampleT'
            for nsrc = 1:length(srcName)
                scans = cell(Nsubject,Ncon);
                for nsub = 1:Nsubject
                    datapath = datafolder{contains(datafolder,char(subject(nsub)))};
                    if srcName(1) ~= ""
                        infilepat = [in.infilepat,'*',allsrc(nsub).(srcName(nsrc)),'*'];
                        idx = find(infilepat == '*');
                        infilepat(idx(diff(idx) == 1)) = [];
                    end
                    confile = string(ls(fullfile(datapath,infilepat)));
                    for ncon = 1:Ncon
                        switch in.SPMorCONN
                            case "SPM"
                                confile = confile(contains(confile,char(conName(ncon))));
                            case "CONN"
                                confile = confile(contains(confile,allcond(nsub).(conName(ncon))));
                        end
                        scans{nsub,ncon} = [fullfile(datapath,char(confile)),',1'];
                    end
                end

                procs_Num = in.procsNum;
                for ncon = 1:Ncon
                    outpath = char(fullfile(outputpath,round,conName(ncon),in.voi,srcName(nsrc)));
                    if ~exist(outpath,'dir'), mkdir(outpath); end
                    % 2nd level design
                    matlabbatch{end+1}.spm.stats.factorial_design.dir = {outpath};
                    matlabbatch{end}.spm.stats.factorial_design.des.t1.scans = scans(:,ncon);
                    matlabbatch{end}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                    matlabbatch{end}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
                    matlabbatch{end}.spm.stats.factorial_design.masking.tm.tm_none = 1;
                    matlabbatch{end}.spm.stats.factorial_design.masking.im = 1;
                    matlabbatch{end}.spm.stats.factorial_design.masking.em = {''};
                    matlabbatch{end}.spm.stats.factorial_design.globalc.g_omit = 1;
                    matlabbatch{end}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
                    matlabbatch{end}.spm.stats.factorial_design.globalm.glonorm = 1;

                    % estimate
                    matlabbatch{end+1}.spm.stats.fmri_est.spmmat = {fullfile(outpath,'SPM.mat')};
                    matlabbatch{end}.spm.stats.fmri_est.write_residuals = 0;
                    matlabbatch{end}.spm.stats.fmri_est.method.Classical = 1;

                    % contrast
                    matlabbatch{end+1}.spm.stats.con.spmmat = {fullfile(outpath,'SPM.mat')};
                    matlabbatch{end}.spm.stats.con.consess{1}.tcon.name = 'GroupLevel';
                    matlabbatch{end}.spm.stats.con.consess{1}.tcon.weights = 1;
                    matlabbatch{end}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                    procs_Num = procs_Num+3;
                end
                
                if ~isempty(matlabbatch)
                    if f.gfunc.sndL
                        % run job
                        % ==================================================
                        spm('defaults', 'FMRI');
                        spm_jobman('initcfg')
                        spm_jobman('run',matlabbatch);
                        % ==================================================
        
                        % save batch
                        % ===================================================
                        in.batchfolder = fullfile(outpath,'batch');
                        if ~exist(in.batchfolder,'dir'), mkdir(in.batchfolder); end
                        if ~isempty(in.batchfolder)
                            if ~exist(in.batchfolder,"dir"), mkdir(in.batchfolder); end
                            save(fullfile(in.batchfolder,'snd_level.mat'),"matlabbatch");
                        end
                        % ===================================================
                    end
                    % plot flapmap
                    % ===================================================
                    if f.gfunc.flapmapPlot
                        Steps.plotflapmat(outpath,'spmT*.nii')
                    end
                    % ===================================================  
                end
                matlabbatch = [];
            end
        case 'flexANOVA' 
            procs_Num = in.procsNum;
            foldername = ['FANOVA_',round,'_'];
            for i = 1:length(in.conName)
                foldername = append(foldername,[in.conName{i}]);
                if i ~= length(in.conName), foldername = append(foldername,'_'); end
            end
            outpath = char(fullfile(outputpath,foldername,in.voi));
            if ~exist(outpath,"dir"), mkdir(outpath); end
    
            % second level design
            matlabbatch{end+1}.spm.stats.factorial_design.dir = {outpath};
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(1).name = 'sub';
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(2).name = 'homophone';
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(3).name = 'frequency';
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(3).dept = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
            for nsub = 1:Nsubject
                scans = cell(Ncon,1);
                for ncon = 1:Ncon
                    indipath = fullfile(filepath,char(subject(nsub)),round,'1st_level',infilepat);
                    filename = {dir(indipath).name};
                    filename = char(filename(contains(filename,conName(ncon))));
                    scans{ncon} = fullfile(filepath,char(subject(nsub)),round,'1st_level',filename);
                end
                subconds = ones(Ncon,1)*nsub;
                conconds = in.conds;
                conds = cat(2,subconds,conconds);
                matlabbatch{end}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(nsub).scans = scans;
                matlabbatch{end}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(nsub).conds = conds;
            end   
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
            matlabbatch{end}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [2 3];
            matlabbatch{end}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{end}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{end}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{end}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{end}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{end}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{end}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{end}.spm.stats.factorial_design.globalm.glonorm = 1;
    
            % estimate second level model
            matlabbatch{end+1}.spm.stats.fmri_est.spmmat = {fullfile(outpath,'SPM.mat')};
            matlabbatch{end}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{end}.spm.stats.fmri_est.method.Classical = 1;
            writecell(in.conName,fullfile(outpath,'contrast_name.txt'));
    
            % define contrast
            matlabbatch{end+1}.spm.stats.con.spmmat = {fullfile(outpath,'SPM.mat')};
            if in.allcon
                for ncon = 1:length(conName)
                    cw = zeros(1,Ncon);
                    cw(conName==(conName(ncon))) = 1;
                    weight = [ones(1,Nsubject)/Nsubject,cw];
                    matlabbatch{end}.spm.stats.con.consess{ncon}.tcon.name = char(conName(ncon));
                    matlabbatch{end}.spm.stats.con.consess{ncon}.tcon.weights = weight;
                    matlabbatch{end}.spm.stats.con.consess{ncon}.tcon.sessrep = 'none';
                end
                conN = length(conName);
            else
                conN = 0;
            end
            if ~isempty(in.contrast)
                if length(in.contrast)~=length(in.weight),error('contrast and weight length is not same'); end
                for ncon = 1:length(in.contrast)
                    weightPos = arrayfun(@(x) find(conName == x),string(in.contrast{ncon}));
                    cw = zeros(size(conName));
                    for i = 1:length(weightPos)
                        cw(weightPos(i)) = in.weight{ncon}(i);
                    end
    
                    weight = [zeros(1,Nsubject),cw];
      
                    matlabbatch{end}.spm.stats.con.consess{ncon+conN}.tcon.name = char(in.conName_2nd(ncon));
                    matlabbatch{end}.spm.stats.con.consess{ncon+conN}.tcon.weights = weight;
                    matlabbatch{end}.spm.stats.con.consess{ncon+conN}.tcon.sessrep = 'none';
                end
            end
            if ~isempty(matlabbatch)
                if f.gfunc.sndL
                    % run job
                    % ==================================================
                    spm('defaults', 'FMRI');
                    spm_jobman('initcfg')
                    spm_jobman('run',matlabbatch);
                    % ==================================================
    
                    % save batch
                    % ===================================================
                    in.batchfolder = fullfile(outpath,'batch');
                    if ~exist(in.batchfolder,'dir'), mkdir(in.batchfolder); end
                    if ~isempty(in.batchfolder)
                        if ~exist(in.batchfolder,"dir"), mkdir(in.batchfolder); end
                        save(fullfile(in.batchfolder,'snd_level.mat'),"matlabbatch");
                    end
                    % ===================================================
                end
            end
    end

  

    % plot flapmap
    % ===================================================
    if f.gfunc.flapmapPlot && in.space ~= "MNI"
        Steps.plotflapmat(outpath,'spmT*.nii')
    end
    % ===================================================  
end
