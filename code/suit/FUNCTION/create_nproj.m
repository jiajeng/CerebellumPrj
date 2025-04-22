function create_nproj(procestep,sliceorder,varargin)
    in = finputcheck(varargin, ...
    {'sub'              'cell'      []  [];
     'sess'             'string'    []  [];
     'round'            'cell'      []  {};
     'condition'        'cell'      []  {};
     'conn_proj_name'   'string'    []  [];
     'conn_proj_path'   'string'    []  [];
     'TR'               'real'      []  [];
     'smooth_kernel'    'real'      []  [];
     'StrucVres'        'real'      []  [];
     'funcVres'         'real'      []  [];
     'filter_band'      'real'      []  [];
     'analysisName'     'string'    []  [];
     'contrast'         'cell'      []  {};
     'outsubpath'       'string'    []  [];
     'mrifilepath'      'string'    []  [];
     'rawdata'          'boolean'   []  [];
     'prepsteps'        'string'    []  'default_mni';
     'roiPath'          'string'    []  '.\..\ROI';
     'subidntfr'        'string'    []  'SUB';
     'AnaSource'        'string'    []  [];
     });
    if isempty(in.sub)
        subject = {dir(in.mrifilepath).name};
        subject = subject(contains(subject,'SUB'));
    else
        subject = in.sub;
    end
    sess = in.sess;
    round = in.round;
    conn_proj_name = in.conn_proj_name;
    conn_proj_path = in.conn_proj_path;
    if ~exist(conn_proj_path,'dir'), mkdir(conn_proj_path); end
    condition = in.condition;
    filepath = in.mrifilepath;
    
    % round is conn toolbox "session"
    if ~exist(filepath,'dir')
        mkdir(filepath);
    end
    
    % ----------------------------define condition-----------------------------
    % sess = 'ses-01';
    % round = {'REST'};
    % conn_proj_name = ['conn_',cell2mat(round)];
    % conn_proj_folder = ['conn_',cell2mat(round)];
    % conn_proj_p = pwd;
    % conn_proj_path = fullfile(conn_proj_p,conn_proj_folder);
    % -----------------------------------------------------------------------
 
    % -----------------------------define process---------------------------
    
    % procestep = 'Add_roi';

    steps = dictionary('proces_use_procD',struct("Setup",1, ...
                                          "Preprocessing",0, ...
                                          "Denoising",1, ...
                                          "Add_Roi",0, ...
                                          "fst_Analysis",1, ...
                                          "snd_Analysis",1, ...
                                          "Add_RoiResult",0) ...
                  ,'Denois_use_procD',struct("Setup",1, ...
                                          "Preprocessing",0, ...
                                          "Denoising",1, ...
                                          "fst_Analysis",0, ...
                                          "Add_Roi",0, ...
                                          "snd_Analysis",0, ...
                                          "Add_RoiResult",0) ...
                  ,'proces_new',struct("Setup",1, ...
                                      "Preprocessing",1, ...
                                      "Denoising",1, ...
                                      "Add_Roi",0, ...
                                      "fst_Analysis",1, ...
                                      "snd_Analysis",1, ...
                                      "Add_RoiResult",0) ...
                  ,'Denois_new',struct("Setup",1, ...
                                      "Preprocessing",1, ...
                                      "Denoising",1, ...
                                      "Add_Roi",0, ...
                                      "fst_Analysis",0, ...
                                      "snd_Analysis",0, ...
                                      "Add_RoiResult",0) ...
                 ,'Add_roi_Reslt',struct("Setup",0, ...
                                       "Preprocessing",0, ...
                                       "Denoising",0, ...
                                       "Add_Roi",0, ...
                                       "fst_Analysis",0, ...
                                       "snd_Analysis",0, ...
                                       "Add_RoiResult",1) ...
                 ,'Add_roi',struct("Setup",0, ...
                                   "Preprocessing",0, ...
                                   "Denoising",0, ...
                                   "Add_Roi",1, ...
                                   "fst_Analysis",0, ...
                                   "snd_Analysis",0, ...
                                   "Add_RoiResult",0) ...
                                   );
    % using step by pre-define steps
    switch class(procestep)
        case "string"
            procesStep = steps(procestep);
        case "char"
            procesStep = steps(procestep);
        case "struct"
            procesStep = procestep;
    end
    procesStepName = dictionary(fieldnames(procesStep),["Setup","Preprocessing","Denoising","Add_Roi","1st Analysis","2nd Analysis","Add_RoiResult"]');
    % ----------------------------------------------------------------------
    
    
    % -----------------------------define filename--------------------------
    if procestep.Preprocessing
        func_filename = '*_4D.nii';
        struct_filename = '*_T1.nii';
    else
        func_filename = 'wwau*_4D.nii';
        struct_filename = 'c0SUB*.nii';
    end
    % -----------------------------------------------------------------------
    
    % -----------------------------defind roi file---------------------------
    roifile = {fullfile(fileparts(which('conn')),'rois','atlas.nii'),...
           fullfile(fileparts(which('conn')),'rois','networks.nii'),...
           };
    adroifile = {dir(in.roiPath).name};
    adroifile = adroifile(contains(adroifile,in.subidntfr)&contains(adroifile,'.nii'));
    adroifile = strcat([in.roiPath,filesep],adroifile);
    roifile = cat(2,roifile,adroifile);
    % -----------------------------------------------------------------------
    
    % ------------------------------define parameter---------------------------
    % setup
    onsetname = 'WordOnset';
    dur = 0; % event duration
    
    % preprocessing
    voxelres = 1; % 1:2mm(default SPM), 2:same as structural, 3:same as functional, 4:surface-based template(Freesurfer)
    Analysis_unit = 1;% 1:PSC uinits   2:raw units
    Disunit = dictionary(1,'PSC units',2,'raw units');
    if isfile(in.prepsteps)
        steps = load(in.prepsteps);
        prepsteps = steps.STEPS;
    else
        prepsteps = in.prepsteps;
    end
    
    % denoising 
    % analysis
    Betw_SUB_eff = {'AllSubjects'};
    Betw_SUB_con = 1;
    Betw_cond_eff = condition;
    Betw_cond_con = cell(1,length(in.contrast));
    if length(condition) == 1
        Betw_cond_con = {1};
    else
        for i = 1:length(in.contrast)
            Betw_cond_con{i} = diag(in.contrast{i});
        end
    end
    
    % ----------------------------------------------------------------------
    
    nsubject = length(subject);
    nround = length(round);
    
    % ------------------------------batch-------------------------------------
    if procesStep.Setup
        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
        projinfo.TR = in.TR;
        strcond = [];
        for i = 1:length(condition)
            if i~=length(condition)
                strcond = [strcond,condition{i},', '];
            else
                strcond = [strcond,condition{i}];
            end
        end
        projinfo.condition = string(strcond);
        %-------------------------
    
        % ------------get all subject func and struct .nii filepath--------------
        % Selects functional / anatomical volumes
        % get all func and anat file filepath --> C:/filepath/filename.nii
        FUNCTIONAL_FILE = cell(nsubject,nround); % nsub * nsess(REST and task)
        STRUCTURAL_FILE = cell(nsubject,1); % nsub * 1(T1)
        for nsub = 1:size(FUNCTIONAL_FILE,1)
            for nrd = 1:nround
                STRUCTURAL_file=dir(fullfile(filepath,char(subject(nsub)),'**',struct_filename));
                % FUNCTIONAL_file=dir(fullfile(filepath,char(subject(nsub)),'niifile',sess,'mri',session{nsess},func_filename));
                FUNCTIONAL_file=dir(fullfile(filepath,char(subject(nsub)),'**',round{nrd},'**',func_filename));
                FUNCTIONAL_FILE(nsub,nrd) = {fullfile(char(FUNCTIONAL_file(1).folder),char(FUNCTIONAL_file(1).name))};
                STRUCTURAL_FILE(nsub,nrd) = {fullfile(char(STRUCTURAL_file(1).folder),char(STRUCTURAL_file(1).name))};
            
              
                % ----------------------log 
                tmp = niftiinfo(FUNCTIONAL_FILE{nsub,nrd});
                fdatainfo(nsub).Subject = subject(nsub);
                fdatainfo(nsub).(['Fpath_',num2str(nrd)]) = FUNCTIONAL_FILE{nsub,nrd};
                fdatainfo(nsub).(['ImagSize_',num2str(nrd)]) = num2str(tmp.ImageSize);
                fdatainfo(nsub).VoxlSize = num2str(tmp.PixelDimensions);
        
                tmp = niftiinfo(STRUCTURAL_FILE{nsub,nrd});
                sdatainfo(nsub).Subject = subject(nsub);
                sdatainfo(nsub).Fpath = STRUCTURAL_FILE{nsub,1};
                sdatainfo(nsub).ImagSize = num2str(tmp.ImageSize);
                sdatainfo(nsub).VoxlSize = num2str(tmp.PixelDimensions);
        
                % ---------------------------
            end
        end
        f = 0;
        tmp = checkImgSize(sdatainfo,'ImagSize');f = f || tmp;
        tmp = checkImgSize(fdatainfo,'ImagSize_1');f = f || tmp;
        % tmp = checkImgSize(fdatainfo,'ImagSize_2');f = f || tmp;
        if f, keyboard; end
        % ------------------------------------------------------------------------
    
        % ---------------batch
        batch = [];
        batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
        % define functional and structual filepath
        batch.Setup.nsubjects = nsubject;
        % set structural, functional file and condition variable
        batch.Setup.conditions.names= condition;
        batch.Setup.sessions.names = round;
        for ncond=1:length(condition)
            for nsub=1:nsubject
                batch.Setup.structurals{nsub}=STRUCTURAL_FILE{nsub};
                for nses=1:nround  
                    batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; 
                    batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;
                    batch.Setup.functionals{nsub}{nses}=FUNCTIONAL_FILE{nsub,nses};
                end
            end
        end
        
        % set task condition
        if string(round) ~= "REST"
            for nsub = 1:nsubject
                for nses = 1:nround
                    behpath = fullfile(filepath,char(subject(nsub)),'BEHAV',round{nses});
                    load(fullfile(behpath,'behD.mat'))
                    for ncond = 1:length(condition)
                        try
                            onsettmp = behD.(condition{ncond}).(onsetname);
                        catch ME
                            onsettmp = behD.(condition{ncond}).Onset;
                        end
                        batch.Setup.conditions.onsets{ncond}{nsub}{nses} = cell2mat(onsettmp)+0.75;
                        batch.Setup.conditions.durations{ncond}{nsub}{nses} = ones(1,length(onsettmp))*dur;
                    end
                end
            end
        end
        
        % define roi
        for nroi = 1:length(roifile)
            for nsub = 1:nsubject
                for nses = 1:nround
                    batch.Setup.rois.files{nroi}{nsub}{nses}=roifile(nroi);
                    if ispc
                        roiname = split(roifile{nroi},'\');
                    elseif isunix
                        roiname = split(roifile{nroi},'/');
                    end
                    roiname = split(roiname{end},'.');
        
                    roiname = roiname{1};
                    batch.Setup.rois.names{nroi}=roiname;
                end
            end
        end
        
        if ~in.rawdata
            % define covariate
            % corvName = {'realignment','QC_timeseries','scrubbing'};
            % corvfile = {'rp_*.txt','art_regression_timeseries*.mat','art_regression_outliers_au*4D.mat'};
            % for ncov = 1:length(corvName)
            %     for nsub = 1:nsubject
            %         for nrd = 1:nround
            %             corvpath = fullfile(filepath,char(subject(nsub)),'niifile',sess,'mri',round{nrd});
            %             corvfname = ls(fullfile(corvpath,corvfile{ncov}));
            %             batch.Setup.covariates.files{ncov}{nsub}{nrd} = {fullfile(corvpath,corvfname)};
            %         end
            %     end
            %     batch.Setup.covariates.names{ncov} = corvName{ncov};
            % end
        end
        
        batch.Setup.RT=in.TR;
        % batch.Setup.conditions.names=condition;       
        batch.Setup.isnew=1; 
        batch.Setup.voxelresolution = voxelres;% same as functionals
        batch.Setup.outputfiles = [0,1,0,0,0,0]; % creates:[confound beta-map,
        %                                                   confound-correlated timeseries,
        %                                                   seed-to-voxel r-maps,
        %                                                   seed-to-voxel p-maps,
        %                                                   seed-to-voxel FDR-p-maps,
        %                                                   ROI-extraction REX files]
        
        batch.Setup.conditions.missingdata = 1;
        if ~in.rawdata
            batch.Setup.done=1;
        end
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'Setup'}));
        conn_batch(batch);
    end
    
    if procesStep.Preprocessing
        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
        projinfo.Struct_Voxel_Resolution = in.StrucVres;
        projinfo.Func_Voxel_Resolution = in.funcVres;
        projinfo.smooth_kernel_fwhm = in.smooth_kernel;
        projinfo.sliceorder = string(sliceorder);
        projinfo.Analysis_units = Disunit(Analysis_unit);
        %--------------------------
        batch = [];
        batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
        %% CONN preprocessing
        batch.Setup.analysisunits = Analysis_unit;
        batch.Setup.preprocessing.sliceorder = sliceorder;
    
        batch.Setup.preprocessing.steps = prepsteps;
        projinfo.PrepStep = string(prepsteps)';
    
        batch.Setup.preprocessing.voxelsize_anat = in.StrucVres;
        batch.Setup.preprocessing.voxelsize_func = in.funcVres;
        batch.Setup.preprocessing.fwhm = in.smooth_kernel;
        batch.Setup.overwrite = 1;% not overwrite data can run faster 
        batch.Setup.conditions.missingdata = 1;
        batch.Setup.done=1;
    
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'Preprocessing'}));
        conn_batch(batch);
    end
    
    if procesStep.Denoising
        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
    
        projinfo.filter = string(mat2str(in.filter_band));
        % projinfo.confound = string(confound)';
        % projinfo.confound_dim = cell2mat(confoundNum');
        %--------------------------
        batch = [];
        batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
        batch.Denoising.filter=in.filter_band;          % frequency filter (band-pass values, in Hz)
        % batch.Denoising.confounds.names = confound;
        % for i = 1:length(confound)
        %     batch.Denoising.confounds.dimensions{i} = confoundNum{i};
        % end
        batch.Denoising.overwrite=1;
        batch.Denoising.done=1;
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'Denoising'}));
        conn_batch(batch);
    end

    if procesStep.Add_Roi
        roifilepath = in.roiPath;
        roifilename = string(ls(fullfile(roifilepath,'*.nii')));
        roiname = split(roifilename,'.');
        roiname = roiname(:,1);
    
        % --------------------check rois name is in project already or not
        CONN_x = load(fullfile(conn_proj_path,[conn_proj_name,'.mat']));
        CONN_x = CONN_x.CONN_x;
        projRoi = CONN_x.Setup.rois.names;
        projRoi = string(projRoi)';
        rep_roi = [];
        for i = 1:length(roiname)
            r = roiname(i);
            if any(projRoi==r)
                rep_roi = cat(2,rep_roi,find(roiname==r));
            end
        end
        roiname(rep_roi) = [];roifilename(rep_roi) = [];
        % ----------------------------------------------------------------
    
        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
        tmp = cellfun(@(x) string(num2str(diag(x))),Betw_cond_con,'UniformOutput',false);
        
        projinfo.Roi = roiname(:,1);
        roicor = [];
        for nfile = 1:length(roifilename)
            roi = load(fullfile(roifilepath,[char(roiname(nfile)),'.mat']));
            roi = roi.roi;
            roicor = cat(1,roicor,string(num2str(c_o_m(roi))));
        end
        projinfo.RoiCoordi = roicor;
        %--------------------------
    
        roifilepath = strcat(repmat({roifilepath},length(roifilename),1),'\',char(roifilename));
        roifilepath = convertStringsToChars(roifilepath);
    
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'Add_Roi'}));
    
        conn_batch( 'filename',fullfile(conn_proj_path,[conn_proj_name,'.mat']), ...
                    'Setup.rois.add',1, ...
                    'Setup.rois.files',roifilepath, ...
                    'Setup.rois.names',convertStringsToChars(roiname), ...
                    'Setup.overwrite',0,...
                    'Setup.done',1, ...
                    'Denoising.overwrite',0,...
                    'Denoising.done',1)
    end
    
    if procesStep.fst_Analysis
        % define Sources
        Analysis_name = in.analysisName;
        if ~isempty(in.AnaSource)
            if in.AnaSource == "SUBroi"
                load(fullfile(conn_proj_path,[conn_proj_name,'.mat']))
                Analysis_Source = CONN_x.Analysis_variables.names(cellfun(@(x) x == "roi",CONN_x.Analysis_variables.types));
                Analysis_Source = Analysis_Source(contains(Analysis_Source,in.subidntfr));
            else
                Analysis_Source = in.AnaSource;
            end
        else
            % all source name 
            load(fullfile(conn_proj_path,[conn_proj_name,'.mat']))
            Analysis_Source = CONN_x.Analysis_variables.names(cellfun(@(x) x == "roi",CONN_x.Analysis_variables.types));
        end

        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
        projinfo.fst_Analysis = string(Analysis_name);
        projinfo.Source = string(Analysis_Source);

        %--------------------------
        batch = [];
        batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
        %% CONN 1st level Analysis
        batch.Analysis.name = Analysis_name;
        batch.Analysis.done = true;
        batch.Analysis.sources = Analysis_Source;
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'fst_Analysis'}));
        conn_batch(batch);
    end
    if procesStep.snd_Analysis
        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
        tmp = cellfun(@(x) string(num2str(diag(x))),Betw_cond_con,'UniformOutput',false);
        for neft = 1:length(Betw_cond_con)
            projinfo.(['snd_Analysis_con_',num2str(neft)]) = cat(1,string(Betw_cond_eff),tmp{neft}');
        end
        %--------------------------
    
        %% CONN 2nd level
        roipath = fullfile(conn_proj_folder,conn_proj_name,'results','firstlevel',Analysis_name);
        load(fullfile(roipath,'_list_sources.mat'));
        com_sourcenames = sourcenames;
        projinfo.roi = com_sourcenames;
        for nroi = 1:length(com_sourcenames)
            for neft = 1:length(Betw_cond_con)
                batch = [];
                batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
                batch.Results.analysis_number = Analysis_name;
                batch.Results.between_subjects.effect_names = Betw_SUB_eff;
                batch.Results.between_subjects.contrast = Betw_SUB_con;
                batch.Results.between_conditions.effect_names = Betw_cond_eff;
                batch.Results.between_conditions.contrast = Betw_cond_con{neft};
                batch.Results.between_sources.effect_names = com_sourcenames(nroi);
                batch.Results.between_sources.contrast = 1;
                batch.Results.display = false;
                conn_batch(batch);
            end
        end
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'snd_Analysis'}));
        conn_batch(batch);
    end
    
    if procesStep.Add_RoiResult
        roifilepath = in.roiPath;
        roifilename = string(ls(fullfile(roifilepath,'*.nii')));
        roiname = split(roifilename,'.');
        roiname = roiname(:,1);
    
        % --------------------check rois name is in project already or not
        CONN_x = load(fullfile(conn_proj_path,[conn_proj_name,'.mat']));
        CONN_x = CONN_x.CONN_x;
        projRoi = CONN_x.Setup.rois.names;
        projRoi = string(projRoi)';
        rep_roi = [];
        for i = 1:length(roiname)
            r = roiname(i);
            if any(projRoi==r)
                rep_roi = cat(2,rep_roi,find(roiname==r));
            end
        end
        roiname(rep_roi) = [];roifilename(rep_roi) = [];
        % ----------------------------------------------------------------
    
        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
        tmp = cellfun(@(x) string(num2str(diag(x))),Betw_cond_con,'UniformOutput',false);
        
        projinfo.Roi = roiname(:,1);
        roicor = [];
        for nfile = 1:length(roifilename)
            roi = load(fullfile(roifilepath,[char(roiname(nfile)),'.mat']));
            roi = roi.roi;
            roicor = cat(1,roicor,string(num2str(c_o_m(roi))));
        end
        projinfo.RoiCoordi = roicor;
        %--------------------------
    
        roifilepath = strcat(repmat({roifilepath},length(roifilename),1),'\',char(roifilename));
        roifilepath = convertStringsToChars(roifilepath);
    
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'Add_Roi'}));
    
        conn_batch( 'filename',fullfile(conn_proj_path,[conn_proj_name,'.mat']), ...
                    'Setup.rois.add',1, ...
                    'Setup.rois.files',roifilepath, ...
                    'Setup.rois.names',convertStringsToChars(roiname), ...
                    'Setup.overwrite',0,...
                    'Setup.done',1,...
                    'Denoising.overwrite',0,...
                    'Denoising.done',1,...
                    'Analysis.Source',convertStringsToChars(roiname), ...
                    'Analysis.overwrite',0,...
                    'Analysis.done',1)
        for nfile = 1:length(roifilename)
            for neft = 1:length(Betw_cond_con)
                conn_batch( 'Results.analysis_number',Analysis_name, ...
                            'Results.between_subjects.effect_names', Betw_SUB_eff, ...
                            'Results.between_subjects.contrast',Betw_SUB_con, ...
                            'Results.between_conditions.effect_names', Betw_cond_eff, ...
                            'Results.between_conditions.contrast',Betw_cond_con{neft}, ...
                            'Results.between_sources.effect_names',{roiname(nfile)}, ...
                            'Results.between_sources.contrast',1, ...
                            'Results.display',false)
                
            end
        end
    end

     if procesStep.snd_Analysis
        % ----------------------log 
        projinfo = struct();
        sdatainfo = struct();
        fdatainfo = struct();
    
        tmp = cellfun(@(x) string(num2str(diag(x))),Betw_cond_con,'UniformOutput',false);
        for neft = 1:length(Betw_cond_con)
            projinfo.(['snd_Analysis_con_',num2str(neft)]) = cat(1,string(Betw_cond_eff),tmp{neft}');
        end
        %--------------------------
    
        %% CONN 2nd level
        roipath = fullfile(conn_proj_folder,conn_proj_name,'results','firstlevel',Analysis_name);
        load(fullfile(roipath,'_list_sources.mat'));
        com_sourcenames = sourcenames;
        projinfo.roi = com_sourcenames;
        for nroi = 1:length(com_sourcenames)
            for neft = 1:length(Betw_cond_con)
                batch = [];
                batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
                batch.Results.analysis_number = Analysis_name;
                batch.Results.between_subjects.effect_names = Betw_SUB_eff;
                batch.Results.between_subjects.contrast = Betw_SUB_con;
                batch.Results.between_conditions.effect_names = Betw_cond_eff;
                batch.Results.between_conditions.contrast = Betw_cond_con{neft};
                batch.Results.between_sources.effect_names = com_sourcenames(nroi);
                batch.Results.between_sources.contrast = 1;
                batch.Results.display = false;
                conn_batch(batch);
            end
        end
        saveLog(projinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,procesStepName({'snd_Analysis'}));
        conn_batch(batch);
    end
    
    % s = [];
    % stpN = fieldnames(procesStep);
    % for i = 1:length(stpN)
    %     if procesStep.(stpN{i})
    %         s = cat(1,s,procesStepName(stpN(i)));
    %     end
    % end
    % projinfo.Steps = s;

    %% function define
    function f = checkImgSize(datainfo,fieldName)
        IN = inputname(1);
        if IN(1) == 's'
            warningmeg = ['Structual Data',newline];
        elseif IN(1) == 'f'
            buf = split(fieldName,'_');
            buf = char(buf(end));
            switch buf
                case '1'
                    warningmeg = ['Functional Data round 1',newline];
                case '2'
                    warningmeg = ['Functional Data round 2',newline];
            end
        end
    
        sIamgSize = string(cell2mat({datainfo.(fieldName)}'));
        usIamgSize = unique(sIamgSize);
        usIamgSizeN = [];
        ssub = {datainfo.Subject};
    
        if size(unique(usIamgSize),1)>1
            for ix = 1:length(usIamgSize)
                usIamgSizeN = cat(1,usIamgSizeN,sum(sIamgSize==usIamgSize(ix)));
            end
            warningmeg = cat(2,warningmeg,['almost IamgeSize is ',char(usIamgSize(usIamgSizeN==max(usIamgSizeN))),newline]);
            usIamgSize(usIamgSizeN == max(usIamgSizeN)) = [];
            for ix = 1:length(usIamgSize)
                wsub = ssub(sIamgSize==usIamgSize(ix));
                warningmeg = cat(2,warningmeg,[cell2mat(strcat(wsub{:}," ")'),'has ImageSize ',char(usIamgSize(ix)),newline]);
            end
            disp(warningmeg)
            f = 1;
        else
            f = 0;
        end
    end
    
    function saveLog(oldinfo,sdatainfo,fdatainfo,conn_proj_path,conn_proj_name,Step)
        oldinfo.Steps = Step;
        % define log info table
        ML = max(cellfun(@(x) length(oldinfo.(x)),fieldnames(oldinfo)));
        fN = fieldnames(oldinfo);
        newinfo = struct();
    
        for ix = 1:length(fN)
            for jx = 1:ML+2
                try
                    newinfo(jx).(fN{ix}) = oldinfo.(fN{ix})(jx);
                catch
                    if jx == ML+2 && ix == length(fN)
                        newinfo(jx).(fN{ix}) = string(datetime("now"));
                    else
                        if class(oldinfo.(fN{ix})) == "double"
                            newinfo(jx).(fN{ix}) = NaN;
                        elseif class(oldinfo.(fN{ix})) == "string"
                            newinfo(jx).(fN{ix}) = "";
                        end
                    end
                end
            end
        end
        
        shtnum = 1;
        if exist(fullfile(conn_proj_path,[conn_proj_name,'_info.xlsx']),'file')
            shtNam = sheetnames(fullfile(conn_proj_path,[conn_proj_name,'_info.xlsx']));
            buf = split(shtNam(contains(shtNam,'PROJECT')),'_');
            shtnum = max(str2double(buf(:,end)))+1;
        end
        if isnan(shtnum), shtnum = 1; end
        if ~exist(fullfile(conn_proj_path),'dir'), mkdir(conn_proj_path); end
        
        if ~isempty(fieldnames(sdatainfo))
            writetable(struct2table(sdatainfo),fullfile(conn_proj_path,[conn_proj_name,'_info.xlsx']),'Sheet','STRUCT_DATA_info');
            writetable(struct2table(fdatainfo),fullfile(conn_proj_path,[conn_proj_name,'_info.xlsx']),'Sheet','FUNCTION_DATA_info');
        end
        
        writetable(struct2table(newinfo),fullfile(conn_proj_path,[conn_proj_name,'_info.xlsx']),'Sheet',['PROJECT_info_',sprintf('%03d',shtnum)]);
    end
end
