%SET VARIABLE
[o,f] = Setup.SetVariable;
UIFlag = false;
SetupUI(f,o);
while 1
    pause(1);
    if UIFlag
        break;
    end
end
ofname = fieldnames(o);
for i = 1:length(ofname)
    eval([ofname{i},'=o.',ofname{i},';']);
end

clear o ofname i;
f.func.resliceMNI = true;
noSUB = [];

% =============== get behave file
if f.beh
    % get folder name in nas
    % Setup.GetBehFolder(BEH,sub,sess,ftpServer)
    
    % download behave file to local
    % Setup.GetBehfile(BEH,ftpServer,subidntfr,subpath)
end
% =================================

if exist("ow","var")
    save('overwrite.mat',"ow")
end

% GET DATA
if f.raw
    convert2nii(sub,sess,round, ...
        'folder_nest',FolderNest, ...
        'subpath',subpath, ...
        'outsubpath',subpath, ...
        'ftpServer',ftpServer, ...
        'DataType','cere')
    if ~f.local
        niiFile2local(sess,round,ftpServer,FolderNest,subpath)
    end
end
if f.GetData
    if ~f.local
        niiFile2local(sess,round,ftpServer,FolderNest,subpath)
    end
end

if ~f.conn
    % GET MRI INFO
    for nround = 1:length(round)
        MRinfo.(round{nround}) = readtable(fullfile(subpath,'MRinfo.xlsx'),'ReadRowNames',true,'Sheet',round{nround});
    end
    
    if f.beh
        % organized behave data and run simple statistic
        Steps.cereb_stat_beh(subpath,fullfile('BEHAV','Word'),ftpServer)
    end
    %% SPM analysis
    % step 1 anatomical analysis to normalize to SUIT space
    % step 2 functional preprocessing (realign and sliceTiming)
    % step 3 functional 1st level   
    % step 4 functional 2nd level

    spm fmri
    % MAIN
    if f.ana.ana || f.func.func
        for nsub = 1:length(sub)
            try
            for nsess = 1:length(sess)
                % for mri data
                indipath = fullfile(subpath,sub{nsub},sess{nsess});
                fuc_round = round(~contains(round,'T1'));
                if f.ana.ana
                    fprintf('===========================================\n')
                    fprintf('anatomical Data Step1 isolate and normalize\n')
                    fprintf('===========================================\n')
                    % structure
                    % ===========================================
                    % step1 isolate the cerebellum and brainstem from structure data
                    % step2 normalize structure data to template space
                    Steps.ana_steps(fullfile(indipath,'T1'),'SUB*T1.nii',f);
                    % ==========================================
                end
        
        
                if f.func.func
                    for nround = 1:length(fuc_round)
                        batchfolder = fullfile(indipath,round{nround},'batchfile');
                        info = MRinfo.(fuc_round{nround})(sub{nsub},:);
                        r = fuc_round{nround};
                        if r == "Rest" || r == "REST"
                            r = "rest";
                        elseif contains(r,'WORD')
                            r = R; 
                        end
                        % set condition
                        fst_con = Setup.Getcondition(r,'behpath',fullfile(indipath,'BEHAV',fuc_round{nround}));
        
        
                        fprintf('===========================================\n')
                        fprintf('functional Data Step1 "%s" preprosessing\n',fuc_round{nround})
                        fprintf('===========================================\n')
                        % functional
                        % ============================================
                        % step1 - preprocessing 
                        %           - realignment
                        %           - slice timing
                        %           - corigister to structurnal space
                        %       - output : auNAME.nii
                        Steps.func_step1(fullfile(indipath,'T1'), ...
                              'SUB*.nii', ...
                              fullfile(indipath,fuc_round{nround}), ...
                              'SUB*.nii', ...
                              str2double(info.RT), ...
                              str2num(info.sliceorder{1}), ...
                              batchfolder,f);
        
                        fprintf('===========================================\n')
                        fprintf('functional Data Step2 "%s" 1st level\n',fuc_round{nround})
                        fprintf('===========================================\n')
                        % step2 - 1st level
                        %           -model design and estimate
                        %           -define contrast
                        %       - output : T_CONTRAST.nii
                        outputFolder = fullfile(fullfile(indipath,fuc_round{nround}),'1st_level');
                        Steps.func_step2(fullfile(indipath,fuc_round{nround}), ...
                            'au*.nii', ...
                            outputFolder, ...
                            str2double(info.RT), ...
                            'secs', ...
                            batchfolder, ...
                            fst_con.condition, ...
                            fst_con.con, ...
                            fst_con.conWt, ...
                            fst_con.conName, ...
                            fst_con.condi, ...
                            f);
        
                        fprintf('===========================================\n')
                        fprintf('plot T1 and epi image cheg_img\n')
                        fprintf('outputpath : %s\n',opath);
                        fprintf('===========================================\n')
                        % check mri image is same using beta map and raw T1
                        oname = ['SUB',sprintf('%02d',nsub),'_',sess{nsess},'_',round{nround},'.jpeg'];
                        % Steps.checkmrifile(fullfile(indipath,'T1','SUB*T1.nii') ...
                        %           ,fullfile(indipath,fuc_round{nround},'1st_level','beta_0001.nii'), ...
                        %           oname, ...
                        %           opath)


                        fprintf('===========================================\n')
                        fprintf('functional Data Step3 "%s" reslice and normalize to SUIT space\n',round{nround})
                        fprintf('===========================================\n')
                        % step 3 reslice and normalize to SUIT space
                        %        - get functional contrast data
                        %           - reslice to SUIT (mask --> T1 'c_name_pcereb.nii')
                        %           - normalize to SUIT 
                        %        - output : wsuit_NAME.nii
                        T1filepath = fullfile(indipath,'T1');
                        T1filename = 'c_*_pcereb.nii';
                        funcfilepath = fullfile(indipath,fuc_round{nround},'1st_level');
                        funcfilename = [];
        
                        Steps.func_step3(T1filepath,T1filename,funcfilepath,funcfilename,batchfolder,f);
                        % ============================================
                    end % function round loop
                end % run functional data or not
            end % session loop
            catch ME
                noSUB = cat(2,noSUB,string(nsub));
                continue;
            end
        end % sub loop
    end
    
    if f.gfunc.func
        fuc_round = round(~contains(round,'T1'));
        rmSUB = {'SUB001','SUB050'};
        msg = 'check mri image from \n';
        for nround = 1:length(fuc_round)
            tmp = fullfile(opath,fuc_round{nround});
            tmp = split(tmp,filesep);
            tmp = strjoin(tmp,string([filesep,filesep]));
            msg = cat(2,[msg,tmp,' \n']);
        end
        msg = cat(2,[msg,'Enter subject name with cell array that want to remove from 2nd level \n(e.g. {''SUB001'',''SUB022''})\n: ']);
        if ~exist("rmSUB","var")
            rmSUB = input(msg);
        end

        if isempty(rmSUB)
            rmSUB = {};
        end
        snd_round = round(~contains(round,'T1'));
        for nround = 1:length(snd_round)
            fst_con = Setup.Getcondition(R,'behpath',fullfile(subpath,sub{1},sess{1},'BEHAV',snd_round{nround}));
            % 2nd level
            fprintf('===========================================\n')
            fprintf('functional Data Step4 second level analysis\n')
            fprintf('===========================================\n')
            
            Steps.second_level(subpath,snd_round{nround},fullfile(subpath,'2nd_level','SUIT'),f, ...
                'method','flexANOVA', ...
                'infilepat','swcmT*.nii', ...
                'subidntfr',subidntfr, ...
                'conName',fst_con.conName, ...
                'conds',fst_con.m, ...
                'contrast',fst_con.contrast_2nd, ...
                'weight',fst_con.weight_2nd, ...
                'conName_2nd',fst_con.condName_2nd, ...
                'rmSUB',rmSUB)

            Steps.second_level(subpath,snd_round{nround},fullfile(subpath,'2nd_level','MNI'),f, ...
                'method','flexANOVA', ...
                'infilepat','swT*.nii', ...
                'subidntfr',subidntfr, ...
                'conName',fst_con.conName, ...
                'conds',fst_con.m, ...
                'contrast',fst_con.contrast_2nd, ...
                'weight',fst_con.weight_2nd, ...
                'conName_2nd',fst_con.condName_2nd, ...
                'rmSUB',rmSUB)
        end
    end
    spm quit
else
    %% create single subject conn project from old project
    spm fmri
    for nsub = 1:length(conn_sub)
        % ROIcor and ROIname will change in somewhere, so store it to
        % another variable
        OROIcor = ROIcor;
        OROIname = ROIname;

        if f.old_connPrj
            % using old project to create single subject project 
            conn_NprjPath = fullfile(subpath,conn_sub{nsub},'connPrj');
            conn_NprjName = conn_prjName;
            if ~exist(fullfile(subpath,conn_sub{nsub},'connPrj',conn_prjName),'dir')
                Steps.Create_From_Oprj(conn_prjPath,conn_prjName,conn_NprjPath,conn_NprjName,nsub)
            end
        else
            conn_NprjPath = fullfile(subpath,conn_sub{nsub},'connPrj');
            conn_NprjName = conn_prjName;
            if ~exist(fullfile(subpath,conn_sub{nsub},'connPrj',conn_prjName),'dir')
               conn_steps =  struct("Setup",1, ...
               "Preprocessing",1, ...
               "Denoising",1, ...
               "Add_Roi",0, ...
               "fst_Analysis",0, ...
               "snd_Analysis",0, ...
               "Add_RoiResult",0);
               Steps.create_nproj(conn_steps,conn_slice_order,...
                'sub',conn_sub(nsub), ...
                'sess',conn_ses, ...
                'round',conn_round, ...
                'condition',conn_conditon, ...
                'TR',conn_TR,...
                'conn_proj_name',conn_NprjName, ...
                'conn_proj_path',conn_NprjPath, ...
                'analysisName',conn_AnalysisName, ...
                'mrifilepath',subpath,...
                'rawdata',conn_rawdata, ...
                'prepsteps',conn_prepsteps, ...
                'StrucVres',conn_strVres,...
                'funcVres',conn_funcVres,...
                'smooth_kernel',conn_fwhm)
            end
        end

        % save ROI coordinate information
        % initial roicor
        ROItxtfile = fullfile(subpath,conn_sub{nsub},'ROI','ROI_cor.txt');
        if ~exist(ROItxtfile,'file')
            roicor = {};
        else
            roicor = readcell(ROItxtfile,Delimiter=' ');
            roicor(cellfun(@(x) any(ismissing(x)),roicor)) = {''};
            tmp = [];
            for i = 1:size(roicor,1)
                if class(roicor{i,1}) == "char" && i~=1
                    tmp = cat(1,tmp,{''});
                end
                tmp = cat(1,tmp,{horzcat(roicor{i,:})});
            end
            tmp{end+1} = '';
            roicor = tmp;
        end

        % check rois are repeat or not 
        nmidx = find(cellfun(@(x) class(x) == "char" & ~isempty(x),roicor));
        nam = roicor(nmidx);
        sidx = nmidx+1;
        eidx = find(cellfun(@isempty,roicor))-1;
        % cat new roi and check is already exist or not
        nROIcor = ROIcor;
        for i = 1:length(ROIname)
            % check ROI name is same
            idx = contains(nam,ROIname{i});
            if any(idx)
                % check ROI cor is same
                cor = roicor(sidx(idx):eidx(idx));
                for j = 1:size(ROIcor{i},1)
                    for k = 1:length(cor)
                        if all(cor{k} == ROIcor{i}(j,:))
                            nROIcor{i}(j,:) = nan(1,3);
                        end
                    end
                end
                tmp = nROIcor{i};
                tmp(isnan(tmp)) = [];
                nROIcor{i} = tmp;
            end
        end
        if any(idx)
            ROIcor(cellfun(@isempty,nROIcor)) = [];
            ROIname(cellfun(@isempty,nROIcor)) = [];
            for i = 1:length(ROIcor)
                c = cellfun(@(x) x==nROIcor{i}, ROIcor(i),'UniformOutput',false);
                c = find(sum(c{:},2)~=3);
                ROIcor{i}(c,:) = [];
            end
        end

        % add new cor in same name
        nroicor = [];
        for i = 1:length(ROIname)
            % check ROI name is same
            idx = contains(nam,ROIname{i});
            if any(idx)
                % check ROI cor is same
                cor = roicor(sidx(idx):eidx(idx));
                cor = cat(1,cor,mat2cell(ROIcor{i},ones(1,size(ROIcor{i},1)),3));
                tmp = nROIcor{i};
                tmp(isnan(tmp)) = [];
                nROIcor{i} = tmp;
                nroicor = cat(1,nroicor,ROIname(i),cor,{''});
            else
                % cat new roi
                tmp = [ROIname(i);mat2cell(ROIcor{i},ones(1,size(ROIcor{i},1)),3);{''}];
                nroicor = cat(1,nroicor,tmp);
            end
        end
        if isempty(nroicor), nroicor = roicor; end
        if ~exist(ROItxtfile,"file"), mkdir(fullfile(subpath,conn_sub{nsub},'ROI')); end
        writecell(nroicor,ROItxtfile,Delimiter=' ');

        % create Native ROI to subject
        subject_roi_path = fullfile(subpath,conn_sub{nsub},'ROI');
        for nsess = 1:length(sess)
            for nROI = 1:length(OROIcor)
                T1filepath = fullfile(subpath,conn_sub{nsub},sess{nsess},'T1');
                Steps.GetNatROI(OROIcor{nROI},T1filepath,subject_roi_path,OROIname{nROI})
            end
        end

        % create new project to single subject only need to add roi and
        % 1st level
        conn_steps =  struct("Setup",0, ...
                       "Preprocessing",0, ...
                       "Denoising",0, ...
                       "Add_Roi",1, ...
                       "fst_Analysis",1, ...
                       "snd_Analysis",0, ...
                       "Add_RoiResult",0);

        Nroi = {dir(fullfile(subject_roi_path,'*.nii')).name}';
        Nroi = split(Nroi,'.nii');
        Nroi(cellfun(@isempty,Nroi)) = [];

        Steps.create_nproj(conn_steps,conn_slice_order,...
            'sub',conn_sub(nsub), ...
            'sess',conn_ses, ...
            'round',conn_round, ...
            'condition',conn_conditon, ...
            'conn_proj_name',conn_NprjName, ...
            'conn_proj_path',conn_NprjPath, ...
            'analysisName',conn_AnalysisName, ...
            'mrifilepath',subpath,...
            'rawdata',conn_rawdata, ...
            'prepsteps',conn_prepsteps, ...
            'roiPath',subject_roi_path, ...
            'AnaSource',char(Nroi))

        % reslice 1st level result to SUIT space and MNI space
        f.func.resliceSUIT = true;
        f.func.resliceMNI = true;
        T1filepath = fullfile(subpath,conn_sub{nsub},'T1');
        T1filename = 'c_*_pcereb.nii';
        funcfilepath = fullfile(conn_NprjPath,conn_NprjName,'results','firstlevel','SBC_01');
        funcfilename = 'BETA*';
        Steps.func_step3(T1filepath,T1filename,funcfilepath,funcfilename,[],f);

        % copy condition and scource list to result folder
        FILE = dir(fullfile(funcfilepath,'_list*'));
        destPath = fullfile(subpath,conn_sub{nsub},conn_ResultF);
        if ~exist(destPath,'dir'),mkdir(destPath); end
        for nfile = 1:length(FILE)
            copyfile(fullfile(funcfilepath,FILE(nfile).name),fullfile(destPath,FILE(nfile).name));
        end

        % remove SUIT reslice and smooth data to specific result folder
        MFILE = dir(fullfile(funcfilepath,'m*'));
        SFILE = dir(fullfile(funcfilepath,'swc*'));
        WFILE = dir(fullfile(funcfilepath,'wc*'));
        FILE = cat(1,MFILE,SFILE);
        FILE = cat(1,FILE,WFILE);
        destPath = fullfile(subpath,conn_sub{nsub},conn_ResultF,'SUIT');
        if ~exist(destPath,'dir'),mkdir(destPath); end
        for nfile = 1:length(FILE)
            subid = conn_sub{nsub};
            subid = strrep(conn_sub{nsub},subidntfr,'Subject');
            Nname = strrep(FILE(nfile).name,'Subject001',subid);
            movefile(fullfile(funcfilepath,FILE(nfile).name),destPath);
            if any(FILE(nfile).name~=Nname)
                movefile(fullfile(destPath,FILE(nfile).name),fullfile(destPath,Nname));
            end
        end

        % remove MNI normalize and smooth data to result folder
        SFILE = dir(fullfile(funcfilepath,'sw*'));
        WFILE = dir(fullfile(funcfilepath,'w*'));
        FILE = cat(1,SFILE,WFILE);
        destPath = fullfile(subpath,conn_sub{nsub},conn_ResultF,'MNI');
        if ~exist(destPath,'dir'),mkdir(destPath); end
        for nfile = 1:length(FILE)
            subid = strrep(conn_sub{nsub},subidntfr,'Subject');
            Nname = strrep(FILE(nfile).name,'Subject001',subid);
            movefile(fullfile(funcfilepath,FILE(nfile).name),fullfile(destPath,FILE(nfile).name));
            if any(FILE(nfile).name~=Nname)
                movefile(fullfile(destPath,FILE(nfile).name),fullfile(destPath,Nname));
            end
        end



        fprintf('===========================================\n')
        fprintf('plot T1 and CONN REST BETA image cheg_img\n')
        fprintf('outputpath : %s\n',opath);
        fprintf('===========================================\n')
        % check mri image is same using beta map and raw T1
        oname = [conn_sub{nsub},'_REST_CONN','.jpeg'];
        Steps.checkmrifile(fullfile(subpath,conn_sub{nsub},'T1','c0SUB*T1.nii') ...
                  ,fullfile(subpath,conn_sub{nsub},'REST','dwwau*.nii'), ...
                  oname, ...
                  fullfile(opath,'REST'))

        fprintf('===========================================\n')
        fprintf('plot norm T1 and norm CONN REST BETA image cheg_img\n')
        fprintf('outputpath : %s\n',opath);
        fprintf('===========================================\n')
        % check mri image is same using norm beta map and norm raw T1
        oname = [conn_sub{nsub},'_REST_CONN_norm','.jpeg'];
        Steps.checkmrifile(fullfile(subpath,conn_sub{nsub},'T1','wc0SUB*T1.nii') ...
                  ,fullfile(subpath,conn_sub{nsub},conn_ResultF,'MNI','w*Source001.nii'), ...
                  oname, ...
                  fullfile(opath,'REST_norm'))
        ROIcor = OROIcor;
        ROIname = OROIname;

    end


    tmp = split(opath,filesep);
    tmp = strjoin(tmp,string([filesep,filesep]));
    rmSUB = {'SUB001','SUB036'};
    if ~exist("rmSUB","var")
        rmSUB = input(sprintf('check mri image from \n%s \n%s \nEnter subject name with cell array that want to remove from 2nd level \n(e.g. {''SUB001'',''SUB022''})\n: ', ...
        [tmp,filesep,filesep,'REST_conn'],[tmp,filesep,filesep,'REST_conn_norm']));
    end
    fst_con = Setup.Getcondition('rest');

    %% SUIT 2nd level pairT
    rnd = 'REST_suit';
    f.gfunc.flapmapPlot = true;
    Steps.second_level(subpath,rnd,fullfile(subpath,'2nd_level'),f, ...
        'method','pairT', ...
        'infilepat','sw*', ...
        'subidntfr',subidntfr, ...
        'conName',fst_con.conName, ...
        'rmSUB',rmSUB, ...
        'SPMorCONN','CONN', ...
        'datafolder',conn_ResultF, ...
        'space','SUIT')

    % move flapmap to same folder
    srcpath = fullfile(subpath,'2nd_level',['PairT_',rnd]);
    src = fullfile(srcpath,'**','spmT_*.jpeg');
    destPath = fullfile(srcpath,'flapmap');
    if ~exist(destPath,'dir'), mkdir(destPath); end
    srcpath = dir(src);
    for i = 1:length(srcpath)
        srcfile = fullfile(srcpath(i).folder,srcpath(i).name);
        copyfile(srcfile,destPath);
        % movefile(fullfile(destPath,'spmT_GroupLevel.jpeg'),fullfile(destPath,Nname))
    end

    %% SUIT 2nd level one sample T
    f.gfunc.flapmapPlot = true;
    Steps.second_level(subpath,rnd,fullfile(subpath,'2nd_level'),f, ...
        'method','onesampleT', ...
        'infilepat','sw*', ...
        'subidntfr',subidntfr, ...
        'conName',fst_con.conName, ...
        'rmSUB',rmSUB, ...
        'SPMorCONN','CONN', ...
        'datafolder',conn_ResultF, ...
        'space','SUIT')

    % move flapmap to same folder
    srcpath = fullfile(subpath,'2nd_level',['OnesmpT_',rnd]);
    src = fullfile(srcpath,'**','spmT_*.jpeg');
    destPath = fullfile(srcpath,'flapmap');
    if ~exist(destPath,'dir'), mkdir(destPath); end
    srcpath = dir(src);
    for i = 1:length(srcpath)
        srcfile = fullfile(srcpath(i).folder,srcpath(i).name);
        copyfile(srcfile,destPath);
        Nname = split(srcpath(i).folder,filesep);
        Nname = [Nname{end-1},'.jpeg'];
        movefile(fullfile(destPath,'spmT_GroupLevel.jpeg'),fullfile(destPath,Nname))
    end

    % get xor and & in R and L cerebellum
    srcpath = fullfile(subpath,'2nd_level',['OnesmpT_',rnd]);
    srcdir = dir(srcpath);
    srcname = {dir(srcpath).name};
    idx = contains(srcname,ROIname{1});
    for j = 2:length(ROIname)
        idx = contains(srcname,ROIname{j}) | idx;
    end
    srcdir = srcdir(idx);
    srcname = srcname(idx);
    for j = 1:length(srcname)
        tmp = split(srcname{j},'L_');
        tmp = split(tmp{end},'R_');
        srcname{j} = tmp{end};
    end
    srcname = unique(srcname);
    for i = 1:length(srcname)
        idx = find(cellfun(@(x) contains(x,srcname{i}),{srcdir.name}));
        if length(idx) > 2, error(sprintf('find more than two roi in this path %s, suppose find one Right roi and Left roi',srcpath)); end
        img1path = fullfile(srcdir(idx(1)).folder,srcdir(idx(1)).name);
        img2path = fullfile(srcdir(idx(2)).folder,srcdir(idx(2)).name);
        img1pat = 'spmT*.nii';
        img2pat = img1pat;
        outpath = fullfile(subpath,'2nd_level',['OnesmpT_',rnd],'flapmap');
        Steps.imgsetcal_suit(img1path,img1pat,img2path,img2pat,outpath)
    end
    
    %% brain 2nd level pairT
    rnd = 'REST';
    f.gfunc.flapmapPlot = false;
    Steps.second_level(subpath,rnd,fullfile(subpath,'2nd_level'),f, ...
        'method','pairT', ...
        'infilepat','sw*', ...
        'subidntfr',subidntfr, ...
        'conName',fst_con.conName, ...
        'rmSUB',rmSUB, ...
        'SPMorCONN','CONN', ...
        'datafolder',conn_ResultF, ...
        'space','MNI')

    % mask out cerebellum
    funcdir = dir(fullfile(subpath,'2nd_level',['PairT_',rnd],'**','spmT*'));
    maskfile = 'C:\Users\user\Documents\MATLAB\spm12\toolbox\suit\templates\srmaskMNI.nii';
    Steps.maskoutCere(maskfile,funcdir,'');

    %% brain 2nd level onesampleT
    f.gfunc.flapmapPlot = false;
    Steps.second_level(subpath,rnd,fullfile(subpath,'2nd_level'),f, ...
        'method','onesampleT', ...
        'infilepat','sw*', ...
        'subidntfr',subidntfr, ...
        'conName',fst_con.conName, ...
        'rmSUB',rmSUB, ...
        'SPMorCONN','CONN', ...
        'datafolder',conn_ResultF, ...
        'space','MNI')

    % mask out cerebellum
    funcdir = dir(fullfile(subpath,'2nd_level',['OnesmpT_',rnd],'**','spmT*'));
    maskfile = 'C:\Users\user\Documents\MATLAB\spm12\toolbox\suit\templates\srmaskMNI.nii';
    Steps.maskoutCere(maskfile,funcdir,'');
    spm quit

end


