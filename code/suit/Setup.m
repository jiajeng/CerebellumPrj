classdef Setup
    methods(Static)
        function [o,f] = SetVariable
            % =========flags 
            % local --> local data or nas data
            % raw --> need to transform to nifti file 
            % ana --> run anatomical file 
            % func --> run functional file 
            % ow --> overwrite spm.mat file 
            f.local = false; % local data or nas data
            f.raw = true; % need to transform to nifti file
            f.beh = true; % get behave data and organized data
            f.GetData = false;
            f.conn = true;
            f.old_connPrj = true;
            % anatomical
            % -------------------------------
            f.ana.SegIso = true; % run anatomical segmentation and isolation
            f.ana.Norm = true; % run anatomical normalize (to SUIT space)
            % ------------ no need for now -----------------
            f.ana.wholemask = false; % get whole brain mask 
            % ----------------------------------------------
            fn = fieldnames(f.ana);
            fn = cellfun(@(x) f.ana.(x),fn);

            ana = true; % run anatomical file
            f.ana.ana = any(fn) && ana;
            % ------------------------------

            % idividual level 
            % -----------------------------------------
            f.func.relign = true; % run functinoal prep realignment
            f.func.sliTim = true; % run functional prep slice Timing
            f.func.corig = true; % run functional perp corigister to structure 
            f.func.fstL = true; % run functional 1st level analysis
            f.func.fstcon = true; % run functional 1st level contrast
            f.func.fstconP = true; % positive contrast [1]
            f.func.fstconN = false; % negative contrast [-1]
            f.func.resliceSUIT = true; % run functional contrast reslice to structure
            f.func.resliceMNI = true; % run functional normalized 
            f.func.smooth = true; % run functional data smooth

            fn = fieldnames(f.func);
            fn = cellfun(@(x) f.func.(x),fn);

            func = true; % run functional file
            f.func.func = any(fn) && func;

            % group level analysis
            % ---------------------------------------------
            f.gfunc.sndL = true;
            f.gfunc.flapmapPlot = true;
    
            fn = fieldnames(f.gfunc);
            fn = cellfun(@(x) f.gfunc.(x),fn);
            
            gfunc = false; % run functional file
            f.gfunc.func = any(fn) && gfunc;
            % ---------------------------------------------

            % -------------------------------------------

            f.ow = 1; % for spm.mat overwrite
            % ===============================
       
            % ========================= nas server parameter 
            % ip --> nas ip
            % account --> nas account
            % password --> nas password
            % infolder --> in which folder to get data
            % outfolder --> push nifti file to which folder 
            o.ftpServer.ip = 'ftp://xxxxxxxxxxx/';
            o.ftpServer.account = 'xxxxxxxxxxx';
            o.ftpServer.password = 'xxxxxxxxxxx';
            o.ftpServer.infolder = 'LabData/cerebellum_NYCU/rawdata'; % Nas Raw Data folder 
                                     %   e.x.LabData/廣達人腦健康資料庫
            o.ftpServer.outfolder = 'LabData/jeng/Cerebellum_suit/rawfile'; % Nas put nii file folder
                                      %   e.x.LabData/jeng/test/Data
            % ==============================================
        
            % =================define subject name and sess name
            % sess --> session, like sess01 ,sess02, 
            % round --> round, like T1 REST WORD ... 
            % folderNest --> in subject folder the folderNest order
            %                e.x. 'sess/mri/round' --> sess01/mri/T1 or sess02/mri/REST
            % subpath --> subject direction
            % sub --> subejct name
            % outputfolder --> save data local folder
            % R --> condition case name
            % subidntfr --> indentefier subject 
            %               e.x sub001 identefier = 'sub'
            o.subidntfr = 'SUB';
            o.R = 'WORD_cere';
            o.sess = {''};
            o.round = {'WORD'};
            o.FolderNest = {'round'};
            o.subpath = 'E:\Cerebellum\Data_co\';
            if f.local
                sub = {dir(o.subpath).name};
                o.sub = sub(contains(sub,o.subidntfr));
            else
                ftpobj = ftp(o.ftpServer.ip,o.ftpServer.account,o.ftpServer.password);
                cd(ftpobj,o.ftpServer.infolder);
                sub = {dir(ftpobj).name};
                sub = sub(contains(sub,o.subidntfr));
                sub = split(sub','_');
                sub = sub(contains(sub,o.subidntfr));
                o.sub = unique(sub);
            end
            o.opath = fullfile(o.subpath,'checkMRIimage');
            % ==============================================

            % ========================behave file variable
            % BEH --> "struct", fieldname is Behave name, e.x. Word
            %                   contains subpath, "string", subject folder path
            %                            folder, "string", behave folder path under sub folder
            %                            ext, "string", file extension
            %                     option idn, "cell"--"string", file contains certain string
            o.BEH.Word.subpath = '/LabData/cerebellum_NYCU/rawdata';% subpath
            o.BEH.Word.folder = '/BEHAV/';% folder under the subject folder
            o.BEH.Word.ext = '.mat';% file extension

            % o.BEH.Physio.subpath = '/LabData/cerebellum_NYCU/prep'; % subpath
            % o.BEH.Physio.folder = 'mri_physio'; % folder under the subject folder 
            % o.BEH.Physio.ext = '.txt'; % file extension
            % o.BEH.Physio.idn = {'Word','Rest'}; % any other identifier or identify file using file extension
            % ===================================================

            % ======================== ROI define
            o.ROIpath = ['E:\Cerebellum\ROI\'];
            o.ROIcor = {[20 -68 -47
                        8 -78 -31
                        24 -72 -48
                        10 -80 -24
                        32 -64 -26
                        22 -68 -26
                        42 -50 -30
                        36 -56 -24
                        38 -70 -28
                        24 -64 -28
                        6 -80 -24
                        24 -74 -50
                        12 -67 -30
                        17 -65 -35], ...
                        [-20 -68 -47
                        -8 -78 -31
                        -24 -72 -48
                        -10 -80 -24
                        -32 -64 -26
                        -22 -68 -26
                        -42 -50 -30
                        -36 -56 -24
                        -38 -70 -28
                        -24 -64 -28
                        -6 -80 -24
                        -24 -74 -50
                        -12 -67 -30
                        -17 -65 -35]};
            o.ROIname = {'R_cere','L_cere'};
            if length(o.ROIcor) ~= length(o.ROIname)
                warning('ROIcor must has same length with ROIName, set ROIname to g1, g2, ...');
                o.ROIname = string(strcat('g',num2str([1:length(o.ROIcor)]')));
                o.ROIname = convertStringsToChars(o.ROIname)';
            end
            % ===========================================

            % ===================== conn project variable

            o.conn_sub = o.sub;
            o.conn_ses = '';
            o.conn_round = {'WORD'};
            o.conn_prjName = ['conn_',cell2mat(o.conn_round)];
            o.conn_prjPath = 'E:\Cerebellum\connPrj';

            % set conn steps
            o.conn_steps =  struct("Setup",0, ...
                           "Preprocessing",0, ...
                           "Denoising",0, ...
                           "Add_Roi",0, ...
                           "fst_Analysis",0, ...
                           "snd_Analysis",0, ...
                           "Add_RoiResult",0);
            
             if exist(fullfile(o.subpath,'MRinfo.xlsx'),'file')
                round = [o.conn_round,{'T1'}];
                empsub = cell(1,length(round));

                for roundi = 1:length(round)
                    info = readtable(fullfile(o.subpath,'MRinfo.xlsx'),'Sheet',round{roundi}); 
                    empsub{roundi} = string(info.subject)=="";
                end
                if any(diff(cellfun(@length,empsub))~=0)
                    nsub = max(cellfun(@length,empsub));
                    for i = 1:length(round)
                        if length(empsub{i})<nsub
                            empsub{i} = cat(1,empsub{i},true);
                        end
                    end
                end
                empsub = any(cell2mat(empsub),2);
                info = readtable(fullfile(o.subpath,'MRinfo.xlsx'),'Sheet',round{1}); 
                info = info(~empsub,:);

                o.conn_sub = info.subject;
                o.conn_sub = o.sub;
                o.conn_TR = unique(str2double(info.RT));
                o.conn_slice_order = str2num(unique(string(info.sliceorder))); % use str2num instead of str2double, using str2double will get nan

                if size(o.conn_slice_order,1) > 1, o.conn_slice_order = o.conn_slice_order'; end
                if size(o.conn_TR,1) > 1
                    error('infoData "TR" has different between subjects');
                end
                if size(o.conn_slice_order,1) > 1
                    warning('infoData "sliceorder" has different between subject');
                    o.conn_slice_order = '';
                end
            else
                o.conn_TR = 2.4;
                o.conn_slice_order = 'interleaved (bottom-up)';
            end

            o.conn_fwhm = 8;
            o.conn_strVres = 1;
            o.conn_funcVres = 2;
            o.conn_filtBand = [0.008, 0.09];
            o.conn_AnalysisName = 'SBC_01';
            conn_AnaSource = cell(1,length(o.ROIname)*length(o.ROIcor));
            for i = 1:length(o.ROIname)
                tmp = o.ROIname{i};
                for j = 1:length(o.ROIcor)
                    conn_AnaSource{(i-1)*i+j} = [tmp,'_',num2str(j)];
                end
            end
            o.conn_AnaSource = conn_AnaSource; % == SUBroi --> select roi contains "SUB" char
            o.conn_contrast = {1};
            o.conn_ROIpath = o.ROIpath;
            o.conn_conditon = o.conn_round;
            o.conn_rawdata = false;
            o.conn_prepsteps = 'conn_prepSteps.mat';
            o.conn_ResultF = ['connWORD_1stresult'];
 
            % ===========================================



        end

        function o = Getcondition(r,varargin)
            % defined condition in 1st level analysis
            % condition --> "struct", set condition in 1st level analysis
            %
            % defined contrast parameter
            % con --> "cell", every cell contains a contrast, name is same
            %                 as condition, 
            %                 e.g. condition A - condition B and condition C - condition D -->
            %                 {{'A','B'},{'C','D'}}
            % conWt --> "cell", every cell contains contrast weight, same 
            %                   length as con
            %                  e.g. condition A - condition B and condition C - condition D -->
            %                 {{1,-1},{1,-1}}
            % condName --> "cell", every cell means each contrast name,
            %                      same length as con and conWt
            %                  e.g. {'A-B','C-D'}
            % condi --> "real", 1 or 0, means set all condition is a
            %                   contrast or not
            %                  e.g. {{'A'},{'B'},{'C'},...} --> {{1},{1},{1},...}
            switch r
                case 'rest'
                    o.condition = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                    % o.con = {{'constant'}};
                    % o.conWt = {{1}};
                    o.con = {};
                    o.conWt = {};
                    o.conName = {'Rest'};
                    o.condi = 1;
                case 'WORD_quanta'
                    inName = varargin(1:2:end);
                    inVar = varargin(2:2:end);
                    try
                        behpath = inVar(string(inName)=="behpath");
                    catch
                        error('need input "behpath" in WORD task condition')
                    end

                    % 1st level
                    o.conName = {'HFHC','HFLC','LFHC','LFLC'};
                    o.condition = Setup.BEHcond(behpath,'behD.mat','ResponseTime');

                    % 2nd level 
                    o.contrast_2nd = {{'HFHC','HFLC','LFHC','LFLC'},{'HFHC','HFLC','LFHC','LFLC'}};
                    o.weight_2nd = {[1,1,-1,-1],[1,-1,1,-1],[1,1,-1,-1]*(-1),[1,-1,1,-1]*(-1)};
                    o.condName_2nd = {'HF-LF','HC-LC','LF-HF','LC-HC'};
                    
                    % define flexible ANOVA coditions matrix
                    ANOVA_feature = {'HF','LF';
                                     'HC','LC'};
                    feaCod = [1,2;
                              1,2];
                    
                    o.m = zeros(length(o.conName),size(ANOVA_feature,2));
                    for i = 1:length(o.conName)
                        for j = 1:size(ANOVA_feature,1)
                            c = cellfun(@(x) contains(o.conName{i},x),ANOVA_feature(j,:));
                            o.m(i,j) = feaCod(j,c);
                        end
                    end
                    o.condi = 1;
                case 'WORD_cere'
                    inName = varargin(1:2:end);
                    inVar = varargin(2:2:end);
                    
                    try
                        behpath = inVar{string(inName)=="behpath"};
                    catch
                        error('need input "behpath" in WORD task condition')
                    end

                    % 1st level
                    o.conName = {'HDHF','HDLF','LDHF','LDLF'};
                    % BEHcond(behpath,behfilename,parName,onsetname,dur,timModul)
                    o.condition = BEHcond(behpath,'behD.mat',{'ResponseTime'},'SoundOnset',0,0);
                    o.con = {};
                    o.conWt = {};
                    % 2nd level
                    o.contrast_2nd = {{'HDHF','HDLF','LDHF','LDLF'},{'HDHF','HDLF','LDHF','LDLF'},{'HDHF','HDLF','LDHF','LDLF'},{'HDHF','HDLF','LDHF','LDLF'}};
                    o.weight_2nd = {[1,-1,1,-1],[1,1,-1,-1],[1,-1,1,-1]*(-1),[1,1,-1,-1]*(-1)};
                    o.condName_2nd = {'HF-LF','HD-LD','LF-HF','LD-HD'};
                    
                    % define flexible ANOVA coditions matrix
                    ANOVA_feature = {'HD','LD';
                                     'HF','LF'};
                    feaCod = [1,2;
                              1,2];
                    
                    o.m = zeros(length(o.conName),size(ANOVA_feature,2));
                    for i = 1:length(o.conName)
                        for j = 1:size(ANOVA_feature,1)
                            c = cellfun(@(x) contains(o.conName{i},x),ANOVA_feature(j,:));
                            o.m(i,j) = feaCod(j,c);
                        end
                    end
                    o.condi = 1;
                otherwise
        
            end

            function con = BEHcond(behpath,behfilename,parName,onsetname,dur,timModul)
                % find file
                behfile = {dir(behpath).name};
                behfile = behfile{contains(behfile,behfilename)};
                behD = load(fullfile(behpath,behfile));
                fn = char(fieldnames(behD));
                behD = behD.(fn);
                condition = fieldnames(behD);
                condition = condition(~(contains(condition,'OUTLIER')| ...
                                        contains(condition,'OutLier')| ...
                                        contains(condition,'outlier')));
                con = repmat(struct(),4,1);
                for ncond = 1:length(condition)
                    for npar = 1:length(parName)
                        con(ncond).name = condition{ncond};
                        con(ncond).onset = cell2mat([behD.(condition{ncond}).(onsetname)]);
                        con(ncond).duration = dur;
                        con(ncond).tmod = timModul;
                        con(ncond).orth = 1;
                
                        param = cell2mat([behD.(condition{ncond}).(parName{npar})]);
                        if condition{ncond} == "ERROR"
                            con(ncond).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        else
                            con(ncond).pmod(npar).name = parName{npar};
                            con(ncond).pmod(npar).param = param;
                            con(ncond).pmod(npar).poly = 1;
                        end
                    end
                end
            end
        end

        function GetBehFolder(BEH,sub,sess,ftpServer)
            BEHname = fieldnames(BEH);
            for nsess = 1:length(sess)
                if ~exist('BEHFold_log.mat','file')
                    log = cell2table(cell(length(unique(sub)),length(BEHname))); 
                    log.Properties.VariableNames = strcat(BEHname,'folder');
                    RowNames = 1:length(unique(sub));
                    log.Properties.RowNames = convertStringsToChars(string(num2str(RowNames')));
                    existF = 0;
                else
                    load('BEHFold_log.mat');
                    existF = 1;
                end
                for nsub = 1:length(sub)
                    for nBEH = 1:length(BEHname)
                        if ~isempty(fieldnames(ftpServer))
                            subpath = BEH.(BEHname{nBEH}).subpath;
                            ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                            cd(ftpobj,BEH.(BEHname{nBEH}).subpath);
                            SUBname = string({dir(ftpobj).name}');
                            SUBname = SUBname(contains(SUBname,sub{nsub}));
                            close(ftpobj);
                            % check if has same subject name 
                            if length(unique(SUBname)) > 1
                                flag = ones(1,length(SUBname));
                                % check subject name does have wanted folder
                                for nsubName = 1:length(SUBname)
                                    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                                    try
                                        BEH.(BEHname{nBEH}).folder(BEH.(BEHname{nBEH}).folder=='\') = '/';
                                        cd(ftpobj,[BEH.(BEHname{nBEH}).subpath,'/',char(SUBname(nsubName)),'/',BEH.(BEHname{nBEH}).folder]);
                                    catch
                                        flag(nsubName) = 0;
                                    end
                                    close(ftpobj);
                                end
                                % if all subject name has same folder, then get the latest one(_n)
                               if sum(flag) > 1
                                    tmp = split(SUBname,'_');
                                    tmp = str2double(tmp(~(contains(tmp,'SUB')|contains(tmp,'sub')|contains(tmp,'Sub'))));
                                    
                                    SUBname = SUBname(tmp == max(tmp));
                                else    
                                    SUBname = SUBname(logical(flag));
                                end
                            end
                        else
                            wd = pwd;
                            subpath = BEH.(BEHname{nBEH}).subpath;
                            ftpobj = BEH.(BEHname{nBEH}).subpath;
                            SUBname = string({dir(ftpobj).name}');
                            SUBname = SUBname(contains(SUBname,sub));
                            % check if has same subject name 
                            if length(unique(SUBname)) > 1
                                flag = ones(1,length(SUBname));
                                % check subject name does have wanted folder
                                for nsubName = 1:length(SUBname)
                                    try
                                        cd(fullfile(ftpobj,char(SUBname(nsubName)),BEH.(BEHname{nBEH}).folder));
                                    catch
                                        flag(nsubName) = 0;
                                    end
                                    cd(wd);
                                end
                                % if all subject name has same folder, then get the latest one
                                if sum(flag) > 1
                                    tmp = split(SUBname,'_');
                                    tmp = str2double(tmp(~(contains(tmp,'SUB')|contains(tmp,'sub')|contains(tmp,'Sub'))));
                                    SUBname = SUBname(tmp == max(tmp));
                                else    
                                    SUBname = SUBname(logical(flag));
                                end
                            end
                            cd(wd);
                        end
                        if existF, subNum = size(log,1)+1; else subNum = nsub; end
                        log(subNum,[BEHname{nBEH},'folder']) = array2table({fullfile(subpath,SUBname,BEH.(BEHname{nBEH}).folder)});
                    end
                end
                save(['BEHFold_log',sess{nsess},'.mat'],"log");
            end
        end
    
        function GetBehfile(BEH,ftpServer,subidntfr,localpath)
           
            % get folder log 
            foldlog = load("BEHFold_log.mat");
            BEHname = foldlog.log.Properties.VariableNames;
            BEHname = split(BEHname','folder');
            BEHname(cellfun(@isempty,BEHname)) = [];
            
            % get subject name 
            for i = 1:size(foldlog.log,2)
                sub = table2array(foldlog.log(:,i));
                sub = split(string(sub),filesep);
                sub = sub(contains(sub,subidntfr));
                sub = split(sub,'_');
                sub = sub(contains(sub,subidntfr));
                osub = sub;
                if osub ~= sub
                    error('subject folder Name is not same, check GetBehFolder function get correct folder');
                end
            end

            % save log file
            log = struct();

            for nBEH = 1:length(BEHname)
                if ~isempty(fieldnames(ftpServer))
                    % get file in Nas
                    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                    for nsub = 1:length(sub)
                        try
                            TarFolder = char(string(table2array(foldlog.log(nsub,[BEHname{nBEH},'folder']))));
                            TarFolder(TarFolder=='\') = '/';
                            cd(ftpobj,TarFolder);
                            tmp = {dir(ftpobj).name};
                            flag = contains(tmp,lower(BEHname{nBEH})) | ...
                                   contains(tmp,BEHname{nBEH}) | ...
                                   contains(tmp,upper(BEHname{nBEH}));
                            file = tmp(contains(tmp,BEH.(BEHname{nBEH}).ext) & flag);
                            flag = flag & ~(contains(tmp,BEH.(BEHname{nBEH}).ext));
                            if isempty(file)
                                if any(flag)
                                    cd(ftpobj,tmp(flag));
                                    tmp = {dir(ftpobj).name};
                                    flag = contains(tmp,lower(BEHname{nBEH})) | ...
                                       contains(tmp,BEHname{nBEH}) | ...
                                       contains(tmp,upper(BEHname{nBEH}));
                                    file = tmp(contains(tmp,BEH.(BEHname{nBEH}).ext) & flag);
                                end     
                            end

                            if length(file)>1
                                flag = false(size(file));
                                for i = 1:length(BEH.(BEHname{nBEH}).idn)
                                    idn = BEH.(BEHname{nBEH}).idn{i};
                                    flag = flag | (contains(file,idn)| ...
                                                    contains(file,upper(idn))| ...
                                                    contains(file,lower(idn)));
                                end
                                file = file(flag);
                            end
                            
                            for nfile = 1:length(file)
                                localfolder = fullfile(localpath,char(sub(nsub)),'BEHAV',BEHname{nBEH});
                                if ~exist(localfolder,'dir'), mkdir(localfolder);  end
                                mget(ftpobj,char(file{nfile}),localfolder)
                                try
                                    for i = 1:length(BEH.(BEHname{nBEH}).idn)
                                        idn = BEH.(BEHname{nBEH}).idn{i};
                                        flag = contains(file{nfile},idn)| ...
                                               contains(file{nfile},upper(idn))| ...
                                               contains(file{nfile},lower(idn));
                                        if flag
                                            idn = BEH.(BEHname{nBEH}).idn{i};
                                            break;
                                        end
                                    end
                                    s = fullfile(localfolder,char(file{nfile}));
                                    d = fullfile(localfolder,[char(sub(nsub)),'_',BEHname{nBEH},'_',idn,BEH.(BEHname{nBEH}).ext]);
                                    if isempty(fieldnames(log))
                                        log.([BEHname{nBEH}]) = cell(length(sub),1);
                                    end
                                    if ~contains(fieldnames(log),[BEHname{nBEH},'_',idn])
                                        log.([BEHname{nBEH},'_',idn]) = cell(length(sub),1);
                                    end
                                    log.([BEHname{nBEH},'_',idn]){nsub} = d;
                                catch
                                    s = fullfile(localfolder,char(file{nfile}));
                                    d = fullfile(localfolder,[char(sub(nsub)),'_',BEHname{nBEH},BEH.(BEHname{nBEH}).ext]);
                                    if isempty(fieldnames(log))
                                        log.([BEHname{nBEH}]) = cell(length(sub),1);
                                    end
                                    if ~contains(fieldnames(log),[BEHname{nBEH}]) 
                                        log.([BEHname{nBEH}]) = cell(length(sub),1);
                                    end
                                    log.([BEHname{nBEH}]){nsub} = d;
                                end
                                movefile(s,d);
                            end
                        catch ME
                            continue;
                        end
                    end
                else
                    % check folder has file or not
                    for nBEH = 1:length(BEHname)
                        
                    end
                end
            end
            log = struct2table(log);
            save("BEHFile_log.mat","log");
        end
    end
end