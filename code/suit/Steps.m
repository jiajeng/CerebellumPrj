classdef Steps
    methods(Static)
        function cereb_stat_beh(subpath,outFolder,ftpServer)
            % subpath --> "string", subject folder path
            % outFolder --> "string", save behave file path under subject folder
            % ftpServer --> "struct", back up result to nas, if not empty(ftpServer = struct())

            subject = string({dir(subpath).name});
            subject = subject(contains(subject,"SUB"));
            feature = {'accuracy','Resp_time'};
            res = struct();
            for nsub = 1:length(subject)
                behpath = fullfile(subpath,char(subject(nsub)),outFolder);
                if ~exist(fullfile(behpath),'dir'), continue; end
                if ~exist(fullfile(behpath,'behD.mat'),'file')

                    % process behave data
                    behfilename = string({dir(behpath).name});
                    behfilename = behfilename(contains(behfilename,'.mat') & ~contains(behfilename,'behD.mat'));
                    BEH = load(fullfile(behpath,behfilename));
                    fn = char(fieldnames(BEH));
                    Result = BEH.(fn);
                    behD = struct();

                    % sort data(only for look better, not necessery)
                    sdata = [Result.Beh.YesTrial];
                    [~,i] = sort(sdata,'descend');
                    Result.Beh = Result.Beh(i);

                    % get condition name 
                    condGr = {Result.Beh.WordGroup};
                    condGr = string(cellfun(@(x) x(1:end-1),condGr,'UniformOutput',false));
                    Gr = unique(condGr);

                    % get condition idex
                    idx = arrayfun(@(x) contains(condGr,x), Gr,'UniformOutput',false)';
                    idx = cell2mat(idx);

                    % get correct trial and put others to idx new row(Error)
                    idx = cat(1,idx,false(2,size(idx,2)));
                    for i = 1:size(idx,1)
                        corR = [Result.Beh.CorrectResponse];
                        aRT = [Result.Beh.ResponseTime];
                        RT = aRT(idx(i,:));
                        RT = RT(~isnan(RT));
                        outnum = [mean(RT)-2*std(RT),mean(RT)+2*std(RT)];
                        idx(end,:) = idx(end,:) | ((aRT<outnum(1) | aRT>outnum(2)) & idx(i,:)); % outlier
                        idx(i,:) = idx(i,:) & corR==1 & aRT>outnum(1) & aRT<outnum(2); % condition remove outlier
                        idx(end-1,:) = idx(end-1,:) | corR~=1; % error trial
                    end

                    % condition --> response error or no responese are category to error
                    Gr = cat(2,Gr,["ERROR","OUTLIER"]);

                    % store result(behD) data
                    field = fieldnames(Result.Beh);
                    for nf = 1:length(field)
                        tmp = {Result.Beh.(field{nf})};
                        for i = 1:length(Gr)
                            behD.(Gr(i)).(field{nf}) = tmp(idx(i,:))';
                        end
                    end 

                    % save data to xlsx
                    outputfolder = fullfile(subpath,char(subject(nsub)),outFolder);
                    if exist(fullfile(outputfolder,'behave.xlsx'),'file')
                        delete(fullfile(outputfolder,'behave.xlsx'));
                    end

                    for nf = 1:length(Gr)
                        behD.(Gr(nf)) = struct2table(behD.(Gr(nf)));
                        behD.(Gr(nf)).Properties.VariableNames = field';
                        
                        writetable(behD.(Gr(nf)),fullfile(outputfolder,'behave.xlsx'),'Sheet',char(Gr(nf)));
                    end

                    % save .mat data
                    save(fullfile(outputfolder,'behD.mat'),"behD");

                    % save to nas
                    if ~isempty(ftpServer)
                        ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                        Targetfolder = [ftpServer.outfolder,'/',char(subject(nsub))];
                        cd(ftpobj,Targetfolder);
                        mput(ftpobj,fullfile(outputfolder,'..'))
                        close(ftpobj)
                    end
                else
                    load(fullfile(behpath,'behD.mat'))
                end

                cond = fieldnames(behD);
                cond = cond(1:4);

                % initial res struct
                if nsub == 1
                    for ncond = 1:length(cond)
                        tab = cell(length(subject),length(feature)); % acc and RT
                        res.(cond{ncond}) = tab;
                    end
                end
                % create table for all subject for accuracy and RT
                for ncond = 1:length(cond)
                    res.(cond{ncond}){nsub,1} = size(behD.(cond{ncond}),1)/30; % accuracy
                    res.(cond{ncond}){nsub,2} = mean(cell2mat(behD.(cond{ncond}).ResponseTime)); % RT
                end
            end

            % save all subject behave data
            output_folder = fullfile(subpath,outFolder);
            if ~exist(output_folder,'dir'), mkdir(output_folder); end
            if exist(fullfile(output_folder,'behave.xlsx'),'file'), delete(fullfile(output_folder,'behave.xlsx'));end
            tab = [];
            for ncond = 1:length(cond)
                tmp = res.(cond{ncond});
                tmp = cat(1,tmp,num2cell(mean(cell2mat(tmp),1)));
                tmp = cat(1,cellfun(@(x) [cond{ncond},'_',x],feature,'UniformOutput',false),tmp);
                
                tab = cat(2,tab,tmp);
                tab = cat(2,tab,cell(size(tmp,1),1));
            end
            Rowname = [{''};convertStringsToChars(subject)';{'average'}];
            tab = cat(2,Rowname,tab);
            writecell(tab,fullfile(output_folder,'behave.xlsx'));
            save(fullfile(output_folder,'res.mat'),'res');
            % save to nas
            if ~isempty(ftpServer)
                ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                ftpfolder = ftpServer.outfolder;
                cd(ftpobj,ftpfolder);
                mput(ftpobj,fullfile(output_folder,'..'))
                close(ftpobj)
            end
        end

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

        function ana_steps(filepath,filename,f)
            % segment anatomical data
            % isolate cerebullum and brainstem
            % cropping? to suit space 
            % normalize to suit templete

            % SUIT segmentation
            % =====================================================
            if f.ana.SegIso
                file = ls(fullfile(filepath,filename));
                Steps.checkfile(file,filepath,filename)
                suit_isolate_seg({fullfile(filepath,file)});
            end
            % =====================================================

            % SUIT normalize
            % ========================================================
            if f.ana.Norm
                fn = split(file,'.');
                fn = fn{1};
                mfile = ['c_',fn,'_pcereb.nii'];
                file = ['c_',fn,'.nii'];
                suit_normalize(fullfile(filepath,file),'mask',fullfile(filepath,mfile));
            end
            % ========================================================

            % get white matter and grey matter in whole brain
            % ========================================================
            if f.ana.wholemask
                mi = [1,2,7,8];
                fl = 0;
                for i = mi
                    % get c1 c2 c7 c8 
                    file = char({dir(fullfile(filepath,['c',num2str(i),'*.nii'])).name});
                    mask = spm_read_vols(spm_vol(fullfile(filepath,file)));
                    if ~fl
                        MASK = zeros(size(mask));
                        fl = 1;
                    end
                    MASK = MASK+mask;
                end
                V = spm_vol(fullfile(filepath,file));
                V.fname = fullfile(filepath,'wholebrain.nii');
                spm_write_vol(V,MASK);
            end
            % ========================================================
            
        end
        
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
                matlabbatch = [];
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

        function checkfile(file,filename,filepath)
            if size(file,1)>1
                msg = 'find more than one file: ';
                for i = 1:size(file,1)-1
                    msg = [msg,file(i),', '];
                end
                msg = [msg,file(end)];
                error(msg); 
            end
            if isempty(file), error(['can not find ',filename,' in ',filepath]); end
        end     

        function checkmrifile(T1file,funcfile,oname,opath)
            T1file = dir(T1file);
            funcfile = dir(funcfile);
            T1file = fullfile(T1file.folder,T1file.name);
            funcfile = [fullfile(funcfile.folder,funcfile.name),',1'];
            spm_check_registration(T1file,funcfile);
            if ~exist(opath,'dir'), mkdir(opath); end
            saveas(gcf,fullfile(opath,oname));
            close(gcf);
        end

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

        function second_level(filepath,round,outputpath,f,varargin)
            in = finputcheck(varargin, ...
                   {'procsNum'      'real'      []  1;
                    'conName'       'cell'      []  {};
                    'method'        'string'    {'onesampleT','flexANOVA','pairT'}  'flexANOVA';
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
                case 'pairT'
                    foldername = ['PairT_',round];
                    for ncon = 1:Ncon
                        if Ncon > 1
                            confoldname = conName{i};
                        else
                            confoldname = '';
                        end
                        pair = cat(2,convertStringsToChars(srcName(contains(srcName,'R'))),convertStringsToChars(srcName(contains(srcName,'L'))));
                        Npair = size(pair,1);
                        SCANS = struct();
                        for npair = 1:Npair
                            pairName = pair{npair,1};
                            pairName = pairName(char(pair{npair,1}) == char(pair{npair,2}));
                            idx = find(pairName =='_');
                            pairName(idx(diff(idx) == 1)) = [];
                            if pairName(1) == "_"
                                pairName(1) = [];
                            end
                            for nsub = 1:Nsubject
                                datapath = datafolder{contains(datafolder,char(subject(nsub)))};
                                pairscans = cell(2,1);
                                for i = 1:2
                                    infilepat = [in.infilepat,'*',allsrc(nsub).(pair{npair,i}),'*'];
                                    idx = find(infilepat == '*');
                                    infilepat(idx(diff(idx) == 1)) = [];
                                    confile = ls(fullfile(datapath,infilepat));
                                    pairscans{i} = fullfile(datapath,[confile,',1']);
                                end
                                SCANS(nsub).scans = pairscans;
                            end
                            outpath = char(fullfile(outputpath,foldername,confoldname,in.voi,pairName));
                            if ~exist(outpath,'dir'), mkdir(outpath); end
                            % pair T test
                            matlabbatch{end+1}.spm.stats.factorial_design.dir = {outpath};
                            matlabbatch{end}.spm.stats.factorial_design.des.pt.pair = SCANS;
                            matlabbatch{end}.spm.stats.factorial_design.des.pt.gmsca = 0;
                            matlabbatch{end}.spm.stats.factorial_design.des.pt.ancova = 0;
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
                            cname = {[pair{npair,1},'-',pair{npair,2}],[pair{npair,2},'-',pair{npair,1}],[pair{npair,1}],[pair{npair,2}]};
                            cw = {[1,-1],[-1,1],[1,0],[0,1]};
                            for i = 1:2
                                matlabbatch{end+1}.spm.stats.con.spmmat = {fullfile(outpath,'SPM.mat')};
                                matlabbatch{end}.spm.stats.con.consess{1}.tcon.name = cname{i};
                                matlabbatch{end}.spm.stats.con.consess{1}.tcon.weights = cw{i};
                                matlabbatch{end}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
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

                                    % plot flapmap
                                    % ===================================================
                                    if f.gfunc.flapmapPlot && in.space ~= "MNI"
                                        Steps.plotflapmat(outpath,'spmT*.nii')
                                    end
                                    % ===================================================  
                                end
                            end
                            matlabbatch = [];
                        end
                    end

                    
                case 'onesampleT'
                    foldername = ['OnesmpT_',round];
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
                                if Ncon > 1
                                    confoldname = conName{i};
                                else
                                    confoldname = '';
                                end
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
                            outpath = char(fullfile(outputpath,foldername,confoldname,in.voi,srcName(nsrc)));
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
                            if f.gfunc.flapmapPlot && in.space ~= "MNI"
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
                        % plot flapmap
                        % ===================================================
                        if f.gfunc.flapmapPlot && in.space ~= "MNI"
                            Steps.plotflapmat(outpath,'spmT*.nii')
                        end
                        % ===================================================  
                    end
            end

        end

        function plotflapmat(filepath,fpat)
            % filepath --> "string", filepath
            % fpat --> "string", filename pattern, e.g. 'spmT*.nii'
            FILE = dir(fullfile(filepath,'**',fpat));
            disp('start running these file ...')
            for nfiel = 1:length(FILE)
                disp(fullfile(FILE.folder,FILE.name));
            end
            for nfile = 1:length(FILE)
                file = fullfile(FILE(nfile).folder,FILE(nfile).name);
                load(fullfile(FILE(1).folder,'SPM.mat'));
                df = SPM.xX.erdf;
                fname = split(FILE(nfile).name,'.');
                if size(fname,2) == 1
                    fname = fname(1);
                else
                    fname = fname(:,1);
                end
                
                conname = split(fname,'spmT_');
                conname = conname(~cellfun(@isempty, conname));
                conname = split(conname,'_');
                conname = strjoin(conname,"\_");
                thres = 0.001;
                % file = strcat(filepath,filesep,file);
                outpath = fullfile(FILE(nfile).folder,'flapmap');
                if ~exist(outpath,'dir'), mkdir(outpath); end

                figure;
                Data = suit_map2surf(file,'space','SUIT','stats',@minORmax);
                Tthres = tinv(1-thres,df);
                % pData = tcdf(Data,df);
                suit_plotflatmap(Data,'threshold',Tthres);
                axis('off')
                stitle = sprintf('%s\nthreshold: %s',conname,num2str(thres));
                title(stitle,'FontSize',13)
                saveas(gcf,fullfile(outpath,[char(fname),'.jpeg']));
                close(gcf);
            end
        end
    
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
                    if string(in.AnaSource) == "SUBroi"
                        load(fullfile(conn_proj_path,[conn_proj_name,'.mat']))
                        Analysis_Source = CONN_x.Analysis_variables.names(cellfun(@(x) x == "roi",CONN_x.Analysis_variables.types));
                        Analysis_Source = Analysis_Source(contains(Analysis_Source,in.subidntfr));
                    else
                        Analysis_Source = convertStringsToChars(string(in.AnaSource))';
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

        function GetNatROI(MNIcor,T1filepath,ROIpath,varargin)
            % input : MNIcor --> array(nx3), MNI Space corordinate, row is
            %                                all transform corordinate 
            %         T1filepath --> string, store T1 file directory(can use **)
            %         ROIpath --> string, directory save ROI file
            %         ROI group Name --> string, in these MNIcor add ROI
            %                                    Name ,"L_cere"
            if nargin == 4
                ROINam = [varargin{1},'_'];
            elseif nargin > 4
                error('not allow enter more than two varargin input')
            else
                ROINam = '';
            end

            % get Affine variable file and load it 
            AfinMatf = 'mc_*_snc.mat';
            AfinMatf = dir(fullfile(T1filepath,AfinMatf));
            load(fullfile(AfinMatf.folder,AfinMatf.name));
            subName = split(AfinMatf.folder,filesep);
            subName = subName(contains(subName,'SUB')|contains(subName,'sub')|contains(subName,'Sub'));

            % get Affine matrix
            Q = VG.mat*inv(Affine)/VF.mat; % Native Space to MNI Space
            
            % MNI Space to Native Space corordinate
            Natcor = zeros(size(MNIcor));
            for i = 1:size(MNIcor,1) 
                tmp = [MNIcor(i,:), 1]*inv(Q');
                Natcor(i,:) = tmp(1:3);
            end
            mid = 0;
            % create Native Space ROI to ROI folder
            % check name is not duplicate
            % old ROI name 
            oroiname = {dir(fullfile(ROIpath,'*.nii')).name};
            if ~isempty(oroiname)
                oroiid = split(oroiname','_');
                if size(oroiid,2) == 1
                    oroiid = oroiid(end);
                else
                    oroiid = oroiid(:,end);
                end
                oroiid = split(oroiid,'.nii');
                oroiid(cellfun(@isempty,oroiid)) = [];
                oroiid = cellfun(@str2double,oroiid);
                mid = max(oroiid(contains(oroiname,ROINam)));
            end
            if isempty(mid), mid = 0; end
            for i = 1:size(Natcor,1)
                d = [ROINam,num2str(mid+i)]; % 'L_cere_1'

                roi = maroi_box(struct('centre',Natcor(i,:)','widths',[6;6;6]));
                roi = descrip(roi,d);
                roi = label(roi,d);
                if ~exist(ROIpath,'dir'), mkdir(ROIpath); end
                save(fullfile(ROIpath,[d,'.mat']),'roi');
                mars_rois2img(fullfile(ROIpath,[d,'.mat']), fullfile(ROIpath,[d,'.nii']), '', 'c')
            end
        end
       
        function Create_From_Oprj(OldprjPath,OldprjName,NprjPath,NprjName,Osubid)
            % OldprjName, 'char' --> old conn project name
            % OldprjPath, 'char' --> old conn project Path
            % NprjName, 'char' --> new conn project name
            % NprjPath, 'char' --> new conn project Path
            % subject, 'int' --> which subjects id in old project get to new project
        
            % Set CONN_x from old project
            CONN_x = load(fullfile(OldprjPath,[OldprjName,'.mat']));
            OldCONN_x = CONN_x.CONN_x;
            NewCONN_x = OldCONN_x;
        
            NewCONN_x.Setup.RT = OldCONN_x.Setup.RT(Osubid);
            NewCONN_x.Setup.nsubjects = length(Osubid);
            NewCONN_x.Setup.nsessions = OldCONN_x.Setup.nsessions(Osubid);
            NewCONN_x.Setup.functional = OldCONN_x.Setup.functional(Osubid);
            NewCONN_x.Setup.structural = OldCONN_x.Setup.structural(Osubid);
            NewCONN_x.Setup.spm = OldCONN_x.Setup.spm(Osubid);
            NewCONN_x.Setup.dicom = OldCONN_x.Setup.dicom(Osubid);
            NewCONN_x.Setup.nscans = OldCONN_x.Setup.nscans(Osubid);
            NewCONN_x.Setup.rois.files = OldCONN_x.Setup.rois.files(Osubid);
            NewCONN_x.Setup.conditions.values = OldCONN_x.Setup.conditions.values(Osubid);
            NewCONN_x.Setup.l1covariates.files = OldCONN_x.Setup.l1covariates.files(Osubid);
            NewCONN_x.Setup.l2covariates.values = OldCONN_x.Setup.l2covariates.values(Osubid);
            NewCONN_x.Setup.spatialresolutionvolume = OldCONN_x.Setup.spatialresolutionvolume(Osubid);
            f = fieldnames(NewCONN_x.folders);
            for nf = 1:length(f)
                NewCONN_x.folders.(f{nf}) = strrep(NewCONN_x.folders.(f{nf}),fullfile(OldprjPath,OldprjName),fullfile(NprjPath,NprjName));
                % create new conn data folders
                if ~exist(NewCONN_x.folders.(f{nf}),'dir'), mkdir(NewCONN_x.folders.(f{nf})); end
            end
            CONN_x = NewCONN_x;
            save(fullfile(NprjPath,[NprjName,'.mat']),"CONN_x");
        
            % copy file from old project
            for nsub = 1:length(Osubid)
                OsubidName = sprintf('Subject%03d',Osubid(nsub));
                file = dir(fullfile(OldprjPath,OldprjName,'**',['*',OsubidName,'*']));
                for nfile = 1:length(file)
                    srcfile = fullfile(file(nfile).folder,file(nfile).name);
                    detfpath = strrep(file(nfile).folder,fullfile(OldprjPath,OldprjName),fullfile(NprjPath,NprjName));
                    detfile = fullfile(detfpath,file(nfile).name);
                    copyfile(srcfile,detfile);
                    Nfname = strrep(file(nfile).name,OsubidName,sprintf('Subject%03d',nsub));
                    if fullfile(detfile,"")~=fullfile(detfpath,Nfname,"")
                        movefile(fullfile(detfile),fullfile(detfpath,Nfname));
                    end
                end
            end

        end
        
        function maskoutCere(maskfile,funcF,prefix)
            % mask out cerebellum
            % maskfile --> "string", mask file, absolute direction, 
            %                        e.g. 'C:\Users\user\Documents\MATLAB\spm12\toolbox\suit\templates\srmaskMNI.nii'
            % funcF --> "cell", every cell contains file that need to mask,
            %                        e.g. {'E:\Cerebellum\Data\2nd_level\REST\REST\L_cere_1\spmT_GroupLevel.nii,1'
            %                              'E:\Cerebellum\Data\2nd_level\REST\REST\L_cere_2\spmT_GroupLevel.nii,1'}
            %       --> or "struct", dir struct of all file that need to
            %                        mask
            %                        e.g. dir(fullfile(subpath,'2nd_level',rnd,'**','spmT*'));
            % prefix --> "string", outputfile prefix, if empty('') then
            %                      rename original file to ['org',filename]
            F = true;
            MASK = spm_read_vols(spm_vol(maskfile));
            for Nfile = 1:length(funcF)
                switch class(funcF)
                    case "cell"
                        FUNCFILE = funcF{Nfile};
                    case "struct"
                        FUNCFILE = [fullfile(funcF(Nfile).folder,funcF(Nfile).name),',1'];
                end
                funcV = spm_vol(FUNCFILE);
                funcData = spm_read_vols(funcV);
                % reslice template file(maskMNI.nii)
                if any(size(funcData) ~= size(MASK)) && F
                    batch{1}.spm.spatial.coreg.write.ref = {FUNCFILE};
                    batch{1}.spm.spatial.coreg.write.source = {[maskfile,',1']};
                    batch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                    batch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                    batch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                    batch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

                    batch{2}.spm.spatial.smooth.data(1) = cfg_dep('Coregister: Reslice: Resliced Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
                    batch{2}.spm.spatial.smooth.fwhm = [8 8 8];
                    batch{2}.spm.spatial.smooth.dtype = 0;
                    batch{2}.spm.spatial.smooth.im = 0;
                    batch{2}.spm.spatial.smooth.prefix = 's';
                    % run job
                    % ==================================================
                    spm('defaults', 'FMRI');
                    spm_jobman('initcfg')
                    spm_jobman('run',batch);
                    % ==================================================
                    maskfile = split(maskfile,filesep);
                    if size(maskfile,1) ~=1, maskfile = maskfile'; end
                    maskfile = cat(2,maskfile(1:end-1),{['sr',maskfile{end}]});
                    maskfile = strjoin(maskfile,filesep);
                    F = false;
                    MASK = spm_read_vols(spm_vol(maskfile));
                end
                if isempty(prefix)
                    filename = split(funcV.fname,filesep);
                    if size(filename,1) ~= 1, filename = filename'; end
                    filename = cat(2,filename(1:end-1),{['org',filename{end}]});
                    filename = strjoin(filename,filesep);
                    movefile(funcV.fname,filename);
                end
                % mask
                if Nfile == 1, MASK = MASK<0.001; end
                filename = split(funcV.fname,filesep);
                if size(filename,1) ~= 1, filename = filename'; end
                filename = cat(2,filename(1:end-1),{[prefix,filename{end}]});
                filename = strjoin(filename,filesep);
                funcV.fname =  filename;
                spm_write_vol(funcV,MASK.*funcData);
            end
        end
    
        function imgsetcal_suit(img1path,img1pat,img2path,img2pat,outpath)
            % Get suit spmT map to get intersection or difference of two set
            % get img1 file
            FILE = dir(fullfile(img1path,'**',img1pat));
            % check get one img1 file
            if length(FILE)>1
                imgemsg = 'Get more than one file for img1.';
                for i = 1:length(FILE)
                    imgemsg = cat(2,imgemsg,sprintf('\n %s ',fullfile(FILE(i).folder,FILE(i).name)));
                end
                error(imgemsg);
            end

            img1file = fullfile(FILE.folder,FILE.name);
            if exist(fullfile(FILE.folder,'SPM.mat'),"file")
                SPM = load(fullfile(FILE.folder,'SPM.mat'));
                df = SPM.SPM.xX.erdf;
                conname1 = split(FILE.folder,filesep);
                conname1 = conname1{end};
                conname1 = split(conname1,'SUB001_');
                conname1 = conname1{end};
                figname1 = conname1;
                conname1 = split(conname1,'_');
                conname1 = char(strjoin(conname1,"\_"));
            else
                df = input('enter degree of freedom : ');
                conname1 = input('enter data1 name : ');
                conname1 = split(conname1,'_');
                conname1 = char(strjoin(conname1,"\_"));
            end
            % conname = split(conname,"L_");
            % conname = split(conname,"R_");
            % conname = conname{end};
            Data1 = suit_map2surf(img1file,'space','SUIT','stats',@minORmax);
            
            % get img2 file
            FILE = dir(fullfile(img2path,'**',img2pat));
            % check get one img1 file
            if length(FILE)>1
                imgemsg = 'Get more than one file for img1.';
                for i = 1:length(FILE)
                    imgemsg = cat(2,imgemsg,sprintf('\n %s ',fullfile(FILE(i).folder,FILE(i).name)));
                end
                error(imgemsg);
            end

            img2file = fullfile(FILE.folder,FILE.name);
            if exist(fullfile(FILE.folder,'SPM.mat'),"file")
                SPM = load(fullfile(FILE.folder,'SPM.mat'));
                df = SPM.SPM.xX.erdf;
                conname2 = split(FILE.folder,filesep);
                conname2 = conname2{end};
                conname2 = split(conname2,'SUB001_');
                conname2 = conname2{end};
                figname2 = conname2;
                conname2 = split(conname2,'_');
                conname2 = char(strjoin(conname2,"\_"));
            else
                df = input('enter degree of freedom : ');
                conname2 = input('enter data2 name : ');
                conname2 = split(conname2,'_');
                conname2 = char(strjoin(conname2,"\_"));
            end
            Data2 = suit_map2surf(img2file,'space','SUIT','stats',@minORmax);
            thres = 0.001;
            Tthres = tinv(1-thres,df);
            setcal = {'union','1-2','2-1'};
            setcalname1 = {[conname1,' & ',conname2],[conname1,' & not ',conname2],['not ',conname1,' & ',conname2]};
            setcalname2 = {[conname2,' & ',conname1],[conname2,' & not ',conname1],['not ',conname2,' & ',conname1]};
            
            for i = 1:length(setcal)
                % suit_plotflatmap(Data1,'threshold',Tthres);
                % hold on;
                figure(1);Steps.suit_plotflatmap_cmp(Data1,Data2,'setcal',setcal{i},'threshold',Tthres)                
                axis('off')
                % start from here title name output path and name ...
                stitle = sprintf('%s\nthreshold: %s',setcalname1{i},num2str(thres));
                title(stitle,'FontSize',13)
                outpath1 = fullfile(outpath,figname1);
                if ~exist(outpath1,'dir'), mkdir(outpath1); end
                saveas(gcf,fullfile(outpath1,[char([figname1,'_',setcal{i}]),'.jpeg']));
                close(gcf);
                figure(2);Steps.suit_plotflatmap_cmp(Data2,Data1,'setcal',setcal{i},'threshold',Tthres)
                axis('off')
                % start from here title name output path and name ...
                stitle = sprintf('%s\nthreshold: %s',setcalname2{i},num2str(thres));
                title(stitle,'FontSize',13)
                outpath2 = fullfile(outpath,figname2);
                if ~exist(outpath2,'dir'), mkdir(outpath2); end
                saveas(gcf,fullfile(outpath2,[char([figname2,'_',setcal{i}]),'.jpeg']));
                close(gcf);
            end
        end
        
        function suit_plotflatmap_cmp(data,data2,varargin)
            % suit_plotflatmap('overlay',file,thresholds);
            % Visualised cerebellar cortical acitivty on a flatmap in a matlab window
            % INPUT:
            %    data                Vector(s) to be plotted. If overlay is NaN,
            %                        underlay is shown. Needs to be a 28935x1 vector
            % VARARGIN:
            %  'underlay',file        Specific a metric or surface_shape file that
            %                         dictates the shape of the grey underlay
            %  'underscale',[min max] Color scale to determine the value to color mapping
            %                         for the underlay
            %  'undermap',map         Color map for the underlay,a N x 3 Matrix specifying RGB values [0..1]
            %  'type':                Data type to be shown, either
            %                           'func':  Functional data
            %                           'label': discrete labels - integer values
            %                           'rgb':   A Nx3 matrix of specifying RGB values [0..1] directly
            %  'cmap',C             Color map for overlap,a N x 3 Matrix specifying RGB values [0..1]
            %  'cscale',[min max]   Color scale: determines the mapping of values to
            %                        color map
            %  'border',borderfile  Specifies a borderfile to plot. Use [] for no borders
            %  'borderstyle','k.'   Color and marker style string for the borders
            %  'bordersize',pt      Border point size in pt (default 8)
            %  'threshold',val      Shows only values above a certain threshold
            %  'xlims',[xmin xmax]  Limits of the x-coordinate
            %  'ylims',[ymin ymax]  Limits of the y-coordinate
            %  'alpha',val          Alpha-value (Opacity) of the overlay 
            %
            % (c) joern.diedrichsen@googlemail.com, 2014, 2019
            % EXAMPLES:
            % 1. Plot a functional volume at a certain treshold + and color scale
            %   Data = suit_map2surf('name.nii','space','SUIT');
            %   suit_plotflatmap(Data,'threshold',0.01,'cscale',[0.01 0.5]);
            % 2. Plot flatmap overlayed with the lobules
            %   Data = suit_map2surf('Cerebellum-lobules.paint','stats',@(x)nanmedian(x,2));
            %   suit_plotflatmap(Data,'interpol',0,'cmap','Cerebellum-lobules.color.txt');
            % 3. Plot flatmap with the 17 resting state networks by Buckner
            %   Data = suit_map2surf('Buckner_17Networks.paint');
            %   suit_plotflatmap('overlay',Data,'cmap','Buckner_17Networks.color.txt');
            % -----------------------------------------------------------------
            
            % -----------------------------------------------------------------
            %   Set Default values and deal with user options
            % -----------------------------------------------------------------
            global defaults;
            if (isempty(defaults))
                spm('Defaults','fmri');
                global defaults;
            end;
            flat_dir   = [];                    %
            surf       = 'FLAT.surf.gii';       % Surface file for flat map
            underlay   = 'SUIT.shape.gii';      % File determining colring of underlay
            underscale = [-1 0.5];              % Color scale [min max] for the underlay
            undermap   = gray;                    % Color map for underlay
            type       = 'func';                % 'func': funtional activation 'label': categories 'rgb': RGB values
            threshold  = [];                    % Threshold for functional overlay
            cscale     = [];                    % Color scale [min max] for the overlay
            cmap       =  colormap;             % Use current colormap by default
            border     = 'fissures_flat.mat';   % File containing fissure information
            bordersize = 8;                     % Size of the border points
            borderstyle = 'k.';                 % Color and symbol for border points 
            xlims      =  [-100 100];           % X-limits for Figure
            ylims      =  [-100 100];           % Y-limits for Figure
            alpha      = 1;                     % Opacity of the overlay 
            setcal     = 'union';               % 'union' : union of two set, '1-2' : data1 contain but data2 not contain, '2-1' : data2 contains but data1 not contains
            vararginoptions(varargin,{'coord','topo','underlay','underscale','undermap',...
                'type','threshold','cscale','cmap','border','bordersize','borderstyle','xlims','ylims',...
                'flat_dir','alpha','setcal'});
            
            SCCSid   = '3.1';
            SPMid    = spm('FnBanner',mfilename,SCCSid);
            
            % -----------------------------------------------------------------
            %  Determine directory for the flatmap files
            % -----------------------------------------------------------------
            spmVer=spm('Ver');
            if ~strcmp(spmVer,'SPM12')
                error('plot flatmap requires SPM12 (for gifti support)');
            end;
            
            if (isempty(flat_dir))
                flat_dir=fullfile(defaults.tbx.dir{1},'suit','flatmap');
            end;
            
            % -----------------------------------------------------------------
            %   Load the surface and determine X,Y coordinates for all tiles
            % -----------------------------------------------------------------
            C=gifti(fullfile(flat_dir,surf));
            P=size(C.vertices,1);
            X(1,:)=C.vertices(C.faces(:,1),1);
            X(2,:)=C.vertices(C.faces(:,2),1);
            X(3,:)=C.vertices(C.faces(:,3),1);
            Y(1,:)=C.vertices(C.faces(:,1),2);
            Y(2,:)=C.vertices(C.faces(:,2),2);
            Y(3,:)=C.vertices(C.faces(:,3),2);
            
            % Find all tiles that have a single vertex (or more) in the image
            k=find(any(X>xlims(1) & X<xlims(2),1) & any(Y>ylims(1) & Y<ylims(2),1));
            X=X(:,k);
            Y=Y(:,k);
            
            % -----------------------------------------------------------------
            %  Determine the underlay and assign color
            % -----------------------------------------------------------------
            UN=gifti(fullfile(flat_dir,underlay));
            
            % Determine the shading of the faces by the vertices and scale the color assignment
            d=[UN.cdata(C.faces(k(:),1),1) UN.cdata(C.faces(k(:),2),1) UN.cdata(C.faces(k(:),3),1)]';
            M=size(undermap,1);
            dindx=ceil((d-underscale(1))/(underscale(2)-underscale(1))*M);
            dindx(dindx<1)=1;
            dindx(dindx>M)=M;
            
            % Now assign the RGB color to each face
            for i=1:3 % Color channel
                for j=1:size(dindx,1) % 1-3rd vertex
                    COL(j,:,i) = undermap(dindx(j,:),i)';
                end;
            end;
            
            % -----------------------------------------------------------------
            %  determine the overlay and assign color
            % -----------------------------------------------------------------
            
            % If input data is empty, make vector of NaNs;
            if (isempty(data))
                data=nan(P,1);
            end;
            if (isempty(data2))
                data=nan(P,1);
            end;
            
            % Check input data
            if (size(data,1))~=P
                error(sprintf('Input data must be a numVert (%d) x 1 vector',P));
            end;
            if (size(data2,1))~=P
                error(sprintf('Input data must be a numVert (%d) x 1 vector',P));
            end;
            
            % Determine colormap
            if (ischar(cmap))
                CM=load(cmap);
                cmap=[CM(:,2) CM(:,3) CM(:,4)]/255;
            end;
            colormax = size(cmap,1);
            
            
            % Determine the color of the overlay
            OCOL = nan(size(COL));    % set all values to NaN
            switch (type)
                case 'func'
                    % Otherwise take the mean value
                    d=[data(C.faces(k(:),1),1) data(C.faces(k(:),2),1) data(C.faces(k(:),3),1)]';
                    d2=[data2(C.faces(k(:),1),1) data2(C.faces(k(:),2),1) data2(C.faces(k(:),3),1)]';
                    if (isempty(cscale) || any(isnan(cscale)))
                        cscale=[min(d(:)) max(d(:))];
                    end;
                    dindx=ceil((d-cscale(1))/(cscale(2)-cscale(1))*colormax);
                    dindx(dindx<=0)=1;
                    dindx(dindx>colormax)=colormax;
                    d2indx=ceil((d-cscale(1))/(cscale(2)-cscale(1))*colormax);
                    d2indx(d2indx<=0)=1;
                    d2indx(d2indx>colormax)=colormax;
                    for i=1:3 % Color channnel
                        for j=1:size(dindx,1)
                            % Apply map thresholding
                            if (isempty(threshold) || isnan(threshold))
                                indx = find(dindx(j,:)>0 & dindx(j,:)<=colormax);
                            else
                                switch setcal
                                    case 'union'
                                        indx = find(d(j,:)>threshold & d2(j,:)>threshold);  
                                    case '1-2'
                                        indx = (d(j,:)>threshold) - (d2(j,:)>threshold);
                                        indx(indx~=1) = 0;
                                        indx = find(indx);
                                    case '2-1'
                                        indx = (d2(j,:)>threshold) - (d(j,:)>threshold);
                                        indx(indx~=1) = 0;
                                        indx = find(indx);
                                    otherwise
                                        indx = find(dindx(j,:)>0 & dindx(j,:)<=colormax);
                                end
                            end;
                            OCOL(j,indx,i) = cmap(dindx(j,indx),i)';
                        end;
                    end;
                otherwise
                    error('unknown data type');
            end;

            % Now blend the overlay with the underlay: 
            indx = ~isnan(OCOL); 
            if (isscalar(alpha))
                COL(indx) = alpha.*OCOL(indx)+(1-alpha)*COL(indx); % Set opacacy of overlay 
            else 
                for i=1:3 
                    ALPHA(:,:,i) = [alpha(C.faces(k(:),1),1) alpha(C.faces(k(:),2),1) alpha(C.faces(k(:),3),1)]';
                end; 
                COL(indx) = ALPHA(indx).*OCOL(indx)+(1-ALPHA(indx)).*COL(indx); 
            end; 
            
            % -----------------------------------------------------------------
            %   Draw the actual patches
            % -----------------------------------------------------------------
            drawmode=get(gca,'NextPlot');
            if (~strcmp(drawmode,'add'))
                cla;
            end;
            p=patch(X,Y,COL);
            set(gca,'XLim',xlims,'YLim',ylims,'XTick',[],'YTick',[]);
            set(p,'LineStyle','none');
            set(p,'EdgeAlpha',0);
            axis equal;
            
            % -----------------------------------------------------------------
            %   Plot the border
            % -----------------------------------------------------------------
            if (~isempty(border))
                hold on;
                load(fullfile(flat_dir,border));
                for i=1:length(Border)
                    xB=Border(i).data(:,1);
                    yB=Border(i).data(:,2);
                    indx=find(xB>xlims(1) & xB<xlims(2) & yB>ylims(1) & yB<ylims(2));
                    p=plot(xB(indx),yB(indx),borderstyle);
                    set(p,'MarkerSize',bordersize);
                end;
                if (~strcmp(drawmode,'add'))
                    hold off;
                end;
            end;
        end

    end
    
end