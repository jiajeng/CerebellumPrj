function [matlabbatch] = second_level(filepath,outputpath,matlabbatch,varargin)
in = finputcheck(varargin, ...
       {'procsNum'      'real'      []  [];
        'conName'       'cell'      []  {};
        'method'        'string'    {'onesampeT','flexANOVA'}  'flexANOVA';
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
       });

subject = string({dir(filepath).name});
subject = subject(contains(subject,'SUB'));

Nsubject = length(subject);


datapath = fullfile(filepath,char(subject(1)),in.analysis,in.voi,in.diff);
load(fullfile(datapath,'SPM.mat'));
oldcondName = string({SPM.xCon.name});

for nsub = 1:Nsubject
    datapath = fullfile(filepath,char(subject(nsub)),in.analysis,in.voi,in.diff);
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
Ncon = length(conName);

switch in.method
    case 'onesampeT'
        scans = cell(Nsubject,Ncon);
        
        for nsub = 1:Nsubject
            datapath = fullfile(filepath,char(subject(nsub)),in.analysis,in.voi,in.diff);
            for ncon = 1:Ncon
                con = [in.TorF,'_',char(conName(ncon)),'.nii'];
                scans{nsub,ncon} = fullfile(datapath,con);
            end
        end

        procs_Num = in.procsNum;
        for ncon = 1:Ncon
            outpath = char(fullfile(outputpath,conName(ncon),in.voi));
            if ~exist(outpath,'dir'), mkdir(outpath); end
            matlabbatch{procs_Num}.spm.stats.factorial_design.dir = {outpath};
            matlabbatch{procs_Num}.spm.stats.factorial_design.des.t1.scans = scans(:,ncon);
            % matlabbatch{in.procsNum+nfile}.spm.stats.factorial_design.cov.c = [1
            %                                                    2
            %                                                    3
            %                                                    4];
            % matlabbatch{in.procsNum+nfile}.spm.stats.factorial_design.cov.cname = 'aCC';
            % matlabbatch{in.procsNum+nfile}.spm.stats.factorial_design.cov.iCFI = 1;
            % matlabbatch{in.procsNum+nfile}.spm.stats.factorial_design.cov.iCC = 1;
            matlabbatch{procs_Num}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{procs_Num}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{procs_Num}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{procs_Num}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{procs_Num}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{procs_Num}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{procs_Num}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{procs_Num}.spm.stats.factorial_design.globalm.glonorm = 1;
            matlabbatch = Estimate(outpath,matlabbatch,'procsNum',procs_Num+1);
            
            matlabbatch{procs_Num+2}.spm.stats.con.spmmat = {fullfile(outpath,'SPM.mat')};
            matlabbatch{procs_Num+2}.spm.stats.con.consess{1}.tcon.name = 'GroupLevel';
            matlabbatch{procs_Num+2}.spm.stats.con.consess{1}.tcon.weights = 1;
            matlabbatch{procs_Num+2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            procs_Num = procs_Num+3;
        end
    case 'flexANOVA' 
        procs_Num = in.procsNum;
        foldername = 'FANOVA_';
        for i = 1:length(in.conName)
            foldername = append(foldername,[in.conName{i}]);
            if i ~= length(in.conName), foldername = append(foldername,'_'); end
        end
        outpath = char(fullfile(outputpath,foldername,in.voi));
        if ~exist(outpath,"dir"), mkdir(outpath); end

        % second level design
        matlabbatch{procs_Num}.spm.stats.factorial_design.dir = {outpath};
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(1).name = 'sub';
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(2).name = 'homophone';
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(3).name = 'frequency';
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(3).dept = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
        for nsub = 1:Nsubject
            scans = cell(Ncon,1);
            for ncon = 1:Ncon
                scans{ncon} = fullfile(filepath,char(subject(nsub)),in.analysis,in.voi,[in.TorF,'_',char(conName(ncon)),'.nii']);
            end
            subconds = ones(Ncon,1)*nsub;
            conconds = in.conds;
            conds = cat(2,subconds,conconds);
            matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(nsub).scans = scans;
            matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(nsub).conds = conds;
        end   
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
        matlabbatch{procs_Num}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [2
                                                                                          3];
        matlabbatch{procs_Num}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{procs_Num}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{procs_Num}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{procs_Num}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{procs_Num}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{procs_Num}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{procs_Num}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{procs_Num}.spm.stats.factorial_design.globalm.glonorm = 1;

        % estimate second level model
        matlabbatch = Estimate(outpath,matlabbatch,'procsNum',procs_Num+1);
        writecell(in.conName,fullfile(outpath,'contrast_name.txt'));

        % define contrast
        matlabbatch{procs_Num+2}.spm.stats.con.spmmat = {fullfile(outpath,'SPM.mat')};
        if in.allcon
            for ncon = 1:length(conName)
                cw = zeros(1,Ncon);
                cw(conName==(conName(ncon))) = 1;
                weight = [ones(1,Nsubject)/Nsubject,cw];
                matlabbatch{procs_Num+2}.spm.stats.con.consess{ncon}.tcon.name = char(conName(ncon));
                matlabbatch{procs_Num+2}.spm.stats.con.consess{ncon}.tcon.weights = weight;
                matlabbatch{procs_Num+2}.spm.stats.con.consess{ncon}.tcon.sessrep = 'none';
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
  
                matlabbatch{procs_Num+2}.spm.stats.con.consess{ncon+conN}.tcon.name = char(in.conName_2nd(ncon));
                matlabbatch{procs_Num+2}.spm.stats.con.consess{ncon+conN}.tcon.weights = weight;
                matlabbatch{procs_Num+2}.spm.stats.con.consess{ncon+conN}.tcon.sessrep = 'none';
            end
        end

end
fileid = 0;
if exist(fullfile(outpath,'batch_secondLevel_0.mat'),'file')
    file = string({dir(outpath).name});
    file = file(contains(file,'batch_secondLevel'));
    file = file(end);
    fileid = split(file,'_');
    fileid = fileid(end);
    fileid = split(fileid,'.');
    fileid = fileid(1);
    fileid = str2double(fileid)+1;
end
save(fullfile(outpath,['batch_secondLevel_',num2str(fileid),'.mat']),"matlabbatch");

end