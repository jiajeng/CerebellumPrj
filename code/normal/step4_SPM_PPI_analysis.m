% select data in server workstation 
% store data in local pc
clear all; close all;

niifilepath = 'D:\jeng\CerebellumPrj\data_task';
ppi_con = 'ppi_allcon';
outputbasepath = ['C:\Users\user\Desktop\jeng\CerebellumPrj\result\',ppi_con];
if ~exist(outputbasepath,'dir'), mkdir(outputbasepath); end

wdir = pwd;

ow = 1;
if exist("ow","var")
    save('overwrite.mat',"ow")
end

if ~exist(outputbasepath,'dir')
    mkdir(outputbasepath);
end

% get subject folder name from Nas
subject = string({dir(niifilepath).name});
subject = subject(contains(subject,'SUB'));
% define parameter
TR = 1.5;

conditionName = {{'HDHF_inter'},{'HDLF_inter'},{'LDHF_inter'},{'LDLF_inter'}};
conditionw = {{1},{1},{1},{1}};
% roi definition
R_cere_cor = [24 -72 -48; ...
              10 -80 -24; ...
              32 -64 -26; ...
              22 -68 -26; ...
              42 -50 -30; ...
              36 -56 -24; ...
              12 -57 -30; ...
              17 -65 -35];

L_cere_cor = [-24 -72 -48; ...
              -10 -80 -24; ...
              -32 -64 -26; ...
              -22 -68 -26; ...
              -42 -50 -30; ...
              -36 -56 -24; ...
              -12 -57 -30; ...
              -17 -65 -35];

auditory = [-62 -6 -2;...
            -58 -10 2];

L_FL = [-42 6 32];

key = [];
value = [];
cor = {R_cere_cor,L_cere_cor,auditory,L_FL};
roiname = {'R_cere_cor','L_cere_cor','auditory','L_FL'};
% give every roi coordinary one name
for ncor = 1:length(cor)
    for nroi = 1:size(cor{ncor},1)
        rname = [roiname{ncor},'_',num2str(nroi)];
        value = cat(1,value,string(rname));
        key = cat(1,key,{cor{ncor}(nroi,:)});
    end
end
roi = dictionary(key,value);

% ------------------------------main------------------------------------
cor = [mat2cell(R_cere_cor,ones(1,size(R_cere_cor,1)),3); ...
       mat2cell(L_cere_cor,ones(1,size(L_cere_cor,1)),3); ...
       mat2cell(auditory,ones(1,size(auditory,1)),3); ...
       mat2cell(L_FL,ones(1,size(L_FL,1)),3)];

cor_ow = false;

% ----------------create 1st level(regressor--condition)
% for nsub = 1:length(subject)
%     % create 1st level(regressor--condition)
%     procsNum = 0;
%     matlabbatch = [];
%     Ifolder = fullfile(niifilepath,char(subject(nsub)),'niifile','WORD');
%     Ofolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)));
%     batchfolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)),'batch');
%     if ~exist(batchfolder,'dir'), mkdir(batchfolder); end
%     [matlabbatch,~] = create_1stMod(Ifolder,Ofolder,matlabbatch,'procsNum',procsNum+1);
%     save(fullfile(batchfolder,'batch1_1st_level.mat'),"matlabbatch");
%     spm('defaults','FMRI');
%     spm_jobman('run',matlabbatch);
% end


% ---------------create VOI.mat file (roi file to SPM)
% ------------------create PPI variable
% for nsub = 1:length(subject)
%     eroi = false;
%     Ofolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)));
%     batchfolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)),'batch');
%     if ~exist(batchfolder,'dir'), mkdir(batchfolder); end
%     if exist(fullfile(Ofolder,'PPI',rname),'dir')
%         exist_roiName = {dir(fullfile(Ofolder,'PPI')).name};
%         exist_roiName = exist_roiName(~contains(exist_roiName,'.'));
%         eroi = true;
%     end
%     matlabbatch = [];
%     procsNum = 0;
%     for ncor = 1:length(cor)
%         rname = char(roi(cor(ncor)));
%         foldfile = {dir(fullfile(Ofolder,'PPI',rname)).name};
%         if eroi
%             if ~cor_ow && any(string(exist_roiName)==string(rname))
%                 continue;
%             end
%         end
%         % ---------------create VOI.mat file (roi file to SPM)
%         [matlabbatch,procsNum] = create_VOI(Ofolder,matlabbatch,'roi_name',rname,'cen_cord',cor{ncor},'procsNum',procsNum+1);
%         % ----------------create PPI variable
%         if string(ppi_con) == "ppi_diff"
%             conName = {'LF_HF_diff'};
%             % HDHF HDLF LDHF LDLF
%             con = [1,1,-1;2,1,1;3,1,-1;4,1,1];
%             [matlabbatch,procsNum] = get_ppi_var(Ofolder,rname,matlabbatch,'conName',conName,'con',con,'procsNum',procsNum+1);
%             PPI_Ofolder = fullfile(Ofolder,'PPI',rname,char(conName));
%             if ~exist(PPI_Ofolder,'dir'), mkdir(PPI_Ofolder); end
%         elseif string(ppi_con) == "ppi_allcon"
%             [matlabbatch,procsNum] = get_ppi_var(Ofolder,rname,matlabbatch,'conName',{'HDHF','HDLF','LDHF','LDLF'},'procsNum',procsNum+1);
%             PPI_Ofolder = fullfile(Ofolder,'PPI',rname);
%             if ~exist(PPI_Ofolder,'dir'), mkdir(PPI_Ofolder); end
%         end
%     end
%     if ~isempty(matlabbatch)
%         save(fullfile(batchfolder,'batch2_PPIvar.mat'),"matlabbatch");
%         spm('defaults','FMRI');
%         spm_jobman('run',matlabbatch);
%     end
%     for ncor = 1:length(cor)
%         % move file to PPI folder
%         rname = char(roi(cor(ncor)));
%         if string(ppi_con) == "ppi_diff"
%             PPI_Ofolder = fullfile(Ofolder,'PPI',rname,char(conName));
%         elseif string(ppi_con) == "ppi_allcon"
%             PPI_Ofolder = fullfile(Ofolder,'PPI',rname);
%         end
% 
%         % remove VOI file
%         fileN = {dir(Ofolder).name};
%         fileN = fileN(contains(fileN,'VOI'));
%         fileN = fileN(contains(fileN,rname));
%         for i = 1:length(fileN)
%             movefile(fullfile(Ofolder,fileN{i}),PPI_Ofolder)
%         end
%         % remove PPI file
%         fileN = {dir(Ofolder).name};
%         fileN = fileN(contains(fileN,'PPI'));
%         fileN = fileN(contains(fileN,rname));
%         fileN = fileN(contains(fileN,'.mat'));
%         for i = 1:length(fileN)
%             movefile(fullfile(Ofolder,fileN{i}),PPI_Ofolder)
%         end
%     end
% end


% create PPI 1st level model
% for nsub = 1:length(subject)
%     procsNum = 0;
%     matlabbatch = [];
%     Ifolder = fullfile(niifilepath,char(subject(nsub)),'niifile','WORD');
%     Ofolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)));
%     batchfolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)),'batch');
%     if ~exist(batchfolder,'dir'), mkdir(batchfolder); end
%     for ncor = 1:length(cor)
%         rname = char(roi(cor(ncor)));
%         % create PPI 1st level model
%         [matlabbatch,procsNum] = create_1stMod(Ifolder,fullfile(Ofolder,'PPI',rname),matlabbatch,'procsNum',procsNum+1,'analysis','ppi','voiName',rname);
%     end
%     save(fullfile(batchfolder,'batch3_PPI_1st_level.mat'),"matlabbatch");
%     spm('defaults','FMRI');
%     spm_jobman('run',matlabbatch);
% end


% ---------------------------define contrast------------------------------
% contrast = {{'HDLF','LDLF','HDHF','LDHF'},{'HDHF','HDLF','LDHF','LDLF'},{'HDHF'},{'HDLF'},{'LDHF'},{'LDLF'}};
% weight = {{1,1,-1,-1},{1,1,-1,-1},{1},{1},{1},{1}};
contrast = {{'HDLF_inter','LDLF_inter','HDHF_inter','LDHF_inter'},...
            {'HDLF_inter','LDLF_inter','HDHF_inter','LDHF_inter'}
            };
contrast = cat(2,contrast,conditionName);
weight = {{1,-1,1,-1}, ...
          {-1,-1,1,1}};
weight = cat(2,weight,conditionw);
% user-define contrast name
conName =[{'HD-LD'},{'LF-HF'},conditionName{:}];
% contrast = {{'HDHF_inter'},{'HDLF_inter'},{'LDHF_inter'},{'LDLF_inter'}};
% weight = {{1},{1},{1},{1}};


% define contrast name
% auto set contrast name or user define contrast name
if isempty(conName)
    % auto-define contrast Name 
    conName = cell(size(contrast));
    for ncon = 1:length(contrast)
        wName = string(cellfun(@num2str, weight{ncon}, 'UniformOutput',false));
        cName = string(contrast{ncon});
        cn = reshape([wName;cName],1,2*length(cName));
        tmp = [];
        for i = 1:2:length(cn)
            tmp = cat(2,tmp,append(cn(i),cn(i+1)));
        end
        connm = tmp(1);
        for i = 2:length(tmp)
            connm = append(connm,'_',tmp(i));
        end 
        conName{ncon} = connm;
    end
end
%-------------------------------------------------------------------------

% -----------------------define first level contrast data
% for nsub = 1:length(subject)
%     matlabbatch = [];
%     for ncor = 1:length(cor)
%         rname = char(roi(cor(ncor)));
%         Ofolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)),'PPI',rname);
%         batchfolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)),'batch');
%         if ~exist(batchfolder,'dir'), mkdir(batchfolder); end
%         matlabbatch = define_contrast(Ofolder,matlabbatch,'procsNum',1,'contrast',contrast,'weight',weight,'conName',conName);
%         save(fullfile(batchfolder,'batch4_defContrast.mat'),"matlabbatch");
%         spm('defaults','FMRI');
%         spm_jobman('run',matlabbatch);
%     end
% end


% second level model
matlabbatch = [];
grouplevelm = 'fANOVA';

if string(grouplevelm) == "fANOVA"
% --------------define contrast(for ANOVA)--------------------
conName = {'HDHF_inter','HDLF_inter','LDHF_inter','LDLF_inter'};
contrast_2nd = {{'HDHF_inter','HDLF_inter','LDHF_inter','LDLF_inter'},{'HDHF_inter','HDLF_inter','LDHF_inter','LDLF_inter'}};
weight_2nd = {[-1,1,-1,1],[1,1,-1,-1]};
condName_2nd = {'LF-HF','HD-LD'};
% ------------------------------------------------------------

% -----------define flexible ANOVA coditions matrix---------------------
ANOVA_feature = {'HD','LD';
                 'HF','LF'};
feaCod = [1,2;
          1,2];

m = zeros(length(conName),size(ANOVA_feature,2));
for i = 1:length(conName)
    for j = 1:size(ANOVA_feature,1)
        c = cellfun(@(x) contains(conName{i},x),ANOVA_feature(j,:));
        m(i,j) = feaCod(j,c);
    end
end
end


% ---------------------------------2nd level----------------------------
% for ncor = 1:length(cor)
%     rname = char(roi(cor(ncor)));
%     Ifolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel');
%     Ofolder = fullfile(outputbasepath,'ALL','RESULT','secondLevel');
%     switch grouplevelm
%         case 'onesampT'
%             % one sample T
%             matlabbatch = second_level(Ifolder,Ofolder,matlabbatch,'method','onesampeT','procsNum',1,'conName',conName,'analysis','PPI','voi',rname);
%         case 'fANOVA'
%             % flexible ANOVA
%             matlabbatch = second_level(Ifolder,Ofolder,matlabbatch,'method','flexANOVA','procsNum',1,'conName',conName, ...
%                                 'conds',m,'contrast',contrast_2nd,'weight',weight_2nd,'conName_2nd',condName_2nd,'analysis','PPI','voi',rname);
%     end
%     spm('defaults','FMRI');
%     spm_jobman('run',matlabbatch);
% end

% save result figure
for ncor = 1:length(cor)
    matlabbatch = [];
    rname = char(roi(cor(ncor)));
    threstype = 'none';
    thres = 0.005;
    Ofolder = fullfile(outputbasepath,'ALL','RESULT','secondLevel','FANOVA_HDHF_inter_HDLF_inter_LDHF_inter_LDLF_inter',rname);
    [matlabbatch,procsNum] = secndLevlresult(Ofolder,matlabbatch,'procsNum',1,'threstype',threstype,'thresh',thres);
    spm('defaults','FMRI');
    spm_jobman('run',matlabbatch);
    % move figure file to ./result folder
    fileN = {dir(Ofolder).name};
    fileN = fileN(contains(fileN,".jpg"));
    for nf = 1:length(fileN)
        cond = extractBefore(extractAfter(fileN{nf},'spm_'),'_0');
        rfolder = fullfile(Ofolder,'..',['REPORT_',threstype,'_',num2str(thres)],cond);
        renames = [extractBefore(fileN{nf},cond),rname,extractAfter(fileN{nf},cond)];
        if ~exist(rfolder,'dir'),mkdir(rfolder); end
        movefile(fullfile(Ofolder,fileN{nf}),fullfile(rfolder,renames));
    end
end



delete('overwrite.mat')


% -------------------------------------------------------------------------

%% function define
function [matlabbatch,procsNum] = get_ppi_var(filepath,voiname,matlabbatch,varargin)
in = finputcheck(varargin, ...
       {'conName'   'cell'     []  [];
        'con'       'real'      []  [];
        'procsNum'  'real'      []  [];
       });
procsNum = in.procsNum;
load(fullfile(filepath,'SPM.mat'))
conName = {SPM.Sess.U.name};
tmp = cell(size(conName));
for i = 1:length(conName)
    if length(conName{i}) < 1
        tmp(i) = conName{i};
    else
        tmp{i} = conName{1,i}{1};
    end
end
conName = string(tmp);
if isempty(in.con)
    for ncon = 1:length(in.conName)
        if ncon~=1, procsNum = procsNum+1; end
        con_pos = find(string(conName)==string(in.conName{ncon}));
        matlabbatch{procsNum}.spm.stats.ppi.spmmat = {fullfile(filepath,'SPM.mat')};
        matlabbatch{procsNum}.spm.stats.ppi.type.ppi.voi = {fullfile(filepath,['VOI_',voiname,'_1.mat'])};
        matlabbatch{procsNum}.spm.stats.ppi.type.ppi.u = [con_pos, 1, 1];
        matlabbatch{procsNum}.spm.stats.ppi.name = [in.conName{ncon},'x',voiname];
        matlabbatch{procsNum}.spm.stats.ppi.disp = 0;
    end
else
    matlabbatch{procsNum}.spm.stats.ppi.spmmat = {fullfile(filepath,'SPM.mat')};
    matlabbatch{procsNum}.spm.stats.ppi.type.ppi.voi = {fullfile(filepath,['VOI_',voiname,'_1.mat'])};
    matlabbatch{procsNum}.spm.stats.ppi.type.ppi.u = in.con;
    matlabbatch{procsNum}.spm.stats.ppi.name = [in.conName{1},'x',voiname];
    matlabbatch{procsNum}.spm.stats.ppi.disp = 0;
end
end

function [matlabbatch,procsNum] = create_VOI(filepath,matlabbatch,varargin)
in = finputcheck(varargin, ...
       {'TR'        'real'      []  1.5;    
        'roi_name'  'string'    []  []; % roi name 
        'roi_type'  'string'    {'sphere','box'} 'box';
        'cen_cord'  'real'      []  []; % roi center coordinate 
        'dim'       'real'      []  [6 6 6]; % roi size for box
        'rad'       'real'      []  []; % roi size for sphere radius
        'procsNum'  'real'      []  1;
       });
procsNum = in.procsNum;
matlabbatch{procsNum}.spm.util.voi.spmmat = {fullfile(filepath,'SPM.mat')};
matlabbatch{procsNum}.spm.util.voi.adjust = 0;
matlabbatch{procsNum}.spm.util.voi.session = 1;
matlabbatch{procsNum}.spm.util.voi.name = in.roi_name;
switch in.roi_type
    case 'box'
        matlabbatch{procsNum}.spm.util.voi.roi{1}.box.centre = in.cen_cord;
        matlabbatch{procsNum}.spm.util.voi.roi{1}.box.dim = in.dim;
        matlabbatch{procsNum}.spm.util.voi.roi{1}.box.move.fixed = 1;
    case 'sphere'
        matlabbatch{procsNum}.spm.util.voi.roi{1}.sphere.centre = in.cen_cord;
        matlabbatch{procsNum}.spm.util.voi.roi{1}.sphere.radius = in.rad;
        matlabbatch{procsNum}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
end
matlabbatch{procsNum}.spm.util.voi.expression = 'i1';

end




