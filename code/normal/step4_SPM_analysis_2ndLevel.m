clear all; close all;
niifilepath = 'C:\Users\user\Desktop\jeng\CerebellumPrj\data_task\ALL\RESULT\firstLevel';
outputbasepath = 'C:\Users\user\Desktop\jeng\CerebellumPrj\data_task';
if ~exist(outputbasepath,'dir'), mkdir(outputbasepath); end

wdir = pwd;
ftpServer.ip = 'ftp://xxxxxxxxxxxxxx/';
ftpServer.account = 'xxxxxxxxxxxxxxx';
ftpServer.password = 'xxxxxxxxxxxxxxx';
ftpServer.folder = 'LabData/cerebellum_NYCU/re_process/file_prep';
ftpServer.behfolder = 'LabData/cerebellum_NYCU/rawdata';

ow = 1;
if exist("ow","var")
    save('overwrite.mat',"ow")
end

if ~exist(niifilepath,'dir')
    mkdir(niifilepath);
end

% get subject folder name from Nas
subject = string({dir(niifilepath).name});
subject = subject(contains(subject,'SUB'));

% define parameter
TR = 1.5;
grouplevelm = 'fANOVA';%onesampT fANOVA

% ---------------------------define contrast------------------------------
% contrast = {{'HDLF','LDLF','HDHF','LDHF'},{'HDHF','HDLF','LDHF','LDLF'},{'HDHF'},{'HDLF'},{'LDHF'},{'LDLF'}};
% weight = {{1,1,-1,-1},{1,1,-1,-1},{1},{1},{1},{1}};

contrast = {{'HDHF'},{'HDLF'},{'LDHF'},{'LDLF'}};
weight = {{1},{1},{1},{1}};


% define contrast name
% auto set contrast name or user define contrast name
if 0 
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
else 
    % user-define contrast name
    conName = {'HDHF','HDLF','LDHF','LDLF'};
end
%-------------------------------------------------------------------------

% --------------define contrast(for ANOVA)--------------------
contrast_2nd = {{'HDHF','HDLF','LDHF','LDLF'},{'HDHF','HDLF','LDHF','LDLF'}};
weight_2nd = {[1,-1,1,-1],[1,1,-1,-1]};
condName_2nd = {'HF-LF','HD-LD'};
% ------------------------------------------------------------

% define first level contrast data
% for nsub = 1:length(subject)
%     matlabbatch = [];
%     Ofolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)));
%     matlabbatch = define_contrast(Ofolder,matlabbatch,'procsNum',1,'contrast',contrast,'weight',weight,'conName',conName);
%     save(fullfile(Ofolder,'batch_defContrast.mat'),"matlabbatch");
%     spm('defaults','FMRI');
%     spm_jobman('run',matlabbatch);
% end

% define flexible ANOVA coditions matrix
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

matlabbatch = [];
% 2nd level
Ifolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel');
Ofolder = fullfile(outputbasepath,'ALL','RESULT','secondLevel');
switch grouplevelm
    case 'onesampT'
        % one sample T
        matlabbatch = second_level(Ifolder,Ofolder,matlabbatch,'method','onesampeT','procsNum',1,'conName',conName);
    case 'fANOVA'
        % flexible ANOVA
        matlabbatch = second_level(Ifolder,Ofolder,matlabbatch,'method','flexANOVA','procsNum',1,'conName',conName, ...
                            'conds',m,'contrast',contrast_2nd,'weight',weight_2nd,'conName_2nd',condName_2nd);
end

spm('defaults','FMRI');
spm_jobman('run',matlabbatch);

matlabbatch = [];

delete('overwrite.mat')

%% function define
