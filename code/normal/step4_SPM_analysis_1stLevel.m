% select data in server workstation 
% store data in local pc
clear all; close all;

niifilepath = '/media/tcnl/SSD1/jeng/cerebellum/data';
outputbasepath = '/home/tcnl/Desktop/C on Player (NoMachine)/Users/user/Desktop/jeng/CerebellumPrj/data';
if ~exist(outputbasepath,'dir'), mkdir(outputbasepath); end

wdir = pwd;
ftpServer.ip = 'ftp://xxxxxxxxxxxxx/';
ftpServer.account = 'xxxxxxxxxxxxx';
ftpServer.password = 'xxxxxxxxxxxxx';
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

matlabbatch = [];
% create 1st level
for nsub = 1:length(subject)
    Ifolder = fullfile(niifilepath,char(subject(nsub)),'niifile','WORD');
    Ofolder = fullfile(outputbasepath,'ALL','RESULT','firstLevel',char(subject(nsub)));
    matlabbatch = create_1stMod(Ifolder,Ofolder,matlabbatch,'procsNum',1);
    spm('defaults','FMRI');
    spm_jobman('run',matlabbatch);
end

delete('overwrite.mat')

%% function define

