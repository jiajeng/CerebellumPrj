clear all; close all;


ftpServer.ip = 'ftp://xxxxxxxxxxxxx/';
ftpServer.account = 'xxxxxxxxxxxxx';
ftpServer.password = 'xxxxxxxxxxxxx';
ftpServer.connfolder = 'LabData/cerebellum_NYCU/process/connfolder';  
ftpServer.prepfolder = '/LabData/cerebellum_NYCU/process/connprepdata';

filepath = 'C:\Users\user\Desktop\jeng\CerebellumPrj\data_task';% '/media/tcnl/SSD1/jeng/cerebellum/data'

if ~exist(filepath,'dir')
    mkdir(filepath);
end
func_filename = '*_4D.nii';
struct_filename = 's*-01.nii';

%--------------------- get subject name from Nas--------------------
subject = string({dir(filepath).name});
subject = subject(contains(subject,'SUB'));
% ------------------------------------------------------------------

% ----------------------------set condition------------------------
session = {'WORD'}; %'WORD' or 'REST'
if session{1} == "WORD"
    condition = {'HDHF','HDLF','LDHF','LDLF','ERROR'};
elseif session{1} == "REST"
    condition = {'REST'};
end
conn_proj_name = ['conn_cerebellum_',session{1}];
conn_proj_path = fullfile(pwd,conn_proj_name);
% -----------------------------------------------------------------


% ----------------------define parameter---------------------------
TR = 1.5;
sliceorder = load('sliceorder.mat'); % 
sliceorder = sliceorder.sliceorder;
onsetname = 'SoundOnset';
dur = 0;

% add user define roi
roifilepath = 'C:\Users\user\Desktop\jeng\CerebellumPrj\roi';
roifilename = string({dir(roifilepath).name});
roifilename = roifilename(contains(roifilename,'.nii'));

% add default roi
roifile = {fullfile(fileparts(which('conn')),'rois','atlas.nii'),...
       fullfile(fileparts(which('conn')),'rois','networks.nii'),...
       };

for nfile = 1:length(roifilename)
    roifile = cat(2,roifile,{[fullfile(roifilepath,char(roifilename(nfile)))]});
end
confound = {'White Matter','CSF','realignment','scrubbing'};
% -----------------------------------------------------------------

nsubject = length(subject);
nsession = length(session);
ncondition = length(condition);


%% CONN experiment
batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
% if ~isempty(dir(batch.filename))
%     Ransw=questdlg([conn_proj_name,' project already exists, Overwrite?'],'warning','Yes','No','No');
%     if ~strcmp(Ransw,'Yes'), return; end; 
% end

%% CONN analysis
batch.Analysis.overwirte = 0;
batch.Analysis.name = 'gPPI';
batch.Analysis.measure = 3;
batch.Analysis.modulation = 1;


conn_batch(batch);

% % put conn project to nas
ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
localfolder = fullfile(pwd,conn_proj_name);
try mkdir(ftpobj,ftpServer.connfolder); catch ME, keyboard; end
cd(ftpobj,ftpServer.connfolder)
mput(ftpobj,localfolder);
close(ftpobj);

% % put data to nas
% ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
% localfolder = fullfile(filepath);
% try mkdir(ftpobj,ftpServer.prepfolder); catch ME, keyboard; end
% cd(ftpobj,ftpServer.prepfolder)
% mput(ftpobj,localfolder);
% close(ftpobj);

%% function define