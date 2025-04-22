clear all; close all;


ftpServer.ip = 'ftp://xxxxxxxxxxx/';
ftpServer.account = 'xxxxxxxxxxx';
ftpServer.password = 'xxxxxxxxxxxx';
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


% ------------get all subject func and struct .nii filepath--------------
% Selects functional / anatomical volumes
% get all func and anat file filepath --> C:/filepath/filename.nii
FUNCTIONAL_FILE = cell(nsubject,nsession); % nsub * nsess(REST and task)
STRUCTURAL_FILE = cell(nsubject,1); % nsub * 1(T1)

for nsub = 1:size(FUNCTIONAL_FILE,1)
    STRUCTURAL_file=dir(fullfile(filepath,char(subject(nsub)),'**',struct_filename));
    for nsess = 1:length(session)
        FUNCTIONAL_file=dir(fullfile(filepath,char(subject(nsub)),'niifile',session{nsess},func_filename));
        FUNCTIONAL_FILE(nsub,nsess) = {fullfile(char(FUNCTIONAL_file(nsess).folder),char(FUNCTIONAL_file(nsess).name))};
    end
    STRUCTURAL_FILE(nsub,1) = {fullfile(char(STRUCTURAL_file(1).folder),char(STRUCTURAL_file(1).name))};
end
% ------------------------------------------------------------------------


%% CONN New experiment
batch.filename=fullfile(conn_proj_path,[conn_proj_name,'.mat']);
% if ~isempty(dir(batch.filename))
%     Ransw=questdlg([conn_proj_name,' project already exists, Overwrite?'],'warning','Yes','No','No');
%     if ~strcmp(Ransw,'Yes'), return; end; 
% end
%% CONN Setup
batch.Setup.nsubjects = nsubject;
% set structural, functional file and condition variable
batch.Setup.conditions.names= condition;
batch.Setup.sessions.names = session;
for ncond=1:length(condition)
    for nsub=1:nsubject
        batch.Setup.structurals{nsub}=STRUCTURAL_FILE{nsub};
        for nses=1:nsession  
            batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; 
            batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;
            batch.Setup.functionals{nsub}{nses}=FUNCTIONAL_FILE{nsub,nses};
        end
    end
end

% set task condition
if session{1} == "WORD"
    for nsub = 1:nsubject
        for nses = 1
            behpath = fullfile(filepath,char(subject(nsub)),'BEHAV');
            load(fullfile(behpath,'behD.mat'))
            for ncond = 1:length(condition)
                onsettmp = behD.(condition{ncond}).(onsetname);
                batch.Setup.conditions.onsets{ncond}{nsub}{nses} = cell2mat(onsettmp);
                batch.Setup.conditions.durations{ncond}{nsub}{nses} = ones(1,length(onsettmp))*dur;
            end
        end
    end
end

% define roi
for nroi = 1:length(roifile)
    for nsub = 1:nsubject
        for nses = 1:nsession
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

batch.Setup.RT=TR;
batch.Setup.conditions.names=condition;       
batch.Setup.isnew=0; 
batch.Setup.voxelresolution = 3;% same as functionals
batch.Setup.outputfiles = [0,1,0,0,0,0]; % creates:[confound beta-map,
%                                                   confound-correlated timeseries,
%                                                   seed-to-voxel r-maps,
%                                                   seed-to-voxel p-maps,
%                                                   seed-to-voxel FDR-p-maps,
%                                                   ROI-extraction REX files]



%% CONN preprocessing
batch.Setup.analysisunits = 1;% 1:PSC uinits   2:raw units
batch.Setup.preprocessing.sliceorder = sliceorder;
batch.Setup.preprocessing.steps = 'default_mni';
% batch.Setup.preprocessing.steps = {'functional_slicetime','functional_art','functional_segment&normalize_direct','functional_label_as_mnispace','structural_center','structural_segment&normalize','functional_smooth','functional_label_as_smoothed'};
batch.Setup.preprocessing.voxelsize_anat = 1;
batch.Setup.preprocessing.voxelsize_func = 2;
batch.Setup.preprocessing.fwhm = 8;
batch.Setup.overwrite = 1;% not overwrite data can run faster 
batch.Setup.done=1;

%% CONN Denoising
batch.Denoising.filter=[0.008, 0.09];          % frequency filter (band-pass values, in Hz)
batch.Denoising.confounds.names = confound;
for i = 1:length(confound)
    c = confound{i};
    if c == "White Matter" || c == "CSF"
        batch.Denoising.confounds.dimensions{i} = 5;
    elseif c == "realignment"
        batch.Denoising.confounds.dimensions{i} = 12;
    elseif c == "scrubbing"
        batch.Denoising.confounds.dimensions{i} = 2;
    end
end
batch.Denoising.overwrite=1;
batch.Denoising.done=1;

if session{1} == "REST"
    batch.Analysis.name = 'SBC';
    batch.Analysis.done = true;
end

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