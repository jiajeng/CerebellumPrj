% load data from nas(/rawdata) and store it to nas(/process/file_prep)

clear all; close all;
wkdir = pwd;

%% nas server parameter (change)
ftpServer.ip = 'ftp://xxxxxxxxxxxxx/';
ftpServer.account = 'xxxxxxxxxxxxxx';
ftpServer.password = 'xxxxxxxxxxxxx';
ftpServer.folder = 'LabData/cerebellum_NYCU/rawdata';
ftpServer.outfolder = 'LabData/cerebellum_NYCU/re_process/file_prep';

% find subject name 
ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
cd(ftpobj,ftpServer.folder);
SUBname = string({dir(ftpobj).name}');
SUBname = split(SUBname(contains(SUBname,'SUB')),'_');
close(ftpobj);
sub = unique(SUBname(:,2));
sess = {'REST','WORD'};

if ~exist('log.mat','file')
    log = cell2table(cell(length(unique(sub)),3)); 
    log.Properties.VariableNames = {'RESTfolder','T1folder','WORDfolder'};
    log.Properties.RowNames = unique(sub);
else
    load('log.mat');
end

%%%%%%%%%%%%%%%%%%
sub = sub(1:2);
%%%%%%%%%%%%%%%%%%
try
for nsub = 1:length(sub)
    for nsess = 1:length(sess)
        Ifolder = fullfile('');
        Ofolder = fullfile('niifile');
        log = func2nii(char(sub(nsub)),char(sess(nsess)),Ifolder,Ofolder,log,'ftpServer',ftpServer);
    end
    Ifolder = fullfile('');
    Ofolder = fullfile('niifile');
    log = anat2nii(char(sub(nsub)),'T1',Ifolder,Ofolder,log,'ftpServer',ftpServer);
end
save('niifldr_log.mat',"log");
catch ME
    rethrow(ME)
    % save('log.mat',"log");
end


%% define function
function log = anat2nii(sub,sess,input_folder,output_folder,log,varargin)
    in = finputcheck(varargin, ...
    {'sub'              'string'    []              sub; ...
     'sess'             'string'    []              sess; ...
     'input_folder'     'string'    []              input_folder; ...
     'output_folder'    'string'    []              output_folder; ...
     'ftpServer'        'struct'    []              struct();
     'bufferfolder'     'string'    []              'buffer';
     });
    wd = pwd;
    ftpServer = in.ftpServer;
    if isempty(fieldnames(ftpServer))
        in.bufferfolder = in.subpath;
    end
  
    % output file 
    outputpath = fullfile(in.bufferfolder,ftpServer.outfolder,in.sub,in.output_folder,in.sess);
    if ~exist(outputpath,'dir')
        mkdir(outputpath);
    end
    rename = in.sub;
    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    cd(ftpobj,ftpServer.folder);
    SUBname = string({dir(ftpobj).name}');
    SUBname = SUBname(contains(SUBname,in.sub));
    cd(ftpobj,'/..');
    if length(unique(SUBname)) > 1
        flag = zeros(1,length(SUBname));
        for nsub = 1:length(SUBname)
            cd(ftpobj,[ftpServer.folder,'/',char(SUBname(nsub))]);
            folder =  string({dir(ftpobj).name}');
            if any(folder == string(in.sess))
                flag(nsub) = 1;
            end
            cd(ftpobj,'/..');
        end
        if sum(flag) > 1
            tmp = split(SUBname,'_');
            tmp = str2double(tmp(:,1));
            SUBname = SUBname(tmp == max(tmp));
        else    
            SUBname = SUBname(logical(flag));
        end
    end
    close(ftpobj);
    in.sub = char(SUBname);
    datapath = fullfile(in.bufferfolder,ftpServer.folder,in.sub,in.sess,in.input_folder);
    if ~exist(datapath,"dir")
        datapath = ftpgetfile(ftpServer,in,in.input_folder);
    end
    dataname = {dir(datapath).name}';
    dataname = string(dataname(3:end));
    dataname = dataname(contains(dataname,'.IMA'));
    dataname = convertStringsToChars(dataname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    log{rename,[in.sess,'folder']} = {in.sub};

    matlabbatch{1}.spm.util.import.dicom.data = fullfile(datapath,dataname);
    %%
    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
    matlabbatch{1}.spm.util.import.dicom.outdir = {outputpath};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch);
    cd(outputpath)
    gzip('*.nii');
    cd(wd);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in.sub = rename;
    % save data to nas
    if ~isempty(fieldnames(ftpServer))
        ftpputfile(ftpServer,in)
        if exist(in.bufferfolder,"dir")
            rmdir(in.bufferfolder,'s')
        end
    end
end

function log = func2nii(sub,sess,input_folder,output_folder,log,varargin)
    in = finputcheck(varargin, ...
    {'sub'              'string'    []              sub; ...
     'sess'             'string'    []              sess; ...
     'input_folder'     'string'    []              input_folder; ...
     'output_folder'    'string'    []              output_folder; ...
     'ftpServer'        'struct'    []              struct();
     'bufferfolder'     'string'    []              'buffer';
     });

    ftpServer = in.ftpServer;
    if isempty(fieldnames(ftpServer))
        in.bufferfolder = in.subpath;
    end
  
    % output file 
    outputpath = fullfile(in.bufferfolder,ftpServer.outfolder,in.sub,in.output_folder,in.sess);
    if ~exist(outputpath,'dir')
        mkdir(outputpath);
    end
    rename = in.sub;
    
    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    cd(ftpobj,ftpServer.folder);
    SUBname = string({dir(ftpobj).name}');
    SUBname = SUBname(contains(SUBname,in.sub));
    cd(ftpobj,'/..');
    if length(unique(SUBname)) > 1
        flag = zeros(1,length(SUBname));
        for nsub = 1:length(SUBname)
            cd(ftpobj,[ftpServer.folder,'/',char(SUBname(nsub))]);
            folder =  string({dir(ftpobj).name}');
            if any(folder == string(in.sess))
                flag(nsub) = 1;
            end
            cd(ftpobj,'/..');
        end
        if sum(flag) > 1
            tmp = split(SUBname,'_');
            tmp = str2double(tmp(:,1));
            SUBname = SUBname(tmp == max(tmp));
        else    
            SUBname = SUBname(logical(flag));
        end
    end
    close(ftpobj);
    in.sub = char(SUBname);
    datapath = fullfile(in.bufferfolder,ftpServer.folder,in.sub,in.sess,in.input_folder);
    if ~exist(datapath,"dir")
        datapath = ftpgetfile(ftpServer,in,in.input_folder);
    end
    dataname = {dir(datapath).name}';
    dataname = string(dataname(3:end));
    dataname = dataname(contains(dataname,'.IMA'));
    dataname = convertStringsToChars(dataname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    log{rename,[in.sess,'folder']} = {in.sub};
    
    matlabbatch{1}.spm.util.import.dicom.data = fullfile(datapath,dataname);
    %%
    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
    matlabbatch{1}.spm.util.import.dicom.outdir = {outputpath};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
    matlabbatch{2}.spm.util.cat.vols(1) = cfg_dep('DICOM Import: Converted Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{2}.spm.util.cat.name = [rename,'_4D.nii'];
    matlabbatch{2}.spm.util.cat.dtype = 4;
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in.sub = rename;
    % save data to nas
    if ~isempty(fieldnames(ftpServer))
        ftpputfile(ftpServer,in)
        if exist(in.bufferfolder,"dir")
            rmdir(in.bufferfolder,'s')
        end
    end
end


function path = ftpgetfile(ftpServer,in,folder)
    if ~isempty(fieldnames(ftpServer))
        ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
        targetfolder = fullfile(ftpServer.folder,in.sub,in.sess,folder);
        targetfolder(targetfolder=='\') = '/';
        mget(ftpobj,targetfolder,in.bufferfolder);
        close(ftpobj)
        path = char(fullfile(in.bufferfolder,ftpServer.folder,in.sub,in.sess,folder));
    else
        path = char(fullfile(in.bufferfolder,in.sub,in.sess,folder));
    end
end


function ftpputfile(ftpServer,in)
    ftpObj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    folder = split(ftpServer.outfolder,'/');
    folder = fullfile(folder{1:3});
    folder(folder=='\') = '/';
    localfolder = fullfile(in.bufferfolder,ftpServer.outfolder);
    try
        cd(ftpObj,folder);
    catch
        mkdir(ftpObj,folder);
        cd(ftpObj,folder);
    end
    mput(ftpObj,localfolder)
    close(ftpObj);
end
