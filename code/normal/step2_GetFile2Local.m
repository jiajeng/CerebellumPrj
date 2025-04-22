% get prep file(.nii) and behave to local spcace
% and Organize behave data to .xlsx and .mat(behD.mat)
clear all; close all;

localfolder = 'C:\Users\user\Desktop\jeng\CerebellumPrj\data_task'; % /media/tcnl/SSD1/jeng/data
ftpServer.ip = 'ftp://xxxxxxxxxxxxx/';
ftpServer.account = 'xxxxxxxxxxxxx';
ftpServer.password = 'xxxxxxxxxxxxx';
ftpServer.folder = 'LabData/cerebellum_NYCU/process/file_prep';
ftpServer.behfolder = 'LabData/cerebellum_NYCU/rawdata';  
ftpServer.behoutfolder = '/LabData/cerebellum_NYCU/process/ALLBEH';
session = {'WORD'};


% get subject folder name from Nas
ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
cd(ftpobj,ftpServer.behfolder);
subject = string({dir(ftpobj).name}');
close(ftpobj);
subject = split(subject(contains(subject,'SUB')),'_');
subject = unique(subject(:,2));


if exist('logbev.mat','file')
    load('logbev.mat'); 
else 
    logbev = cell2table(cell(length(subject),1)); 
    logbev.Properties.VariableNames = {'behave_folder'};
    logbev.Properties.RowNames = subject;
end
for nsub = 1:length(subject)
    localfolder = fullfile(localfolder);
    %% get behave data
    if ~exist(fullfile(fullfile(localfolder,subject(nsub),'BEHAV')),'dir')
        try
        logbev = getbehFile(ftpServer,localfolder,char(subject(nsub)),logbev);
        catch ME, continue; 
        end
    end
    folderpath = fullfile(localfolder,subject(nsub),'BEHAV');
    folder = string({dir(folderpath).name});
    if ~sum(contains(folder,'.mat'))
        if sum(contains(folder,'word')) || sum(contains(folder,'Word')) ||sum(contains(folder,'WORD'))
            folder = folder(contains(folder,'word') | contains(folder,'Word') | contains(folder,'WORD'));
            n = string({dir(fullfile(folderpath,folder)).name})';
            n(n=="."|n=="..") = [];
            f = repmat(fullfile(folderpath,folder),1,2)';
            for i = 1:size(n,1)
                s = fullfile(f(i,:),n(i,:));
                movefile(s,folderpath);
            end
        end
        rmdir(fullfile(folderpath,folder))
    end
    %% get functional file
    subj = char(subject(nsub));
    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    cd(ftpobj,ftpServer.folder);
    for nsess = 1:length(session)
        targetfile = fullfile(subj,'niifile',session{nsess},[subj,'_4D.nii']);
        localfolder = fullfile(localfolder);
        if ~exist(fullfile(localfolder,targetfile),"file")
            targetfile(targetfile=='\') = '/';
            mget(ftpobj,targetfile,localfolder)
        end
    end
    %% get structrure file
    targetfolder = [subj,'/niifile/T1'];
    cd(ftpobj,targetfolder);
    targetfile = {dir(ftpobj).name}';
    targetfile = char(targetfile(cellfun(@(x) x(end-3:end) == ".nii",targetfile)));
    for i = 1:3
        cd(ftpobj,'..');
    end
    if ~exist(fullfile(localfolder,targetfolder,targetfile),"file")
        mget(ftpobj,[targetfolder,'/',targetfile],localfolder)
    end
    close(ftpobj)
end
save('logbev',"logbev");

%% organize behave data
Ofolder = fullfile(localfolder,'all','BEHAV');
cereb_stat_beh(localfolder,Ofolder,ftpServer)




%% define function
function log = getbehFile(ftpServer,localfolder,sub,log)
    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    cd(ftpobj,ftpServer.behfolder);
    SUBname = string({dir(ftpobj).name}');
    SUBname = SUBname(contains(SUBname,sub));
    cd(ftpobj,'../../..');
    if length(unique(SUBname)) > 1
        flag = zeros(1,length(SUBname));
        for nsub = 1:length(SUBname)
            cd(ftpobj,[ftpServer.behfolder,'/',char(SUBname(nsub))]);
            folder =  string({dir(ftpobj).name}');
            if sum(contains(folder,'BEHAV'))
                flag(nsub) = 1;
            end
            cd(ftpobj,'../../../..');
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
    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    targetfolder = fullfile(ftpServer.behfolder,char(SUBname));
    targetfolder(targetfolder=='\') = '/';
    cd(ftpobj,targetfolder);
    mget(ftpobj,'BEHAV',fullfile(localfolder,sub))
    close(ftpobj);
    log{sub,'behave_folder'} = {char(SUBname)};
end


function cereb_stat_beh(niifilepath,output_folder,ftpServer)
    subject = string({dir(niifilepath).name});
    subject = subject(contains(subject,"SUB"));
    feature = {'accuracy','Resp_time'};
    res = struct();
    for nsub = 1:length(subject)
        behpath = fullfile(niifilepath,char(subject(nsub)),'BEHAV');
        if ~exist(fullfile(behpath,'behD.mat'),'file')
            % process behave data
            behfilename = string({dir(behpath).name});
            behfilename = behfilename(contains(behfilename,'.mat') & ~contains(behfilename,'behD.mat'));
            load(fullfile(behpath,behfilename)); % Result
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
            outputfolder = fullfile(niifilepath,char(subject(nsub)),'BEHAV');
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
                cd(ftpobj,ftpServer.folder);
                mput(ftpobj,outputfolder);
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
        ftpfolder = ftpServer.behoutfolder ;
        try
            mkdir(ftpobj,ftpfolder);
        catch ME

        end
        cd(ftpobj,ftpfolder);
        mput(ftpobj,output_folder)
        close(ftpobj)
    end
end
