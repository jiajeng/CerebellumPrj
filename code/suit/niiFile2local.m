% get prep file(.nii)
function niiFile2local(sess,round,ftpServer,folder_nest,localfolder)
    % input : sess --> "cell", all session name `
    %         round --> "cell", all round name (task name or T1 is in here) 
    %         ftpServer --> "struct", store ftp information, include account, ip, ... 
    %         folder_nest --> "cell", the mri data folder nest under the subject, 
    %                         e.x. mri data under subject folder /sess-01/mri/REST/
    %                              folder_nest = {'sess','mri','round'}  
    %         localfolder --> "string", local path to put nii file


    % get subject folder name from Nas
    ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
    cd(ftpobj,ftpServer.outfolder);
    subject = string({dir(ftpobj).name}');
    subject = subject(contains(subject,'Sub')|contains(subject,'SUB')|contains(subject,'sub'));
    close(ftpobj);
    for nsub = 1:length(subject)
        %% get functional file
        subj = char(subject(nsub));
        for nsess = 1:length(sess)
            for nround = 1:length(round)
                try
                    Fnest = folder_nest;
                    Fnest(cellfun(@(x) x=="sess",Fnest)) = sess(nsess);
                    Fnest(cellfun(@(x) x=="round",Fnest)) = round(nround);
                    Fnest = char(strjoin(string(Fnest),filesep));
                    if round{nround} == "T1"
                        %% get structrure file
                        ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                        cd(ftpobj,ftpServer.outfolder);
                        targetfolder = fullfile(subj,Fnest);
                        targetfolder(targetfolder=='\') = '/';
                        cd(ftpobj,targetfolder);
                        targetfile = {dir(ftpobj).name}';
                        targetfile = fullfile(subj,Fnest,char(targetfile(cellfun(@(x) x(end-3:end) == ".nii" & x(1) == "s",targetfile))));
                        close(ftpobj);
                        ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                        cd(ftpobj,ftpServer.outfolder);
                        if ~exist(fullfile(localfolder,targetfolder,targetfile),"file")
                            targetfile(targetfile=='\') = '/';
                            mget(ftpobj,targetfile,localfolder)
                        end
                        close(ftpobj)
                    else
                        ftpobj = ftp(ftpServer.ip,ftpServer.account,ftpServer.password);
                        cd(ftpobj,ftpServer.outfolder);
                        targetfile = fullfile(subj,Fnest,[subj,'_4D.nii']);
                        if ~exist(fullfile(localfolder,targetfile),"file")
                            targetfile(targetfile=='\') = '/';
                            mget(ftpobj,targetfile,localfolder)
                        end
                        close(ftpobj)
                    end
                catch ME
                    try rmdir(fullfile(localfolder,subj),'s');catch,end
                    continue;
                end
            end
        end   
    end
end