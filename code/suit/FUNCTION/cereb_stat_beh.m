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
