function plotflapmat(filepath,fpat)
    file = {dir(fullfile(filepath,fpat)).name}';
    load(fullfile(filepath,'SPM.mat'));
    df = SPM.xX.erdf;
    fname = split(file,'.');
    if size(fname,2) == 1
        fname = fname(1);
    else
        fname = fname(:,1);
    end
    
    conname = split(fname,'_');

    if size(conname,2) == 1
        conname = conname(2);
    else
        conname = conname(:,2);
    end
    thres = 0.001;
    file = strcat(filepath,filesep,file);
    outpath = fullfile(filepath,'flapmap');
    if ~exist(outpath,'dir'), mkdir(outpath); end
    for nfile = 1:length(file)
        figure;
        Data = suit_map2surf(file{nfile},'space','SUIT','stats',@minORmax);
        Tthres = tinv(1-thres,df);
        % pData = tcdf(Data,df);
        suit_plotflatmap(Data,'threshold',Tthres);
        axis('off')
        stitle = sprintf('%s\nthreshold: %s',conname{nfile},num2str(thres));
        title(stitle,'FontSize',13)
        saveas(gcf,fullfile(outpath,[fname{nfile},'.jpeg']));
        close(gcf);
    end
end
