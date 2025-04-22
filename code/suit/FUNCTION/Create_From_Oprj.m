function Create_From_Oprj(OldprjPath,OldprjName,NprjPath,NprjName,Osubid)
    % OldprjName, 'char' --> old conn project name
    % OldprjPath, 'char' --> old conn project Path
    % NprjName, 'char' --> new conn project name
    % NprjPath, 'char' --> new conn project Path
    % subject, 'int' --> which subjects id in old project get to new project

    % Set CONN_x from old project
    CONN_x = load(fullfile(OldprjPath,[OldprjName,'.mat']));
    OldCONN_x = CONN_x.CONN_x;
    NewCONN_x = OldCONN_x;

    NewCONN_x.Setup.RT = OldCONN_x.Setup.RT(Osubid);
    NewCONN_x.Setup.nsubjects = length(Osubid);
    NewCONN_x.Setup.nsessions = OldCONN_x.Setup.nsessions(Osubid);
    NewCONN_x.Setup.functional = OldCONN_x.Setup.functional(Osubid);
    NewCONN_x.Setup.structural = OldCONN_x.Setup.structural(Osubid);
    NewCONN_x.Setup.spm = OldCONN_x.Setup.spm(Osubid);
    NewCONN_x.Setup.dicom = OldCONN_x.Setup.dicom(Osubid);
    NewCONN_x.Setup.nscans = OldCONN_x.Setup.nscans(Osubid);
    NewCONN_x.Setup.rois.files = OldCONN_x.Setup.rois.files(Osubid);
    NewCONN_x.Setup.conditions.values = OldCONN_x.Setup.conditions.values(Osubid);
    NewCONN_x.Setup.l1covariates.files = OldCONN_x.Setup.l1covariates.files(Osubid);
    NewCONN_x.Setup.l2covariates.values = OldCONN_x.Setup.l2covariates.values(Osubid);
    NewCONN_x.Setup.spatialresolutionvolume = OldCONN_x.Setup.spatialresolutionvolume(Osubid);
    f = fieldnames(NewCONN_x.folders);
    for nf = 1:length(f)
        NewCONN_x.folders.(f{nf}) = strrep(NewCONN_x.folders.(f{nf}),fullfile(OldprjPath,OldprjName),fullfile(NprjPath,NprjName));
        % create new conn data folders
        if ~exist(NewCONN_x.folders.(f{nf}),'dir'), mkdir(NewCONN_x.folders.(f{nf})); end
    end
    CONN_x = NewCONN_x;
    save(fullfile(NprjPath,[NprjName,'.mat']),"CONN_x");

    % copy file from old project
    for nsub = 1:length(Osubid)
        OsubidName = sprintf('Subject%03d',Osubid(nsub));
        file = dir(fullfile(OldprjPath,OldprjName,'**',['*',OsubidName,'*']));
        for nfile = 1:length(file)
            srcfile = fullfile(file(nfile).folder,file(nfile).name);
            detfpath = strrep(file(nfile).folder,fullfile(OldprjPath,OldprjName),fullfile(NprjPath,NprjName));
            detfile = fullfile(detfpath,file(nfile).name);
            copyfile(srcfile,detfile);
            Nfname = strrep(file(nfile).name,OsubidName,sprintf('Subject%03d',nsub));
            movefile(fullfile(detfile),fullfile(detfpath,Nfname));
        end
    end

end
