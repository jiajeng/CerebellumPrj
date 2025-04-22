function checkmrifile(T1file,funcfile,oname,opath)
    T1file = dir(T1file);
    funcfile = dir(funcfile);
    T1file = fullfile(T1file.folder,T1file.name);
    funcfile = [fullfile(funcfile.folder,funcfile.name),',1'];
    spm_check_registration(T1file,funcfile);
    if ~exist(opath,'dir'), mkdir(opath); end
    saveas(gcf,fullfile(opath,oname));
    close(gcf);
end
