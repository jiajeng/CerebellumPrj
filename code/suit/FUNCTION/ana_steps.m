function ana_steps(filepath,filename,f)
    % segment anatomical data
    % isolate cerebullum and brainstem
    % cropping? to suit space 
    % normalize to suit templete

    % SUIT segmentation
    % =====================================================
    if f.ana.SegIso
        file = ls(fullfile(filepath,filename));
        Steps.checkfile(file,filepath,filename)
        suit_isolate_seg({fullfile(filepath,file)});
    end
    % =====================================================

    % SUIT normalize
    % ========================================================
    if f.ana.Norm
        fn = split(file,'.');
        fn = fn{1};
        mfile = ['c_',fn,'_pcereb.nii'];
        file = ['c_',fn,'.nii'];
        suit_normalize(fullfile(filepath,file),'mask',fullfile(filepath,mfile));
    end
    % ========================================================

    % get white matter and grey matter in whole brain
    % ========================================================
    if f.ana.wholemask
        mi = [1,2,7,8];
        fl = 0;
        for i = mi
            % get c1 c2 c7 c8 
            file = char({dir(fullfile(filepath,['c',num2str(i),'*.nii'])).name});
            mask = spm_read_vols(spm_vol(fullfile(filepath,file)));
            if ~fl
                MASK = zeros(size(mask));
                fl = 1;
            end
            MASK = MASK+mask;
        end
        V = spm_vol(fullfile(filepath,file));
        V.fname = fullfile(filepath,'wholebrain.nii');
        spm_write_vol(V,MASK);
    end
    % ========================================================
    
end
