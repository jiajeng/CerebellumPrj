function GetNatROI(MNIcor,T1filepath,ROIpath,varargin)
    % input : MNIcor --> array(nx3), MNI Space corordinate, row is
    %                                all transform corordinate 
    %         T1filepath --> string, store T1 file directory(can use **)
    %         ROIpath --> string, directory save ROI file
    %         ROI group Name --> string, in thess MNIcor add ROI
    %                                    Name ,"L_cere"
    if nargin == 4
        ROINam = [varargin{1},'_'];
    elseif nargin > 4
        error('not allow enter more than two varargin input')
    else
        ROINam = '';
    end

    % get Affine variable file and load it 
    AfinMatf = 'mc_*_snc.mat';
    AfinMatf = dir(fullfile(T1filepath,AfinMatf));
    load(fullfile(AfinMatf.folder,AfinMatf.name));
    subName = split(AfinMatf.folder,filesep);
    subName = subName(contains(subName,'SUB')|contains(subName,'sub')|contains(subName,'Sub'));

    % get Affine matrix
    Q = VG.mat*inv(Affine)/VF.mat; % Native Space to MNI Space
    
    % MNI Space to Native Space corordinate
    Natcor = zeros(size(MNIcor));
    for i = 1:size(MNIcor,1) 
        tmp = [MNIcor(i,:), 1]*inv(Q');
        Natcor(i,:) = tmp(1:3);
    end

    % create Native Space ROI to ROI folder
    for i = 1:size(Natcor,1)
        d = [ROINam,num2str(i)]; % 'SUB002_L_cere_1'
        roi = maroi_box(struct('centre',Natcor(i,:)','widths',[6;6;6]));
        roi = descrip(roi,d);
        roi = label(roi,d);
        if ~exist(ROIpath,'dir'), mkdir(ROIpath); end
        save(fullfile(ROIpath,[d,'.mat']),'roi');
        mars_rois2img(fullfile(ROIpath,[d,'.mat']), fullfile(ROIpath,[d,'.nii']), '', 'c')
    end
end
