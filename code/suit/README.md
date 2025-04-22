# Work steps
![image](https://github.com/user-attachments/assets/45875cc7-2b3f-4b29-b9d6-5242cf521a4d)  

- open conn first to load some enviroment variable from CONN  

- step 1. Set variables [[Setup.m](Setup.m)]  
  Set some needed variables **[now using UI to set variables(20241115)]**  
  condition read from *behD.mat*  
  ```matlab
  [o,f] = Setup.SetVariable;
  SetupUI(f,o);
  ```
- step 2. BEHAVE data folder in NAS (task condition) [[Setup.m]((Setup.m) -- GetBehFolder]
  Get behave data folder path and save in a .mat file -- ("*BEHFold_log.mat*")
  ```matlab
  % input :  
  % BEH --> "struct", fieldname is Behave name, e.x. Word
  %                   contains subpath, "string", subject folder path
  %                            folder, "string", behave folder path under sub folder
  %                            ext, "string", file extension
  %                     option idn, "cell"--"string", file contains certain string
  % sub --> "cell", all subject name, e.x.{'SUB01',SUB02',SUB03'}
  % sess --> "cell", all sess name, e.x.{'sess-01'}, if no sess {''}
  % ftpServer --> "struct", ip --> nas ip
  %                         account --> nas account
  %                         password --> nas password
  %                         infolder --> in which folder to get data
  %                         outfolder --> push nifti file to which folder
  % output : BEHFold_log.mat
  Setup.GetBehFolder(BEH,sub,sess,ftpServer)
  ```
- step 3. BEHAVE data filepath (physio and task condition) [[Setup.m]((Setup.m) -- GetBehFolder]
  if file in nas then read ("*BEHFold_log.mat*") info to get file
  store all local path that store behave file in a .mat file -- ("*BEHFile_log.mat*")
  ```matlab
  % input :  
  % BEH --> "struct", fieldname is Behave name, e.x. Word
  %                   contains subpath, "string", subject folder path
  %                            folder, "string", behave folder path under sub folder
  %                            ext, "string", file extension
  %                     option idn, "cell"--"string", file contains certain string
  % ftpServer --> "struct", ip --> nas ip
  %                         account --> nas account
  %                         password --> nas password
  %                         infolder --> in which folder to get data
  %                         outfolder --> push nifti file to which folder
  % subidntfr --> "string", indentefier subject, e.x sub001 identefier = 'sub'
  % subpath --> "string", subject folder path
  % output : BEHFile_log.mat
  Setup.GetBehfile(BEH,ftpServer,subidntfr,subpath)
  ```
  
- step 4. convert IMA to NII [[conver2nii.m](convert2nii.m)]
  get IMA file in nas or local and convert to NII file
  ```matlab
  % input : sub --> "cell", all subject name   
  %         sess --> "cell", all session name    
  %         round --> "cell", all round name (task name or T1 is in here)  
  %         folder_nest --> "cell", the mri data folder nest under the subject,  
  %              e.x. mri data under subject folder /sess-01/mri/REST/, folder_nest = {'sess','mri','round'}  
  % local input : subpath --> "string", the direction that store all subject folder 
  %               outsubpath --> "string", the direction that put subject nii file folder 
  % nas input : ftpServer --> "struct", store ftp information, include
  %                                      ip
  %                                      account
  %                                      password
  %                                      infolder
  %                                      outfolder
  % option input : DataType --> "string", rename the niifile set method in line 124(in function DICOM2NII), 
  %                                       if empty then do not rename
  % outputfile : .nii file in setting outputpath  
  %               MRinfo.xlsx --> store all process MRI information, include IamgeSize, slice order, filepath ...
  
  convert2nii(sub,sess,round, ...
          'folder_nest',FolderNest, ...
          'subpath',subpath, ...
          'outsubpath',subpath, ...
          'ftpServer',ftpServer, ...
          'DataType','cere')
  ```
> [!Note]
> MRinfo.xlsx store in this run data infomation  
> so run all data in one run is better, or need to manully fix MRinfo.xlsx file.

- step 5. Get NII file to local *if file in nas* [[niiFile2local.m](niiFile2local.m)]
  get NII file from nas to local  
  ```matlab
  % input : sess --> "cell", all session name `
  %         round --> "cell", all round name (task name or T1 is in here) 
  %         ftpServer --> "struct", store ftp information, include account, ip, ... 
  %         folder_nest --> "cell", the mri data folder nest under the subject, 
  %                         e.x. mri data under subject folder /sess-01/mri/REST/
  %                              folder_nest = {'sess','mri','round'}  
  %         localfolder --> "string", local path to put nii file
  % output : get need mri file to local path
  niiFile2local(sess,round,ftpServer,FolderNest,subpath)
  ```
- step 6. Get MRinfo
  read MRinfo.xlsx to MRinfo struct
  ```matlab
  for nround = 1:length(round)
    MRinfo.(round{nround}) = readtable(fullfile(subpath,'MRinfo.xlsx'),'ReadRowNames',true,'Sheet',round{nround});
  end
  ```

- step 7. Get BEHAVE data
  Get behave data (*from* "*BEHFold_log.mat*") folder info  
  Get individual subject behave data and save file in individual folder -- ("*behD.mat and behave.xlsx*")
  ```matlab
  % subpath --> "string", subject folder path
  % outFolder --> "string", save behave file path under subject folder
  % ftpServer --> "struct", back up result to nas, if not empty(ftpServer = struct())
  % Get behD.mat and behave.xlsx under subject folder 

  Steps.cereb_stat_beh(subpath,fullfile('BEHAV','Word'),ftpServer)
  ```
- step 8. fmri analysis
  - 1.  structure data
    ```matlab
    % structure
    % ===========================================
    % step1 isolate the cerebellum and brainstem from structure data
    % step2 normalize structure data to template space
    
    % input : filepath,"string" T1 file folder path
    %         filename,"string" T1 file name 
    %         f, "struct", flags for steps
    Steps.ana_steps(fullfile(indipath,'T1'),'SUB*T1.nii',f);
    % ==========================================
    ```
  - 2. get conditon
    ``` matlab
    % get condition variable
    fst_con = Setup.Getcondition(r,'behpath',fullfile(indipath,'BEHAV',fuc_round{nround}));
    ```
  - 3. functional data preprocessing
    ```matlab
    % functional
    % ============================================
    % step1 - preprocessing 
    %           - realignment
    %           - slice timing
    %           - corigister to structurnal space
    %       - output : auNAME.nii
    Steps.func_step1(fullfile(indipath,'T1'), ...
          'SUB*.nii', ...
          fullfile(indipath,fuc_round{nround}), ...
          'SUB*.nii', ...
          str2double(info.RT), ...
          str2num(info.sliceorder{1}), ...
          batchfolder,f);
    ```
  - 4. function data 1st level
    ```matlab
    % step2 - 1st level
    %           -model design and estimate
    %           -define contrast
    %       - output : T_CONTRAST.nii
    outputFolder = fullfile(fullfile(indipath,fuc_round{nround}),'1st_level');
    Steps.func_step2(fullfile(indipath,fuc_round{nround}), ...
        'au*.nii', ...
        outputFolder, ...
        str2double(info.RT), ...
        'secs', ...
        batchfolder, ...
        fst_con.condition, ...
        fst_con.con, ...
        fst_con.conWt, ...
        fst_con.conName, ...
        fst_con.condi, ...
        f);
    ```
  - 5. function data reslice to SUIT
    ```matlab
    % step 3 reslice and normalize to SUIT space
    %        - get functional contrast data
    %           - reslice to SUIT (mask --> T1 'c_name_pcereb.nii')
    %           - normalize to SUIT 
    %        - output : wsuit_NAME.nii
    T1filepath = fullfile(indipath,'T1');
    T1filename = 'c_*_pcereb.nii';
    funcfilepath = fullfile(indipath,fuc_round{nround},'1st_level');
    funcfilename = [];
  
    Steps.func_step3(T1filepath,T1filename,funcfilepath,funcfilename,batchfolder,f);
    ```
- step 6 add Native ROI for individual subject
  ```matlab
  % from ROI name and coordinate
  % Get ROI .nii file in subject dir "ROI" folder 
  subject_roi_path = fullfile(subpath,conn_sub{nsub},'ROI');
  for nsess = 1:length(sess)
      for nROI = 1:length(ROIcor)
          T1filepath = fullfile(subpath,sub{nsub},sess{nsess},'T1');
          Steps.GetNatROI(ROIcor{nROI},T1filepath,subject_roi_path,ROIname{nROI})
      end
  end
  ```
  
- step 7 create conn project to get connectivity(for one subject)
  ```matlab
  % set steps and run conn batch
  %   get conn project
  %   get conn project <prjNAME>_info.xlsx represent run steps in this project
  Steps.create_nproj(conn_steps,conn_slice_order,...
  'sub',conn_sub(nsub), ...
  'sess',conn_ses, ...
  'round',conn_round, ...
  'condition',conn_conditon, ...
  'conn_proj_name',conn_prjName, ...
  'conn_proj_path',conn_prjPath, ...
  'TR',conn_TR, ...
  'smooth_kernel',conn_fwhm, ...
  'StrucVres',conn_strVres, ...
  'funcVres',conn_funcVres, ...
  'filter_band',conn_filtBand, ...
  'analysisName',conn_AnalysisName, ...
  'contrast',conn_contrast,...
  'mrifilepath',subpath,...
  'rawdata',conn_rawdata, ...
  'prepsteps',conn_prepsteps, ...
  'roiPath',subject_roi_path, ...
  'AnaSource',conn_AnaSource)
  ```

- step 8 Get conn 1st level result, reslice to SUIT and normalize to MNI then move to 1st level result folder
```matalb
line 377 to 449
```

- step 9 run 2nd level for 1st level result
```matlab
Steps.second_level(subpath,rnd,fullfile(subpath,'2nd_level'),f, ...
        'method','pairT', ...
        'infilepat','sw*', ...
        'subidntfr',subidntfr, ...
        'conName',fst_con.conName, ...
        'rmSUB',rmSUB, ...
        'SPMorCONN','CONN', ...
        'datafolder',conn_ResultF, ...
        'space','MNI')
```


    
# SUIT-Toolbox   
#### SUIT_website: [[SUIT website](https://www.diedrichsenlab.org/imaging/suit.htm)]  
#### toolbox Downloads: [[author GitHub](https://github.com/jdiedrichsen/suit)]

## for fmri analysis
### functional data
- step 1 : preprocessing for functional data using SPM
  - realignment -- `input : <NAME>.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;`output : u<NAME>.nii`
  ``` matlab
  SPM batch
  ```  
  - slice timing -- `input : u<NAME>.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;`output : au<NAME>.nii`  
  ``` matlab
  SPM batch
  ```  
  - coregister functional data to structure data  -- `input : au<NAME>.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  `output : au<NAME>.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`au<NAME>.mat`
  ``` matlab
  SPM batch
  ```  
- step 2 : 1st level analyis (define contrast) -- `input : au<NAME>.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`output : beta*.nii or con*.nii`  
  ``` matlab
  SPM batch
  ```  
- step 3 : reslice functional data to SUIT space -- `input : beta*.nii or con*.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `mc_*_snc.mat (structure file)`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `c_*_pcereb.nii (structure mask file)`    
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`output : wc*.nii`  
  ```matlab
  suit_reslice(beta*.nii or con*.nii, mc_*_snc.mat, 'mask', c_*_pcereb.nii)
  ```  
- step 4 : smooth output data -- `input : wc*.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`output : swc*.nii`
  ``` matlab
  SPM batch
  ```  
- step 5 : group level analysis using -- `input : swc*.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`output : swc*.nii`
  ``` matlab
  SPM batch
  ```  

### structure data
- step 1 : segmentation using SUIT -- `input : <NAME>.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  `output : c1_<NAME>.nii (cereb grey matter)`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`c2_<NAME>.nii (cereb white matter)`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`c3_<NAME>.nii (cereb CSF)`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`c7_<NAME>.nii (brain grey matter)`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`c8_<NAME>.nii (brain white matter)`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`c_<NAME>.nii (cereb bound)`   
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`c_N<NAME>_pcereb.nii (cereb mask)`
  ```matlab
  suit_isolate_seg({<NAME>.nii});
  ```  
- step 2 : normalize -- `input : c_<NAME>.nii`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`c_<NAME>_pcereb.nii (mask)`  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`output : wsuit_mc_<NAME>.nii`
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`output : mc_<NAME>.nii`
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`output : mc_<NAME>.mat`

  
  ```matlab
  suit_normalize(c_<NAME>.nii,'mask',c_<NAME>_pcereb.nii (mask));
  ```  





