
## step 1 convert .IMA file to .nii file [[conver2nii.m]](https://github.com/jiajeng/Cerebellumprj/blob/main/code/step1_convert2nii.m)
- 從nas(/rawdata)讀取資料
- 將轉完的檔案放在nas(/process/file_prep)上

## step 2 get file from Nas and Organize behave data [[Getfile.m]](https://github.com/jiajeng/Cerebellumprj/blob/main/code/step2_GetFile2Local.m)
- 從nas(/process/file_prep/SUBxxx/niifile/xxx)上，讀取需要的.nii檔案，每個session跟結構各一個檔案  
-- *T1* --> s*.nii  
-- *Rest* --> SUBxxx_4D.nii  
-- *WORD* --> SUBxxx_4D.nii  
- 將行為資料重新整理，把WordGroup分成六個組別
  - **OUTLIER** -- 任何組別的outlier( ResponseTime mean+-2*std )
  - **HDHF** -- 屬於HDHF的回答正確的trial以及去掉outlier的
  - **HDLF** -- 屬於HDLF的回答正確的trial以及去掉outlier的
  - **LDLF** -- 屬於LDLF的回答正確的trial以及去掉outlier的
  - **LDHF** -- 屬於LDHF的回答正確的trial以及去掉outlier的
  - **ERROR** -- 回答錯誤的trial

  | WordGroup | Word | WordOnset | Sound | SoundOnset | YesTrial | Response | ResponseTime | CorrectResponse |
  | --------- | ---- | --------- | ----- | ---------- | -------- | -------- | ------------ | --------------- |
  | HDHFY | .. | .. | .. | .. | .. | .. | .. | 1 |
  | HDHFN | .. | .. | .. | .. | .. | .. | .. | 1 |
  
  | WordGroup | Word | WordOnset | Sound | SoundOnset | YesTrial | Response | ResponseTime | CorrectResponse |
  | --------- | ---- | --------- | ----- | ---------- | -------- | -------- | ------------ | --------------- |
  | HDLFY | .. | .. | .. | .. | .. | .. | .. | 1 |
  | HDLFN | .. | .. | .. | .. | .. | .. | .. | 1 |

  | WordGroup | Word | WordOnset | Sound | SoundOnset | YesTrial | Response | ResponseTime | CorrectResponse |
  | --------- | ---- | --------- | ----- | ---------- | -------- | -------- | ------------ | --------------- |
  | LDHFY | .. | .. | .. | .. | .. | .. | .. | 1 |
  | LDHFN | .. | .. | .. | .. | .. | .. | .. | 1 |

  | WordGroup | Word | WordOnset | Sound | SoundOnset | YesTrial | Response | ResponseTime | CorrectResponse |
  | --------- | ---- | --------- | ----- | ---------- | -------- | -------- | ------------ | --------------- |
  | LDLFY | .. | .. | .. | .. | .. | .. | .. | 1 |
  | LDLFN | .. | .. | .. | .. | .. | .. | .. | 1 |

  | WordGroup | Word | WordOnset | Sound | SoundOnset | YesTrial | Response | ResponseTime | CorrectResponse |
  | --------- | ---- | --------- | ----- | ---------- | -------- | -------- | ------------ | --------------- |
  | LDHFY | .. | .. | .. | .. | .. | .. | .. | 0 |
  | LDLFY | .. | .. | .. | .. | .. | .. | .. | 0 |


*p.s. run this step in 804 ubuntu pc*  
*file --> /media/tcnl/SSD1/jeng/cerebellum*

## step 3 preprocessing using conn toolbox [[conn_perp.m]](https://github.com/jiajeng/Cerebellumprj/blob/main/code/step3_conn_prep.m)
### - define parameter
  - TR = 1.5
  - voxelresolution = same as functionals(2mm)
  - sliceorder = [sliceorder.mat]
  - codition = using SoundOnset time to condition onset time
  - duration = 0(event)
  - FWHM = 8
  - voxelsize_anat = 1;
  - voxelsize_func = 2;
### - preprocessing[[conn_toolbox]](https://web.conn-toolbox.org/fmri-methods/preprocessing-pipeline)
![image](https://github.com/jiajeng/Cerebellumprj/assets/173427049/9dcc15d2-381b-4980-a02e-69c9f1b6b75d)
### - denoising
  - confound
    -  white Matter(5P,PCA top 5 component)
    -  CSF(5P,PCA top 5 component)
    -  realignment(12P,PCA top 2 component) --> 6+6(x y z-shifting and rotation)
    -  scrubbing(2P,PCA top 2 component)
  -  filter = [0.008 0.09]
  -  linear detrend = 1st order
### - fileDir
  - 在原本fmri資料夾內會產生處理完的資料，要拿dswauSUBxxx.nii的檔案出來處理
  - T1 拿wc0csSUBxxxxx.nii的檔案

### -first level alalysis (REST session)
 - [seed base connectivity](https://web.conn-toolbox.org/fmri-methods/connectivity-measures/seed-based)

*p.s. run this step in 804 ubuntu pc*  
*file --> /media/tcnl/SSD1/jeng/cerebellum*

> [!warning]
> **rest 跟 task 不要放在一起跑**

## step 4 connectivity analysis using SPM toolbox(WORD task analysis) [[SPM_analysis.m]](./../code/step4_SPM_analysis_1stLevel.m)
### define 1st level module
- 5 condition 4 condition parameter 1 constant
![image](https://github.com/user-attachments/assets/733992d9-6c45-4b10-aeba-d3e413601578)

| HDHF | HDHF_RT | HDLF | HDLF_RT | LDHF | LDHF_RT | LDLF | LDLF_RT | ERROR | constant |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

***p.s.  run this step in 804 ubuntu pc***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ***store file in 811 pc --> C:\Users\user\Desktop\jeng\CerebellumPrj\data\all\RESULT***  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ***run below step in 811 pc***  

### define contrast T
every condition is 1 (for anova only need define this)

| HDHF | HDLF | LDHF | LDLF | 
| -- | -- | -- | -- |
| 1 | 1 | -1 | -1 |

| HDHF | LDHF | 
| -- | -- |
| 1 | -1 |

| HDLF | LDLF | 
| -- | -- |
| 1 | -1 |

13 contrast


### 2nd level statitistical

- one sample t-test  
![image](https://github.com/jiajeng/Cerebellumprj/assets/173427049/29568a09-179d-4aa5-b4a1-88aafb2e82a4)  
select one first level contrast for all subject(only one contrast for 2nd level contrast)
  
- flexible ANOVA  
![image](https://github.com/jiajeng/Cerebellumprj/assets/173427049/6b6a59f8-9d3b-4b8a-8728-30b0adf77b40)  
select all wanted first level contrast for all subject (top 47 is subject and other is conditions)

condition matrix(code condition)
|   | D | F |
|--|--|--|
| H | 1 | 1 |
| L | 2 | 2 |
  
HDHF --> 1 1  
HDLF --> 1 2  
LDHF --> 2 1  
LDLF --> 2 2  

| subject_id | condition code |
|--|--|
|1| 1 1|
|2| 2 1|
|..|..|


--> subject 1 condition HDHF
--> subject 2 condition LDHF

see this paper [Contrast weights in flexible factorial design with multiple groups of subjects](https://www.researchgate.net/publication/267779738)  

> [!Note]
> - replace below file in ./spm12/
>   - [spm_contrasts.m](./../rep_code/spm_contrasts.m)
>     - automatically overwriting exist constrast define and rename .nii file to contrast name instead of append to con_00xx.nii
>   - [spm_run_fmri_spec.m](./../rep_code/spm_run_fmri_spec.m)
>     - automatically overwriting SPM.m file(read overwrite.m file store in code work dir)
>   - [spm_run_factorial_design](./../rep_code/spm_run_factorial_design.m)
>     - automatically overwriting SPM.m file(read overwrite.m file store in code work dir)


## step 4 psychopysiological interactions(PPI) analysis using SPM toolbox(WORD task analysis) [[SPM_analysis.m]](./../code/step4_SPM_PPI_analysis.m)
### create 1st level
- regressor  

| HDHF | HDLF | LDHF | LDLF | ERROR |
|--|--|--|--|--|

### create VOI file (roi file to SPM analysis)  

| id |right_cerebellum|left_cerebellum| auditory|
| -- | -- | -- | -- |
| 1 | 24 -72 -48 | -24 -72 -48 |  62 -6 -2 |
| 2 | 10 -80 -24 | -10 -80 -24 | -58 -10 2 | 
| 3 | 32 -64 -26 | -32 -64 -26 | .. |
| 4 | 22 -68 -26 | -22 -68 -26 | .. |
| 5 | 42 -50 -30 | -42 -50 -30 | .. |
| 6 | 36 -56 -24 | -36 -56 -24 | .. |
| 7 | 12 -57 -30 | -12 -57 -30 | .. |
| 8 | 17 -65 -35 | -17 -65 -35 | .. |

### create PPI variable
- compute every interaction term of all condition
- or define contrast first

### create PPI 1st level model
- use all condition(2*4+1 regressor)
  
| HDHF_inter | HDHF | HDLF_inter | HDLF | LDHF_inter | LDHF | LDLF_inter | LDLF | seed BOLD |
|--|--|--|--|--|--|--|--|--|

- use contrast data(3 regressor)
  
| HD-LD_inter | HD-LD | seed BOLD |
|--|--|--|

### create 1st level contrast

### 2nd level analysis

> [!Note]
> - replace below file in ./spm12/
>   - [spm_contrasts.m](./../rep_code/spm_contrasts.m)
>     - automatically overwriting exist constrast define and rename .nii file to contrast name instead of append to con_00xx.nii
>   - [spm_run_fmri_spec.m](./../rep_code/spm_run_fmri_spec.m)
>     - automatically overwriting SPM.m file(read overwrite.m file store in code work dir)
>   - [spm_run_factorial_design](./../rep_code/spm_run_factorial_design.m)
>     - automatically overwriting SPM.m file(read overwrite.m file store in code work dir)
>   - [spm_peb_ppi.m](./../rep_code/spm_peb_ppi.m)
>     - change save file dir(SPM.mat dir create PPI folder ./PPI/ppi*.mat)


### step 5 conn second level[[step5_conn_2nd_level]](./code/step5_conn_2nd_level.m)

- REST condition compute second level analysis using code
- select ROI in first level folder list_source.txt
- compute all ROI in 2nd level analysis
