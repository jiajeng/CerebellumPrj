# Work steps
- step 1. Set variables [[Setup.m](Setup.m)]
  Set some needed variables
  condition read from *behD.mat*
  ![Screenshot 2024-10-12 004900](https://github.com/user-attachments/assets/b813a64a-2adc-4fd2-881e-805f7f17577c)

- step 2. BEHAVE data folder in NAS (task condition) [[Setup.m]((Setup.m) -- GetBehFolder]
  Get behave data folder path and save in a .mat file -- ("*BEHFold_log.mat*")  
  ![Screenshot 2024-10-12 004900](https://github.com/user-attachments/assets/ae98aec9-d995-43f4-a728-7fe440db568b)


- step 3. BEHAVE data filepath (physio and task condition) [[Setup.m]((Setup.m) -- GetBehFolder]
  if file in nas then read ("*BEHFold_log.mat*") info to get file
  store all local path that store behave file in a .mat file -- ("*BEHFile_log.mat*")
  ![Screenshot 2024-10-12 004900](https://github.com/user-attachments/assets/a1ccc783-adff-4c4f-818c-7f68144a2807)

- step 4. convert IMA to NII [[conver2nii.m](convert2nii.m)]
  get IMA file in nas or local and convert to NII file and get MRinfo.xlsx
![Screenshot 2024-10-12 005450](https://github.com/user-attachments/assets/285fd194-6852-40e2-8d2e-4e33bfdbf359)
MRinfo.xlsx  
![Screenshot 2024-10-12 005545](https://github.com/user-attachments/assets/d24066dc-afe2-4d62-ac5f-97da7e1b844b)

> [!Note]
> MRinfo.xlsx store in this run data infomation  
> so run all data in one run is better, or need to manully fix MRinfo.xlsx file.

- step 5. Get NII file to local *if file in nas* [[niiFile2local.m](niiFile2local.m)]
  get NII file from nas to local  
  ![Screenshot 2024-10-12 005450](https://github.com/user-attachments/assets/c3e3c16c-dc85-40f1-a649-2cc394057ef9)

- step 6. Get MRinfo
  read MRinfo.xlsx to MRinfo struct

- step 7. Get BEHAVE data
  Get behave data (*from* "*BEHFold_log.mat*") folder info  
  Get individual subject behave data and save file in individual folder -- ("*behD.mat and behave.xlsx*")  
  ![Screenshot 2024-10-12 005545](https://github.com/user-attachments/assets/77045e22-f6f1-4cc1-b8c4-d4a1787fedf6)  
  ![Screenshot 2024-10-12 004900](https://github.com/user-attachments/assets/4b97ca8f-bcb4-45d9-a2fc-4adeeab11095)

- step 8. fmri analysis
 - 1.  structure data
        ![image](https://github.com/user-attachments/assets/e65472b8-f184-4a24-b1a6-292e26cd03b9)
 - 2. get conditon  
    ![image](https://github.com/user-attachments/assets/52cb4c16-1767-4533-b4f3-bcd930d2b5f1)
 - 3. functional data preprocessing  
    ![image](https://github.com/user-attachments/assets/7c8eea2d-5be1-40d3-9317-caa5966a962e)
 - 4. function data 1st level  
    ![image](https://github.com/user-attachments/assets/4e644f6c-0a5d-40e5-a384-17e71f0843ad)
    ![image](https://github.com/user-attachments/assets/4a18365a-b897-456d-8b12-96a6081c5aaf)
  > [!Note]
  > for original spm it should be con0001.nii and SPMT0001.nii  
  > need to change SPM code for get new name  
 - 5. function data reslice to SUIT
![image](https://github.com/user-attachments/assets/10398657-c852-4080-a164-425901169a03)


    
