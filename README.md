# Cerebellum Project
## step 0 check data from NAS(120.xxx.xxx.xxx)
#### **行為資料路徑** --> /LabData/cerebellum_NYCU/rawdata/date_SUBxxx/
- *行為資料 --> /BEHAV/word/ or /BEHAV*     
- *檔案  --> sub-xxx_task_wordread.csv(.mat)*
  
  | WordGroup | Word | WordOnset | Sound | SoundOnset | YesTrial | Response | ResponseTime | CorrectResponse |
  | --------- | ---- | --------- | ----- | ---------- | -------- | -------- | ------------ | --------------- |
  | LDLFY | 拯 | 6.9605 |	ㄓㄥˇ |	6.5064552 |	1	|	NaN |	NaN | NaN |
  | HDHFY |	節 | 16.8851757 | ㄐㄧㄝˊ | 16.2312326 |	1 |	2@ |	0.587837 | 1 |
  | .. | .. | .. | .. | .. | .. | .. | .. | .. |

- **total trial : 120**  
  -xDxFY : 22  
  -xDxFN : 8  
  -total : 30  

  | homophone  | Frequency | sound is match to the word or not |
  | -- | -- | -- |
  | Low | Low | Yes |
  | High | High | no |


#### **mri資料路徑** --> /LabData/cerebellum_NYCU/rawdata/date_SUBxxx/
- *fmri --> /T1  (結構資料.IMA, 192 file)* 
- *fmri --> /REST  (REST資料.IMA, 240 file)* 
- *fmri --> /WORD  (task資料.IMA, 401 file)*
