# Cerebellum Project
## content   
[check data](#checkdata)  
[DataProcess](#dataprocess)  
[Result Folder](#Result)
     
## <a name="checkdata"></a> check data from NAS(120.xxx.xxx.xxx)
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

     
  
## <a name="dataprocess"></a> Data Process(SUIT) [function step info](./code/suit/README.md)  
![image](https://github.com/user-attachments/assets/e4bd7c81-f994-4a26-a9b5-87acabb8e039)  

### preprocessing  
- structure : segmentation  
- functional : realign, slice Timing, outlier detection, corigester to structure, denoise (conn toolbox)  
     
> [!Note]  
> 不做normalize and smooth，為了不讓visual cortex的資料影響到小腦的位置 
  
### first level  
- task : define condition for images (PPI analysis?)  
![image](https://github.com/user-attachments/assets/f09f7d03-1faa-47ff-b703-df6ca2a31f13)  
    
- rest : connectivity (ROI define by task)   
  ROI define in MNI space, so need to get every subject native space ROI coordinate    
  how? get every subject affine matrix(native space to MNI space), saved in T1 preprocessing file mc_*_snc.mat, variable `Affine`   
  
### normalize and smooth
- normalize to MNI(for whole brain) and SUIT(only for cerebellum)
    ###### normalize to MNI 的結果有把小腦的部分給切掉(mask out), 所以在結果上是不會有小腦的結果的。
    
### 2nd level
##### 46 subject, remove sub001, sub007, sub024, sub035, sub036, sub50 (behave outlier, no rest or word, bad data)
task : flexible ANOVA for smooth file
![image](https://github.com/user-attachments/assets/b3d79ddd-614a-4c7d-8b53-306c47e603d8)
    
rest : one sample T and pair T
![image](https://github.com/user-attachments/assets/0a615284-d042-4dd7-8d3e-c3d74c3322e5)
![image](https://github.com/user-attachments/assets/0117f515-22e5-4561-a1eb-bc496659a268)
  
  
## <a name="Result"></a> Result Folder
#### Nas
- basedir : `/LabData/jeng/Cerebellum_suit/ProcessData/2nd_level`
     
- ROI_coordinate : `./ROI_cor.txt`
     
- task : `./sub46/WORD/MNI` --> whole brain
- &ensp;&ensp;&ensp;&ensp;&ensp;`./sub48/WORD/SUIT` --> cerebellum
   
- Rest(one sample t) : `./sub48/OnesmpT_REST` --> brain
- Rest(one sample t) : `./sub48/OnesmpT_REST_suit` --> cerebellum
- Rest(one sample t) : `./sub48/OnesmpT_REST_wholeBrain` --> include brain and cerebellum
  
  

