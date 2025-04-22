## PPI(gPPI)
- psychophysiological interaction  
![image](https://github.com/user-attachments/assets/18aa77dc-81b8-4e8a-b6df-bf87f0e67561)  
由task(condition)的contrast開始，figure 1最上面的圖
- 在一整個任務階段中屬於conditionA的標為1 其餘的標為 -1(for conditionA PPI)

- 在一個PPI的model中，有2*N+1個regressor
  1. condition A~N Main effect
  2. seed region time Series
  3. condition A~N PPI (ROI(seed) Bold signal x condition)  
       
  ###### *p.s* 因為Bold signal跟刺激會有延遲(haemodynamic lag)，在跟task condition convolve前要先把時間對好，在CONN以及FSL中是先把task跟HRF做convolution，但在SPM裡則是把Bold訓號跟HRF做deconvolution。




> - ### [Understanding psychophysiological interaction and its relations to beta series correlation](https://www.biorxiv.org/content/10.1101/322073v3.full)
> #### 1.1 Modeling of task main effect
>  假設有condition A and B使用regression model 來表示有兩種方式  
> 1. regressor1 屬於A的設為1 其餘為0 | regressor2 屬於B的設為1 其餘為0
> 2. regressor1 interceptrion | regressor2 屬於A的設為1 其餘為0    
> ![image](https://github.com/user-attachments/assets/06dd793d-df29-4216-8087-23749e312d98)   
> 這兩種方式在數學上是相等的(模型2進行線性組合可以得到模型1)  
> 但在解釋上就不一樣，方法1的解釋是每個condition的effect，方法2的解釋則是mean effect跟2個condition的差異。  
> ![image](https://github.com/user-attachments/assets/fde1e1e1-09af-4d1b-ac4e-911f0cbfa096)
> 方法2的模型可以表達成這樣
> x_Psych 由 1 0 組成
> 其中x_Psych表示的就是condtion A跟B的差異(x_Psych可以不用進行centering，不會影響B1的estimate，代表的都是兩個codntion的差異)
> y表示大腦的某區域的訊號
>
> #### 1.2 functional connectivity
> 要找到某個seed region跟其他region的connectivity(相關性)，常見的模型如下
> ![image](https://github.com/user-attachments/assets/0c1c8b5b-74dc-4f4a-aa78-9d2a212cd3bc)
> x_physio 代表Seed region的訊號，y一樣代表大腦某個區域(或是voxel)的訊號
> 把公式1跟2合起來，然後加上相互項，再整理一下  
> ![image](https://github.com/user-attachments/assets/a5b9c618-eaa2-4398-9f69-abd1802d614c)
> 所以大腦某個區域跟seed區域的effect = (B2+B3*x_psycho), 如果B3有顯著代表顯著這兩個區域(x_physio, y)的connectivity是跟任務是有相關的(受到任務的調節??)
>   
> 因為把psycho跟physio的資料要算interaction，所以Psych的數值要進行demean --> 假設coded成0 1 --> demean以後就會變成-0.5 0.5(讓physio跟psycho的資料正交)，如果不demean的話，可能會產生A是正相關B沒有相關或是A沒有相關B是負相關，但是都可以得到A跟B的差異就是了。
> ![image](https://github.com/user-attachments/assets/b4a3c2de-79b6-45e7-871d-59af48a6de1b)
> 模型的樣子長這樣，3 4是interaction項(1 2乘上x_physio)
> ![image](https://github.com/user-attachments/assets/0c6662e8-2989-49f5-99dc-121a561d1590)
> 或是把模型設計成這樣，3 4是interaction項，那這樣的解釋就是每個condition自己的Effect
> ![image](https://github.com/user-attachments/assets/89b790d9-8624-4110-9340-f806dfef07bc)
>
> #### more than two condition
> 
> #### convolution and deconvolution 
> - 通常來說fmri分成
>   - neuronal activity --> hypothetical
>   - blood-oxygen-level dependent(BOLD) --> observed
> - 當有一個event產生的時候，神經的活動我們視為一個脈衝波，那因為在BOLD訊號上，會在脈衝波上乘上gamma distribution，通稱為HRF(hemodynamic response function)
> - 在fmri中，會使用定義好event的位置(box-car function)，跟gamma分布做捲積，產生理想(假設)的BOLD訊號。這是在BOLD訊號上去檢驗相關性。  
> ![image](https://github.com/user-attachments/assets/009dce6e-3452-4689-9d7b-645da6ab5536)  
> z_Psych : box-car function  
> h : hrf   
> x_Physio : BOLD signal  
> - 在大腦不同的區域中，interaction會被認定是在neuroal level上，而不是在BOLD level上。  
> - 所以另一個方式是先對BOLD signal先做deconvolution再乘上box-car的序列，最後再跟HRF convolve  
> - ![image](https://github.com/user-attachments/assets/2d1d98fd-e16a-4315-95c4-78b4ae05438a)  
> z_Psych : box-car function  
> z_Physio : bold signal deconvolution  
> h : hrf  
> - 這兩種方式在Block design的設計上結果會類似，但在event design中可能還是會需要



