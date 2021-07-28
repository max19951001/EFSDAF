

-----------------------------------------------------------------------------------------

​																				EFSDAF Pogram

​										Developed by Chenlie Shi, email:max1995@stumail.nwu.edu.cn

-----------------------------------------------------------------------------------------

## 2021.7.28  Update 
-----------------------------------------------------------------------------------------
### Conetent
   #### 
   (1) modify the exist bug for EFSDAF
   
   (2) add the instrution for EFSDAF Called instrution.txt
   
   (3) add a demo video for EFSDAF
   
   (4) provide the test data
   




-----------------------------------------------------------------------------------------
### Cite: 
 ##### Shi, C.; Wang, X.; Zhang, M.; Liang, X.; Niu, L.; Han, H.; Zhu, X. A Comprehensive and Automated Fusion Method: The Enhanced Flexible Spatiotemporal DAta Fusion Model for Monitoring Dynamic Changes of Land Surface. Appl. Sci. 2019, 9, 3693.




### Files

1.  SVD_Endmembers.csv 

   ​	the Global SVD Endmembers of Landsat images 

2. abundancecaculatemodule.exe

   ​	the fully constrained least squares (FCLS) program

3.  efsdaf.pro

   ​    main code of EFSDAF



### Usage

1.  The program is written in IDL; IDL and ENVI are required,
2.  Before running the program, build a new  IDL project , copy all the files to the project folder;
3.  Compile and execute (enter the file as prompted)



### Note

1. input images must be reflectance image. the range of dn values is 0<dn<1. data type is float

2. the resolution of MODIS image should be sampled to the resolution of Landsat image.



