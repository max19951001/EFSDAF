

-----------------------------------------------------------------------------------------

​																				EFSDAF PROGRAM

​										Developed by Chenlie Shi, email:max1995@stumail.nwu.edu.cn

-----------------------------------------------------------------------------------------



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



### note

1. input images must be reflectance image. the range of dn values is 0<dn<1. data type is float

2. the resolution of MODIS image should be sampled to the resolution of Landsat image.



