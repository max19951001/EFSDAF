;+
; :Author: max
;date:2019/5/5


;****************************改进的FSDAF算法*************************
;**********************************************************************************************
;
;改进一：主要是针对FSDAF中landsat数据未考虑混合像元的问题，特别是在城市等异质性区域。
;改进二：针对相似像素的搜寻，采用地理学第一定律获取邻近邻近像元当作相似像素
;改进三：针对Modis和landsat像素之间的传感器差异进行调整。
;
;总体来说，结合istrum和fsdaf的算法。
;
;特点：  (1)考虑了以往基于像素融合算法未考虑高分辨率混合像元分解的问题。克服了类内可变性以及更符合实际情况
;     (2)将传感器差异考虑进去，提高了预测的准确性
;     (3)引入tps插值对突变的事件有更好的预测
;     (4)结合starfm的距离加权，平滑预测结果。
;更新：
;     第二版      2019/6/3 更新，加上判别landsat基准日期会出现异常值的问题
;     第三版      2019/6/3 更新，在解混前是否判断选择地物类型未发生变化的像元
;     第三版      2019/6/4 更新，薄板样条函数插值的参数对最终结果的影响
;     第五版     2019/6/5 更新，关于最终权重函数的建立，考虑相似像素阈值的设定
;     第六版     2019/7/1 更新，改进输入和fsdaf一样，保证tps预测一致
;     第七版     2019/7/2 更新，改进同质系数的确定，tb和tp的tps插值
; 其他可能改进的地方：(1)残差分配函数的选择、
;               (2)距离加权的确认

;代价函数
function cost_fun, g
  common v_pub1,ind_v
  common v_pub2,dep_v
  l=total((ind_v ## g-dep_v)^2)
  return, l
end

pro   efsdaf


  compile_opt idl2
  envi,/restore_base_save_files
  envi_batch_init

  ;************************************************输入参数***************************************************
  ;*********************************************************************************************************
  fbase_fp=dialog_pickfile(title='打开基准landsat数据')                                ;基准landsat
  cbase_fp=dialog_pickfile(title='打开基准modis数据')                               ;基准modis
  cpre_fp=dialog_pickfile(title='打开预测modis数据')                                 ;预测modis                          ;预测结果
  em=dialog_pickfile(title='打开csv数据')                                    ;端元数据
  fc_ratio=20                                                          ;modis和landsat像素比例
  half_win_size=2                                                      ;移动窗口半径，值为modis像素，移动窗口大小要大于端元数
  ab_threshold=0.05                                                    ;端元丰度阈值，默认0.05
  ;**************************************************************************************************************
  ;******************************************************打开数据**************************************************
  routine_dir=file_dirname(routine_filepath('efsdaf'))+'\'  ;获取程序目录
  print,'routine_dir=========',routine_dir
  cd,routine_dir

  ;打开数据为重采样为30m的landsat和modis数据
  envi_open_file,fbase_fp,r_fid=fb_fid
  if fb_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  envi_open_file,cbase_fp,r_fid=cb_fid
  if cb_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  envi_open_file,cpre_fp,r_fid=cp_fid
  if cp_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  ;文件查询
  envi_file_query,fb_fid,dims=fb_dims,nl=fb_lines,ns=fb_samples,nb=fb_bands,data_type=fb_dt
  envi_file_query,cb_fid,dims=cb_dims,nl=cb_lines,ns=cb_samples,nb=cb_bands,data_type=cb_dt
  envi_file_query,cp_fid,dims=cp_dims,nl=cp_lines,ns=cp_samples,nb=cp_bands,data_type=cp_dt
  ;获取空间参考
  fb_mapinfo =envi_get_map_info(fid=fb_fid)
  cb_mapinfo =envi_get_map_info(fid=cb_fid)

  cb_img   =make_array(cb_samples,cb_lines,cb_bands,type=cb_dt)
  ;fb_img,cb_img,cp_img分别为输入的landsat,基准日期modis，预测日期modis
  cp_img   =make_array(cp_samples,cp_lines,cp_bands,type=cp_dt)
  fb_img   =make_array(fb_samples,fb_lines,fb_bands,type=fb_dt)
  ;将读取的影像存在数组中
  for nb=0,cb_bands-1 do cb_img[*,*,nb]=envi_get_data(fid=cb_fid,dims=cb_dims,pos=nb)
  for nb=0,cp_bands-1 do cp_img[*,*,nb]=envi_get_data(fid=cp_fid,dims=cp_dims,pos=nb)
  for nb=0,fb_bands-1 do fb_img[*,*,nb]=envi_get_data(fid=fb_fid,dims=fb_dims,pos=nb)
  envi_file_mng,id=cb_fid,/remove
  envi_file_mng,id=cp_fid,/remove
  envi_file_mng,id=fb_fid,/remove
  ;******************************以上输入数据都为landat 30m分辨率*********************************************************


  ;校正输入landsat fb_img异常值，
  for ib=0,fb_bands-1, 1 do begin
    sortindex = sort(fb_img[*,*,ib]);返回从小到大排序的索引
    percentiles=[0.0001, 0.9999]
    f_sortindices = (findgen(float(fb_samples)*fb_lines+1))/(float(fb_samples)*fb_lines);创建一个百分比数组
    f_dataindices = value_locate(f_sortindices, percentiles);f_dataindices返回0.0001和0.9999值在f_sortindices的位置，
    data_1_4= (fb_img[*,*,ib])[sortindex[f_dataindices]];返回一个两个元素数组，存储允许的最小和最大值
    ;查找fb_img中最小允许值和最大允许值，f_sortindex[f_dataindices]返回最小值所在fine1中的索引
    ;print,data_1_4
    ;校正过于小的值
    ind_small=where(fb_img[*,*,ib] le data_1_4[0] or fb_img[*,*,ib] lt 0.0)
    ;print,'ind_small',ind_small
    temp=fb_img[*,*,ib]
    temp[ind_small]=min((fb_img[*,*,ib])[where(fb_img[*,*,ib] gt data_1_4[0] and fb_img[*,*,ib] ge 0.0)])
    fb_img[*,*,ib]=temp
    ;校正过于大的值
    ind_large=where(fb_img[*,*,ib] ge data_1_4[1] or fb_img[*,*,ib] gt 1.0)
    temp=fb_img[*,*,ib]
    temp[ind_large]=max((fb_img[*,*,ib])[where(fb_img[*,*,ib] lt data_1_4[1] and fb_img[*,*,ib] le 1.0)])
    fb_img[*,*,ib]=temp
  endfor
  ;(2)判断base modis数据
  for ib=0,nb-1, 1 do begin
    ;获取大于0.0001小于0.9999的索引位置
    percentiles=[0.0001, 0.9999]
    sortindices = (findgen(float(fb_samples)*fb_lines+1))/(float(fb_samples)*fb_lines);创建一个百分比数组
    dataindices = value_locate(sortindices, percentiles);dataindices返回0.0001和0.9999值在sortindices的位置
    cb_sortindex = sort(cb_img[*,*,ib]);返回从小到大排序的索引
    cb_data= (cb_img[*,*,ib])[cb_sortindex[dataindices]];返回一个两个元素数组，存储允许的最小和最大值
    ;校正过于小的值
    ind_small=where(cb_img[*,*,ib] le cb_data[0] or cb_img[*,*,ib] lt 0.0)
    temp=cb_img[*,*,ib]
    temp[ind_small]=min((cb_img[*,*,ib])[where(cb_img[*,*,ib] gt cb_data[0] and cb_img[*,*,ib] ge 0.0)])
    cb_img[*,*,ib]=temp
    ;校正过于大的值
    ind_large=where(cb_img[*,*,ib] ge cb_data[1] or cb_img[*,*,ib] gt 1.0)
    temp=cb_img[*,*,ib]
    temp[ind_large]=max((cb_img[*,*,ib])[where(cb_img[*,*,ib] lt cb_data[1] and cb_img[*,*,ib] le 1.0)])
    cb_img[*,*,ib]=temp
  endfor
  ;(3)判断pre modis数据
  for ib=0,nb-1, 1 do begin
    ;获取大于0.0001小于0.9999的索引位置
    percentiles=[0.0001, 0.9999]
    sortindices = (findgen(float(fb_samples)*fb_lines+1))/(float(fb_samples)*fb_lines);创建一个百分比数组
    dataindices = value_locate(sortindices, percentiles);dataindices返回0.0001和0.9999值在sortindices的位置
    cp_sortindex = sort(cp_img[*,*,ib]);返回从小到大排序的索引
    cp_data= (cp_img[*,*,ib])[cp_sortindex[dataindices]];返回一个两个元素数组，存储允许的最小和最大值
    ;校正过于小的值
    ind_small=where(cp_img[*,*,ib] le cp_data[0] or cp_img[*,*,ib] lt 0.0)
    temp=cp_img[*,*,ib]
    temp[ind_small]=min((cp_img[*,*,ib])[where(cp_img[*,*,ib] gt cp_data[0] and cp_img[*,*,ib] ge 0.0)])
    cp_img[*,*,ib]=temp
    ;校正过于大的值
    ind_large=where(cp_img[*,*,ib] ge cp_data[1] or cp_img[*,*,ib] gt 1.0)
    temp=cp_img[*,*,ib]
    temp[ind_large]=max((cp_img[*,*,ib])[where(cp_img[*,*,ib] lt cp_data[1] and cp_img[*,*,ib] le 1.0)])
    cp_img[*,*,ib]=temp
  endfor

  ;初始化index_f和index_c索引数组,目的是后面进行选择
  ns_c=ceil(float(fb_samples)/fc_ratio)
  nl_c=ceil(float(fb_lines)/fc_ratio) ;ns_c和nl_c分别为粗分辨率像素的个数
  print,ns_c,nl_c
  ii=0
  index_f=intarr(fb_samples,fb_lines)
  index_c=intarr(ns_c,nl_c)
  for i=0, ns_c-1, 1 do begin
    for j=0,nl_c-1,1 do begin
      index_f[i*fc_ratio:(i+1)*fc_ratio-1, j*fc_ratio:(j+1)*fc_ratio-1]=ii
      index_c[i,j]=ii
      ii=ii+1.0
    endfor
  endfor


  ;landsat行列索引值
  row_ind=intarr(fb_samples,fb_lines)
  col_ind=intarr(fb_samples,fb_lines)
  for i=0,fb_samples-1,1 do begin
    col_ind[i,*]=i ;列索引
  endfor
  for i=0,fb_lines-1,1 do begin
    row_ind[*,i]=i
  endfor

  ;采样modis像素到modis分辨率，此处是480米
  fine_c1=fltarr(ns_c,nl_c,fb_bands);输入的基准的landsat采样到480米，具体是取一个modis像素内的平均值
  coarse_c1=fltarr(ns_c,nl_c,fb_bands)
  coarse_c2=fltarr(ns_c,nl_c,fb_bands);coarse_c2和coarse_c1为480米modis像素值
  row_c=fltarr(ns_c,nl_c)
  col_c=fltarr(ns_c,nl_c);modis像素行列索引
  for ic=0,ns_c-1, 1 do begin
    for jc=0,nl_c-1, 1 do begin
      ind_c=where(index_f eq index_c[ic,jc])
      row_c[ic,jc]= mean(row_ind[ind_c])
      col_c[ic,jc]= mean(col_ind[ind_c])
      for ib=0,fb_bands-1,1 do begin
        fine_c1[ic,jc,ib]=mean((fb_img[*,*,ib])[ind_c])
        coarse_c1[ic,jc,ib]=mean((cb_img[*,*,ib])[ind_c])
        coarse_c2[ic,jc,ib]=mean((cp_img[*,*,ib])[ind_c])
      endfor
    endfor
  endfor
  print,"预处理完成"

  ;*************************************************预处理包括校正输入影像异常值，对输入的landsat像素采样到modis像素，即30m-480m.
  ;step 1  时间预测，此过程只能预测土地类型未发生变化的情况
  ;1 打开端元光谱文件
  ;*****************************************************************************
  em_spec=read_csv(em,header=em_name)
  em_samples=n_elements(em_name) ;返回列数，也就是端元数
  em_lines  =n_elements(em_spec.(0));em_lines等于波段数
  temp=fltarr(em_samples,em_lines)
  for i=0,em_samples-1 do begin
    temp[i,*]=float(em_spec.(i))
  endfor
  em_spec=temporary(temp);em_spec存储各个波段光谱值

  ;2 完全最小二乘解混
  ;*******************************************************************
  ;移动窗口大小
  win_size =2*half_win_size+1
  pixelcnts=win_size*win_size
  if pixelcnts le (em_samples) then begin ;移动窗口的平方大于端元数
    print,'移动窗口太小'
    envi_batch_exit
    return
  endif
  ;fcls
  print,"ok"
  fbaseabd_fp=file_dirname(fbase_fp)+'\'+file_basename(fbase_fp)+'_abundance.tif'
  cmdstr='abundancecaculatemodule.exe '+fbase_fp+' '+em+' '+fbaseabd_fp
  spawn,cmdstr,/hide
  envi_open_file,fbaseabd_fp,r_fid=fabd_fid
  if fabd_fid eq -1 then begin
    envi_batch_exit
    print,'光谱解混失败'
    return
  endif
  print,"约束最小二乘结束"

  ;3 丰度聚集，求取modis像素的丰度
  ;*********************************************************
  ;求解光谱角，判断相似端元
  lut=uintarr(em_samples,em_samples);存放光谱相似度的索引值，其值大小一次递减
  for i=0,em_samples-1 do begin
    cossam=fltarr(em_samples)
    for j=0,em_samples-1 do begin
      if j eq i then continue ;若i=j,跳出本次循环，进入下一次循环
      cossam[j]=total(em_spec[i,*]*em_spec[j,*])/(sqrt(total(em_spec[i,*]^2))*sqrt(total(em_spec[j,*]^2)))
    endfor
    lut[*,i]=reverse(sort(cossam)) ;lun存放数据值由大到小的索引值，也就是第一个值为相似性最大的索引
    ;print,"reverse",lut[*,i]
  endfor
  print,"光谱角计算"

  ;丰度聚集
  fabd_img=fltarr(fb_samples,fb_lines,em_samples)
  for ib=0,em_samples-1 do begin
    fabd_img[*,*,ib]=envi_get_data(fid=fabd_fid,pos=ib,dims=fb_dims)
    index_fabd=where(fabd_img[*,*,ib] lt 0.0,counts)
    if counts gt 0 then begin
      (fabd_img[*,*,ib])[index_fabd]=0.0 ;若存在异常值，赋值为0
    endif
  endfor
  envi_file_mng,id=fabd_fid,/remove,/delete
  cabd_img=fltarr(ns_c,nl_c,em_samples)   ;modis像素内的丰度值
  for cabd_i=0,nl_c-1 do begin
    for cabd_j=0,ns_c-1 do begin
      win_data=fabd_img[(cabd_j*fc_ratio):((cabd_j+1)*fc_ratio-1),(cabd_i*fc_ratio):((cabd_i+1)*fc_ratio-1),*]
      cabd=fltarr(em_samples)
      for num_em=0,em_samples-1 do begin
        cabd[num_em]=mean(win_data[*,*,num_em])
      endfor
      print,'cabd',cabd;输出端元的丰度值
      ;merge endmembers to their similar types if the abundance is lt ab_threshold
      index=where((cabd gt 0),cnts)
      print,cnts,'cnts'        ;endmembers with abundance gt 0
      abd_min=min(cabd[index],min_ind)
      while abd_min lt ab_threshold do begin
        min_ind=index[min_ind]
        abdtemp=cabd[min_ind]
        cabd[min_ind]=0.0
        for sc=0,em_samples-1 do begin
          if cabd[lut[sc,min_ind]] ne 0 then begin
            cabd[lut[sc,min_ind]]+=abdtemp
            break
          endif
        endfor
        index=where((cabd gt 0),cnts)
        abd_min=min(cabd[index],min_ind)
      endwhile
      ;判断modis像素丰度异常值
      if total(cabd) lt 0.999 then begin ;将丰度小于1的残差赋值给端元丰度最大的那个
        redi=1.0-total(cabd);残差值
        sortCabd=sort(cabd)
        maxAbundance=cabd[sortCabd[em_samples-1]]
        maxAbundance=maxAbundance+redi
        cabd[sortCabd[em_samples-1]]=maxAbundance
      endif
      cabd_img[cabd_j,cabd_i,*]=cabd
    endfor
  endfor
  print,"开始进行最小二乘解混"

  ;4 约束最小二乘解混，求取modis像素内每一类的变化
  ;*************************************************************
  change_f =make_array(fb_samples,fb_lines,cb_bands,type=cb_dt)
  change_c=coarse_c2-coarse_c1
  ;cf_img=make_array(cb_samples,cb_lines,cb_bands,type=cb_dt);存储的是基准日期一个modis像素内所有landsat的平均值
  he_index=fltarr(fb_samples,fb_lines,fb_bands)

  ;设置端元的变化阈值
  min_allow=fltarr(em_samples,cb_bands)
  max_allow=fltarr(em_samples,cb_bands)
  for ib=0,cb_bands-1,1 do begin
    min_allow[*,ib]=min(change_c[*,*,ib])-stddev(change_c[*,*,ib])
    max_allow[*,ib]=max(change_c[*,*,ib])+stddev(change_c[*,*,ib])
  endfor
  ;约束最小二乘参数设置

  common v_pub1
  common v_pub2
  gbnd    =[0,100]
  nobj    = 0
  lcomp   = 'cost_fun'
  for ci=0,nl_c-1 do begin
    for cj=0,ns_c-1 do begin
      ai=max([0,ci-half_win_size])       ; 移动窗口直径
      bi=min([nl_c-1,ci+half_win_size])
      aj=max([0,cj-half_win_size])
      bj=min([ns_c-1,cj+half_win_size])
      ;**********************窗口的大小直接影响着解混的结果，但是最终的端元变化还是针对一个modis像素而言的**************************************
      fai=ci*fc_ratio
      fbi=(ci+1)*fc_ratio-1
      faj=cj*fc_ratio
      fbj=(cj+1)*fc_ratio-1 ;landsat像素窗口位置，为一个modis像素
      c_win_pixels=(bi-ai+1)*(bj-aj+1)
      ;    print,c_win_pixels
      if (c_win_pixels gt em_samples) then begin   ;窗口大小大于端元数，继续执行
        fabd_temp=fabd_img[fai:fbi,faj:fbj,*]        ;窗口内landsat像素存在fabd_temp
        ;      cabd_temp=cabd_img[ci,cj,*]
        ind_v   = transpose(reform(cabd_img[aj:bj,ai:bi,*],c_win_pixels,em_samples))
        for nb=0,cb_bands-1 do begin
          dep_v = double(reform(change_c[aj:bj,ai:bi,nb],c_win_pixels))
          x     = fltarr(1,em_samples);x为modis像素内端元的变化值
          xbnd  = [[min_allow[*,nb]], [max_allow[*,nb]]]
          constrained_min, x, xbnd, gbnd, nobj, lcomp,inform, nstop = 5 ;此处x返回的是每个modis像素内的端元变化范围
          ;print,'xxxxx:      ',x
          ds_change=fabd_temp*rebin(reform(x,1,1,em_samples),fc_ratio,fc_ratio,em_samples,/sample);最邻近采样为了保证值不变的情况下进行矩阵相乘
          change_f[faj:fbj,fai:fbi,nb]=total(ds_change,3);change_f为每一landsat的变化值
          ;cf_img[cj,ci,nb]=mean(fb_img[faj:fbj,fai:fbi,nb]);cf_img为一landsat像素聚集为modis像素的值，便于调整传感器差异
        endfor
      endif
    endfor
  endfor
  print,'约束最小二乘结束'

  ;注意：*********************************此处求得landsat像素变化值存在一个modis块效应，即在一个modis像素中，landat的变化值是一样的所以，此处，应用一个移动窗口进行加权计算*************
  ;5 传感器差异调整，求取系数变化
  ;*********************************************************************
  fp_img_temp=make_array(fb_samples,fb_lines,fb_bands,type=fb_dt)
  for ib=0,fb_bands-1  do begin
    x =reform(fine_c1[*,*,ib],ns_c*nl_c)
    y =reform(coarse_c1[*,*,ib],ns_c*nl_c)
    ; print,"x",x
    ;  print,"y",y
    coefs=linfit(x,y)
    fp_img_temp[*,*,ib]=fb_img[*,*,ib]+(coefs[1]*change_f[*,*,ib])
  endfor

  ;相似性计算
  similar_th=fltarr(fb_bands)
  for iband=0,fb_bands-1,1 do begin
    similar_th[iband]=stddev(fb_img[*,*,iband])*2.0/em_samples
  endfor

  ;时间预测值进行聚合

  cp_img_temp=fltarr(ns_c,nl_c,cb_bands)
  for ci=0,nl_c-1 do begin
    for cj=0,ns_c-1 do begin
      for ib=0,fb_bands-1 do begin
        cp_img_temp[cj,ci,ib]=mean(fp_img_temp[cj:((cj+1)*fc_ratio-1),ci:((ci+1)*fc_ratio-1),ib])
      endfor
    endfor
  endfor
  ;
  ;
  print,'时间预测结束'
  ;6 时间预测，求取最终预测值
  ;***********************************************************************
  fpre_fp=dialog_pickfile(title="输出临时文件")
  openw,lun,fpre_fp,/get_lun
  writeu,lun,fp_img_temp
  free_lun,lun
  envi_setup_head,fname=fpre_fp,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write

  ;step 2  空间预测


  ;7 TPS薄板函数插值，求取土地覆盖发生变化信息
  ;************************************************************************
  ;tps插值
  
  ;基准日期modis插值
  tps1=fltarr(fb_samples,fb_lines,fb_bands)
  for ib=0,fb_bands-1,1 do begin
    tps1[*,*,ib] = min_curve_surf(coarse_c1[*,*,ib], col_c, row_c,/tps, xpout=col_ind,  ypout=row_ind)
  endfor
  ;预测日期modis插值
  tps2=fltarr(fb_samples,fb_lines,fb_bands)
  for ib=0,fb_bands-1,1 do begin
    tps2[*,*,ib] = min_curve_surf(coarse_c2[*,*,ib], col_c, row_c,/tps, xpout=col_ind,  ypout=row_ind)
  endfor
  similar_tps=fltarr(fb_bands)
  ;进行筛选地物类型发生变化的信息，假设地物类型变化发生在一个较小的区域。那么对于一个移动窗口内，若每一像素的tps差值大于similar_tps,则我们认为是地物类型发生变化的区域
  for ib=0,fb_bands-1 do begin
    similar_tps[ib]=abs(mean(tps2[*,*,ib])-mean(tps1[*,*,ib]))
  endfor
  print,'similar========',similar_tps
  fn=dialog_pickfile(title="tps结果")
  openw,lun,fn,/get_lun
  writeu,lun,tps1
  free_lun,lun
  envi_setup_head,fname=fn,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write
  print,'tps插值结束'
  ;*******************************************************************************************
  ;**********************************************************************************************
  he_index=fltarr(fb_samples,fb_lines,fb_bands)
  for fi=0,fb_lines-1 do begin
    for fj=0,fb_samples-1 do begin
      ai=max([0,fi-fc_ratio])
      bi=min([fb_lines-1,fi+fc_ratio]) ;搜索窗口，[fj,fi]为中心像素
      aj=max([0,fj-fc_ratio])
      bj=min([fb_samples-1,fj+fc_ratio])
      ;    print,"lajissssssssssssssssssssssssssss",ai,bi,aj,bj
      temp=fltarr(bj-aj+1,bi-ai+1,fb_bands)
      for ib=0,fb_bands-1 do begin
        temp[*,*,ib]=abs(tps2[aj:bj,ai:bi,ib]-tps1[aj:bj,ai:bi,ib])
        index=where(temp[*,*,ib] gt (similar_tps[ib]+0.001),nums)
        print,'nums======',nums
        if nums gt 0 then begin
          he_index[fj,fi,ib]=nums/((bj-aj+1)*(bi-ai+1)) ;此处he_index表示的是异质性区域的系数，越大的值表示赋值更多给异质性区域
        endif else begin
          he_index[fj,fi,ib]=0.0
        endelse
      endfor
    endfor
  endfor
  ;;*********************************************************************************************
  ;********************************************************************************************
  ;8 残差分配
  ;**************************************************************************
  predict_change_c=cp_img_temp-fine_c1
  real_change_c=coarse_c2-coarse_c1
  change_r=real_change_c-predict_change_c ;change_r为残差值
  change_21_c=fltarr(ns_c,nl_c,cb_bands) ;残差分配
  change_21=fltarr(fb_samples,fb_lines,fb_bands);landsat每一像素残差分配结果
  for ci=0,nl_c-1 do begin
    for cj=0,ns_c-1 do begin
      fai=ci*fc_ratio
      fbi=(ci+1)*fc_ratio-1
      faj=cj*fc_ratio
      fbj=(cj+1)*fc_ratio-1 ;modis像素内landat的索引
      fb_nums=float(fc_ratio)*fc_ratio ;一个modis内多少个lansat像素
      for ib=0,cb_bands-1 do begin
        diff_change=change_r[cj,ci,ib]
        w_change_tps=(tps2[*,*,ib])[faj:fbj,fai:fbi]-(fp_img_temp[*,*,ib])[faj:fbj,fai:fbi];文中Eh0
        if (diff_change le 0) then begin ;diff-change 小于0，则w_change_tps也要小于 0
          ind_noc=where(w_change_tps gt 0, num_noc)
          if (num_noc gt 0) then begin
            w_change_tps[ind_noc]=0
          endif
        endif else begin
          ind_noc=where(w_change_tps lt 0, num_noc)
          if (num_noc gt 0) then begin
            w_change_tps[ind_noc]=0
          endif
        endelse
        w_change_tps=abs(w_change_tps)
        w_unform=fltarr(fb_nums)     ;对于landsat像素级别
        w_unform[*]=abs(diff_change)
        ;w_change=w_change_tps*he_index[faj:fbj,fai:fbi,ib]+w_unform*(1.0-he_index[faj:fbj,fai:fbi,ib])+0.000001  ;combine these two weights
        w_change=w_unform*he_index[faj:fbj,fai:fbi,ib]+w_change_tps*(1.0-he_index[faj:fbj,fai:fbi,ib])+0.000001
        w_change=w_change/(mean(w_change)) ;归一化
        ;去除异常值
        ind_extrem=where(w_change gt 10, num_extrem)
        if (num_extrem gt 0) then begin
          w_change[ind_extrem]=mean(w_change)
        endif
        w_change=w_change/(mean(w_change)) ;w_change为加权值
        change_21[faj:fbj,fai:fbi,ib]=w_change*diff_change
      endfor
    endfor
  endfor
  ;将change_21加到fp_img_temp上
  ;允许t2变化范围
  min_allow=fltarr(fb_bands)
  max_allow=fltarr(fb_bands)
  for ib=0,fb_bands-1 do begin
    min_allow[ib]=min([min(coarse_c2[*,*,ib]),min(fp_img_temp[*,*,ib])])
    max_allow[ib]=max([max(coarse_c2[*,*,ib]),max(fp_img_temp[*,*,ib])])
  endfor

  fp_img_redis=fltarr(fb_samples,fb_lines,fb_bands)
  fp_img_redis=fP_img_temp+change_21
  for ib=0,fb_bands-1 do begin
    temp=fP_img_redis[*,*,ib]
    index_min=where(temp lt min_allow[ib],min_nums)
    if min_nums gt 0 then begin
      temp[index_min]=min_allow[ib]
    endif
    index_max=where(temp gt max_allow[ib],max_nums)
    if max_nums gt 0 then begin
      temp[index_max]=max_allow[ib]
    endif
    fp_img_redis[*,*,ib]=temp
    change_21[*,*,ib]=fp_img_redis[*,*,ib]-fp_img_temp[*,*,ib] ;残差值
  endfor
  ;fpre_fp1="C:\Users\15845\Desktop\ifsdaf.dat"
  ;openw,lun,fpre_fp1,/get_lun
  ;writeu,lun,fp_img_redis
  ;free_lun,lun
  ;envi_setup_head,fname=fpre_fp1,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write
  print,'残差分配结束'

  ;9 加权函数进行最终预测
  ;d_d_all为窗口内到中心像素的距离
  fp_img=fltarr(fb_samples,fb_lines,fb_bands);最终预测结果
  w=half_win_size*fc_ratio;landsat像素
  d_d_all=((w-indgen(w*2+1)#(intarr(1,w*2+1)+1))^2+(w-(intarr(w*2+1)+1)#indgen(1,w*2+1))^2)^0.5 ;d_d_all为完整窗口中心像素到周围其他像素的距离
  d_d_all=reform(d_d_all,(w*2+1)*(w*2+1))



  ;第二次写的
  for fi=0,fb_lines-1 do begin
    for fj=0,fb_samples-1 do begin
      fai=max([0,fi-w])
      fbi=min([fb_lines-1,fi+w])
      faj=max([0,fj-w])
      fbj=min([fb_samples-1,fj+w]) ;窗口范围
      i_tar=fi-fai
      j_tar=fj-faj ;目标像素位置
      col_wind=indgen(fbi-fai+1)#(intarr(1,fbj-faj+1)+1)*1.0 ;列的值都一样
      row_wind=(intarr(fbi-fai+1)+1)#indgen(1,fbj-faj+1)*1.0 ;行的值都一样
      ;search similar pixels within window
      similar_cand=fltarr((fbi-fai+1)*(fbj-faj+1)) ;pleace the similarity measure between each pixel and the target pixel
      position_cand=intarr((fbi-fai+1)*(fbj-faj+1))+1  ;place the location of each similar pixel
      for ib=0,fb_bands-1,1 do begin
        cand_band=intarr((fbi-fai+1)*(fbj-faj+1))
        wind_fine=fb_img[faj:fbj,fai:fbi,ib]
        s_s=abs(wind_fine-wind_fine[j_tar,i_tar])
        similar_cand=similar_cand+s_s/(wind_fine[j_tar,i_tar]+0.00000001)
        ind_cand=where(s_s lt similar_th[ib])
        cand_band[ind_cand]=1
        position_cand=position_cand*cand_band
      endfor
      indcand=where(position_cand ne 0,number_cand)  ;select similar pixel initially
      order_dis=sort(similar_cand[indcand])
      ;number_cand=min([number_cand0,num_similar_pixel])
      ind_same_class=indcand[order_dis[0:number_cand-1]]           ; select the n most similar samples
      ;compute weight for each simialr pixel
      ;spatial distance
      if ((fbi-fai+1)*(fbj-faj+1) lt (w*2.0+1)*(w*2.0+1)) then begin   ;not an integrate window
        d_d_cand=((i_tar-col_wind[ind_same_class])^2+(j_tar-row_wind[ind_same_class])^2)^0.5+0.00001
      endif else begin
        d_d_cand=d_d_all[ind_same_class]      ;integrate window
      endelse

      ;normalize these distances
      d_d_cand=(1.0+d_d_cand/w)*(similar_cand[ind_same_class]+1.0)
      c_d=1.0/d_d_cand
      weight=c_d/total(c_d)
      w_pixel=float(fc_ratio)*fc_ratio
      for ib=0,fb_bands-1,1 do begin
        change_cand=(change_21[faj:fbj,fai:fbi,ib])[ind_same_class]
        fp_img[fj,fi,ib]=fb_img[fj,fi,ib]+total(weight*change_cand)
        x=fb_img[fj,fi,ib]+cp_img[fj,fi,ib]-cb_img[fj,fi,ib]
        if fp_img[fj,fi,ib] lt 0.0 then begin
          another_predict=max([0,x])
          fp_img[fj,fi,ib]=min([1.0,another_predict])
        endif
        if fp_img[fj,fi,ib] gt 1.0 then begin
          another_predict=min([1.0,x])
          fp_img[fj,fi,ib]=max([0,another_predict])
        endif
      endfor
    endfor
  endfor
  print,"加权结束"

  fpre_fp1="C:\Users\15845\Desktop\ifsdaf_result2.dat"
  openw,lun,fpre_fp1,/get_lun
  writeu,lun,fp_img
  free_lun,lun
  envi_setup_head,fname=fpre_fp1,nb=fb_bands, ns=fb_samples,nl=fb_lines, interleave=0,data_type=fb_dt,map_info=fb_mapinfo,/write


end