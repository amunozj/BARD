pro amj_pnr_dt, CRD_in, mdi_i, PRs, NRs, instr, seg_const=seg_const, display=display, detdisp=detdisp, prt=prt, info=info, not_merge = not_merge, pnr_lbl = pnr_lbl
;+
; NAME:
;      amj_pnr_dt
;
; PURPOSE:
;       For a given MDI magnetogram, determine the areas of interest
;       using a multilayered approach to detection, giving priority to objects with more flux.
;       It operates by calling amj_pnr_dt_core sequentially       
;
; CALLING SEQUENCE:
;       amj_pnr_dt, CRD_in, mdi_i, PRs, NRs, instr
;
;
; INPUTS:
;       CRD_in: a structure holding the original and corrected images, as well as the calculated
;                heliographic coordinates
;            1.  im_raw: original magnetogram
;            2.  im_crr: magnetogram corrected using LOS component and assuming al field is radial
;            3.  hdr: header file of the input magnetogram
;            4.  mgnt_ar: elements of area corresponding to each pixel
;            5.  mgnt_flx: flux elements corresponding to each pixel
;            6.  Xar: x coordinate in a heliographic coordinate system 
;            7.  Yar: y coordinate in a heliographic coordinate system 
;            8.  Zar: z coordinate in a heliographic coordinate system
;            9.  Lath: heliographic latitude
;           10.  Lonh: heliographic longitude
;
;        mdi_i:  Mission day used to tag regions  
;
;       PRs and NRs: Structures holding positive and negative regions:
;           1. seg_const: used parameters in the region detection (see above)
;           2. hdr: header of the MDI magnetogram
;           3. mask: mask holding the position of each detected region
;           4. num_p: Number of positive regions
;           5. pregions: structure holding details from all detected positive regions 
;              (see below for details on the rgn structure) 
;           6. num_p: Number of negative regions
;           7. nregions: rgn structure holding details from all detected negative regions 
;              (see below for details on the rgn structure)
;           8. r_sw: Switch that indicates if there is at least one region of each sign
;
;   	instr: 	Variable that indicates which instrument is used:
;   			1. KPVT-512 files
;   			2. SPMG files
;   			3. MDI files
;   			4. HMI files
;
; OPTIONAL INPUTS:
;    seg_const: struct holding control parameters for segmentation used 
;               in the detection of positive and negative regions
;            1. ker_th:  kernel pixel upper threshold, default 150 G; used in segmentation
;            2. ker_th2: kernel pixel lower threshold, default 60 G; used in segmentation
;            3. ar_th: region pixel threshold, default 50 G;
;               used in region growth
;            4. ar_th2: post detection pixel threshold, default 40 G;
;               used in region growth
;            5. qr_th: quiet Sun region threshold, default 60 G
;            6. dila_size: dilation structural size, used to merge
;               the neibouring AR elements into a single AR used in morphological closing operation
;            7. eros_size: erosion structural size used to remove
;               small AR kernel pixels after the kernel pixel
;               segmentation; used in morphological opending operation
;            8. npssu: Number of levels used in the multi-layer detection
;               thresholds are choosen between ker_th and ker_th2
;            9. dis_lim:  Maximum distance in pixel radius inside which post-detection merging 
;               of PNRs is allowed
;           10. ovr_lim:  Minimum overlap required between regions in order to allow post-detection merging
;           11. mxB: Maximum field strength allowed for region storage
;           12. k_sig (alternate to 1): kernel pixel threshold = k_sig times the background
;               standard deviation; default 15
;           13. ar_grow_sig (alternate to 2): active region pixel
;               threshold = ar_grow_sig times background standard deviation; default 5           
;           14. valid_range: value range of valid magnetogram pixels,
;               default=[-20000,20000]; used to deal bad pixels in the iamges
;           15. deg_lim: maximum angle from disk center taken from
;               the snapshot
;        
;        pnr_lbl: Region label within a mission day.  Used in pre-fragmentation to associate
;                 a set of fragmented regions with their mother.
;
; OUTPUTS: 
;        PRs: Array of previusly detected positive regions. Positive regions
;             detected by this program will be appended to it.  Elements of this array
;             contain the rgn structure
;
;        NRs: Array of previusly detected negative regions. Negative regions
;             detected by this program will be appended to it.  Elements of this array
;             contain the rgn structure
;             
;       rgn structure:
;           1. lnk_sw: Flag that indicates if a region is part of a bipolar magnetic region
;           2. mdi_i:  MDI mission day used to tag regions
;           3. num: region index in that day
;           4. indx: indexes of all region pixels. Retrieve via long(strsplit(rgn.indx,/extract))
;           5. flux: magnetic flux (Mx)
;           6. area: region area (cm^2)
;           7. dm:   Dipolar moment 
;           8. fcn_lt: latitude of the flux weighted geographical center (deg)   
;           9. fcn_ln: longitude of the flux weighted geographical center (deg)   
;          10. dcen: average distance to geographical center (deg) (flux weighted)
;          11. fcenxp: x pixel coordinate of geographical center 
;          12. fcenyp: y pixel coordinate of geographical center 
;          13. dcenp: average distance to geographical center in pixels (flux weighted)
;          14. dm: Dipolar Moment
;          15. qm: Quadrupolar Moment  
;          
;                        
; KEYWORDS:
;      display: plots detected regions
;      detdisp: plots detected regions and all intermediate steps
;      prt: saves detected regions as an EPS figure
;      info: prints process information on screen
;      not_merge:  Prevents the code
;
; 
; MODIFICATION HISTORY:
;   2014/05/30:  Andres Munoz-Jaramillo:   initiated.
;   2015/06/09:  Andres Munoz-Jaramillo:   Adapted for KPVT data, including pre-fragmentation
;   2015/09/16:  Andres Munoz-Jaramillo:   Adapted to work with different instruments
;-


;define the usage
if n_params() lt 5 then begin
   print,'USAGE: amj_pnr_dt, CRD_in, mdi_i, PRs, NRs, instr'
   return
endif


;get the input image and dimension
im=CRD_in.im_raw
imsel=CRD_in.im_crr
sz=size(im)

;
;define the segmentation control structure and the default value
;
seg_const_def={ker_th:500.0, ker_th2:250.0, ar_th:250.0, ar_th2:50.0, eros_size:9.0, npssu: 3, dis_lim:2.0, ovr_lim: 0.4, qr_th:60.0, dila_size:2, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}


if not keyword_set(seg_const) then begin
    seg_const=seg_const_def
endif

;
;set the display window size and display threshold
;
display_zoom=0.25
display_xsize=sz[1]*display_zoom & display_ysize=sz[2]*display_zoom  

print_zoom=1.0
print_xsize=sz[1]*print_zoom & print_ysize=sz[2]*print_zoom 

display_thresholdu = 350.0
display_thresholdd = - display_thresholdu       

print, 'Detecting positive and negative regions'

;-----------------------------------------------------------------------------------------------------------
;First pass. Detectiong positive and negative regions using multiple sets of detecting criteria to enhance detection--------------------------------------------------------

;Decrement for lowering the detection threshold
d_thrs = 0;
if (seg_const.npssu gt 1) then begin
   d_thrs = (seg_const.ker_th - seg_const.ker_th2)/(seg_const.npssu-1)
endif

;Number of detected positive and negative regions
n_regp = 0
n_regn = 0

for n=0,seg_const.npssu-1 do begin

    ;New detection threshold
    n_thrs = seg_const.ker_th - d_thrs*n
;    n_thrs = seg_const.ker_th2 + d_thrs/(3.0^n)
;    n_thrs = seg_const.ker_th/(n+1)
;    n_thrs = seg_const.ker_th/(2.0^n)
;     n_thrs = seg_const.ker_th2
    
;    print, n_thrs
    
    ;Creating new detection structure and input magnetogram
    seg_const_tmp = seg_const
    seg_const_tmp.ker_th = n_thrs
;    seg_const_tmp.ar_th = seg_const.ar_th/(n+1)
;    seg_const_tmp.ar_th = seg_const.ar_th/(2.0^n)
    ;seg_const_tmp.eros_size = round(seg_const.eros_size*(seg_const.npssu-n)/seg_const.npssu )
    seg_const_tmp.eros_size = seg_const.eros_size-n

    
    
    CRD_t = {im_raw: CRD_in.im_raw, im_crr: imsel, hdr: CRD_in.hdr, mgnt_ar: CRD_in.mgnt_ar, mgnt_flx: CRD_in.mgnt_flx, Xar: CRD_in.Xar, Yar: CRD_in.Yar, Zar: CRD_in.Zar, Lath: CRD_in.Lath, Lonh: CRD_in.Lonh}

    
    ;New detection of positive and negative regions
;    amj_pnr_dt_core, CRD_t,  mdi_i, PNR_t, instr, seg_const = seg_const_tmp, pltn = (6*n), /All_ng, /detdisp;, /display, /info 
    amj_pnr_dt_core, CRD_t,  mdi_i, PNR_t, instr, seg_const = seg_const_tmp, pltn = (8*n);, /detdisp, /prt;, /display, /info 
    
    ;stop 
    ;removal of detected positive regions from the magnetogram
    if (PNR_t.num_p ge 1) then begin
      
        if ( n_regp eq 0 ) then begin
            pregions = PNR_t.pregions
        endif else begin
            pregions = [pregions,PNR_t.pregions]
        endelse
        
        n_regp = n_regp + PNR_t.num_p
        
        ;Removing pixels from magnetogram used for detecting regions with lower thresholds
        for i=0,PNR_t.num_p-1 do begin

            indx_pr = long(strsplit(PNR_t.pregions[i].indx,/extract))
            imsel[indx_pr] = 0.0
        
        endfor

    endif
    
    ;removal of detected negative regions from the magnetogram
    if (PNR_t.num_n ge 1) then begin
      
        if ( n_regn eq 0 ) then begin
            nregions = PNR_t.nregions
        endif else begin
            nregions = [nregions,PNR_t.nregions]
        endelse
        
        n_regn = n_regn + PNR_t.num_n
        
        ;Removing pixels from magnetogram used for detecting regions with lower thresholds
        for i=0,PNR_t.num_n-1 do begin
            indx_nr = long(strsplit(PNR_t.nregions[i].indx,/extract))
            imsel[indx_nr] = 0.0        
        endfor

    endif

endfor


; Creating reference image


im_mask = im*0.0
if (n_regp ge 1) then begin
    for i=0,n_regp-1 do begin
        inp = long(strsplit(pregions[i].indx,/extract))
        im_mask[inp] = 1 
    endfor
endif

if (n_regn ge 1) then begin
    for i=0,n_regn-1 do begin
        inn= long(strsplit(nregions[i].indx,/extract))
        im_mask[inn] = -1 
    endfor
endif

if keyword_set(detdisp) then begin
    window,1,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    loadct, 0, /silent
    ;tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
    tv,bytscl(congrid(im_mask,display_xsize,display_ysize),min=-1,max=1)
    loadct, 0, /silent

    if (n_regp ge 1) then begin
        for i=0,n_regp-1 do begin
            ;str='P'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom,str,charsize=1.5, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*display_zoom, pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom, color=255,/device
        endfor
    endif

    if (n_regn ge 1) then begin
        for i=0,n_regn-1 do begin
            ;str='N'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,nregions[i].fcenxp*display_zoom,nregions[i].fcenyp*display_zoom,str,charsize=1.5, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*display_zoom, nregions[i].fcenxp*display_zoom, nregions[i].fcenyp*display_zoom, color=0,/device
        endfor
    endif
    
endif

if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename='Mgntgram_reg_all.eps'
    device,/color,bits_per_pixel=8,/portr,/inches,xsize=6.0,ysize=6.0,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    !x.window=[0.0,1.0]
    !y.window=[0.0,1.0]
    px = !x.window * !d.x_vsize ;Get size of window in device units
    py = !y.window * !d.y_vsize
    
    loadct, 0, /silent
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(im_mask,display_xsize,display_ysize),min=-1,max=1)
    ;tv,bytscl(congrid(imgs0,print_xsize,print_ysize),min=display_thresholdd,max=display_thresholdu)
   
    if (n_regp gt 0) then begin
        for i=0,n_regp-1 do begin
            ;str='P'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*print_zoom/sz[1]*px[1], pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1], color=255, thick=3,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            ;str='N'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,nregions[i].fcenxp*print_zoom/sz[1]*px[1],nregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*print_zoom/sz[1]*px[1], nregions[i].fcenxp*print_zoom/sz[1]*px[1], nregions[i].fcenyp*print_zoom/sz[2]*py[1], color=0, thick=3,/device
        endfor
    endif
    
    device, /close
    set_plot,'X'
    loadct, 0, /silent           
   
endif


print, 'Enlarging positive and negative regions'

;-----------------------------------------------------------------------------------------------------------
;Second pass. Enlarging Detected PNRs-------------------------------------------------------------------------------------------
Th_max = max(abs(seg_const.valid_range))

;Decrement for lowering the detection threshold
d_thrs = (seg_const.ar_th - seg_const.ar_th2)/(seg_const.npssu)
;d_thrs = 0.0


for n=0,seg_const.npssu-1 do begin
  
    Th_min = seg_const.ar_th - d_thrs*(n+1)
;    Th_min = seg_const.ar_th2
;    print, Th_min 
    
    ;Positive regions
    if (n_regp gt 0) then begin
    
        Flx = pregions.flux
        
        for i=0,n_regp-1 do begin
        
            ;Finding the smallest region
            minrg = min(abs(Flx),minrgn,/nan)
            
            ;Extracting indices
            inp = long(strsplit(pregions[minrgn].indx,/extract))
            
            ;Adding field to the unused magnetogram in order to perform region growth
            imsel[inp] =   seg_const.ar_th
                    
            ;Enlargement
            ind_gn=idl_region_grow(imsel,inp,thresh=[Th_min,Th_max])
            
            ;Recalculating parameters----------------------------------------------------
            flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux

            ;Storing positive region
           
            ;Adding region indexes
            pregions[minrgn].indx = strjoin(string(ind_gn))                    
            
            ;Updating relevant quantities
            tmp_r = {flux: flux} 
                           
            prtmp = pregions[minrgn]                        
            STRUCT_ASSIGN, tmp_r, prtmp, /nozero
             
            pregions[minrgn] = prtmp  
                             
            ;Removing region from the unused magnetogram
            imsel[ind_gn] = 0.0
                    
            Flx[minrgn] = !values.f_nan
     
        endfor
     
    endif


    ;Negative regions
    if (n_regn gt 0) then begin
    
        Flx = nregions.flux
        
        for i=0,n_regn-1 do begin
        
            ;Finding the smallest region
            minrg = min(abs(Flx),minrgn,/nan)
            
            ;Extracting indices
            inp = long(strsplit(nregions[minrgn].indx,/extract))
            
            ;Adding field to the unused magnetogram in order to perform region growth
            imsel[inp] = - seg_const.ar_th
                    
            ;Enlargement                
            ind_gn=idl_region_grow(-imsel,inp,thresh=[Th_min,Th_max])
            
            ;Recalculating parameters----------------------------------------------------
            flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
            
            ;Storing positive region
           
            ;Adding region indexes
            nregions[minrgn].indx = strjoin(string(ind_gn))                    
            
            ;Updating relevant quantities
            tmp_r = {flux: flux} 
                           
            prtmp = nregions[minrgn]                        
            STRUCT_ASSIGN, tmp_r, prtmp, /nozero 
            nregions[minrgn] = prtmp  
                             
            ;Removing region from the unused magnetogram
            imsel[ind_gn] = 0.0
                    
            Flx[minrgn] = !values.f_nan
     
        endfor
     
    endif

endfor

;Calculating Region parameters----------------------------------------------------
;Positive regions
if (n_regp gt 0) then begin

    for i=0,n_regp-1 do begin
            
        ;Extracting indices
        ind_gn = long(strsplit(pregions[i].indx,/extract))
        
        ;Calculating parameters----------------------------------------------------
        flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
        area = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)         ;total area
        
        ;Geographical Center
        fcenx = total(CRD_in.Xar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        fceny = total(CRD_in.Yar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        fcenz = total(CRD_in.Zar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        
        fcenlt = atan(fcenz,sqrt( fcenx^2.0 + fceny^2.0 ) )/!dtor    ;Latitude center  
        fcenln = atan(fceny,fcenx)/!dtor                             ;Longitude center
        
        ;Center pixels
        lat1 = fcenlt*!dtor
        lon1 = fcenln*!dtor
        
        lat2 = CRD_in.Lath*!dtor
        lon2 = CRD_in.Lonh*!dtor
        
        dlat = abs(lat1-lat2)
        dlon = abs(lon1-lon2)
    
        dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
        
        mindis = min( dis_rg, minin, /nan )
        tmp_2d = array_indices(im,minin)
        fcenxp = tmp_2d[0,0];    ;Center in x
        fcenyp = tmp_2d[1,0];    ;Center in y
        
        ;Mean and Pixel Radius
        dcen = mean( dis_rg[ind_gn],/double, /nan)
        dcenfw = total( dis_rg[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        
        dcenp  = dcen/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor            
        dcenpfw  = dcenfw/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor
        
        dm = 0
        qm = 0
        
        ;Storing positive region
       
        ;Adding region indexes
        pregions[i].indx = strjoin(string(ind_gn))                    
        
        ;Updating relevant quantities
        tmp_r = {flux: flux, area: area, $
                  fcn_lt: fcenlt, fcn_ln: fcenln, dcen: dcen, $ 
                  fcenxp: fcenxp, fcenyp: fcenyp, dcenp: dcenp, dm: dm, qm: qm} 
                       
        prtmp = pregions[i]                        
        STRUCT_ASSIGN, tmp_r, prtmp, /nozero
         
        pregions[i] = prtmp  
                         
    endfor
 
endif


;Negative regions
if (n_regn gt 0) then begin
    
    for i=0,n_regn-1 do begin
    
        ;Extracting indices
        ind_gn = long(strsplit(nregions[i].indx,/extract))

        ;Calculating parameters----------------------------------------------------
        flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
        area = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)         ;total area
        
        ;Geographical Center
        fcenx = total(CRD_in.Xar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        fceny = total(CRD_in.Yar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        fcenz = total(CRD_in.Zar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        
        fcenlt = atan(fcenz,sqrt( fcenx^2.0 + fceny^2.0 ) )/!dtor    ;Latitude center  
        fcenln = atan(fceny,fcenx)/!dtor                             ;Longitude center
        
        ;Center pixels
        lat1 = fcenlt*!dtor
        lon1 = fcenln*!dtor
        
        lat2 = CRD_in.Lath*!dtor
        lon2 = CRD_in.Lonh*!dtor
        
        dlat = abs(lat1-lat2)
        dlon = abs(lon1-lon2)
    
        dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
        
        mindis = min( dis_rg, minin, /nan )
        tmp_2d = array_indices(im,minin)
        fcenxp = tmp_2d[0,0];    ;Center in x
        fcenyp = tmp_2d[1,0];    ;Center in y
        
        ;Mean and Pixel Radius
        dcen = mean( dis_rg[ind_gn],/double, /nan)
        dcenfw = total( dis_rg[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
        
        dcenp  = dcen/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor            
        dcenpfw  = dcenfw/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor
        
        dm = 0
        qm = 0
        
        ;Storing negative region
       
        ;Adding region indexes
        nregions[i].indx = strjoin(string(ind_gn))                    
        
        ;Updating relevant quantities
        tmp_r = {flux: flux, area: area, $
                  fcn_lt: fcenlt, fcn_ln: fcenln, dcen: dcen, $ 
                  fcenxp: fcenxp, fcenyp: fcenyp, dcenp: dcenp, dm: dm, qm: qm} 
                       
        prtmp = nregions[i]                        
        STRUCT_ASSIGN, tmp_r, prtmp, /nozero 
        nregions[i] = prtmp  
                         
    endfor
 
endif


im_mask = im*0.0
if (n_regp ge 1) then begin
    for i=0,n_regp-1 do begin
        inp = long(strsplit(pregions[i].indx,/extract))
        im_mask[inp] = 1 
    endfor
endif

if (n_regn ge 1) then begin
    for i=0,n_regn-1 do begin
        inn= long(strsplit(nregions[i].indx,/extract))
        im_mask[inn] = -1 
    endfor
endif


if keyword_set(detdisp) then begin
    window,2,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    loadct, 0, /silent
    ;tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
    tv,bytscl(congrid(im_mask,display_xsize,display_ysize),min=-1,max=1)
    loadct, 0, /silent

    if (n_regp ge 1) then begin
        for i=0,n_regp-1 do begin
            ;str='P'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom,str,charsize=1.5, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*display_zoom, pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom, color=255,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            ;str='N'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,nregions[i].fcenxp*display_zoom,nregions[i].fcenyp*display_zoom,str,charsize=1.5, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*display_zoom, nregions[i].fcenxp*display_zoom, nregions[i].fcenyp*display_zoom, color=0,/device
        endfor
    endif
    
endif

if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename='Mgntgram_reg_all_join.eps'
    device,/color,bits_per_pixel=8,/portr,/inches,xsize=6.0,ysize=6.0,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    !x.window=[0.0,1.0]
    !y.window=[0.0,1.0]
    px = !x.window * !d.x_vsize ;Get size of window in device units
    py = !y.window * !d.y_vsize
    
    loadct, 0, /silent
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(im_mask,display_xsize,display_ysize),min=-1,max=1)
    ;tv,bytscl(congrid(imgs0,print_xsize,print_ysize),min=display_thresholdd,max=display_thresholdu)
   
    if (n_regp gt 0) then begin
        for i=0,n_regp-1 do begin
            ;str='P'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*print_zoom/sz[1]*px[1], pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1], color=255, thick=3,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            ;str='N'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,nregions[i].fcenxp*print_zoom/sz[1]*px[1],nregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*print_zoom/sz[1]*px[1], nregions[i].fcenxp*print_zoom/sz[1]*px[1], nregions[i].fcenyp*print_zoom/sz[2]*py[1], color=0, thick=3,/device
        endfor
    endif
    
    device, /close
    set_plot,'X'
    loadct, 0, /silent           
   
endif

;stop

if not keyword_set(not_merge) then begin

    print, 'Merging regions'
    ;-----------------------------------------------------------------------------------------------------------
    ;Third pass. Merging PNRs-------------------------------------------------------------------------------------------
    
    ;Merging distance
    n_dis = seg_const.dis_lim
    
    ;Density limit
    dns_lim = 0.0
    
    ;Positive regions
    if (n_regp gt 1) then begin
      
        ;Creating Overlapping Matrix
        OvM = fltarr(n_regp,n_regp) ; Overlap matrix
        
        for i=0,n_regp-1 do begin    
      
            for j=0,n_regp-1 do begin    
                
                ;Ignoring auto overlap
                if (i eq j) then begin
                   OvM(i,j) = !values.f_nan
                endif else begin
                
                    ;Calculating distances between its centroids
                    lat1 = pregions[i].fcn_lt*!dtor
                    lon1 = pregions[i].fcn_ln*!dtor
                    
                    lat2 = pregions[j].fcn_lt*!dtor
                    lon2 = pregions[j].fcn_ln*!dtor
                    
                    dlat = abs(lat1-lat2)
                    dlon = abs(lon1-lon2)
            
                    dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
                                    
                    ;Finding regions withing prescribed distance
                    if ( dis_rg lt n_dis*(pregions[i].dcen + pregions[j].dcen) ) then begin
                      
                        r1 = n_dis*pregions[i].dcen 
                        r2 = n_dis*pregions[j].dcen 
    
                        ;If the region is contained inside another
                        if ( ( dis_rg + r1 ) lt r2 ) then begin
                          
                            OvM(i,j) = 1.0;
                          
                        endif
                        
                        ;If the region contains another
                        if ( ( dis_rg + r2 ) lt r1 ) then begin
                          
                            OvM(i,j) = (!pi*r2^2.0)/(!pi*r1^2.0);
                          
                        endif
    
                        ;Make sure neither region is contained inside the other                   
                        if ( ( dis_rg + r1 ge r2 ) and ( dis_rg + r2 ge r1 ) ) then begin
                                          
                            Aint =   r1^2.0*acos( ( dis_rg^2.0 + r1^2.0 - r2^2.0 )/(2.0*dis_rg*r1) ) $ 
                                   + r2^2.0*acos( ( dis_rg^2.0 - r1^2.0 + r2^2.0 )/(2.0*dis_rg*r2) ) $
                                   - 0.5*sqrt( ( - dis_rg + r1 + r2 )*( dis_rg - r1 + r2 )*( dis_rg + r1 - r2 )*( dis_rg + r1 + r2 ) )
                            
                            OvM(i,j) = Aint/(!pi*r1^2.0) 
    
                        endif
                        
                        
                        ;Comparing densities to prevent swallowing of dense regions by less dense regions
                        if( abs(pregions[j].flux/pregions[j].area) lt abs(pregions[i].flux/pregions[i].area)*dns_lim ) then begin
                          
                            OvM(i,j) = !values.f_nan
                          
                        end
                        
                      
                    endif else begin
                      
                       OvM(i,j) = !values.f_nan
                       
                    endelse 
                    
                endelse
          
            endfor        
      
        endfor
            
        ;Keeping those above aggregation overlap
        tmp_in = where(OvM lt seg_const.ovr_lim, n_below)
        if (n_below ne 0) then OvM[where(OvM lt seg_const.ovr_lim)] = !values.f_nan       
                      
        while (total(finite(OvM)) ne 0) do begin
          
            mdis = max( max(OvM,mt,dimension=2,/nan) , ii,/nan)
            mdis = max(OvM(ii,*),jj,/nan)
            
            
            ;Extracting indices of Growing region
            inpg = long(strsplit(pregions[jj].indx,/extract))
    
            ;Extracting indices of Killed region
            inpk = long(strsplit(pregions[ii].indx,/extract))
                            
            ;Enlargement
            ind_gn= [inpg, inpk]
            
            ;Recalculating parameters----------------------------------------------------
            flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
            area = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)         ;total area
            
            ;Geographical Center
            fcenx = total(CRD_in.Xar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            fceny = total(CRD_in.Yar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            fcenz = total(CRD_in.Zar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            
            fcenlt = atan(fcenz,sqrt( fcenx^2.0 + fceny^2.0 ) )/!dtor    ;Latitude center  
            fcenln = atan(fceny,fcenx)/!dtor                             ;Longitude center
            
            ;Center pixels
            lat1 = fcenlt*!dtor
            lon1 = fcenln*!dtor
            
            lat2 = CRD_in.Lath*!dtor
            lon2 = CRD_in.Lonh*!dtor
            
            dlat = abs(lat1-lat2)
            dlon = abs(lon1-lon2)
        
            dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
            
            mindis = min( dis_rg, minin, /nan )
            tmp_2d = array_indices(im,minin)
            fcenxp = tmp_2d[0,0];    ;Center in x
            fcenyp = tmp_2d[1,0];    ;Center in y
            
            ;Mean and Pixel Radius
            dcen = mean( dis_rg[ind_gn],/double, /nan)
            dcenfw = total( dis_rg[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            
            dcenp  = dcen/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor            
            dcenpfw  = dcenfw/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor
            
            dm = 0
            qm = 0
            
            ;Storing positive region
           
            ;Adding region indexes
            pregions[jj].indx = strjoin(string(ind_gn))                    
            
            ;Updating relevant quantities
            tmp_r = {flux: flux, area: area, $
                      fcn_lt: fcenlt, fcn_ln: fcenln, dcen: dcen, $ 
                      fcenxp: fcenxp, fcenyp: fcenyp, dcenp: dcenp, dm: dm, qm: qm} 
                           
            prtmp = pregions[jj]                        
            STRUCT_ASSIGN, tmp_r, prtmp, /nozero
             
            pregions[jj] = prtmp       
                             
            ;Removing pair from the matrix
            Keep_inx = findgen(n_regp)
            OvM = OvM[where(Keep_inx ne ii),*]
            OvM = OvM[*,where(Keep_inx ne ii)]
            
            pregions = pregions[where(Keep_inx ne ii)]
            n_regp = n_elements(pregions)
          
        endwhile
            
    endif
    
    
    
    ;Negative regions
    if (n_regn gt 1) then begin
      
        ;Creating Overlapping Matrix
        OvM = fltarr(n_regn,n_regn) ; Overlap matrix
        
        for i=0,n_regn-1 do begin    
      
            for j=0,n_regn-1 do begin    
                
                ;Ignoring auto overlap
                if (i eq j) then begin
                   OvM(i,j) = !values.f_nan
                endif else begin
                
                    ;Calculating distances between its centroids
                    lat1 = nregions[i].fcn_lt*!dtor
                    lon1 = nregions[i].fcn_ln*!dtor
                    
                    lat2 = nregions[j].fcn_lt*!dtor
                    lon2 = nregions[j].fcn_ln*!dtor
                    
                    dlat = abs(lat1-lat2)
                    dlon = abs(lon1-lon2)
            
                    dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
                                    
                    ;Finding regions withing prescribed distance
                    if ( dis_rg lt n_dis*(nregions[i].dcen + nregions[j].dcen) ) then begin
                      
                        r1 = n_dis*nregions[i].dcen 
                        r2 = n_dis*nregions[j].dcen 
    
                        ;If the region is contained inside another
                        if ( ( dis_rg + r1 ) lt r2 ) then begin
                          
                            OvM(i,j) = 1.0;
                          
                        endif
                        
                        ;If the region contains another
                        if ( ( dis_rg + r2 ) lt r1 ) then begin
                          
                            OvM(i,j) = (!pi*r2^2.0)/(!pi*r1^2.0);
                          
                        endif
    
                        ;Make sure neither region is contained inside the other                   
                        if ( ( dis_rg + r1 ge r2 ) and ( dis_rg + r2 ge r1 ) ) then begin
                                          
                            Aint =   r1^2.0*acos( ( dis_rg^2.0 + r1^2.0 - r2^2.0 )/(2.0*dis_rg*r1) ) $ 
                                   + r2^2.0*acos( ( dis_rg^2.0 - r1^2.0 + r2^2.0 )/(2.0*dis_rg*r2) ) $
                                   - 0.5*sqrt( ( - dis_rg + r1 + r2 )*( dis_rg - r1 + r2 )*( dis_rg + r1 - r2 )*( dis_rg + r1 + r2 ) )
                            
                            OvM(i,j) = Aint/(!pi*r1^2.0) 
    
                        endif
                        
                        ;Comparing densities to prevent swallowing of dense regions by less dense regions
                        if( abs(nregions[j].flux/nregions[j].area) lt abs(nregions[i].flux/nregions[i].area)*dns_lim ) then begin
                          
                            OvM(i,j) = !values.f_nan
                          
                        end                        
                        
                      
                    endif else begin
                      
                       OvM(i,j) = !values.f_nan
                       
                    endelse 
                    
                endelse
          
            endfor        
      
        endfor
            
        ;Keeping those above aggregation overlap
        
        
        tmp_in = where(OvM lt seg_const.ovr_lim, n_below)
        if (n_below ne 0) then OvM[where(OvM lt seg_const.ovr_lim)] = !values.f_nan               
        while (total(finite(OvM)) ne 0) do begin
          
            mdis = max( max(OvM,mt,dimension=2,/nan) , ii,/nan)
            mdis = max(OvM(ii,*),jj,/nan)
            
            
            ;Extracting indices of Growing region
            inng = long(strsplit(nregions[jj].indx,/extract))
    
            ;Extracting indices of Killed region
            innk = long(strsplit(nregions[ii].indx,/extract))
                            
            ;Enlargement
            ind_gn= [inng, innk]
            
            ;Recalculating parameters----------------------------------------------------
            flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
            area = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)         ;total area
            
            ;Geographical Center
            fcenx = total(CRD_in.Xar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            fceny = total(CRD_in.Yar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            fcenz = total(CRD_in.Zar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            
            fcenlt = atan(fcenz,sqrt( fcenx^2.0 + fceny^2.0 ) )/!dtor    ;Latitude center  
            fcenln = atan(fceny,fcenx)/!dtor                             ;Longitude center
            
            ;Center pixels
            lat1 = fcenlt*!dtor
            lon1 = fcenln*!dtor
            
            lat2 = CRD_in.Lath*!dtor
            lon2 = CRD_in.Lonh*!dtor
            
            dlat = abs(lat1-lat2)
            dlon = abs(lon1-lon2)
        
            dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
            
            mindis = min( dis_rg, minin, /nan )
            tmp_2d = array_indices(im,minin)
            fcenxp = tmp_2d[0,0];    ;Center in x
            fcenyp = tmp_2d[1,0];    ;Center in y
            
            ;Mean and Pixel Radius
            dcen = mean( dis_rg[ind_gn],/double, /nan)
            dcenfw = total( dis_rg[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
            
            dcenp  = dcen/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor            
            dcenpfw  = dcenfw/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor
            
            dm = 0
            qm = 0
            
            ;Storing Negative region
           
            ;Adding region indexes
            nregions[jj].indx = strjoin(string(ind_gn))                    
            
            ;Updating relevant quantities
            tmp_r = {flux: flux, area: area, $
                      fcn_lt: fcenlt, fcn_ln: fcenln, dcen: dcen, $ 
                      fcenxp: fcenxp, fcenyp: fcenyp, dcenp: dcenp, dm: dm, qm: qm} 
                           
            prtmp = nregions[jj]                        
            STRUCT_ASSIGN, tmp_r, prtmp, /nozero
             
            nregions[jj] = prtmp       
                             
            ;Removing pair from the matrix
            Keep_inx = findgen(n_regn)
            OvM = OvM[where(Keep_inx ne ii),*]
            OvM = OvM[*,where(Keep_inx ne ii)]
            
            nregions = nregions[where(Keep_inx ne ii)]
            n_regn = n_elements(nregions)
          
        endwhile
            
    endif
    
endif


  
;Labeling regions within the magnetogram
for i=0,n_regp-1 do begin  
    if keyword_set(pnr_lbl) then begin
        pregions[i].fr_lbl = pnr_lbl
    endif else begin
        pregions[i].fr_lbl = i+1
    endelse  
endfor

for i=0,n_regn-1 do begin  
    if keyword_set(pnr_lbl) then begin
        nregions[i].fr_lbl = pnr_lbl
    endif else begin
        nregions[i].fr_lbl = -(i+1)
    endelse  
endfor


;Storing positive and negative regions
if (n_regp gt 0) then begin
    PRs = [PRs, pregions]
endif

if (n_regn gt 0) then begin
    NRs = [NRs, nregions]
endif

;
;display the extracted regions
;

im_mask = im*0.0
if (n_regp ge 1) then begin
    for i=0,n_regp-1 do begin
        inp = long(strsplit(pregions[i].indx,/extract))
        im_mask[inp] = 1 
    endfor
endif

if (n_regn ge 1) then begin
    for i=0,n_regn-1 do begin
        inn= long(strsplit(nregions[i].indx,/extract))
        im_mask[inn] = -1 
    endfor
endif


if keyword_set(detdisp) then begin
    window,3,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    loadct, 0, /silent
    ;tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
    tv,bytscl(congrid(im_mask,display_xsize,display_ysize),min=-1,max=1)
    loadct, 0, /silent

    if (n_regp ge 1) then begin
        for i=0,n_regp-1 do begin
            ;str='P'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom,str,charsize=1.5, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*display_zoom, pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom, color=255,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            ;str='N'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,nregions[i].fcenxp*display_zoom,nregions[i].fcenyp*display_zoom,str,charsize=1.5, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*display_zoom, nregions[i].fcenxp*display_zoom, nregions[i].fcenyp*display_zoom, color=0,/device
        endfor
    endif
    
endif


if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename='Mgntgram_reg_all_mrg.eps'
    device,/color,bits_per_pixel=8,/portr,/inches,xsize=6.0,ysize=6.0,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    !x.window=[0.0,1.0]
    !y.window=[0.0,1.0]
    px = !x.window * !d.x_vsize ;Get size of window in device units
    py = !y.window * !d.y_vsize
    
    loadct, 0, /silent
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(im_mask,display_xsize,display_ysize),min=-1,max=1)
    ;tv,bytscl(congrid(imgs0,print_xsize,print_ysize),min=display_thresholdd,max=display_thresholdu)
   
    if (n_regp gt 0) then begin
        for i=0,n_regp-1 do begin
            ;str='P'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*print_zoom/sz[1]*px[1], pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1], color=255, thick=3,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            ;str='N'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,nregions[i].fcenxp*print_zoom/sz[1]*px[1],nregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*print_zoom/sz[1]*px[1], nregions[i].fcenxp*print_zoom/sz[1]*px[1], nregions[i].fcenyp*print_zoom/sz[2]*py[1], color=0, thick=3,/device
        endfor
    endif
    
    device, /close
    set_plot,'X'
    loadct, 0, /silent           
   
endif


if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename='Mgntgram_reg_all_mgn.eps'
    device,/color,bits_per_pixel=8,/portr,/inches,xsize=6.0,ysize=6.0,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    !x.window=[0.0,1.0]
    !y.window=[0.0,1.0]
    px = !x.window * !d.x_vsize ;Get size of window in device units
    py = !y.window * !d.y_vsize
    
    loadct, 0, /silent
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    ;tv,bytscl(congrid(im_mask,display_xsize,display_ysize),min=-1,max=1)
    tv,bytscl(congrid(CRD_in.im_crr,print_xsize,print_ysize),min=display_thresholdd,max=display_thresholdu)
   
    if (n_regp gt 0) then begin
        for i=0,n_regp-1 do begin
            ;str='P'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*print_zoom/sz[1]*px[1], pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1], color=255, thick=3,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            ;str='N'+strtrim(string(i+1,'(I2)'),2)
            ;xyouts,nregions[i].fcenxp*print_zoom/sz[1]*px[1],nregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*print_zoom/sz[1]*px[1], nregions[i].fcenxp*print_zoom/sz[1]*px[1], nregions[i].fcenyp*print_zoom/sz[2]*py[1], color=0, thick=3,/device
        endfor
    endif
    
    device, /close
    set_plot,'X'
    loadct, 0, /silent           
   
endif

;stop
;print, junk

return
end
