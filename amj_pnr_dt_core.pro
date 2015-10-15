pro amj_pnr_dt_core, CRD_in, mdi_i, PNR_out, instr, seg_const=seg_const, display=display, detdisp=detdisp, pltn = pltn, prt=prt, info=info, All_ng = All_ng
;+
; NAME:
;       amj_pnr_dt_core
;
; PURPOSE:
;       For a given MDI magnetogram, determine the areas of interest
;       using thresholding segmentation
;       morphological opening, morphological closing, morphological
;       distance, region growth techniques;   
; CALLING SEQUENCE:
;       amj_pnr_dt_core, CRD_in, mdi_i, PNR_out, instr
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
;        mdi_i: Mission day used to tag regions
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
;         pltn: plotting offset used so that detections performed for different days don't
;               overwrite each other 
;   
;   
; OUTPUTS: 
;       PNR_out: a structure holding the detected regions
;           1. seg_const: used parameters in the region detection (see above)
;           2. mask: mask holding the position of each detected region
;           3. num_p: Number of positive regions
;           4. pregions: structure holding details from all detected positive regions 
;              (see below for details on the rgn structure) 
;           5. num_n: Number of negative regions
;           6. nregions: rgn structure holding details from all detected negative regions 
;              (see below for details on the rgn structure)
;           7. r_sw: Switch that indicates if there is at least one region of each sign
;
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
;       All_ng:  Also uses diagonal neighbors in region growth
;       display: plots detected regions
;       detdisp: plots detected regions and all intermediate steps
;       prt: saves detected regions as an EPS figure
;       info: prints process information on screen
;       
; MODIFICATION HISTORY:
;   2012/07/05:  Andres Munoz-Jaramillo:   initiated.  Based on an algorithm by Jie Zhang
;   2014/05/31:  Andres Munoz-Jaramillo:   Adapted for human interface
;   2015/06/09:  Andres Munoz-Jaramillo:   Adapted for KPVT data
;   2015/09/16:  Andres Munoz-Jaramillo:   Adapted to work with different instruments
;-


;define the usage
if n_params() lt 4 then begin
    print,'USAGE:  amj_pnr_dt_core, CRD_in, mdi_i, PNR_out, instr'
    return
endif

;get the input image and dimension
im=CRD_in.im_crr
sz=size(im)




;
;define the segmentation control structure and the default value-----------------------------
;
seg_const_def={seg_const, ker_th:150.0, ker_th2:60.0, ar_th:60.0, ar_th2:40.0, qr_th:60.0, dila_size:4, eros_size:6.0, npssu: 4, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}


if not keyword_set(seg_const) then begin
    seg_const=seg_const_def
endif

if not keyword_set(pltn) then begin
    pltn=0;
endif

print, 'Ker_th', seg_const.ker_th, 'Ar_th', seg_const.ar_th


;
;set the display window size and threshold----------------------------------------------------------
display_zoom=0.5
display_xsize=sz[1]*display_zoom & display_ysize=sz[2]*display_zoom  

print_zoom=1.0
print_xsize=sz[1]*print_zoom & print_ysize=sz[2]*print_zoom 

display_thresholdu = 1000.00
display_thresholdd = - display_thresholdu       



;Define x and y scales---------------------------------------------------------------------------------------

;HMI uses structures for header values
if instr eq 4 then begin
	date = CRD_in.hdr.DATE_OBS

	;Define center and radius
	hfx = CRD_in.hdr.CRPIX1 ;  Location of the center in x pixels 
	hfy = CRD_in.hdr.CRPIX2 ;    Location of the center in y pixels
	di = CRD_in.hdr.RSUN_OBS/CRD_in.hdr.CDELT1;

	;Load Solar Coordinates
	P0 = 0.0
	RD = CRD_in.hdr.DSUN_OBS/CRD_in.hdr.RSUN_REF
	B0 = CRD_in.hdr.CRLT_OBS
	L0 = CRD_in.hdr.CRLN_OBS

	;Observer Coordinates
	X_scl = CRD_in.hdr.CDELT1/60.0
	Y_scl = CRD_in.hdr.CDELT2/60.0
	
endif else begin

	date = fxpar(CRD_in.hdr, 'DATE_OBS')

	;KPVT-512
	if instr eq 1 then begin

		;Define center and radius
		hfx = fxpar(CRD_in.hdr, 'CRPIX1A');    Location of the center in x pixels 
		hfy = fxpar(CRD_in.hdr, 'CRPIX2A');    Location of the center in y pixels
		di =  fxpar(CRD_in.hdr,'EPH_R0');

		;Observer Coordinates
		dx = fxpar(CRD_in.hdr, 'CDELT1')*fxpar(CRD_in.hdr, 'CRR_SCLX')/60.0
		dy = fxpar(CRD_in.hdr, 'CDELT2')*fxpar(CRD_in.hdr, 'CRR_SCLY')/60.0

	endif

	;MDI
	if instr eq 3 then begin

		;Define center and radius
		hfx = fxpar(CRD_in.hdr, 'X0');  Location of the center in x pixels 
		hfy = fxpar(CRD_in.hdr, 'Y0');  Location of the center in y pixels
		di = fxpar(CRD_in.hdr,'R_SUN');

		;Observer Coordinates
		dx = fxpar( CRD_in.hdr, 'XSCALE' )
		dy = fxpar( CRD_in.hdr, 'YSCALE' )

	endif	
	
endelse

Xg = transpose(double(floor(findgen(sz[1],sz[1])/sz[1]))) - hfx 
Yg = double(floor(findgen(sz[1],sz[1])/sz[1])) - hfy 
R = sqrt( Xg^2 + Yg^2 )




imgs0=im
im = CRD_in.im_crr

im[where( (R gt di*sin(seg_const.deg_lim*!dtor)) and finite(im) )] = 0.0


;
;determine the intensity threshold for segmentation 
g_th_min=seg_const.ar_th
g_th_max=max(abs(seg_const.valid_range))


;display Image------------------------------------------------
if keyword_set(detdisp) then begin
    window,2+pltn,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
endif



;
;choose kernel pixels
k_th=seg_const.ker_th

;kernel pixle index
;Positive polarity
imp=im
ind_kp=where(imp ge k_th,n_kp); le min(abs(seg_const.valid_range)),n_kp)
;Negative polarity
imn=im
ind_kn=where(-imn ge k_th,n_kn); and -imn le min(abs(seg_const.valid_range)),n_kn)

if ( n_kp eq 0 or n_kp eq 0 ) then begin
    if keyword_set(info) then print,'no kernel pixels found, aborted'
endif


;display kernel pixels
;Positive polarity
im_kp=fltarr(sz[1],sz[2])
if n_kp gt 0 then im_kp(ind_kp)=1.0
;Negative polarity
im_kn=fltarr(sz[1],sz[2])
if n_kn gt 0 then im_kn(ind_kn)=1.0

;
;Kernel Display------------------------------------------------

if keyword_set(detdisp) then begin
    window,0+pltn,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(im_kp-im_kn,display_xsize,display_ysize))
endif


;
;To obtain the valid active region kernels
;
;
;morphological opening operation: remove small kernel patches, e.g.,
;small bipoles; done by erosion operation first followed by dilation
;image erosion operation: shrink the features by erosion size 
;        effectively remove features whose size is smaller than
;        erosion size
;image dilation operation: grow the feature by growth size; 
;        effectively grow back to the orignal size (if eroded first);
;        the image may not be fully recovered
;        
;;morphological closing operation: Merge Kernels.
;image dilation operation: grow the feature by dilation size;
;       merging Kernels
;image erosion operation: shrink the features by dilation size 
;        effectively grow back to the orignal size
;--------------------------------------------
eros_pixel=seg_const.eros_size
r=floor(eros_pixel/2)
s=SHIFT(DIST(2*r+1,2*r+1), r, r) LE r

dila_size=seg_const.dila_size
r=floor(dila_size/2)
s2=SHIFT(DIST(2*r+1,2*r+1), r, r) LE r

;Positive Polarity
im_op  = erode(im_kp,s)
im_op  = dilate(im_op,s)   ;image after morphological opening operation
im_op  = dilate(im_op,s2)  ;Second dilation to merge cores
im_op  = erode(im_op,s2)   ;image after morphological closing operation
ind_op = where(im_op gt 0,n_op)

;Negative Polarity
im_on  = erode(im_kn,s)
im_on  = dilate(im_on,s)   ;image after morphological opening operation
im_on  = dilate(im_on,s2)  ;Second dilation to merge cores
im_on  = erode(im_on,s2)   ;image after morphological closing operation
ind_on = where(im_on gt 0,n_on)


if keyword_set(info) then print,';IMAGE: percentage of kernel pixels after morphological opening and closing operation=',string(float(n_op+n_on)/n_valid*100,'(F7.3)')+"%"

if ( n_op eq 0 or n_on eq 0 ) then begin
    regions={num:0}
    if keyword_set(info) then print,'No kernel pixel after erosion' 
endif

;display kernel pixels after erosion and merging------------------------------------------------
if keyword_set(detdisp) then begin
    window,1+pltn,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(double(im_op)-double(im_on),display_xsize,display_ysize))
endif

;---------------------------------------------------------------------------------------------------------
;number of discrete regions
;by iteratively extract the largest closed region from the closed image
;
n_regp=0  ;initiate the variable holding number of positive regions
n_regn=0  ;initiate the variable holding number of negative regions

;calculate a distance map, the distance also indicates the size of the region

;Image for removing visited pixels
imsel = imgs0

;Positive Polarity
im_dp=morph_distance(im_op,neighbor=2)
d_maxp=max(im_dp)       ; maximum distance for the collected pixels

if keyword_set(info) then print,';IMAGE: largest positive distance= ',string(d_maxp,'(I3)') 
srg_d_maxp=where(im_dp eq d_maxp,n_d_maxp)

;iteratively extract the regions from the segmented image
while (d_maxp gt 0) do begin

    ;indices of grown region
    srg_lrp_raw = idl_region_grow(im_dp,srg_d_maxp[0],/All,thresh=[1,d_maxp])
    index = where(im_op[srg_lrp_raw] gt 0. and im[srg_lrp_raw] gt 0., nvldin)
    
    if (nvldin gt 0) then begin
        srg_lrp = srg_lrp_raw[index]
        
        fluxp  = total(CRD_in.mgnt_flx[srg_lrp],/double, /nan)       ;total flux
            
        ; Store region indexes in a string, retrieve via:
        ; indx = long(strsplit(pregion.indx,/extract))
        
        pregion = {rgn,  lnk_sw: 0, mdi_i: long(mdi_i), ar_lbl: 0, fr_lbl: 0L, date: date, indx:strjoin(string(srg_lrp)), flux:fluxp, area:0.0, fcn_lt:0.0, fcn_ln:0.0, dcen:0.0, fcenxp: 0.0, fcenyp: 0.0, dcenp:0.0, dm:0.0, qm:0.0}
        if (n_regp eq 0) then begin
            pregions = pregion         
        endif else begin
            pregions = [pregions, pregion]
        endelse
        n_regp=n_regp+1        
    
        ;display the largest region
        ;    if keyword_set(display) then begin
        ;        if not keyword_set(ps) then window,5,xsize=display_xsize,ysize=display_ysize,retain=2
        ;        im_lrp=fltarr(sz[1],sz[2])
        ;        im_lrp[srg_lrp]=1
        ;        plot,[1,1],/nodata,xstyle=5,ystyle=5
        ;        tv,bytscl(congrid(im_lrp,display_xsize,display_ysize))
        ;    endif
    
        imsel[srg_lrp] = 0.0

    endif

    ;move to the next largest region
    im_dp(srg_lrp_raw)=0.0
    d_maxp=max(im_dp)   ; maximum distance from the centroid to the edge
    if keyword_set(info) then print,';IMAGE: largest positive distance= ',string(d_maxp,'(I3)') 
    srg_d_maxp=where(im_dp eq d_maxp,n_dp)
endwhile

if (n_regp eq 0) then begin
    pregions = {rgn,  lnk_sw: 0,  mdi_i: long(mdi_i), ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm:!values.f_nan, qm:!values.f_nan}
endif


;Negative Polarity
im_dn=morph_distance(im_on,neighbor=2)
d_maxn=max(im_dn)       ; maximum distance for the collected pixels

if keyword_set(info) then print,';IMAGE: largest negative distance= ',string(d_maxn,'(I3)') 
srg_d_maxn=where(im_dn eq d_maxn,n_d_maxn)

;iteratively extract the regions from the segmented image
while (d_maxn gt 0) do begin
  
    ;indices of grown region
    srg_lrn_raw = idl_region_grow(im_dn,srg_d_maxn[0],/All,thresh=[1,d_maxn])
    index = where(im_on[srg_lrn_raw] gt 0. and im[srg_lrn_raw] lt 0.,nvldin)

    if (nvldin gt 0) then begin

        srg_lrn = srg_lrn_raw[index]
        
        fluxn  = total(CRD_in.mgnt_flx[srg_lrn],/double, /nan)       ;total flux    
        
        ; Store region indexes in a string, retrieve via:
        ; indx = long(strsplit(nregion.indx,/extract))
        
        nregion = {rgn,  lnk_sw: 0,  mdi_i: long(mdi_i), ar_lbl: 0, fr_lbl: 0L, date: date, indx:strjoin(string(srg_lrn)), flux:fluxn, area:0.0, fcn_lt:0.0, fcn_ln:0.0, dcen:0.0, fcenxp: 0.0, fcenyp: 0.0, dcenp:0.0, dm:0.0, qm:0.0}
        if (n_regn eq 0) then begin  ;initial assignment        
            nregions = nregion         
        endif else begin
            nregions = [nregions,nregion]             
        endelse
        n_regn=n_regn+1    
    
        ;display the largest region
        ;    if keyword_set(display) then begin
        ;        if not keyword_set(ps) then window,6,xsize=display_xsize,ysize=display_ysize,retain=2
        ;        im_lrn=fltarr(sz[1],sz[2])
        ;        im_lrn[srg_lrn]=1
        ;        plot,[1,1],/nodata,xstyle=5,ystyle=5
        ;        tv,bytscl(congrid(im_lrn,display_xsize,display_ysize))
        ;    endif
    
        imsel[srg_lrn] = 0.0

    endif
    
;move to the next largest region
    im_dn(srg_lrn_raw)=0.0
    d_maxn=max(im_dn)   ; maximum distance from the centroid to the edge
    if keyword_set(info) then print,';IMAGE: largest negative distance= ',string(d_maxn,'(I3)') 
    srg_d_maxn=where(im_dn eq d_maxn,n_dn)
endwhile

if (n_regn eq 0) then begin
    nregions = {rgn,  lnk_sw: 0,  mdi_i: long(mdi_i), ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm:!values.f_nan, qm:!values.f_nan}
endif

if keyword_set(info) then print,'IMAGE: number of positive regions extracted= ',string(n_regp,'(I3)') 
if keyword_set(info) then print,'IMAGE: number of negative regions extracted= ',string(n_regn,'(I3)') 


;To obtain the full size of the magnetic regions
;region growth operation
;recover the region interested: growing from the AR kernel pixels to all
;connected pixels above the threshold seg_const.ar_th


;-----------------------------------------------------------------------------------------------------------
;Second pass. Enlarging Detected PNRs-------------------------------------------------------------------------------------------
Th_max = max(abs(seg_const.valid_range))

;Decrement for lowering the detection threshold
d_thrs = (2.0*seg_const.ker_th - seg_const.ar_th)/(2.0*seg_const.npssu-1)
;d_thrs = 0.0


for n=0,2.0*seg_const.npssu-1 do begin
  
    Th_min = 2.0*seg_const.ker_th - d_thrs*(n)
;    print, Th_min
    
    ;Positive regions
    if (n_regp gt 0) then begin
      
        Flx = pregions.flux      
            
        for i = n_regp-1,0,-1 do begin

            ;Finding the smallest region
            minrg = min(abs(Flx),minrgn,/nan)
        
            ;Extracting indices
            inp = long(strsplit(pregions[minrgn].indx,/extract))
            
            ;Adding field to the unused magnetogram in order to perform region growth
            imsel[inp] =   seg_const.ker_th
                    
            ;Enlargement
            ind_gn=idl_region_grow(imsel,inp,thresh=[Th_min,Th_max])
                        
            ;Storing positive region                        
            if ind_gn[0] ne -1 then begin
              
                flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
                
                ;Updating relevant quantities
                tmp_r = {flux: flux} 
                               
                prtmp = pregions[minrgn]                        
                STRUCT_ASSIGN, tmp_r, prtmp, /nozero
                 
                pregions[minrgn] = prtmp                  
              
                ;Adding region indexes
                pregions[minrgn].indx = strjoin(string(ind_gn))                    
                                             
                ;Removing region from the unused magnetogram
                imsel[ind_gn] = 0.0
                Flx[minrgn] = !values.f_nan
                            
            endif
                         
        endfor
     
    endif


    ;Negative regions
    if (n_regn gt 0) then begin
    
        Flx = nregions.flux

        for i = n_regn-1,0,-1 do begin
        
            ;Finding the smallest region
            minrg = min(abs(Flx),minrgn,/nan)

            ;Extracting indices
            inp = long(strsplit(nregions[minrgn].indx,/extract))
            
            ;Adding field to the unused magnetogram in order to perform region growth
            imsel[inp] = - seg_const.ker_th
                    
            ;Enlargement                
            ind_gn=idl_region_grow(-imsel,inp,thresh=[Th_min,Th_max])
            
            ;Storing negative region
           
            if ind_gn[0] ne -1 then begin
              
                flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
              
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
                          
            endif  
                         
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





;
;display the extracted regions
;
if (keyword_set(display) or keyword_set(detdisp)) then begin
    window,5+pltn,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    loadct, 0, /silent
    tv,bytscl(congrid(imgs0,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
    loadct, 0, /silent

    if (n_regp gt 0) then begin
        for i=0,n_regp-1 do begin
            str='P'+strtrim(string(i+1,'(I2)'),2)
            xyouts,pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom,str,charsize=1.5, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*display_zoom, pregions[i].fcenxp*display_zoom,pregions[i].fcenyp*display_zoom, color=255,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            str='N'+strtrim(string(i+1,'(I2)'),2)
            xyouts,nregions[i].fcenxp*display_zoom,nregions[i].fcenyp*display_zoom,str,charsize=1.5, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*display_zoom, nregions[i].fcenxp*display_zoom, nregions[i].fcenyp*display_zoom, color=0,/device
        endfor
    endif
    
endif

if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename='Mgntgram_reg.eps'
    device,/color,bits_per_pixel=8,/portr,/inches,xsize=6.0,ysize=6.0,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    !x.window=[0.0,1.0]
    !y.window=[0.0,1.0]
    px = !x.window * !d.x_vsize ;Get size of window in device units
    py = !y.window * !d.y_vsize
    
    loadct, 0, /silent
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(imgs0,print_xsize,print_ysize),min=display_thresholdd,max=display_thresholdu)
   
    if (n_regp gt 0) then begin
        for i=0,n_regp-1 do begin
            str='P'+strtrim(string(i+1,'(I2)'),2)
            xyouts,pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=255,/device
            tvcircle,2.0*pregions[i].dcenp*print_zoom/sz[1]*px[1], pregions[i].fcenxp*print_zoom/sz[1]*px[1],pregions[i].fcenyp*print_zoom/sz[2]*py[1], color=255, thick=3,/device
        endfor
    endif

    if (n_regn gt 0) then begin
        for i=0,n_regn-1 do begin
            str='N'+strtrim(string(i+1,'(I2)'),2)
            xyouts,nregions[i].fcenxp*print_zoom/sz[1]*px[1],nregions[i].fcenyp*print_zoom/sz[2]*py[1],str,charsize=2, color=0,/device
            tvcircle,2.0*nregions[i].dcenp*print_zoom/sz[1]*px[1], nregions[i].fcenxp*print_zoom/sz[1]*px[1], nregions[i].fcenyp*print_zoom/sz[2]*py[1], color=0, thick=3,/device
        endfor
    endif
    
    device, /close
    set_plot,'X'
    loadct, 0, /silent           
   
endif

if ( (n_regp gt 0) and (n_regn gt 0) ) then begin
    r_sw = 1
endif else begin
    r_sw = 0
endelse

PNR_out = {seg_const:seg_const, hdr:CRD_in.hdr, num_p: n_regp, pregions:pregions, num_n: n_regn, nregions:nregions, r_sw:r_sw}

return
end
