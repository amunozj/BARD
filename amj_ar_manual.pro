pro amj_ar_manual, CRD_in, mdi_i, inp, inn, lbl, PRs, NRs, ARs, date, ar_cnst=ar_cnst, seg_const=seg_const, track = track, display=display, prt=prt, info=info, sqs_nm = sqs_nm, pltn = pltn
;+
; NAME:
;      amj_ar_dt_track_dr
;
; PURPOSE:
;       Agregates a given set of positive and negative regions into active regions (ARs)
;       and calculate physical AR parameters
;       
;
; CALLING SEQUENCE:
;       amj_ar_dt_track_dr, CRD_in, PNR_p, PNR_c, File, lbl, AR_p, AR_c, AR_lb_Out
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
;		 inp : Index of the positive region that will be used to create the BMR
;
;		 inn : Index of the negative region that will be used to create the BMR
;
;        lbl:   Number to label each individual active region
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
;       File:  Path to the fits file containing the MDI magnetogram
;           
;       rgn structure:
;           1. num: region index
;           2. indx: indexes of all region pixels. Retrieve via long(strsplit(rgn.indx,/extract))
;           3. flux: magnetic flux (Mx)
;           4. area: region area (cm^2)
;           5. fcn_lt: latitude of the flux weighted geographical center    
;           6. fcn_ln: longitude of the flux weighted geographical center   
;           7. dcen: average distance to geographical center (not flux weighted)
;           8. fcenxp: x pixel coordinate of geographical center 
;           9. fcenyp: y pixel coordinate of geographical center 
;          10. dcenp: average distance to geographical center in pixels ;           
;           
;		date: Date of observation used to tag each AR.
;
;
; OPTIONAL INPUTS:
;       ar_cnst: struct holding control parameters for AR pairing
;            1. dis_lim1: multiplier which defines the maximum
;               distance used for pairing up regions in the first pass;
;            2. dis_lim2: multiplier which defines the maximum
;               distance used for pairing up regions in the second pass;
;            3. exp_f: Flux weighting exponent for distance function   
;            4. exp_d: Distance weighting exponent for distance function   
;            5. exp_s: Size weighting exponent for distance function   
;
;
; OUTPUTS: 
;       ARs: a structure holding the detected active regions
;            1. seg_const: used parameters in the positive and negative region detection (see amj_pnr_dt.pro)
;            2. ar_cnst: struct holding control parameters for AR pairing (see above) 
;            3. hdr: heather of the MDI magnetogram
;            4. num_ar: Number of detected active regions
;            5. ARs: ar structure holding details from all detected active regions 
;              (see below for details on the ar structure)
;            6. ar_sw: Switch that indicates if there is at least one active region detected
;             
;       ar structure: 
;            1. labl:  Unique labl to identify each active region
;            2. clr: plotting color
;            3. num: ar index
;            4. pregions: index of positive regions in the AR
;            5. nregions: index of positive regions in the AR
;            6. indxp: indexes of all positive pixels in the AR. Retrieve via long(strsplit(ar.indxp,/extract))   
;            7. indxn: indexes of all negative pixels in the AR. Retrieve via long(strsplit(ar.indxn,/extract))   
;            8. fluxp: positive flux in AR (Mx)
;            9. fluxn: negative flux in AR (Mx)
;           10. areap: area of positive part of the AR (cm^2)
;           11. arean: area of negative part of the AR (cm^2)
;           12. fcn_ltp: latitude of the flux weighted geographical center (deg) (positive part)
;           13. fcn_lnp: longitude of the flux weighted geographical center (deg) (positive part)
;           14. dcenp: average distance to geographical center (deg) (positive part; not flux weighted)
;           15. fcn_ltn: latitude of the flux weighted geographical center (deg) (negative part)
;           16. fcn_lnn: longitude of the flux weighted geographical center (deg) (negative part)
;           17. dcenn: average distance to geographical center (deg) (negative part; not flux weighted)
;           18. fcenxpp: x pixel coordinate of geographical center (positive part)
;           19. fcenypp: y pixel coordinate of geographical center (positive part)
;           20. dcenpp: average distance to geographical center in pixels (positive part)
;           21. dis: Distance between each polarities's centroids (deg)
;           22. tilt: AR tilt (deg)
;           23. lp: sign of the leading polarity
;  
;
; KEYWORDS:
;       display: plots intermediate steps
;       prt: saves intermediate steps to EPS figures
;       info: prints process information on screen
;
; 
; MODIFICATION HISTORY:
;   2012/07/11:  Andres Munoz-Jaramillo:   initiated.
;   2015/09/16:  Andres Munoz-Jaramillo:   Adapted to work with different instruments
;-


;define the usage
;if n_params() lt 6 then begin
;    print,'USAGE: amj_ar_dt_track_dr, CRD_in, File, lbl, AR_p, AR_c, AR_lb_Out'
;    return
;endif


;get the input image and dimension
im=CRD_in.im_raw
imsel=CRD_in.im_crr
sz=size(im)

;
;define the segmentation control structure and the default value
;

;Playing around
ar_cnst_def={dis_lim1:6.0, dis_lim2:3.0, exp_f: 1.0,  exp_d: 4.0, exp_s: 0.5, mlth: 40.0, mxB: 180.0, MxFlxim:3.0, Imb_tol: 0.10, Imb_it: 10, lim_lon: 0.0, k_sig:15.0, npr: 5, nmgnt: 5 , vld_thr: 0.69, valid_range:[-20000.,20000.]}
seg_const_def={seg_const, ker_th:150.0, ker_th2:60.0, ar_th:60.0, ar_th2:40.0, qr_th:60.0, dila_size:4, eros_size:6.0, npssu: 4, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}

if not keyword_set(ar_cnst) then begin
    ar_cnst=ar_cnst_def
endif
if not keyword_set(seg_const) then begin
    seg_const=seg_const_def
endif

if not keyword_set(pltn) then begin
    pltn=0;
endif
;
;set the display window size
;
display_zoom=0.5
display_xsize=sz[1]*display_zoom & display_ysize=sz[2]*display_zoom  

print_zoom=1.0
print_xsize=sz[1]*print_zoom & print_ysize=sz[2]*print_zoom 


;
;image basic information-------------------------------------------------------------
;
ind_valid=where(im ge ar_cnst.valid_range[0] and im le ar_cnst.valid_range[1],n_valid,com=ind_invalid,ncom=n_invalid)
mom=moment(im(ind_valid))
im_mean=mom[0]
im_sigma = sqrt(mom[1]) ; the standard deviation of the whole effective image
im_med=median(im(ind_valid))
display_thresholdu = im_mean+im_sigma*ar_cnst.k_sig
display_thresholdd = - display_thresholdu       



;Defining colors to maximise distance
clrsqs = [0,255,127]

ss = size(clrsqs)
ss = ss[1]-1
n = 2; 
while ss lt 15 do begin
    
    tmp = round(clrsqs - 128/n)
    tmp = tmp[1:ss]

    for i = 0,ss-1 do begin    
        mxd = where( abs(clrsqs[ss+i]-tmp) eq max(abs(clrsqs[ss+i]-tmp)) )     
        clrsqs = [clrsqs, tmp[mxd[0]]]
        
        first = 0
        last = ss-i-1
        
        if last gt 1 then begin
            case mxd[0] of
                first: tmp = tmp[1:*]
                last:  tmp = tmp[first:last-1]
                else:  tmp = [ tmp[first:mxd[0]-1], tmp[mxd[0]+1:last] ]
            endcase
        endif

    endfor
    n = n*2;
    ss = size(clrsqs)
    ss = ss[1]-1
  
endwhile




;-----------------------------------------------------------------------------------------------------------
;First pass. Associating the regions provided by the user--------------------------------------------------------



;Making Sure all visible longitudes are in the same hemisphere
;---------------------------------------------------------

Lonp = PRs[inp].fcn_ln
Lonn = NRs[inn].fcn_ln
        
sun_data = get_sun(date,carr=carr,he_lon=he_lon)
if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
        
  tmpin = where(Lonp lt 0, n_tmp)
  if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
        
  tmpin = where(Lonn lt 0, n_tmp)
  if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
        
endif 


;Finding leading polarity
if ( Lonp ge Lonn ) then begin

    LatL  = PRs[inp].fcn_lt*!dtor
    LongL = Lonp*!dtor

    LatF  = NRs[inn].fcn_lt*!dtor
    LongF = Lonn*!dtor
    
    lp = 1;

endif else begin

    LatF  = PRs[inp].fcn_lt*!dtor
    LongF = Lonp*!dtor

    LatL  = NRs[inn].fcn_lt*!dtor
    LongL = Lonn*!dtor
    
    lp = -1;

endelse

dlat = abs(LatF-LatL)
dlon = abs(LongF-LongL)

Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt

dm = 0
qm = 0

arin  = where(ARs.mdi_i eq mdi_i, nars)
; Store positive and region number in a string, retrieve via:
; nmbr = long(strsplit(ar.pregions,/extract)) or nmbr = long(strsplit(ar.nregions,/extract))
Nwar = {ar, mdi_i: mdi_i, date: date, labl: lbl, clr: clrsqs[lbl mod 17], indxp: PRs[inp].indx, indxn: NRs[inn].indx, fluxp: PRs[inp].flux, fluxn: NRs[inn].flux, areap:PRs[inp].area, arean:NRs[inn].area, $
          fcn_ltp:PRs[inp].fcn_lt, fcn_lnp:PRs[inp].fcn_ln, dcenp:PRs[inp].dcen, $ 
          fcn_ltn:NRs[inn].fcn_lt, fcn_lnn:NRs[inn].fcn_ln, dcenn:NRs[inn].dcen, $
          fcenxpp: PRs[inp].fcenxp, fcenypp: PRs[inp].fcenyp, dcenpp:PRs[inp].dcenp, $
          fcenxpn: NRs[inn].fcenxp, fcenypn: NRs[inn].fcenyp, dcenpn:NRs[inn].dcenp, $               
          dis: Dis, tilt: Tilt, lp: lp, dm: dm, qm: qm}                 

if keyword_set(info) then print,';AR No. ',string(n_ars,'(I3)'),' combines regions P', string(n_inx[inn],'(I3)'), ' and N', string(n_inx[inn],'(I3)')
    
;Removing pixels from magnetogram used for enlarging regions
tmp_in = long(strsplit(PRs[inp].indx,/extract))
imsel[tmp_in] = 0.0

tmp_in = long(strsplit(NRs[inn].indx,/extract))
imsel[tmp_in] = 0.0

;Marking regions as linked and labeling them
PRs[inp].lnk_sw = 1
PRs[inp].ar_lbl = lbl

NRs[inn].lnk_sw = 1  
NRs[inn].ar_lbl = lbl  



lbl = lbl + 1


;;-----------------------------------------------------------------------------------------------------------
;;Second pass. Enlarging new AR if unbalanced-------------------------------------------------------------------------------------------
;
;
;;removing existing ARs from the expansion magnetogram
;if (nars ne 0) then begin
;
;    for i=0,nars-1 do begin 
;
;        ;Removing pixels from magnetogram used for enlarging regions
;        tmp_in = long(strsplit(ARs[arin[i]].indxp,/extract))
;        imsel[tmp_in] = 0.0
;        
;        tmp_in = long(strsplit(ARs[arin[i]].indxn,/extract))
;        imsel[tmp_in] = 0.0
;
;    endfor
;
;endif
;
;
;Flx_im = Nwar.fluxp + Nwar.fluxn
;Flxmin = min(abs([Nwar.fluxp, Nwar.fluxn]))
;
;Flx_imr = Flx_im/Flxmin
;nflxim = Flx_imr
;if keyword_set(info) then print,'FLUX:'
;if keyword_set(info) then print, Flx_im
;if keyword_set(info) then print, Flx_imr
;
;
;if ( abs(Flx_imr) gt ar_cnst.Imb_tol ) then begin
;        
;    ;Extracting indices
;    ingp = long(strsplit(Nwar.indxp,/extract))
;    ingn = long(strsplit(Nwar.indxn,/extract))
;    
;    ;Adding field to the unused magnetogram in order to perform region growth
;    imsel[ingp] =   seg_const.ar_th
;    imsel[ingn] = - seg_const.ar_th
;    
;    ;Storing variables to keep track of the best iteration
;    bstingp = ingp
;    bstingn = ingn
;    bstflxim = Flx_imr
;    
;    ;Iterative enlargement
;    n = 1
;    Th_min = seg_const.ar_th/2;   
;    Th_max = max(abs(ar_cnst.valid_range))  
;    while (abs(nflxim) gt ar_cnst.Imb_tol) and (n le ar_cnst.Imb_it) do begin
;      
;        ;Enlarging the positive region
;        if Flx_imr lt 0.0 then begin
;            
;            ind_gn=idl_region_grow(imsel,ingp,thresh=[Th_min,Th_max])
;            ;Calculating new flux
;            nflux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)
;            nflxim = (Nwar.fluxn + nflux)/nflux 
;            area  = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)
;
;            ;comparing performance and ensuring size is not too large
;            if (abs(nflxim) lt abs(bstflxim)) and (area le 2.0*Nwar.arean) then begin
;                bstflxim = nflxim
;                bstingp = ind_gn
;            endif                
;          
;        ;Enlarging the negative region                
;        endif else begin
;
;            ind_gn=idl_region_grow(-imsel,ingn,thresh=[Th_min,Th_max])
;            ;Calculating new flux
;            nflux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)
;            nflxim = (Nwar.fluxp + nflux)/nflux 
;            area  = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)
;
;            ;comparing performance and ensuring size is not too large
;            if (abs(nflxim) lt abs(bstflxim)) and (area le 2.0*Nwar.areap) then begin
;                bstflxim = nflxim
;                bstingn = ind_gn
;            endif
;                            
;        endelse
;                      
;        ;Performing binary search
;        if abs(nflxim) gt 1.0 then begin
;            Th_min = Th_min - seg_const.ar_th/(2.0^(n+1))
;        endif else begin
;            Th_min = Th_min + seg_const.ar_th/(2.0^(n+1))
;        endelse
;                     
;        n = n+1          
;      
;    endwhile
;    
;    ;recovering best indices
;    if Flx_imr lt 0.0 then begin
;        ind_gn = bstingp
;    endif else begin
;        ind_gn = bstingn          
;    endelse
;    
;    ;Recalculating parameters
;    flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
;    area = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)         ;total area
;    
;    ;Geographical Center
;    fcenx = total(CRD_in.Xar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;    fceny = total(CRD_in.Yar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;    fcenz = total(CRD_in.Zar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;    
;    fcenlt = atan(fcenz,sqrt( fcenx^2.0 + fceny^2.0 ) )/!dtor    ;Latitude center  
;    fcenln = atan(fceny,fcenx)/!dtor                             ;Longitude center
;    
;    ;Center pixels
;    lat1 = fcenlt*!dtor
;    lon1 = fcenln*!dtor
;    
;    lat2 = CRD_in.Lath*!dtor
;    lon2 = CRD_in.Lonh*!dtor
;    
;    dlat = abs(lat1-lat2)
;    dlon = abs(lon1-lon2)
;
;    dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
;    
;    mindis = min( dis_rg, minin, /nan )
;    tmp_2d = array_indices(im,minin)
;    fcenxp = tmp_2d[0,0];    ;Center in x
;    fcenyp = tmp_2d[1,0];    ;Center in y
;    
;    ;Mean and Pixel Radius
;    dcen = mean( dis_rg[ind_gn],/double, /nan)
;    dcenfw = total( dis_rg[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;    
;    dcenp  = dcen/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor            
;    dcenpfw  = dcenfw/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor
;    
;                    
;    ;Storing positive region
;    if Flx_imr lt 0.0 then begin
;                  
;        ;Adding region indexes
;        Nwar.indxp = strjoin(string(ind_gn))                    
;
;        ;Making Sure all visible longitudes are in the same hemisphere
;        ;---------------------------------------------------------
;        sun_data = get_sun(date,carr=carr,he_lon=he_lon)
;    
;        Lonp = fcenln
;        Lonn = Nwar.fcn_lnn
;                                
;        if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
;                
;          tmpin = where(Lonp lt 0, n_tmp)
;          if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
;                
;          tmpin = where(Lonn lt 0, n_tmp)
;          if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
;                
;        endif 
;        
;        
;        ;Finding leading polarity
;        if ( Lonp ge Lonn ) then begin
;        
;            LatL  = fcenlt*!dtor
;            LongL = Lonp*!dtor
;    
;            LatF  = Nwar.fcn_ltn*!dtor
;            LongF = Lonn*!dtor
;            
;            lp = 1;
;        
;        endif else begin
;    
;            LatF  = fcenlt*!dtor
;            LongF = Lonp*!dtor
;    
;            LatL  = Nwar.fcn_ltn*!dtor
;            LongL = Lonn*!dtor
;            
;            lp = -1;
;        
;        endelse
;                
;        dlat = abs(LatF-LatL)
;        dlon = abs(LongF-LongL)
;    
;        Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
;        Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
;        
;        ;Updating relevant quantities
;        tmp_ar = {fluxp: flux, areap: area, $
;                  fcn_ltp: fcenlt, fcn_lnp: fcenln, dcenp: dcen, $ 
;                  fcenxpp: fcenxp, fcenypp: fcenyp, dcenpp: dcenp, $
;                  dis: Dis, tilt: Tilt, lp: lp} 
;                       
;        STRUCT_ASSIGN, tmp_ar, Nwar, /nozero 
;                                                       
;    ;Storing negative region                
;    endif else begin
;      
;        ;Adding region indexes
;        Nwar.indxn = strjoin(string(ind_gn))                    
;
;        ;Making Sure all visible longitudes are in the same hemisphere
;        ;---------------------------------------------------------
;        sun_data = get_sun(date,carr=carr,he_lon=he_lon)
;    
;        Lonp = Nwar.fcn_lnp
;        Lonn = fcenln
;                                
;        if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
;                
;          tmpin = where(Lonp lt 0, n_tmp)
;          if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
;                
;          tmpin = where(Lonn lt 0, n_tmp)
;          if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
;                
;        endif 
;        
;        
;        ;Finding leading polarity
;        if ( Lonp ge Lonn ) then begin
;        
;            LatL  = Nwar.fcn_ltp*!dtor
;            LongL = Lonp*!dtor
;    
;            LatF  = fcenlt*!dtor
;            LongF = Lonn*!dtor
;            
;            lp = 1;
;        
;        endif else begin
;    
;            LatF  = Nwar.fcn_ltp*!dtor
;            LongF = Lonp*!dtor
;    
;            LatL  = fcenlt*!dtor
;            LongL = Lonn*!dtor
;            
;            lp = -1;
;        
;        endelse
;                
;        dlat = abs(LatF-LatL)
;        dlon = abs(LongF-LongL)
;    
;        Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
;        Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
;        
;        ;Updating relevant quantities
;        tmp_ar = {fluxn: flux, arean: area, $
;                  fcn_ltn: fcenlt, fcn_lnn: fcenln, dcenn: dcen, $ 
;                  fcenxpn: fcenxp, fcenypn: fcenyp, dcenpn: dcenp, $
;                  dis: Dis, tilt: Tilt, lp: lp} 
;                       
;        STRUCT_ASSIGN, tmp_ar, Nwar, /nozero 
;                                               
;    endelse
;  
;endif



;Adding region
arin  = where(ARs.mdi_i le mdi_i, nars)
arin2 = where(ARs.mdi_i gt mdi_i, nars2)

if ( (nars eq 0) and (nars2 eq 0) ) then begin
  ARs = Nwar
endif else if (nars eq 0) then begin
  ARs = [Nwar,ARs[arin2]]
endif else if (nars2 eq 0) then begin
  ARs = [ARs[arin],Nwar]
endif else begin
  ARs = [ARs[arin],Nwar,ARs[arin2]]
endelse


if keyword_set(track) then begin

    ;Tracking Backward
    amj_ar_sngl_trck, CRD_in, mdi_i, lbl, Nwar, PRs, NRs, ARs, -1, date, ar_cnst=ar_cnst

    ;Tracking Forward
    amj_ar_sngl_trck, CRD_in, mdi_i, lbl, Nwar, PRs, NRs, ARs, 1, date, ar_cnst=ar_cnst

  
endif
            


;
;display the extracted region
;

if keyword_set(display) then begin
    set_plot,'X'
    window,8+pltn,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    loadct, 0, /silent
    tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
 
    if (n_ars gt 1) then begin
        for i=0,n_ars-2 do begin
             
            loadct, 13, /silent
            
            str='P'+strtrim(string(AR_c.ARs[i].labl,'(I5)'),2)
            xyouts,ARs[i].fcenxpp*display_zoom,ARs[i].fcenypp*display_zoom,str,charsize=1.5,color= long(AR_c.ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpp*display_zoom, ARs[i].fcenxpp*display_zoom,ARs[i].fcenypp*display_zoom, color = strtrim(string(AR_c.ARs[i].clr,'(I3)'),2),/device
            
            str='N'+strtrim(string(AR_c.ARs[i].labl,'(I5)'),2)
            xyouts,ARs[i].fcenxpn*display_zoom,ARs[i].fcenypn*display_zoom,str,charsize=1.5,color= long(AR_c.ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpn*display_zoom, ARs[i].fcenxpn*display_zoom,ARs[i].fcenypn*display_zoom, color = strtrim(string(AR_c.ARs[i].clr,'(I3)'),2),/device
    
            loadct, 0, /silent
        endfor    
    endif 
    
endif

flnm = 'Mgntgram_ARs.eps'

if keyword_set(sqs_nm) then begin
    flnm = '/home/amunoz/AR_Detection/imag/AR_dec' + strtrim(string(sqs_nm/1000.,'(F5.3)'),2) + '.eps'
endif

if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename=flnm,/color, encapsulated = 1, decomposed = 0, bits_per_pixel=8,/portr,/inches,xsize=6.0*print_zoom,ysize=6.0*print_zoom,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    px = !d.x_vsize ;Get size of window in device units
    py = !d.y_vsize

    plot,[1,1],/nodata,xstyle=5,ystyle=5
    ;tv,bytscl(congrid(im_s_gp+im_s_gn,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu) 
    loadct, 0, /silent
    tv,bytscl(congrid(im,print_xsize,print_ysize),min=display_thresholdd,max=display_thresholdu)
    
    if (n_ars gt 1) then begin
        for i=0,n_ars-2 do begin
            str='P'+strtrim(string(AR_c.ARs[i].labl,'(I5)'),2)
            
            loadct, 13, /silent
            
            !P.COLOR=long(AR_c.ARs[i].clr)
            xyouts,ARs[i].fcenxpp*px/sz[1],ARs[i].fcenypp*py/sz[2],str,charsize=1.5,charthick = 4,/device
;            xyouts,ARs[i].fcenxpp*px/sz[1],ARs[i].fcenypp*py/sz[2],str,charsize=1.5,charthick = 2,color = long(AR_c.ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpp*px/sz[1], ARs[i].fcenxpp*px/sz[1],ARs[i].fcenypp*py/sz[2], color = strtrim(string(AR_c.ARs[i].clr,'(I3)'),2), thick=4,/device
            
            str='N'+strtrim(string(AR_c.ARs[i].labl,'(I5)'),2)
            !P.COLOR=long(AR_c.ARs[i].clr)
            xyouts,ARs[i].fcenxpn*px/sz[1],ARs[i].fcenypn*py/sz[2],str,charsize=1.5,charthick = 4,/device
;            xyouts,ARs[i].fcenxpn*px/sz[1],ARs[i].fcenypn*py/sz[2],str,charsize=1.5,charthick = 2,color= long(AR_c.ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpn*px/sz[1], ARs[i].fcenxpn*px/sz[1],ARs[i].fcenypn*py/sz[2], color = strtrim(string(AR_c.ARs[i].clr,'(I3)'),2), thick=4,/device
    
            
        endfor    
    endif
    
    loadct, 0, /silent
    
    xyouts,0.0*px,0.0*py,date,charsize=1.75, charthick = 6,color= 254,/device
    
    device, /close_file
    set_plot,'X'
    loadct, 0, /silent           
    
    
endif

;print, junk

return
end