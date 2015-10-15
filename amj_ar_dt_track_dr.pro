pro amj_ar_dt_track_dr, CRD_in, mdi_i, lbl, PRs, NRs, AR_c, date, ar_cnst=ar_cnst, seg_const=seg_const, bck_trck = bck_trck, display=display, prt=prt, info=info, sqs_nm = sqs_nm, pltn = pltn
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
;       amj_ar_dt_track_dr, CRD_in, mdi_i, lbl, PRs, NRs, AR_c, date
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
;       AR_c: a structure holding the detected active regions
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
ar_cnst_def={ar_cnst, dis_lim1:6.0, dis_lim2:3.0, exp_f: 1.0,  exp_d: 4.0, exp_s: 0.5, mlth: 40.0, mxB: 180.0, MxFlxim:3.0, Imb_tol: 0.10, Imb_it: 10, lim_lon: 0.0, k_sig:15.0, npr: 5, nmgnt: 5 , vld_thr: 0.69, valid_range:[-20000.,20000.]}
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
display_zoom=0.125
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

print, 'Tracking bipolar regions'
;-----------------------------------------------------------------------------------------------------------
;First pass. Tracking positive and negative regions using the past--------------------------------------------------------

n_ars = 1 ;Number of found active regions
vis_in = -1 ;Variable that stores indices of revisited ARs


;finding Available PNRs
pr_ix = where((PRs.mdi_i eq mdi_i) and (PRs.lnk_sw eq 0), n_regp)
if (n_regp gt 0) then begin
    pregions = PRs[pr_ix]
endif

nr_ix = where((NRs.mdi_i eq mdi_i) and (NRs.lnk_sw eq 0), n_regn)
if (n_regn gt 0) then begin
    nregions = NRs[nr_ix]
endif


lat2p = pregions.fcn_lt
lon2p = pregions.fcn_ln

lat2n = nregions.fcn_lt
lon2n = nregions.fcn_ln

sun_data = get_sun(date,carr=carr,he_lon=he_lon)

;Making Sure all visible longitudes are in the same hemisphere
;--------------------------------------------------------- 
if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
        
  tmpin = where(lon2p lt 0, n_tmp)
  if (n_tmp gt 0) then lon2p[tmpin] = lon2p[tmpin]+360.0

  tmpin = where(lon2n lt 0, n_tmp)
  if (n_tmp gt 0) then lon2n[tmpin] = lon2n[tmpin]+360.0
        
endif

lat2p = lat2p*!dtor
lon2p = lon2p*!dtor

lat2n = lat2n*!dtor
lon2n = lon2n*!dtor

arinp  = where(AR_c.mdi_i eq mdi_i, narsp)

for n = 1, ar_cnst.npr do begin

    ;Check for AR in the previous day
    arin  = where(AR_c.mdi_i eq mdi_i-n, nars)
    
    ;Look for regions
    if ((nars ne 0) and (n_regp ne 0) and (n_regn ne 0) ) then begin
            
        tARs = AR_c[arin]
        timep = tARs[0].date
        timec = date    
                
        LatpAR = tARs.fcn_ltp        
        LonpAR = tARs.fcn_lnp
        
        LatnAR = tARs.fcn_ltn        
        LonnAR = tARs.fcn_lnn
        

        
        ;Making Sure all visible longitudes are in the same hemisphere
        ;--------------------------------------------------------- 
    
        if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
                
          tmpin = where(LonpAR lt 0, n_tmp)
          if (n_tmp gt 0) then LonpAR[tmpin] = LonpAR[tmpin]+360.0
    
          tmpin = where(LonnAR lt 0, n_tmp)
          if (n_tmp gt 0) then LonnAR[tmpin] = LonnAR[tmpin]+360.0
                
        endif
        
        ;Addition of differential rotation
        ;-------------------------------------       
        
        dtim = anytim(timec,/utime)-anytim(timep,/utime)
        
        ; Differential rotation profile from Snodgrass (1983), gives rotation
        ; rate vs. lat. in microradians per sec:
        ;
        ;  omega = snod_A + snod_B*sin(latitude)^2 + snod_C*sin(latitude)^4
        ;====================================================================
        snod_A =  0.0367;2.902  ; magnetic rot. coeffs, in microrad.    
        ;Set to 0.0367 because the heliographic coordinates include the carrington rotation
        
        snod_B = -0.464
        snod_C = -0.328   
        
        ;Accounting for differential rotation
        omegap = 1e-6*(snod_A + $ ; Omega at each pixel's lat., in radians
                      snod_B*sin(abs(LatpAR)*!dtor)^2 + $
                      snod_C*sin(abs(LatpAR)*!dtor)^4 )/!dtor
        
        LonpAR = LonpAR + omegap*dtim       
        
        
        omegan = 1e-6*(snod_A + $ ; Omega at each pixel's lat., in radians
                      snod_B*sin(abs(LatnAR)*!dtor)^2 + $
                      snod_C*sin(abs(LatnAR)*!dtor)^4 )/!dtor
        
        LonnAR = LonnAR + omegan*dtim       
        
        ;Checking if past active regions are present in this image
        ;---------------------------------------------------------
        
        
    ;    ;Making Sure all visible longitudes are in the same hemisphere
    ;    ;---------------------------------------------------------
    ;
    ;    ;Finding the minimum and maximum longitudes for all regions
    ;    minlon = min([min(lon2p), min(lon2n), min(LonpAR), min(LonnAR)])
    ;    maxlon = max([max(lon2p), max(lon2n), max(LonpAR), max(LonnAR)])
    ;    
    ;    ;Verifying they are in the same hemisphere
    ;    if (maxlon-minlon gt 180.0) then begin
    ;      
    ;      tmpin = where(lon2p lt 0, n_tmp)
    ;      if (n_tmp gt 0) then lon2p[tmpin] = lon2p[tmpin]+360.0
    ;
    ;      tmpin = where(lon2n lt 0, n_tmp)
    ;      if (n_tmp gt 0) then lon2n[tmpin] = lon2n[tmpin]+360.0
    ;      
    ;      tmpin = where(LonpAR lt 0, n_tmp)
    ;      if (n_tmp gt 0) then LonpAR[tmpin] = LonpAR[tmpin]+360.0
    ;
    ;      tmpin = where(LonnAR lt 0, n_tmp)
    ;      if (n_tmp gt 0) then LonnAR[tmpin] = LonnAR[tmpin]+360.0
    ;            
    ;    endif   
        
        for i = 0, nars - 1 do begin
          
            arinp  = where(AR_c.mdi_i eq mdi_i, narsp)
            pres_sw = 0
            
            if (narsp eq 0) then begin
              pres_sw = 1
            endif else begin
              if (total(tARs[i].labl eq AR_c[arinp].labl) eq 0) then pres_sw = 1             
            endelse
            
            
            ;Making sure the region doesn't exist already in the snapshot
            if ( pres_sw eq 1 ) then begin            
          
                sw_tr = 1; Switch which idicates whereas matching positive and negative regions have been found   
            
                ;Positive Side
                lat1 = LatpAR[i]*!dtor
                lon1 = LonpAR[i]*!dtor
                                    
                dlat = abs(lat1-lat2p)
                dlon = abs(lon1-lon2p)
            
                dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2p)*sin(dlon/2.0)^2.0 ))/!dtor                    
                mindis = min( dis_rg, minin, /nan )
                
                if (mindis le min([pregions[minin].dcen, tARs[i].dcenp])*ar_cnst.dis_lim2) then begin
                    jj = minin
                endif else begin
                    sw_tr = 0; No matching positive regions were found
                endelse
                
            
                ;Negative Side
                lat1 = LatnAR[i]*!dtor
                lon1 = LonnAR[i]*!dtor
                                    
                dlat = abs(lat1-lat2n)
                dlon = abs(lon1-lon2n)
            
                dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2n)*sin(dlon/2.0)^2.0 ))/!dtor                    
                mindis = min( dis_rg, minin, /nan )
                
                if (mindis le min([nregions[minin].dcen, tARs[i].dcenn])*ar_cnst.dis_lim2) then begin
                    ii = minin
                endif else begin
                    sw_tr = 0; No matching negative regions were found
                endelse
                
                ;Making sure previously used regions are not repeated        
                if ( (sw_tr eq 1) and ( total(tARs[i].labl eq vis_in) gt 0 ) ) then begin
                    sw_tr = 0
                endif
                
                ;Agregating found positive and negative regions to make ARs
                
                if (sw_tr eq 1) then begin
                                    
                    ;Making Sure all visible longitudes are in the same hemisphere
                    ;---------------------------------------------------------
        
                    Lonp = pregions[jj].fcn_ln
                    Lonn = nregions[ii].fcn_ln

                    if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
                            
                      tmpin = where(Lonp lt 0, n_tmp)
                      if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
                            
                      tmpin = where(Lonn lt 0, n_tmp)
                      if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
                            
                    endif         
        
                    ;Finding leading polarity
                    if ( Lonp ge Lonn ) then begin
                    
                        LatL  = pregions[jj].fcn_lt*!dtor
                        LongL = Lonp*!dtor
                
                        LatF  = nregions[ii].fcn_lt*!dtor
                        LongF = Lonn*!dtor
                        
                        lp = 1;
                    
                    endif else begin
                
                        LatF  = pregions[jj].fcn_lt*!dtor
                        LongF = Lonp*!dtor
                
                        LatL  = nregions[ii].fcn_lt*!dtor
                        LongL = Lonn*!dtor
                        
                        lp = -1;
                    
                    endelse
        
                    dlat = abs(LatF-LatL)
                    dlon = abs(LongF-LongL)                   
                                               
                    Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
                    Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
                    
                    dm = 0
                    qm = 0
                    
                    ; Store positive and region number in a string, retrieve via:
                    ; nmbr = long(strsplit(ar.pregions,/extract)) or nmbr = long(strsplit(ar.nregions,/extract))
                    ar = {ar, mdi_i: mdi_i, date: date, labl: tARs[i].labl, clr: tARs[i].clr, indxp: pregions[jj].indx, indxn: nregions[ii].indx, fluxp: pregions[jj].flux, fluxn: nregions[ii].flux, areap:pregions[jj].area, arean:nregions[ii].area, $
                              fcn_ltp:pregions[jj].fcn_lt, fcn_lnp:pregions[jj].fcn_ln, dcenp:pregions[jj].dcen, $ 
                              fcn_ltn:nregions[ii].fcn_lt, fcn_lnn:nregions[ii].fcn_ln, dcenn:nregions[ii].dcen, $
                              fcenxpp: pregions[jj].fcenxp, fcenypp: pregions[jj].fcenyp, dcenpp:pregions[jj].dcenp, $
                              fcenxpn: nregions[ii].fcenxp, fcenypn: nregions[ii].fcenyp, dcenpn:nregions[ii].dcenp, $               
                              dis: Dis, tilt: Tilt, lp: lp, dm: dm, qm: qm}    
            
                    ;set_plot,'X'
                    ;window,i,xsize=display_xsize,ysize=display_ysize,retain=2
                    ;plot,[1,1],/nodata,xstyle=5,ystyle=5
                    ;loadct, 0, /silent
                    ;tv,bytscl(congrid(tmp_msk,display_xsize,display_ysize),min=minmsk,max=maxmsk)
            
                    ;Removing pixels from magnetogram used for enlarging regions
                    tmp_in = long(strsplit(pregions[jj].indx,/extract))
                    imsel[tmp_in] = 0.0
                    
                    tmp_in = long(strsplit(nregions[ii].indx,/extract))
                    imsel[tmp_in] = 0.0
            
                    ;Removing paired regions
                    lat2p[jj] = !values.f_nan
                    lon2p[jj] = !values.f_nan
                    PRs[pr_ix[jj]].lnk_sw = 1
                    PRs[pr_ix[jj]].ar_lbl = tARs[i].labl
        
                    
                    lat2n[ii] = !values.f_nan
                    lon2n[ii] = !values.f_nan
                    NRs[nr_ix[ii]].lnk_sw = 1
                    NRs[nr_ix[ii]].ar_lbl = tARs[i].labl
                                            
                    vis_in = [vis_in, tARs[i].labl]                     
                
                    if ( n_ars eq 1 ) then begin
                        ARs = ar
                    endif else begin
                        ARs = [ARs, ar]
                    endelse
                    
                    if keyword_set(info) then print,';AR No. ',string(n_ars,'(I3)'),' combines regions P', string(jj,'(I3)'), ' and N', string(ii,'(I3)')
                    n_ars = n_ars + 1    
                                
                endif
            
            endif
                    
        endfor
    
    endif

endfor
       
        
if keyword_set(display) then begin
    set_plot,'X'
    window,1,xsize=display_xsize,ysize=display_ysize,retain=2
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    loadct, 0, /silent
    tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
 
    if (n_ars gt 1) then begin
        for i=0,n_ars-2 do begin
             
            loadct, 13, /silent
            
            str='P'+strtrim(string(ARs[i].labl,'(I5)'),2)
            xyouts,ARs[i].fcenxpp*display_zoom,ARs[i].fcenypp*display_zoom,str,charsize=1.5,color= long(ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpp*display_zoom, ARs[i].fcenxpp*display_zoom,ARs[i].fcenypp*display_zoom, color = strtrim(string(ARs[i].clr,'(I3)'),2),/device
            
            str='N'+strtrim(string(ARs[i].labl,'(I5)'),2)
            xyouts,ARs[i].fcenxpn*display_zoom,ARs[i].fcenypn*display_zoom,str,charsize=1.5,color= long(ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpn*display_zoom, ARs[i].fcenxpn*display_zoom,ARs[i].fcenypn*display_zoom, color = strtrim(string(ARs[i].clr,'(I3)'),2),/device
    
            loadct, 0, /silent
        endfor    
    endif 
    
endif



print, 'Pairing bipolar regions'
;-----------------------------------------------------------------------------------------------------------
;Second pass. Paring remaining positive and negative regions--------------------------------------------------------


;finding Available PRs
pr_ix = where((PRs.mdi_i eq mdi_i) and (PRs.lnk_sw eq 0), n_regp)
if (n_regp gt 0) then begin
    pregions = PRs[pr_ix]
endif

;finding Available NRs
nr_ix = where((NRs.mdi_i eq mdi_i) and (NRs.lnk_sw eq 0), n_regn)
if (n_regn gt 0) then begin
    nregions = NRs[nr_ix]
endif

;Pairing distance
n_dis = ar_cnst.dis_lim1
if ( (n_regp gt 0) and (n_regn gt 0) ) then begin
  
    ;Creating distance matrix
    DisM = fltarr(n_regn,n_regp) ; Distance matrix
    
    ;Filling Distance matrix up 
    for i=0,n_regn-1 do begin
    
        for j=0,n_regp-1 do begin
          
        
            ;Calculating distances between its centroids
            lat1 = pregions[j].fcn_lt*!dtor
            lon1 = pregions[j].fcn_ln*!dtor
            
            lat2 = nregions[i].fcn_lt*!dtor
            lon2 = nregions[i].fcn_ln*!dtor
            
            dlat = abs(lat1-lat2)
            dlon = abs(lon1-lon2)
    
            dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
        
            ;Calculating net flux
            del_flx = abs(nregions[i].flux + pregions[j].flux)/1.0e22
    
            ;Calculating difference in size
            del_sz = abs(nregions[i].dcen - pregions[j].dcen)
                
            ;Finding regions withing prescribed distance
            if ( ( dis_rg le n_dis*nregions[i].dcen ) or ( dis_rg le n_dis*pregions[j].dcen ) ) then begin
                
                DisM(i,j) = (del_flx^ar_cnst.exp_f)*(dis_rg^ar_cnst.exp_d)*(del_sz^ar_cnst.exp_s)
            
            endif else begin
            
                DisM(i,j) = !values.f_nan
            
            endelse 
                     
            ;Ignoring regions that have too much flux imbalance
            minflx = min([abs(pregions[j].flux), abs(nregions[i].flux)])
            flxim = (abs(nregions[i].flux + pregions[j].flux))/minflx
            
            if ( flxim ge ar_cnst.MxFlxim ) then begin
                DisM(i,j) = !values.f_nan
            endif                   
        
        endfor    
    
    endfor
    
    ;Agregating positive and negative regions to make ARs        
    while (total(finite(DisM)) ne 0) do begin
    
        szM = size(DisM)
        
        if (szM(0) eq 2) then begin
        
            mdis = min( min(DisM,mt,dimension=2,/nan) , ii,/nan)
            mdis = min(DisM(ii,*),jj,/nan)
        
        endif else begin
    
            ii = 0;
            mdis = min(DisM(ii,*),jj,/nan)
        
        endelse
        
        ;Making Sure all visible longitudes are in the same hemisphere
        ;---------------------------------------------------------
        
        Lonp = pregions[jj].fcn_ln
        Lonn = nregions[ii].fcn_ln
                
        sun_data = get_sun(date,carr=carr,he_lon=he_lon)
        
        if( he_lon gt 180.00) then begin
          he_lon2 = he_lon - 360.00          
        endif else begin
          he_lon2 = he_lon
        endelse

        if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
                
          tmpin = where(Lonp lt 0, n_tmp)
          if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
                
          tmpin = where(Lonn lt 0, n_tmp)
          if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
          
          he_lon2 = he_lon
                
        endif 
        
        
        ;Finding leading polarity
        if ( Lonp ge Lonn ) then begin
        
            LatL  = pregions[jj].fcn_lt*!dtor
            LongL = Lonp*!dtor
    
            LatF  = nregions[ii].fcn_lt*!dtor
            LongF = Lonn*!dtor
            
            lp = 1;
        
        endif else begin
    
            LatF  = pregions[jj].fcn_lt*!dtor
            LongF = Lonp*!dtor
    
            LatL  = nregions[ii].fcn_lt*!dtor
            LongL = Lonn*!dtor
            
            lp = -1;
        
        endelse
        
        ;only detecting ahead of limit meridian
        if (LongL/!dtor ge he_lon2 + ar_cnst.lim_lon) then begin
    
            dlat = abs(LatF-LatL)
            dlon = abs(LongF-LongL)
        
            Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
            Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
    
            dm = 0
            qm = 0
         
            ; Store positive and region number in a string, retrieve via:
            ; nmbr = long(strsplit(ar.pregions,/extract)) or nmbr = long(strsplit(ar.nregions,/extract))
            ar = {ar, mdi_i: mdi_i, date: date, labl: lbl, clr: clrsqs[lbl mod 17], indxp: pregions[jj].indx, indxn: nregions[ii].indx, fluxp: pregions[jj].flux, fluxn: nregions[ii].flux, areap:pregions[jj].area, arean:nregions[ii].area, $
                      fcn_ltp:pregions[jj].fcn_lt, fcn_lnp:pregions[jj].fcn_ln, dcenp:pregions[jj].dcen, $ 
                      fcn_ltn:nregions[ii].fcn_lt, fcn_lnn:nregions[ii].fcn_ln, dcenn:nregions[ii].dcen, $
                      fcenxpp: pregions[jj].fcenxp, fcenypp: pregions[jj].fcenyp, dcenpp:pregions[jj].dcenp, $
                      fcenxpn: nregions[ii].fcenxp, fcenypn: nregions[ii].fcenyp, dcenpn:nregions[ii].dcenp, $               
                      dis: Dis, tilt: Tilt, lp: lp, dm: dm, qm: qm}           
                
            if ( n_ars eq 1 ) then begin
                ARs = ar
            endif else begin
                ARs = [ARs, ar]
            endelse
            
            
            if keyword_set(info) then print,';AR No. ',string(n_ars,'(I3)'),' combines regions P', string(n_inx[ii],'(I3)'), ' and N', string(n_inx[ii],'(I3)')
            n_ars = n_ars + 1
     
            ;Removing pixels from magnetogram used for enlarging regions
            tmp_in = long(strsplit(pregions[jj].indx,/extract))
            imsel[tmp_in] = 0.0
            
            tmp_in = long(strsplit(nregions[ii].indx,/extract))
            imsel[tmp_in] = 0.0
                
            ;Removing paired regions
            DisM(ii,*) = !values.f_nan
            DisM(*,jj) = !values.f_nan
            
            PRs[pr_ix[jj]].lnk_sw = 1
            PRs[pr_ix[jj]].ar_lbl = lbl
    
            NRs[nr_ix[ii]].lnk_sw = 1
            NRs[nr_ix[ii]].ar_lbl = lbl

            lbl = lbl + 1
                        
            ;Tracking Backward
            if keyword_set(bck_trck) then amj_ar_sngl_trck, CRD_in, mdi_i, lbl, ar, PRs, NRs, AR_c, -1, date, ar_cnst=ar_cnst
                    
        endif else begin
          
            DisM(ii,jj) = !values.f_nan
                      
        endelse
        
                         
     
    endwhile
      
endif



;;Third pass. Enlarging imbalanced ARs-------------------------------------------------------------------------------------------
;
;if (n_ars ne 1) then begin
;
;    Flx_im = ARs.fluxp + ARs.fluxn
;    Flxmin = min(abs([[ARs.fluxp], [ARs.fluxn]]), dimension = 2)
;    
;    Flx_imr = Flx_im/Flxmin
;    if keyword_set(info) then print,'FLUX:'
;    if keyword_set(info) then print, Flx_im
;    if keyword_set(info) then print, Flx_imr
;    
;    im_in = where( abs(Flx_imr) le ar_cnst.Imb_tol , nbal)    
;    if (nbal gt 0) then Flx_imr[im_in] = !values.f_nan
;    
;    if keyword_set(info) then print, Flx_imr
;    
;    for i=0,n_ars-2-nbal do begin
;    
;        ;Finding the largest imbalance
;        maxim = max(abs(Flx_imr),maximn,/nan)
;        nflxim = Flx_imr[maximn]
;        
;        ;Extracting indices
;        inp = long(strsplit(ARs[maximn].indxp,/extract))
;        inn = long(strsplit(ARs[maximn].indxn,/extract))
;        
;        ;Adding field to the unused magnetogram in order to perform region growth
;        imsel[inp] =   seg_const.ar_th
;        imsel[inn] = - seg_const.ar_th
;        
;        ;Storing variables to keep track of the best iteration
;        bstinp = inp
;        bstinn = inn
;        bstflxim = nflxim
;        
;        ;Iterative enlargement
;;        print, nflxim
;        n = 1
;        Th_min = seg_const.ar_th/2;   
;        Th_max = max(abs(ar_cnst.valid_range))  
;        while (abs(nflxim) gt ar_cnst.Imb_tol) and (n le ar_cnst.Imb_it) do begin
;          
;            ;Enlarging the positive region
;            if Flx_imr[maximn] lt 0.0 then begin
;                
;                ind_gn=idl_region_grow(imsel,inp,thresh=[Th_min,Th_max])
;                ;Calculating new flux
;                nflux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)
;                nflxim = (ARs[maximn].fluxn + nflux)/nflux 
;                area  = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)
;
;                ;comparing performance and ensuring size is not too large
;                if (abs(nflxim) lt abs(bstflxim)) and (area le 2.0*ARs[maximn].arean) then begin
;                    bstflxim = nflxim
;                    bstinp = ind_gn
;                endif                
;              
;            ;Enlarging the negative region                
;            endif else begin
;
;                ind_gn=idl_region_grow(-imsel,inn,thresh=[Th_min,Th_max])
;                ;Calculating new flux
;                nflux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)
;                nflxim = (ARs[maximn].fluxp + nflux)/nflux 
;                area  = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)
;
;                ;comparing performance and ensuring size is not too large
;                if (abs(nflxim) lt abs(bstflxim)) and (area le 2.0*ARs[maximn].areap) then begin
;                    bstflxim = nflxim
;                    bstinn = ind_gn
;                endif
;                                
;            endelse
;            
;;            print, nflxim
;              
;            ;Performing binary search
;            if abs(nflxim) gt 1.0 then begin
;                Th_min = Th_min - seg_const.ar_th/(2.0^(n+1))
;            endif else begin
;                Th_min = Th_min + seg_const.ar_th/(2.0^(n+1))
;            endelse
;                         
;            n = n+1          
;          
;        endwhile
;        
;        ;recovering best indices
;        if Flx_imr[maximn] lt 0.0 then begin
;            ind_gn = bstinp
;        endif else begin
;            ind_gn = bstinn          
;        endelse
;;        print, bstflxim
;        
;        ;Recalculating parameters
;        flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
;        area = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)         ;total area
;        
;        ;Geographical Center
;        fcenx = total(CRD_in.Xar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;        fceny = total(CRD_in.Yar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;        fcenz = total(CRD_in.Zar[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;        
;        fcenlt = atan(fcenz,sqrt( fcenx^2.0 + fceny^2.0 ) )/!dtor    ;Latitude center  
;        fcenln = atan(fceny,fcenx)/!dtor                             ;Longitude center
;        
;        ;Center pixels
;        lat1 = fcenlt*!dtor
;        lon1 = fcenln*!dtor
;        
;        lat2 = CRD_in.Lath*!dtor
;        lon2 = CRD_in.Lonh*!dtor
;        
;        dlat = abs(lat1-lat2)
;        dlon = abs(lon1-lon2)
;    
;        dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
;        
;        mindis = min( dis_rg, minin, /nan )
;        tmp_2d = array_indices(im,minin)
;        fcenxp = tmp_2d[0,0];    ;Center in x
;        fcenyp = tmp_2d[1,0];    ;Center in y
;        
;        ;Mean and Pixel Radius
;        dcen = mean( dis_rg[ind_gn],/double, /nan)
;        dcenfw = total( dis_rg[ind_gn]*CRD_in.mgnt_flx[ind_gn],/double, /nan)/flux
;        
;        dcenp  = dcen/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor            
;        dcenpfw  = dcenfw/mean([ dis_rg[fcenxp-1,fcenyp], dis_rg[fcenxp+1,fcenyp], dis_rg[fcenxp,fcenyp-1], dis_rg[fcenxp,fcenyp+1] ],/double, /nan);/!dtor
;        
;                        
;        ;Storing positive region
;        if Flx_imr[maximn] lt 0.0 then begin
;                      
;            ;Adding region indexes
;            ARs[maximn].indxp = strjoin(string(ind_gn))                    
;
;            ;Making Sure all visible longitudes are in the same hemisphere
;            ;---------------------------------------------------------            
;            Lonp = fcenln
;            Lonn = ARs[maximn].fcn_lnn
;                                    
;            if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
;                    
;              tmpin = where(Lonp lt 0, n_tmp)
;              if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
;                    
;              tmpin = where(Lonn lt 0, n_tmp)
;              if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
;              
;              he_lon2 = he_lon
;                    
;            endif 
;            
;            
;            ;Finding leading polarity
;            if ( Lonp ge Lonn ) then begin
;            
;                LatL  = fcenlt*!dtor
;                LongL = Lonp*!dtor
;        
;                LatF  = ARs[maximn].fcn_ltn*!dtor
;                LongF = Lonn*!dtor
;                
;                lp = 1;
;            
;            endif else begin
;        
;                LatF  = fcenlt*!dtor
;                LongF = Lonp*!dtor
;        
;                LatL  = ARs[maximn].fcn_ltn*!dtor
;                LongL = Lonn*!dtor
;                
;                lp = -1;
;            
;            endelse
;                    
;            dlat = abs(LatF-LatL)
;            dlon = abs(LongF-LongL)
;        
;            Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
;            Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
;            
;            ;Updating relevant quantities
;            tmp_ar = {fluxp: flux, areap: area, $
;                      fcn_ltp: fcenlt, fcn_lnp: fcenln, dcenp: dcen, $ 
;                      fcenxpp: fcenxp, fcenypp: fcenyp, dcenpp: dcenp, $
;                      dis: Dis, tilt: Tilt, lp: lp} 
;                           
;            ar = ARs[maximn]                        
;            STRUCT_ASSIGN, tmp_ar, ar, /nozero 
;            ARs[maximn] = ar  
;                             
;            ;Updating indices to reset unused magnetogram
;            inp = ind_gn
;                              
;        ;Storing negative region                
;        endif else begin
;          
;            ;Adding region indexes
;            ARs[maximn].indxn = strjoin(string(ind_gn))                    
;
;            ;Making Sure all visible longitudes are in the same hemisphere
;            ;---------------------------------------------------------            
;            Lonp = ARs[maximn].fcn_lnp
;            Lonn = fcenlt
;                                    
;            if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
;                    
;              tmpin = where(Lonp lt 0, n_tmp)
;              if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
;                    
;              tmpin = where(Lonn lt 0, n_tmp)
;              if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
;              
;              he_lon2 = he_lon
;                    
;            endif 
;            
;            
;            ;Finding leading polarity
;            if ( Lonp ge Lonn ) then begin
;            
;                LatL  = ARs[maximn].fcn_ltp*!dtor
;                LongL = Lonp*!dtor
;        
;                LatF  = fcenlt*!dtor
;                LongF = Lonn*!dtor
;                
;                lp = 1;
;            
;            endif else begin
;        
;                LatF  = ARs[maximn].fcn_ltp*!dtor
;                LongF = Lonp*!dtor
;        
;                LatL  = fcenlt*!dtor
;                LongL = Lonn*!dtor
;                
;                lp = -1;
;            
;            endelse
;                    
;            dlat = abs(LatF-LatL)
;            dlon = abs(LongF-LongL)
;        
;            Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
;            Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
;            
;            ;Updating relevant quantities
;            tmp_ar = {fluxn: flux, arean: area, $
;                      fcn_ltn: fcenlt, fcn_lnn: fcenln, dcenn: dcen, $ 
;                      fcenxpn: fcenxp, fcenypn: fcenyp, dcenpn: dcenp, $
;                      dis: Dis, tilt: Tilt, lp: lp} 
;                           
;            ar = ARs[maximn]                        
;            STRUCT_ASSIGN, tmp_ar, ar, /nozero 
;            ARs[maximn] = ar
;            
;            ;Updating indices to reset unused magnetogram
;            inn = ind_gn
;                                                   
;        endelse
;      
;        ;Removing region from the unused magnetogram
;        imsel[inp] = 0.0
;        imsel[inn] = 0.0
;        
;        
;        Flx_imr[maximn] = !values.f_nan
;;        print, Flx_imr
;
;    endfor
;    
;endif


;Adding regions
if (n_ars gt 1) then begin
    
    arin  = where(AR_c.mdi_i le mdi_i, nars)
    arin2 = where(AR_c.mdi_i gt mdi_i, nars2)
    
    if ( (nars eq 0) and (nars2 eq 0) ) then begin
      AR_c = ARs
    endif else if (nars eq 0) then begin
      AR_c = [ARs,AR_c[arin2]]
    endif else if (nars2 eq 0) then begin
      AR_c = [AR_c[arin],ARs]
    endif else begin
      AR_c = [AR_c[arin],ARs,AR_c[arin2]]
    endelse 

    
endif



;
;display the extracted regions
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
            
            str='P'+strtrim(string(ARs[i].labl,'(I5)'),2)
            xyouts,ARs[i].fcenxpp*display_zoom,ARs[i].fcenypp*display_zoom,str,charsize=1.5,color= long(ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpp*display_zoom, ARs[i].fcenxpp*display_zoom,ARs[i].fcenypp*display_zoom, color = strtrim(string(ARs[i].clr,'(I3)'),2),/device
            
            str='N'+strtrim(string(ARs[i].labl,'(I5)'),2)
            xyouts,ARs[i].fcenxpn*display_zoom,ARs[i].fcenypn*display_zoom,str,charsize=1.5,color= long(ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpn*display_zoom, ARs[i].fcenxpn*display_zoom,ARs[i].fcenypn*display_zoom, color = strtrim(string(ARs[i].clr,'(I3)'),2),/device
    
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
            str='P'+strtrim(string(ARs[i].labl,'(I5)'),2)
            
            loadct, 13, /silent
            
            !P.COLOR=long(ARs[i].clr)
            xyouts,ARs[i].fcenxpp*px/sz[1],ARs[i].fcenypp*py/sz[2],str,charsize=1.5,charthick = 4,/device
;            xyouts,ARs[i].fcenxpp*px/sz[1],ARs[i].fcenypp*py/sz[2],str,charsize=1.5,charthick = 2,color = long(ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpp*px/sz[1], ARs[i].fcenxpp*px/sz[1],ARs[i].fcenypp*py/sz[2], color = strtrim(string(ARs[i].clr,'(I3)'),2), thick=4,/device
            
            str='N'+strtrim(string(ARs[i].labl,'(I5)'),2)
            !P.COLOR=long(ARs[i].clr)
            xyouts,ARs[i].fcenxpn*px/sz[1],ARs[i].fcenypn*py/sz[2],str,charsize=1.5,charthick = 4,/device
;            xyouts,ARs[i].fcenxpn*px/sz[1],ARs[i].fcenypn*py/sz[2],str,charsize=1.5,charthick = 2,color= long(ARs[i].clr),/device
            tvcircle,2.0*ARs[i].dcenpn*px/sz[1], ARs[i].fcenxpn*px/sz[1],ARs[i].fcenypn*py/sz[2], color = strtrim(string(ARs[i].clr,'(I3)'),2), thick=4,/device
    
            
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