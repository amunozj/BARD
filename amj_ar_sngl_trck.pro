pro amj_ar_sngl_trck, CRD_in, mdi_i, lbl, tARs, PRs, NRs, ARs, trk_sw, date, ar_cnst=ar_cnst
;+
; NAME:
;      amj_ar_sngl_trck
;
; PURPOSE:
;       looks for instances of the same active region into the past, or into the future, and assings it
;       the same label 
;       
;
; CALLING SEQUENCE:
;       amj_ar_dt_track_dr, CRD_in, mdi_i, lbl, tARs, PRs, NRs, ARs, trk_sw, date
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
;       ARs: a structure holding the active region to be paired
;            1. seg_const: used parameters in the positive and negative region detection (see amj_pnr_dt.pro)
;            2. ar_cnst: struct holding control parameters for AR pairing (see above) 
;            3. hdr: heather of the MDI magnetogram
;            4. num_ar: Number of detected active regions
;            5. ARs: ar structure holding details from all detected active regions 
;              (see below for details on the ar structure)
;            6. ar_sw: Switch that indicates if there is at least one active region detected
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

td = trk_sw*(-1); Time direction

;Tracking
n = 1
trk_sw = 1;

sun_data = get_sun(date,carr=carr,he_lon=he_lon)
dtim = td*n*24.0*60.0*60.0
snod_A_raw = 2.8653  ; carrington rot. coeffs, in microrad. 
            
while (trk_sw eq 1) do begin

    ;finding Available PNRs in the previous day
    pr_ix = where((PRs.mdi_i eq mdi_i-td*n) and (PRs.lnk_sw eq 0), n_regp)
    if (n_regp gt 0) then begin
        pregions = PRs[pr_ix]
    endif
    
    nr_ix = where((NRs.mdi_i eq mdi_i-td*n) and (NRs.lnk_sw eq 0), n_regn)
    if (n_regn gt 0) then begin
        nregions = NRs[nr_ix]
    endif
    
    LonpAR = tARs.fcn_lnp
    LonnAR = tARs.fcn_lnn
    
    if ( (he_lon ge 90.0) and (he_lon lt 270.0) ) then begin
                
      tmpin = where(LonpAR lt 0, n_tmp)
      if (n_tmp gt 0) then LonpAR[tmpin] = LonpAR[tmpin]+360.0

      tmpin = where(LonnAR lt 0, n_tmp)
      if (n_tmp gt 0) then LonnAR[tmpin] = LonnAR[tmpin]+360.0
      
      hm_sw = 0
            
    endif
    
    dtim = 24.0*60.0*60.0*n*td;        

    if ((n_regp ne 0) and (n_regn ne 0) ) then begin
  
        timep = pregions[0].date
        timec = date
        
        
        LatpAR = tARs.fcn_ltp        
        LonpAR = tARs.fcn_lnp
        
        LatnAR = tARs.fcn_ltn        
        LonnAR = tARs.fcn_lnn
                        
        lat2p = pregions.fcn_lt
        lon2p = pregions.fcn_ln
        
        lat2n = nregions.fcn_lt
        lon2n = nregions.fcn_ln

        ;Making Sure all visible longitudes are in the same hemisphere
        ;---------------------------------------------------------
    
        if ( (he_lon ge 90.0) and (he_lon lt 270.0) ) then begin
          
          tmpin = where(lon2p lt 0, n_tmp)
          if (n_tmp gt 0) then lon2p[tmpin] = lon2p[tmpin]+360.0
    
          tmpin = where(lon2n lt 0, n_tmp)
          if (n_tmp gt 0) then lon2n[tmpin] = lon2n[tmpin]+360.0
          
          tmpin = where(LonpAR lt 0, n_tmp)
          if (n_tmp gt 0) then LonpAR[tmpin] = LonpAR[tmpin]+360.0
    
          tmpin = where(LonnAR lt 0, n_tmp)
          if (n_tmp gt 0) then LonnAR[tmpin] = LonnAR[tmpin]+360.0
          
          hm_sw = 0
                
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
        
        LonpAR = LonpAR - omegap*dtim       
        
        
        omegan = 1e-6*(snod_A + $ ; Omega at each pixel's lat., in radians
                      snod_B*sin(abs(LatnAR)*!dtor)^2 + $
                      snod_C*sin(abs(LatnAR)*!dtor)^4 )/!dtor
        
        LonnAR = LonnAR - omegan*dtim       
        

        lat2p = lat2p*!dtor
        lon2p = lon2p*!dtor
        
        lat2n = lat2n*!dtor
        lon2n = lon2n*!dtor
        
        sw_tr = 1; Switch which idicates whereas matching positive and negative regions have been found   
    
        ;Positive Side
        lat1 = LatpAR*!dtor
        lon1 = LonpAR*!dtor
                            
        dlat = abs(lat1-lat2p)
        dlon = abs(lon1-lon2p)
    
        dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2p)*sin(dlon/2.0)^2.0 ))/!dtor                    
        mindis = min( dis_rg, minin, /nan )
        
        if (mindis le min([pregions[minin].dcen, tARs.dcenp])*ar_cnst.dis_lim2) then begin
            jj = minin
        endif else begin
            sw_tr = 0; No matching positive regions were found
        endelse
        
    
        ;Negative Side
        lat1 = LatnAR*!dtor
        lon1 = LonnAR*!dtor
                            
        dlat = abs(lat1-lat2n)
        dlon = abs(lon1-lon2n)
    
        dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2n)*sin(dlon/2.0)^2.0 ))/!dtor                    
        mindis = min( dis_rg, minin, /nan )
        
        if (mindis le min([nregions[minin].dcen, tARs.dcenn])*ar_cnst.dis_lim2) then begin
            ii = minin
        endif else begin
            sw_tr = 0; No matching negative regions were found
        endelse
        
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
            ar = {ar, mdi_i: mdi_i-td*n, date: timep, labl: tARs.labl, clr: tARs.clr, indxp: pregions[jj].indx, indxn: nregions[ii].indx, fluxp: pregions[jj].flux, fluxn: nregions[ii].flux, areap:pregions[jj].area, arean:nregions[ii].area, $
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

        
            ;Removing paired regions
            PRs[pr_ix[jj]].lnk_sw = 1
            PRs[pr_ix[jj]].ar_lbl = tARs.labl
            
            NRs[nr_ix[ii]].lnk_sw = 1
            NRs[nr_ix[ii]].ar_lbl = tARs.labl
                                            
            ;Adding region
        
            arin  = where(ARs.mdi_i le mdi_i-td*n, nars)
            arin2 = where(ARs.mdi_i gt mdi_i-td*n, nars2)
           
            if ( (nars eq 0) and (nars2 eq 0) ) then begin
              ARs = ar
            endif else if (nars eq 0) then begin
              ARs = [ar,ARs[arin2]]
            endif else if (nars2 eq 0) then begin
              ARs = [ARs[arin],ar]
            endif else begin
              ARs = [ARs[arin],ar,ARs[arin2]]
            endelse


        
        endif
              
    endif

    if ( ( (he_lon ge  90.0) and (he_lon lt 270.0) ) and ( td*(max([LonnAR, LonpAR]) - 1e-6*dtim*snod_A_raw/!dtor ) lt td*( he_lon - (td*90.0) ) ) ) then begin
        trk_sw = 0        
    endif
    
    if ( ( (he_lon ge 270.0) and (he_lon le 360.0) ) and ( td*(max([LonnAR, LonpAR]) - 1e-6*dtim*snod_A_raw/!dtor ) lt td*( he_lon - (360.0 + td*90 ) ) ) ) then begin
        trk_sw = 0        
    endif

    if ( ( (he_lon ge   0.0) and (he_lon lt  90.0) ) and ( td*(max([LonnAR, LonpAR]) - 1e-6*dtim*snod_A_raw/!dtor ) lt td*( he_lon - (td*90.0) ) ) ) then begin
        trk_sw = 0        
    endif
    
    n = n+1

endwhile


return
end