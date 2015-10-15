PRO amj_pnr_post_merge,  CRD_in, mdi_i, lbl, PRs, NRs, AR_c, seg_const=seg_const


;finding Available PRs
pr_ix = where((PRs.mdi_i eq mdi_i), n_regpt)
if (n_regpt gt 0) then begin
    pregions = PRs[pr_ix]
endif

;finding unlinked PRs
if (n_regpt gt 0) then begin
    pr_ix = where((pregions.lnk_sw eq 0), n_regp)
endif else begin
    n_regp = 0
endelse 

;finding Available NRs
nr_ix = where((NRs.mdi_i eq mdi_i), n_regnt)
if (n_regnt gt 0) then begin
    nregions = NRs[nr_ix]
endif

;finding unlinked NRs
if (n_regnt gt 0) then begin
    nr_ix = where((nregions.lnk_sw eq 0), n_regn)
endif else begin
    n_regn = 0
endelse



ar_ix = where((AR_c.mdi_i eq mdi_i) , n_ars)
if (n_ars gt 0) then begin
    ARs = AR_c[ar_ix]
endif





;Enlarging imbalanced ARs-------------------------------------------------------------------------------------------
if (n_ars gt 0) then begin

    ;Merging distance
    n_dis = seg_const.dis_lim


    ;Creating Overlap Matrices---------------------------------------------------------------------------------------
    ;
    ;
    ;Positive regions
    if (n_regp gt 1) then begin
      
        ;Creating Overlapping Matrix
        OvMp = fltarr(n_regp,n_ars) ; Overlap matrix
        
        for i=0,n_regp-1 do begin    
      
            for j=0,n_ars-1 do begin    
                
                ;Calculating distances between its centroids
                lat1 = pregions[pr_ix[i]].fcn_lt*!dtor
                lon1 = pregions[pr_ix[i]].fcn_ln*!dtor
                
                lat2 = ARs[j].fcn_ltp*!dtor
                lon2 = ARs[j].fcn_lnp*!dtor
                
                dlat = abs(lat1-lat2)
                dlon = abs(lon1-lon2)
        
                dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
                                
                ;Finding regions withing prescribed distance
                if ( dis_rg lt n_dis*(pregions[pr_ix[i]].dcen + ARs[j].dcenp) ) then begin
                  
                    r1 = n_dis*pregions[pr_ix[i]].dcen 
                    r2 = n_dis*ARs[j].dcenp 

                    ;If the region is contained inside another
                    if ( ( dis_rg + r1 ) lt r2 ) then begin
                      
                        OvMp(i,j) = 1.0;
                      
                    endif
                    
                    ;If the region contains another
                    if ( ( dis_rg + r2 ) lt r1 ) then begin
                      
                        OvMp(i,j) = (!pi*r2^2.0)/(!pi*r1^2.0);
                      
                    endif

                    ;Make sure neither region is contained inside the other                   
                    if ( ( dis_rg + r1 ge r2 ) and ( dis_rg + r2 ge r1 ) ) then begin
                                      
                        Aint =   r1^2.0*acos( ( dis_rg^2.0 + r1^2.0 - r2^2.0 )/(2.0*dis_rg*r1) ) $ 
                               + r2^2.0*acos( ( dis_rg^2.0 - r1^2.0 + r2^2.0 )/(2.0*dis_rg*r2) ) $
                               - 0.5*sqrt( ( - dis_rg + r1 + r2 )*( dis_rg - r1 + r2 )*( dis_rg + r1 - r2 )*( dis_rg + r1 + r2 ) )
                        
                        OvMp(i,j) = Aint/(!pi*r1^2.0) 

                    endif                    
                  
                endif else begin
                  
                   OvMp(i,j) = !values.f_nan
                   
                endelse
          
            endfor        
      
        endfor
            
        ;Keeping those above aggregation overlap
        tmp_in = where(OvMp lt seg_const.ovr_lim, n_below)
        if (n_below ne 0) then OvMp[where(OvMp lt seg_const.ovr_lim)] = !values.f_nan       

    endif



    ;Negative regions
    if (n_regn gt 1) then begin
      
        ;Creating Overlapping Matrix
        OvMn = fltarr(n_regn,n_ars) ; Overlap matrix
        
        for i=0,n_regn-1 do begin    
      
            for j=0,n_ars-1 do begin    
                
                ;Calculating distances between its centroids
                lat1 = nregions[nr_ix[i]].fcn_lt*!dtor
                lon1 = nregions[nr_ix[i]].fcn_ln*!dtor
                
                lat2 = ARs[j].fcn_ltn*!dtor
                lon2 = ARs[j].fcn_lnn*!dtor
                
                dlat = abs(lat1-lat2)
                dlon = abs(lon1-lon2)
        
                dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor
                                
                ;Finding regions withing prescribed distance
                if ( dis_rg lt n_dis*(nregions[nr_ix[i]].dcen + ARs[j].dcenn) ) then begin
                  
                    r1 = n_dis*nregions[nr_ix[i]].dcen 
                    r2 = n_dis*ARs[j].dcenn 

                    ;If the region is contained inside another
                    if ( ( dis_rg + r1 ) lt r2 ) then begin
                      
                        OvMn(i,j) = 1.0;
                      
                    endif
                    
                    ;If the region contains another
                    if ( ( dis_rg + r2 ) lt r1 ) then begin
                      
                        OvMn(i,j) = (!pi*r2^2.0)/(!pi*r1^2.0);
                      
                    endif

                    ;Make sure neither region is contained inside the other                   
                    if ( ( dis_rg + r1 ge r2 ) and ( dis_rg + r2 ge r1 ) ) then begin
                                      
                        Aint =   r1^2.0*acos( ( dis_rg^2.0 + r1^2.0 - r2^2.0 )/(2.0*dis_rg*r1) ) $ 
                               + r2^2.0*acos( ( dis_rg^2.0 - r1^2.0 + r2^2.0 )/(2.0*dis_rg*r2) ) $
                               - 0.5*sqrt( ( - dis_rg + r1 + r2 )*( dis_rg - r1 + r2 )*( dis_rg + r1 - r2 )*( dis_rg + r1 + r2 ) )
                        
                        OvMn(i,j) = Aint/(!pi*r1^2.0) 

                    endif                    
                  
                endif else begin
                  
                   OvMn(i,j) = !values.f_nan
                   
                endelse
          
            endfor        
      
        endfor
            
        ;Keeping those above aggregation overlap
        tmp_in = where(OvMn lt seg_const.ovr_lim, n_below)
        if (n_below ne 0) then OvMn[where(OvMn lt seg_const.ovr_lim)] = !values.f_nan       

    endif


    ;Enlarging imbalanced ARs-------------------------------------------------------------------------------------------
  
    Flx_im = ARs.fluxp + ARs.fluxn
    Flxmin = min(abs([[ARs.fluxp], [ARs.fluxn]]), dimension = 2)
    
    Flx_imr = Flx_im/Flxmin
;    lstinx  = !values.f_nan  ;Last visited index in terms of flux imbalance

    while ( (total(finite(Flx_imr)) ne 0) or (total(finite(OvMp)) ne 0) or (total(finite(OvMn)) ne 0)) do begin
      
        ;Finding the largest imbalance
        maxim  = max(abs(Flx_imr),jj,/nan)
        
        ;Negative Excess flux
        if (Flx_imr[jj] lt 0) then begin
        
            ;Finding largest overlap
            mdis = max(OvMp(*,jj),ii,/nan)
            
            ; Checking if there is a valid overlapping region
            if finite(mdis) then begin
                
                ;Finding linked region
                lnk_in = where(ARs[jj].labl eq pregions.ar_lbl)
              
                ;Extracting indices of Growing region
                inpg = long(strsplit(pregions[lnk_in].indx,/extract))
                
                ;Extracting indices of Killed region
                inpk = long(strsplit(pregions[pr_ix[ii]].indx,/extract))
                
                ;Enlargement
                ind_gn= [inpg, inpk]
                  
                ;Recalculating parameters----------------------------------------------------
                flux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)       ;total flux
                area  = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)        ;total area
                
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
                tmp_2d = array_indices(CRD_in.mgnt_flx,minin)
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
                pregions[lnk_in].indx = strjoin(string(ind_gn))                    
                
                ;Updating relevant quantities
                tmp_r = {flux: flux, area: area, $
                          fcn_lt: fcenlt, fcn_ln: fcenln, dcen: dcen, $ 
                          fcenxp: fcenxp, fcenyp: fcenyp, dcenp: dcenp, dm:dm, qm:qm, fr_lbl: -999} 
                               
                prtmp = pregions[lnk_in]                        
                STRUCT_ASSIGN, tmp_r, prtmp, /nozero
                 
                pregions[lnk_in] = prtmp  


                ;Updating AR---------------------------------------------------------------------------------------
                ;Adding region indexes
                ARs[jj].indxp = strjoin(string(ind_gn))                    
            
                ;Making Sure all visible longitudes are in the same hemisphere
                ;---------------------------------------------------------
                sun_data = get_sun(date,carr=carr,he_lon=he_lon)
            
                Lonp = fcenln
                Lonn = ARs[jj].fcn_lnn
                                        
                if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
                        
                  tmpin = where(Lonp lt 0, n_tmp)
                  if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
                        
                  tmpin = where(Lonn lt 0, n_tmp)
                  if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
                  
                  he_lon2 = he_lon
                        
                endif 
                
                
                ;Finding leading polarity
                if ( Lonp ge Lonn ) then begin
                
                    LatL  = fcenlt*!dtor
                    LongL = Lonp*!dtor
            
                    LatF  = ARs[jj].fcn_ltn*!dtor
                    LongF = Lonn*!dtor
                    
                    lp = 1;
                
                endif else begin
            
                    LatF  = fcenlt*!dtor
                    LongF = Lonp*!dtor
            
                    LatL  = ARs[jj].fcn_ltn*!dtor
                    LongL = Lonn*!dtor
                    
                    lp = -1;
                
                endelse
                        
                dlat = abs(LatF-LatL)
                dlon = abs(LongF-LongL)
            
                Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
                Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
                
                dm = 0
                qm = 0
                
                ;Updating relevant quantities
                tmp_ar = {fluxp: flux, areap: area, $
                          fcn_ltp: fcenlt, fcn_lnp: fcenln, dcenp: dcen, $ 
                          fcenxpp: fcenxp, fcenypp: fcenyp, dcenpp: dcenp, $
                          dis: Dis, tilt: Tilt, lp: lp, dm:dm, qm:qm} 
                               
                ar = ARs[jj]                        
                STRUCT_ASSIGN, tmp_ar, ar, /nozero 
                ARs[jj] = ar

                ;Removing pair from the matrix
                Keep_inx = findgen(n_regp)
                OvMp = OvMp[where(Keep_inx ne ii),*]
                OvMp = OvMp[*,where(Keep_inx ne ii)]
                
                Keep_inx  = findgen(n_regpt)
                pregions  = pregions[Keep_inx ne pr_ix[ii])]
                n_regpt   = n_elements(pregions)
                pr_ix     = where((pregions.lnk_sw eq 0), n_regp)
                
                ;Updating flux imbalance
                Flx_im[jj] = ARs[jj].fluxp + ARs[jj].fluxn
                if (Flxmin gt abs(Flx_im[jj])) then begin
                    Flxmin = Flx_im[jj]
                    Flx_imr = Flx_im/Flxmin
                endif                
                
            ;No available regions to merge
            endif else begin
                Flx_imr[jj] = !values.f_nan              
            endelse        
        
        ;Positive Excess flux    
        endif else begin
              
      
        endelse
      
    endwhile

  
endif






if (n_ars gt 0) then begin

    Flx_im = ARs.fluxp + ARs.fluxn
    Flxmin = min(abs([[ARs.fluxp], [ARs.fluxn]]), dimension = 2)
    
    Flx_imr = Flx_im/Flxmin
    if keyword_set(info) then print,'FLUX:'
    if keyword_set(info) then print, Flx_im
    if keyword_set(info) then print, Flx_imr
    
    im_in = where( abs(Flx_imr) le ar_cnst.Imb_tol , nbal)    
    if (nbal gt 0) then Flx_imr[im_in] = !values.f_nan
    
    if keyword_set(info) then print, Flx_imr
    
    for i=0,n_ars-2-nbal do begin
    
        ;Finding the largest imbalance
        maxim = max(abs(Flx_imr),maximn,/nan)
        nflxim = Flx_imr[maximn]
        
        ;Extracting indices
        inp = long(strsplit(ARs[maximn].indxp,/extract))
        inn = long(strsplit(ARs[maximn].indxn,/extract))
        
        ;Adding field to the unused magnetogram in order to perform region growth
        imsel[inp] =   seg_const.ar_th
        imsel[inn] = - seg_const.ar_th
        
        ;Storing variables to keep track of the best iteration
        bstinp = inp
        bstinn = inn
        bstflxim = nflxim
        
        ;Iterative enlargement
;        print, nflxim
        n = 1
        Th_min = seg_const.ar_th/2;   
        Th_max = max(abs(ar_cnst.valid_range))  
        while (abs(nflxim) gt ar_cnst.Imb_tol) and (n le ar_cnst.Imb_it) do begin
          
            ;Enlarging the positive region
            if Flx_imr[maximn] lt 0.0 then begin
                
                ind_gn=idl_region_grow(imsel,inp,thresh=[Th_min,Th_max])
                ;Calculating new flux
                nflux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)
                nflxim = (ARs[maximn].fluxn + nflux)/nflux 
                area  = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)

                ;comparing performance and ensuring size is not too large
                if (abs(nflxim) lt abs(bstflxim)) and (area le 2.0*ARs[maximn].arean) then begin
                    bstflxim = nflxim
                    bstinp = ind_gn
                endif                
              
            ;Enlarging the negative region                
            endif else begin

                ind_gn=idl_region_grow(-imsel,inn,thresh=[Th_min,Th_max])
                ;Calculating new flux
                nflux  = total(CRD_in.mgnt_flx[ind_gn],/double, /nan)
                nflxim = (ARs[maximn].fluxp + nflux)/nflux 
                area  = total(CRD_in.mgnt_ar[ind_gn],/double, /nan)

                ;comparing performance and ensuring size is not too large
                if (abs(nflxim) lt abs(bstflxim)) and (area le 2.0*ARs[maximn].areap) then begin
                    bstflxim = nflxim
                    bstinn = ind_gn
                endif
                                
            endelse
            
;            print, nflxim
              
            ;Performing binary search
            if abs(nflxim) gt 1.0 then begin
                Th_min = Th_min - seg_const.ar_th/(2.0^(n+1))
            endif else begin
                Th_min = Th_min + seg_const.ar_th/(2.0^(n+1))
            endelse
                         
            n = n+1          
          
        endwhile
        
        ;recovering best indices
        if Flx_imr[maximn] lt 0.0 then begin
            ind_gn = bstinp
        endif else begin
            ind_gn = bstinn          
        endelse
;        print, bstflxim
        
        ;Recalculating parameters
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
        
                        
        ;Storing positive region
        if Flx_imr[maximn] lt 0.0 then begin
                      
            ;Adding region indexes
            ARs[maximn].indxp = strjoin(string(ind_gn))                    

            ;Making Sure all visible longitudes are in the same hemisphere
            ;---------------------------------------------------------            
            Lonp = fcenln
            Lonn = ARs[maximn].fcn_lnn
                                    
            if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
                    
              tmpin = where(Lonp lt 0, n_tmp)
              if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
                    
              tmpin = where(Lonn lt 0, n_tmp)
              if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
              
              he_lon2 = he_lon
                    
            endif 
            
            
            ;Finding leading polarity
            if ( Lonp ge Lonn ) then begin
            
                LatL  = fcenlt*!dtor
                LongL = Lonp*!dtor
        
                LatF  = ARs[maximn].fcn_ltn*!dtor
                LongF = Lonn*!dtor
                
                lp = 1;
            
            endif else begin
        
                LatF  = fcenlt*!dtor
                LongF = Lonp*!dtor
        
                LatL  = ARs[maximn].fcn_ltn*!dtor
                LongL = Lonn*!dtor
                
                lp = -1;
            
            endelse
                    
            dlat = abs(LatF-LatL)
            dlon = abs(LongF-LongL)
        
            Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
            Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
            
            ;Updating relevant quantities
            tmp_ar = {fluxp: flux, areap: area, $
                      fcn_ltp: fcenlt, fcn_lnp: fcenln, dcenp: dcen, $ 
                      fcenxpp: fcenxp, fcenypp: fcenyp, dcenpp: dcenp, $
                      dis: Dis, tilt: Tilt, lp: lp} 
                           
            ar = ARs[maximn]                        
            STRUCT_ASSIGN, tmp_ar, ar, /nozero 
            ARs[maximn] = ar  
                             
            ;Updating indices to reset unused magnetogram
            inp = ind_gn
                              
        ;Storing negative region                
        endif else begin
          
            ;Adding region indexes
            ARs[maximn].indxn = strjoin(string(ind_gn))                    

            ;Making Sure all visible longitudes are in the same hemisphere
            ;---------------------------------------------------------            
            Lonp = ARs[maximn].fcn_lnp
            Lonn = fcenlt
                                    
            if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
                    
              tmpin = where(Lonp lt 0, n_tmp)
              if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
                    
              tmpin = where(Lonn lt 0, n_tmp)
              if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
              
              he_lon2 = he_lon
                    
            endif 
            
            
            ;Finding leading polarity
            if ( Lonp ge Lonn ) then begin
            
                LatL  = ARs[maximn].fcn_ltp*!dtor
                LongL = Lonp*!dtor
        
                LatF  = fcenlt*!dtor
                LongF = Lonn*!dtor
                
                lp = 1;
            
            endif else begin
        
                LatF  = ARs[maximn].fcn_ltp*!dtor
                LongF = Lonp*!dtor
        
                LatL  = fcenlt*!dtor
                LongL = Lonn*!dtor
                
                lp = -1;
            
            endelse
                    
            dlat = abs(LatF-LatL)
            dlon = abs(LongF-LongL)
        
            Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
            Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt
            
            ;Updating relevant quantities
            tmp_ar = {fluxn: flux, arean: area, $
                      fcn_ltn: fcenlt, fcn_lnn: fcenln, dcenn: dcen, $ 
                      fcenxpn: fcenxp, fcenypn: fcenyp, dcenpn: dcenp, $
                      dis: Dis, tilt: Tilt, lp: lp} 
                           
            ar = ARs[maximn]                        
            STRUCT_ASSIGN, tmp_ar, ar, /nozero 
            ARs[maximn] = ar
            
            ;Updating indices to reset unused magnetogram
            inn = ind_gn
                                                   
        endelse
      
        ;Removing region from the unused magnetogram
        imsel[inp] = 0.0
        imsel[inn] = 0.0
        
        
        Flx_imr[maximn] = !values.f_nan
;        print, Flx_imr

    endfor
    
endif











return
END
