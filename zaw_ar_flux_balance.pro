;-------------------------------------------------------------------------------------------------------------------------------
; NAME:       zaw_ar_flux_balance
;
; PURPOSE:      Takes the active regions (ARs) as parameters and balances the flux by growing the regions in an
;           iterative fashion.
;
; CALLING SEQUENCE: zaw_ar_flux_balance,  
;
;
;
;-------------------------------------------------------------------------------------------------------------------------------
pro zaw_ar_flux_balance, CRD_in, mdi_i, ARs, n_ars, date, ar_cnst=ar_cnst, seg_const=seg_const, display=display, prt=prt, info=info, sqs_nm = sqs_nm, pltn = pltn

if (n_ars ne 0) then begin
  
  im=CRD_in.im_raw
  imsel=CRD_in.im_crr
  
  ;Clearing magnetogram of detected regions
  for i=0,n_ars-1 do begin

    inp = long(strsplit(ARs[i].indxp,/extract))
    inn = long(strsplit(ARs[i].indxn,/extract))

    imsel[inp] = 0.0
    imsel[inn] = 0.0      

  endfor

  ;Calculating flux imbalance

  Flx_im = ARs.fluxp + ARs.fluxn
  print, Flx_im[0]
  print, "Break"

  Flxmin = min(abs([[ARs.fluxp], [ARs.fluxn]]), dimension = 2)

  ;Relative flux imbalance per active region
  Flx_imr = Flx_im/Flxmin
  if keyword_set(info) then print,'FLUX:'
  if keyword_set(info) then print, Flx_im
  if keyword_set(info) then print, Flx_imr

  ;Fail safe
  im_in = where( abs(Flx_imr) le ar_cnst.Imb_tol , nbal)    
  if (nbal gt 0) then Flx_imr[im_in] = !values.f_nan

  if keyword_set(info) then print, Flx_imr

  for i=0,n_ars-2-nbal do begin
   
    ;Finding the largest imbalance
    maxim = max(abs(Flx_imr),maximn,/nan)
    nflxim = Flx_imr[maximn]

    ;Extracting pixel indices
    inp = long(strsplit(ARs[maximn].indxp,/extract))
    inn = long(strsplit(ARs[maximn].indxn,/extract))

    ;Adding field to the unused magnetogram in order to perform region growth
    imsel[inp] =   seg_const.ar_th
    imsel[inn] = - seg_const.ar_th

    ;Storing variables to keep track of the best iteration
    bstinp = inp
    bstinn = inn
    bstflxim = nflxim

    n = 1
    ;Defining lower and higher thresholds
    Th_minL = seg_const.ar_th*4.0/6.0
    Th_minH = seg_const.ar_th*5.0/6.0    
    Th_max = max(abs(ar_cnst.valid_range))

    ;Iterative enlargement------------------------------------------------------------------------------
    while (abs(nflxim) gt ar_cnst.Imb_tol) and (n le ar_cnst.Imb_it) do begin
     
      ;Trying lowest growth threshold
      ind_gpL = idl_region_grow(imsel,inp,thresh=[Th_minL,Th_max])
      ind_gnL = idl_region_grow(-imsel,inn,thresh=[Th_minL,Th_max])

      nfluxpL  = total(CRD_in.mgnt_flx[ind_gpL],/double, /nan)
      nfluxnL  = total(CRD_in.mgnt_flx[ind_gnL],/double, /nan)

      FlxminL = min(abs([[nfluxpL], [nfluxnL]]), dimension = 2)
      nflximL = (nfluxpL + nfluxnL)/FlxminL

      areapL  = total(CRD_in.mgnt_ar[ind_gpL],/double, /nan) 
      areanL  = total(CRD_in.mgnt_ar[ind_gnL],/double, /nan) 


      ;Trying highest growth threshold
      ind_gpH = idl_region_grow(imsel,inp,thresh=[Th_minH,Th_max])
      ind_gnH = idl_region_grow(-imsel,inn,thresh=[Th_minH,Th_max])

      nfluxpH  = total(CRD_in.mgnt_flx[ind_gpH],/double, /nan)
      nfluxnH  = total(CRD_in.mgnt_flx[ind_gnH],/double, /nan)

      FlxminH = min(abs([[nfluxpH], [nfluxnH]]), dimension = 2)
      nflximH = (nfluxpH + nfluxnH)/FlxminH

      areapH  = total(CRD_in.mgnt_ar[ind_gpH],/double, /nan) 
      areanH  = total(CRD_in.mgnt_ar[ind_gnH],/double, /nan)           

      ;Checking if either of the new thresholds improves flux balance while keeping the areas within reasonable size
      if ((abs(nfluximL) lt bstflxim) and (areapL le 2.0*ARs[maximn].areap) and (areanL le 2.0*ARs[maximn].arean) ) or ( (abs(nfluximH) lt bstflxim) and (areapH le 2.0*ARs[maximn].areap) and (areanH le 2.0*ARs[maximn].arean) )then begin

        ;Checking if the lowest threshold is better than the highest threshold given a certain tolerance
        ; if (abs(nfluximL) lt (1.0-2.0*ar_cnst.Imb_tol)*nfluximH)) then begin
        if (abs(nfluximL) lt abs(nfluximH)) then begin
          
          ;Updating best set
          bstflxim = nflximL
          bstinp = ind_gpL
          bstinn = ind_gnL

          ;Defining new thresholds
          Th_minL = Th_minL - seg_const.ar_th/(6.0*2.0^n)
          Th_minH = Th_minL + seg_const.ar_th/(6.0*2.0^n)               

        ;If not use the highest threshold
        endif else begin

          ;Updating best set
          bstflxim = nflximH
          bstinp = ind_gpH
          bstinn = ind_gnH
          ;Defining new thresholds
          Th_minL = Th_minH - seg_const.ar_th/(6.0*2.0^n)
          Th_minH = Th_minH + seg_const.ar_th/(6.0*2.0^n)  

        endelse

        n = n+1

      ;Leaving search loop because there was no improvement in flux balance
      endif else begin
        n = ar_cnst.Imb_it+1
      endelse

    endwhile
   
    ;Recalculating parameters for positive region------------------------------------------------------

    fluxp  = total(CRD_in.mgnt_flx[bstinp],/double, /nan)       ;total flux
    areap  = total(CRD_in.mgnt_ar[bstinp],/double, /nan)         ;total area

    ;Geographical Center
    fcenxp = total(CRD_in.Xar[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flux
    fcenyp = total(CRD_in.Yar[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flux
    fcenzp = total(CRD_in.Zar[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flux

    fcenltp = atan(fcenzp,sqrt( fcenxp^2.0 + fcenyp^2.0 ) )/!dtor    ;Latitude center  
    fcenlnp = atan(fcenyp,fcenxp)/!dtor                             ;Longitude center

    ;Center pixels
    lat1 = fcenltp*!dtor
    lon1 = fcenlnp*!dtor

    lat2 = CRD_in.Lath*!dtor
    lon2 = CRD_in.Lonh*!dtor

    dlat = abs(lat1-lat2)
    dlon = abs(lon1-lon2)

    dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor

    mindis = min( dis_rg, minin, /nan )
    tmp_2d = array_indices(im,minin)
    fcenxpp = tmp_2d[0,0];    ;Center in x
    fcenypp = tmp_2d[1,0];    ;Center in y

    ;Mean and Pixel Radius
    dcenp = mean( dis_rg[bstinp],/double, /nan)
    dcenfwp = total( dis_rg[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flux

    dcenpp    = dcen/mean([ dis_rg[fcenxpp-1,fcenypp], dis_rg[fcenxpp+1,fcenypp], dis_rg[fcenxpp,fcenypp-1], dis_rg[fcenxpp,fcenypp+1] ],/double, /nan);/!dtor            
    dcenpfwp  = dcenfw/mean([ dis_rg[fcenxpp-1,fcenypp], dis_rg[fcenxpp+1,fcenypp], dis_rg[fcenxpp,fcenypp-1], dis_rg[fcenxpp,fcenypp+1] ],/double, /nan);/!dtor
     

    ;Recalculating parameters for negative region------------------------------------------------------
    fluxn  = total(CRD_in.mgnt_flx[bstinn],/double, /nan)       ;total flux
    arean  = total(CRD_in.mgnt_ar[bstinn],/double, /nan)         ;total area

    ;Geographical Center
    fcenxn = total(CRD_in.Xar[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flux
    fcenyn = total(CRD_in.Yar[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flux
    fcenzn = total(CRD_in.Zar[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flux

    fcenltn = atan(fcenzn,sqrt( fcenxn^2.0 + fcenyn^2.0 ) )/!dtor    ;Latitude center  
    fcenlnn = atan(fcenyn,fcenxn)/!dtor                             ;Longitude center

    ;Center pixels
    lat1 = fcenltn*!dtor
    lon1 = fcenlnn*!dtor

    lat2 = CRD_in.Lath*!dtor
    lon2 = CRD_in.Lonh*!dtor

    dlat = abs(lat1-lat2)
    dlon = abs(lon1-lon2)

    dis_rg = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(lat1)*cos(lat2)*sin(dlon/2.0)^2.0 ))/!dtor

    mindis = min( dis_rg, minin, /nan )
    tmp_2d = array_indices(im,minin)
    fcenxpn = tmp_2d[0,0];    ;Center in x
    fcenypn = tmp_2d[1,0];    ;Center in y

    ;Mean and Pixel Radius
    dcenn = mean( dis_rg[bstinn],/double, /nan)
    dcenfwn = total( dis_rg[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flux

    dcenpn    = dcen/mean([ dis_rg[fcenxpn-1,fcenypn], dis_rg[fcenxpn+1,fcenypn], dis_rg[fcenxpn,fcenypn-1], dis_rg[fcenxpn,fcenypn+1] ],/double, /nan);/!dtor            
    dcenpfwn  = dcenfw/mean([ dis_rg[fcenxpn-1,fcenypn], dis_rg[fcenxpn+1,fcenypn], dis_rg[fcenxpn,fcenypn-1], dis_rg[fcenxpn,fcenypn+1] ],/double, /nan);/!dtor           


    ;Storing active region--------------------------------------------------------------------------

    ;Adding region indexes
    ARs[maximn].indxp = strjoin(string(bstinp)) 
    ARs[maximn].indxn = strjoin(string(bstinn)) 

    ;Making Sure all visible longitudes are in the same hemisphere
    ;---------------------------------------------------------            
    Lonp = fcenlnp
    Lonn = fcenlnn

    sun_data = get_sun(date,carr=carr,he_lon=he_lon)
                           
    if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin

      if (Lonp lt 0) then Lonp = Lonp+360.0
      if (Lonn lt 0) then Lonn = Lonn+360.0
           
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
    tmp_ar = {fluxp: fluxp, areap: areap, $
             fcn_ltp: fcenltp, fcn_lnp: fcenlnp, dcenp: dcenp, $ 
             fcenxpp: fcenxpp, fcenypp: fcenypp, dcenpp: dcenpp, $
             fluxn: fluxn, arean: arean, $
             fcn_ltn: fcenltn, fcn_lnn: fcenlnn, dcenn: dcenn, $ 
             fcenxpn: fcenxpn, fcenypn: fcenypn, dcenpn: dcenpn, $
             dis: Dis, tilt: Tilt, lp: lp} 
                  
    ar = ARs[maximn]                        
    STRUCT_ASSIGN, tmp_ar, ar, /nozero 
    ARs[maximn] = ar  
                    
    ;Updating indices to reset unused magnetogram
    inp = ind_gn

    ;Removing region from the unused magnetogram
    imsel[bstinp] = 0.0
    imsel[bstinn] = 0.0


    Flx_imr[maximn] = !values.f_nan
    tempflxim = tmp_ar.fluxp + tmp_ar.fluxn
    print, tempflxim
  endfor
   
endif else begin

  print, 'No valid ARs in call to zaw_ar_flux_balance.pro'

endelse


end
