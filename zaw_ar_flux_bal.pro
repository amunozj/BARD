;-------------------------------------------------------------------------------------------------------------------------------
; NAME:       zaw_ar_flx_balance
;
; PURPOSE:      Takes the active regions (ARs) as parameters and balances the flx by growing the regions in an
;           iterative fashion.
;
; CALLING SEQUENCE: zaw_ar_flx_balance,  
;
;
;
;-------------------------------------------------------------------------------------------------------------------------------
pro zaw_ar_flux_bal, CRD_in, mdi_i, ARs, ar_in, n_ars, date, ar_cnst=ar_cnst, seg_const=seg_const, display=display, prt=prt, info=info, sqs_nm = sqs_nm, pltn = pltn

tempAR = ARs[ar_in]

im=CRD_in.im_raw
imsel=CRD_in.im_crr

bal_count = 0

;Clearing magnetogram of detected regions
for i=0,n_ars-1 do begin

  inp = long(strsplit(tempAR[i].indxp,/extract))
  inn = long(strsplit(tempAR[i].indxn,/extract))

  imsel[inp] = 0.0
  imsel[inn] = 0.0      

endfor

;Calculating flux imbalance
Flx_im = tempAR.fluxp + tempAR.fluxn

Flxmin = min(abs([[tempAR.fluxp], [tempAR.fluxn]]), dimension = 2)

;Relative flx imbalance per active region
Flx_imr = Flx_im/Flxmin
if keyword_set(info) then print,'flx:'
if keyword_set(info) then print, Flx_im
if keyword_set(info) then print, Flx_imr

;Fail safe
im_in = where( abs(Flx_imr) le ar_cnst.Imb_tol , nbal)    
if (nbal gt 0) then Flx_imr[im_in] = !values.f_nan

if keyword_set(info) then print, Flx_imr

for i=0,n_ars-1-nbal do begin
 
  ;Finding the largest imbalance
  maxim = max(abs(Flx_imr),maxim_ind,/nan)
  nflxim = Flx_imr[maxim_ind]

  ;Extracting pixel indices
  inp = long(strsplit(tempAR[maxim_ind].indxp,/extract))
  inn = long(strsplit(tempAR[maxim_ind].indxn,/extract))

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
  while (abs(bstflxim) gt ar_cnst.Imb_tol) and (n le ar_cnst.Imb_it) do begin
   
    ;Trying lowest growth threshold
    ind_gpL = idl_region_grow(imsel,inp,thresh=[Th_minL,Th_max])
    ind_gnL = idl_region_grow(-imsel,inn,thresh=[Th_minL,Th_max])

    nflxpL  = total(CRD_in.mgnt_flx[ind_gpL],/double, /nan)
    nflxnL  = total(CRD_in.mgnt_flx[ind_gnL],/double, /nan)

    FlxminL = min(abs([[nflxpL], [nflxnL]]), dimension = 2)
    nflximL = (nflxpL + nflxnL)/FlxminL

    areapL  = total(CRD_in.mgnt_ar[ind_gpL],/double, /nan) 
    areanL  = total(CRD_in.mgnt_ar[ind_gnL],/double, /nan) 


    ;Trying highest growth threshold
    ind_gpH = idl_region_grow(imsel,inp,thresh=[Th_minH,Th_max])
    ind_gnH = idl_region_grow(-imsel,inn,thresh=[Th_minH,Th_max])

    nflxpH  = total(CRD_in.mgnt_flx[ind_gpH],/double, /nan)
    nflxnH  = total(CRD_in.mgnt_flx[ind_gnH],/double, /nan)

    FlxminH = min(abs([[nflxpH], [nflxnH]]), dimension = 2)
    nflximH = (nflxpH + nflxnH)/FlxminH

    areapH  = total(CRD_in.mgnt_ar[ind_gpH],/double, /nan) 
    areanH  = total(CRD_in.mgnt_ar[ind_gnH],/double, /nan)           

    ;Checking if either of the new thresholds improves flux balance while keeping the areas within reasonable size
    if ((abs(nflximL) lt abs(bstflxim) ) and (areapL le 2.0*tempAR[maxim_ind].areap) and (areanL le 2.0*tempAR[maxim_ind].arean) ) or ( (abs(nflximH) lt abs(bstflxim) ) and (areapH le 2.0*tempAR[maxim_ind].areap) and (areanH le 2.0*tempAR[maxim_ind].arean) )then begin

      ;Checking if the lowest threshold is better than the highest threshold given a certain tolerance
      ; if (abs(nflximL) lt (1.0-2.0*ar_cnst.Imb_tol)*nflximH)) then begin
      if (abs(nflximL) lt abs(nflximH)) then begin
        
        ;Updating best set
        bstflxim = nflximL
        bstinp = ind_gpL
        bstinn = ind_gnL
        bal_count = bal_count + 1

        ;Defining new thresholds
        Th_minL = Th_minL - seg_const.ar_th/(6.0*2.0^n)
        Th_minH = Th_minL + seg_const.ar_th/(6.0*2.0^n)               

      ;If not use the highest threshold
      endif else begin

        ;Updating best set
        bstflxim = nflximH
        bstinp = ind_gpH
        bstinn = ind_gnH
        bal_count = bal_count + 1

        ;Defining new thresholds
        Th_minL = Th_minH - seg_const.ar_th/(6.0*2.0^n)
        Th_minH = Th_minH + seg_const.ar_th/(6.0*2.0^n)  

      endelse

      n = n+1

    ;Leaving search loop because there was no improvement in flux balance
    endif else begin
      n = ar_cnst.Imb_it+1
    endelse
    ;print, bstflxim
  endwhile
 
  ;Recalculating parameters for positive region------------------------------------------------------

  flxp  = total(CRD_in.mgnt_flx[bstinp],/double, /nan)       ;total flux
  areap  = total(CRD_in.mgnt_ar[bstinp],/double, /nan)         ;total area

  ;Geographical Center
  fcenxp = total(CRD_in.Xar[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flxp
  fcenyp = total(CRD_in.Yar[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flxp
  fcenzp = total(CRD_in.Zar[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flxp

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
  temp_2d = array_indices(im,minin)
  fcenxpp = temp_2d[0,0];    ;Center in x
  fcenypp = temp_2d[1,0];    ;Center in y

  ;Mean and Pixel Radius
  dcenp = mean( dis_rg[bstinp],/double, /nan)
  dcenfwp = total( dis_rg[bstinp]*CRD_in.mgnt_flx[bstinp],/double, /nan)/flxp

  dcenpp    = dcenp/mean([ dis_rg[fcenxpp-1,fcenypp], dis_rg[fcenxpp+1,fcenypp], dis_rg[fcenxpp,fcenypp-1], dis_rg[fcenxpp,fcenypp+1] ],/double, /nan);/!dtor            
  dcenpfwp  = dcenfwp/mean([ dis_rg[fcenxpp-1,fcenypp], dis_rg[fcenxpp+1,fcenypp], dis_rg[fcenxpp,fcenypp-1], dis_rg[fcenxpp,fcenypp+1] ],/double, /nan);/!dtor
   

  ;Recalculating parameters for negative region------------------------------------------------------
  flxn  = total(CRD_in.mgnt_flx[bstinn],/double, /nan)       ;total flx
  arean  = total(CRD_in.mgnt_ar[bstinn],/double, /nan)         ;total area

  ;Geographical Center
  fcenxn = total(CRD_in.Xar[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flxn
  fcenyn = total(CRD_in.Yar[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flxn
  fcenzn = total(CRD_in.Zar[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flxn

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
  temp_2d = array_indices(im,minin)
  fcenxpn = temp_2d[0,0];    ;Center in x
  fcenypn = temp_2d[1,0];    ;Center in y

  ;Mean and Pixel Radius
  dcenn = mean( dis_rg[bstinn],/double, /nan)
  dcenfwn = total( dis_rg[bstinn]*CRD_in.mgnt_flx[bstinn],/double, /nan)/flxn

  dcenpn    = dcenn/mean([ dis_rg[fcenxpn-1,fcenypn], dis_rg[fcenxpn+1,fcenypn], dis_rg[fcenxpn,fcenypn-1], dis_rg[fcenxpn,fcenypn+1] ],/double, /nan);/!dtor            
  dcenpfwn  = dcenfwn/mean([ dis_rg[fcenxpn-1,fcenypn], dis_rg[fcenxpn+1,fcenypn], dis_rg[fcenxpn,fcenypn-1], dis_rg[fcenxpn,fcenypn+1] ],/double, /nan);/!dtor           


  ;Storing active region--------------------------------------------------------------------------

  ;Adding region indexes
  tempAR[maxim_ind].indxp = strjoin(string(bstinp)) 
  tempAR[maxim_ind].indxn = strjoin(string(bstinn)) 

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

     LatL  = fcenltp*!dtor
     LongL = Lonp*!dtor

     LatF  = tempAR[maxim_ind].fcn_ltn*!dtor
     LongF = Lonn*!dtor
     
     lp = 1;

  endif else begin

     LatF  = fcenltn*!dtor
     LongF = Lonp*!dtor

     LatL  = tempAR[maxim_ind].fcn_ltn*!dtor
     LongL = Lonn*!dtor
     
     lp = -1;

  endelse
         
  dlat = abs(LatF-LatL)
  dlon = abs(LongF-LongL)

  Dis = 2.0*asin(sqrt( sin(dlat/2.0)^2.0 + cos(LatL)*cos(LatF)*sin(dlon/2.0)^2.0 ))/!dtor ; Angular distance between polarities
  Tilt  = real_part(asin( sin( (LatF-LatL) )/sin(Dis*!dtor) ))/!dtor; Tilt

  ;Updating relevant quantities
  temp_ar = {indxp: bstinp, indxn: bstinn, fluxp: flxp, areap: areap, $
           fcn_ltp: fcenltp, fcn_lnp: fcenlnp, dcenp: dcenp, $ 
           fcenxpp: fcenxpp, fcenypp: fcenypp, dcenpp: dcenpp, $
           fluxn: flxn, arean: arean, $
           fcn_ltn: fcenltn, fcn_lnn: fcenlnn, dcenn: dcenn, $ 
           fcenxpn: fcenxpn, fcenypn: fcenypn, dcenpn: dcenpn, $
           dis: Dis, tilt: Tilt, lp: lp} 
                
  ar = tempAR[maxim_ind]                        
  STRUCT_ASSIGN, temp_ar, ar, /nozero 
  ARs[ar_in[maxim_ind]] = ar  
                  
  ;Updating indices to reset unused magnetogram
  ;inp = ind_gn
  ;inp = bstinp
  ;inn = bstinn

  ;Removing region from the unused magnetogram
  imsel[bstinp] = 0.0
  imsel[bstinn] = 0.0


  Flx_imr[maxim_ind] = !values.f_nan
endfor
bal_msg = "Regions balanced:" + string(bal_count)
print, bal_msg
end
