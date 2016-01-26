;-------------------------------------------------------------------------------------------------------------------------------
;	NAME:				zaw_ar_flux_balance
;
;	PURPOSE:			Takes the active regions (ARs) as parameters and balances the flux by growing the regions in an
;						iterative fashion.
;
;	CALLING SEQUENCE:	zaw_ar_flux_balance, 	
;
;
;
;-------------------------------------------------------------------------------------------------------------------------------
pro zaw_ar_flux_balance, CRD_in, mdi_i, lbl, PRs, NRs, AR_c, date, ar_cnst=ar_cnst, seg_const=seg_const, bck_trck = bck_trck, display=display, prt=prt, info=info, sqs_nm = sqs_nm, pltn = pltn


if (n_ars ne 1) then begin

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

end
