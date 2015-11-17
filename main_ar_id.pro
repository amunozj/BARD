FUNCTION button_choice, pls, pad, plxy, plxx

ftsz = 14

;Defining control (right) buttons
;labs = ['RE-DETECT', 'RE-LABEL', 'AGGREGATE', 'CREATE', 'DELETE', 'QUIT' ]
labs = ['LON-GUIDE','PRINT SCREEN','CONTRAST -','CONTRAST +','SYNCHRONIZE','RE-DETECT', 'FRAGMENT', 'MERGE', 'RE-LABEL', 'CREATE/TRACK', 'CREATE', 'DELETE ALL', 'DELETE ONE', 'QUIT' ]

nl = n_elements( labs )

x0c = 2.0*pls+4.5*pad
x1c = 2.0*pls+3.5*pad+plxx
ybc = 0.25*pad + float(findgen(nl+1))/float(nl)*(pls+plxy+ 3.0*pad/4.0)

;Plotting control (right) buttons
plots, x0c*[1, 1],  [ybc[0], ybc[nl]], /device
plots, x1c*[1, 1],  [ybc[0], ybc[nl]], /device

plots, [ x1c,  x0c ], ybc[0]*[1, 1], /device
FOR i=0, nl-1 DO BEGIN
  y1 = ybc[i+1]
  y0 = 0.5*( ybc[i] + ybc[i+1] ) - ftsz/2
  xyouts, 0.5*(x1c+x0c), y0, labs[i], /device, align=0.5, charsize=2
  plots, [ x1c,  x0c ], y1*[1, 1], /device
ENDFOR



;Defining revision (bottom left) buttons
labs2 = ['< x10', '<', 'JUMP TO','>',  '> x10']
n2 = n_elements( labs2 )

xbr = 2.0*pad + float(findgen(n2+1))/float(n2)*pls

y0r = pad/4
y1r = plxy - pad*3/4

;Plotting revision (bottom left) buttons
plots, [xbr[0], xbr[n2]], y0r*[1, 1], /device
plots, [xbr[0], xbr[n2]], y1r*[1, 1], /device

plots, xbr[0]*[1, 1], [ y1r,  y0r ], /device
FOR i=0, n2-1 DO BEGIN
  x1b = xbr[i+1]
  x0b = 0.5*( xbr[i] + xbr[i+1] )
  xyouts, x0b, 0.5*(y1r+y0r) - ftsz/2, labs2[i], /device, align=0.5, charsize=2
  plots, x1b*[1, 1], [ y0r, y1r ], /device
ENDFOR

x0b = 0.5*( xbr[0] + xbr[n2] )
xyouts, x0b, plxy - pad/2, 'CLICK MAGNETOGRAM ABOVE TO CYCLE OVERLAY', /device, align=0.5, charsize=2


;Defining detection advancement (bottom right) buttons
labs3 = ['<', 'JUMP TO', '>']
n3 = n_elements( labs3 )

xba = 4.0*pad + pls + float(findgen(n3+1))/float(n3)*pls

y0a = pad/4
y1a = plxy - pad*3/4

;Plotting detection advancement (bottom right) buttons
plots, [xba[0], xba[n3]], y0a*[1, 1], /device
plots, [xba[0], xba[n3]], y1a*[1, 1], /device

plots, xba[0]*[1, 1], [ y1a,  y0a ], /device
FOR i=0, n3-1 DO BEGIN
  x1b = xba[i+1]
  x0b = 0.5*( xba[i] + xba[i+1] )
  xyouts, x0b, 0.5*(y1a+y0a) - ftsz/2, labs3[i], /device, align=0.5, charsize=2
  plots, x1b*[1, 1], [ y0a, y1a ], /device
ENDFOR

x0b = 0.5*( xba[0] + xba[n3] )
xyouts, x0b, plxy - pad/2, 'CLICK MAGNETOGRAM ABOVE TO CYCLE OVERLAY', /device, align=0.5, charsize=2



;Plotting Syncronization button
;plots, [x0c, x1c], y0a*[1, 1], /device
;plots, [x0c, x1c], y1a*[1, 1], /device
;plots, x0c*[1, 1], [y0a, y1a], /device
;plots, x1c*[1, 1], [y0a, y1a], /device

;x0b = 0.5*( x0c + x1c )
;y0b = 0.5*( y0a + y1a )

;xyouts, x0b, y0b, 'SYNCHRONIZE', /device, align=0.5, charsize=2



;Identificaiton of buttons
cursor, u, v, 3, /device

;Switch for outside
out_sw = 1

;Cursor inside the control (right) buttons
if ( ( u ge x0c ) and ( u le x1c ) and ( v ge ybc[0] ) and ( v le ybc[nl] ) ) then begin
  ii = where( v LE ybc)
  ii = nl - ii + 1
  out_sw = 0
endif

;Cursor inside the reference (bottom left) buttons
if ( ( u ge xbr[0] ) and ( u le xbr[n2] ) and ( v ge y0r ) and ( v le y1r ) ) then begin
  ii = where( u LE xbr)
  ii = ii + nl
  out_sw = 0
endif

;Cursor inside the operations (bottom right) buttons
if ( ( u ge xba[0] ) and ( u le xba[n3] ) and ( v ge y0a ) and ( v le y1a ) ) then begin
  ii = where( u LE xba)
  ii = ii + nl + n2
  out_sw = 0
endif

;Cursor inside the left overlay (reference)
if ( ( u ge 2.0*pad ) and ( u le pls+2.0*pad ) and ( v ge pad+plxy ) and ( v le pls+pad+plxy ) ) then begin
  ii = 1 + nl + n2 + n3
  out_sw = 0
endif

;Cursor inside the right overlay (control)
if ( ( u ge pls+4.0*pad ) and ( u le 2.0*pls+4.0*pad ) and ( v ge pad+plxy ) and ( v le pls+pad+plxy ) ) then begin
  ii = 2 + nl + n2 + n3
  out_sw = 0
endif

;Cursor inside synchronize button
;if ( ( u ge x0c ) and ( u le x1c ) and ( v ge y0r ) and ( v le y1r ) ) then begin
;  ii = 3 + nl + n2 + n3
;  out_sw = 0
;endif


;Cursor outside buttons area
if  (out_sw eq 1) then begin
  ii = 0
endif

return, ii[0]
END



;----------------------------------------------------------------------------------------------------------


PRO amj_pick_ar, mdi_ir, mdi_il, ARs, inx, pls, pad, plxy, plxx, d_zoom

print, 'Left click once on AR  ||  Rigth Click to cancel'
scc_sw = 0 ;Switch to indicate success

while ( (scc_sw eq 0) and ( !mouse.button ne 4 )) do begin
  
  cursor, u, v, 3, /device
  
  ;Check Window  
  if (u le pls+3.0*pad) then begin
    pos = [2.0*pad, pad+plxy]
    mdi_t = mdi_il
  endif else begin
    pos = [pls+4.0*pad, pad+plxy]
    mdi_t = mdi_ir
  end
  
  ;Check for ARs
  arin  = where(ARs.mdi_i eq mdi_t, nars)
  
  ;Look for region
  if (nars ne 0) then begin
    
    tARs = ARs[arin]
    
    ;Positions of regions
    xp = tARs.fcenxpp*d_zoom + pos[0]
    yp = tARs.fcenypp*d_zoom + pos[1]
    
    xn = tARs.fcenxpn*d_zoom + pos[0]
    yn = tARs.fcenypn*d_zoom + pos[1]
  
    ;Find distances
    disp = sqrt((u-xp)^2.0 + (v-yp)^2.0)
    disn = sqrt((u-xn)^2.0 + (v-yn)^2.0)
    
    ;Find Minima
    minp = min(disp,inp)
    minn = min(disn,inn)
    
    ;making sure it is close enough
    if ( (minn lt 2.0*tARs[inn].dcenpn*d_zoom) or (minp lt 2.0*tARs[inp].dcenpp*d_zoom) ) then begin
      
      ;Choosing active region
      if (minp lt minn) then begin
         t_in = inp
      endif else begin
         t_in = inn
      endelse
      inx = arin[t_in]
      scc_sw = 1 
      
    endif else begin
      
      print, 'No AR selected'
      print, 'Left click once on AR  ||  Rigth Click to cancel'
      
    endelse
  
  endif else begin
    
    print, 'No ARs on that window'
    print, 'Left click once on AR  ||  Rigth Click to cancel'
  endelse
  
endwhile

if (!mouse.button eq 4) then inx = -2
                
return
END

;----------------------------------------------------------------------------------------------------------



PRO amj_pick_long, hdr_r, hdr_l, ARs, latlr, Lonhl, Lonhr, pls, pad, plxy, plxx, d_zoom, instr

print, 'Left click once on disk  ||  Rigth Click to cancel'
scc_sw = 0 ;Switch to indicate success


;LEFT MAGNETOGRAM-------------------------------------------------------

;HMI uses structures for header values
if instr eq 4 then begin
	datel = hdr_l.DATE_OBS

	;Define center and radius
	hfxl = hdr_l.CRPIX1 ;  Location of the center in x pixels 
	hfyl = hdr_l.CRPIX2 ;    Location of the center in y pixels
	dil = hdr_l.RSUN_OBS/hdr_l.CDELT1;

	;Load Solar Coordinates
	P0l = 0.0
	RDl = hdr_l.DSUN_OBS/hdr_l.RSUN_REF
	B0l = hdr_l.CRLT_OBS
	L0l = hdr_l.CRLN_OBS

	;Observer Coordinates
	X_scll = hdr_l.CDELT1/60.0
	Y_scll = hdr_l.CDELT2/60.0
	
endif else begin


	datel = fxpar(hdr_l, 'DATE_OBS')

	;KPVT-512
	if instr eq 1 then begin
	
		;Define center and radius
		hfxl = fxpar(hdr_l, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
		hfyl = fxpar(hdr_l, 'CRPIX2A');+1.0;    Location of the center in y pixels
		dil = fxpar(hdr_l,'EPH_R0');

		;Load Solar Coordinates
		P0l = 0.0
		RDl = !values.f_nan
		B0l = fxpar(hdr_l, 'EPH_B0')
		L0l = fxpar(hdr_l, 'EPH_L0')

		;Observer Coordinates
		X_scll = fxpar(hdr_l, 'CDELT1')*fxpar(hdr_l, 'CRR_SCLX')/60.0
		Y_scll = fxpar(hdr_l, 'CDELT2')*fxpar(hdr_l, 'CRR_SCLY')/60.0

	endif

	;MDI
	if instr eq 3 then begin
	
		;Define center and radius
		hfxl = fxpar(hdr_l, 'X0');  Location of the center in x pixels 
		hfyl = fxpar(hdr_l, 'Y0');  Location of the center in y pixels
		dil = fxpar(hdr_l,'R_SUN');

		;Load Solar Coordinates
		P0l = fxpar(hdr_l, 'P_ANGLE')
		RDl = fxpar(hdr_l, 'OBS_DIST')/0.0046491
		B0l = fxpar(hdr_l, 'B0')
		L0l = fxpar(hdr_l, 'L0')

		;Observer Coordinates
		X_scll = fxpar(hdr_l, 'XSCALE')
		Y_scll = fxpar(hdr_l, 'YSCALE')	
	
	endif

endelse

;RIGHT MAGNETOGRAM-------------------------------------------------------

;HMI uses structures for header values
if instr eq 4 then begin
	dater = hdr_r.DATE_OBS

	;Define center and radius
	hfxr = hdr_r.CRPIX1 ;  Location of the center in x pixels 
	hfyr = hdr_r.CRPIX2 ;    Location of the center in y pixels
	dir = hdr_r.RSUN_OBS/hdr_r.CDELT1;

	;Load Solar Coordinates
	P0r = 0.0
	RDr = hdr_r.DSUN_OBS/hdr_r.RSUN_REF
	B0r = hdr_r.CRLT_OBS
	L0r = hdr_r.CRLN_OBS

	;Observer Coordinates
	X_sclr = hdr_r.CDELT1/60.0
	Y_sclr = hdr_r.CDELT2/60.0
	
endif else begin


	dater = fxpar(hdr_r, 'DATE_OBS')

	;KPVT-512
	if instr eq 1 then begin
	
		;Define center and radius
		hfxr = fxpar(hdr_r, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
		hfyr = fxpar(hdr_r, 'CRPIX2A');+1.0;    Location of the center in y pixels
		dir = fxpar(hdr_r,'EPH_R0');

		;Load Solar Coordinates
		P0r = 0.0
		RDr = !values.f_nan
		B0r = fxpar(hdr_r, 'EPH_B0')
		L0r = fxpar(hdr_r, 'EPH_L0')

		;Observer Coordinates
		X_sclr = fxpar(hdr_r, 'CDELT1')*fxpar(hdr_r, 'CRR_SCLX')/60.0
		Y_sclr = fxpar(hdr_r, 'CDELT2')*fxpar(hdr_r, 'CRR_SCLY')/60.0

	endif

	;MDI
	if instr eq 3 then begin
	
		;Define center and radius
		hfxr = fxpar(hdr_r, 'X0');  Location of the center in x pixels 
		hfyr = fxpar(hdr_r, 'Y0');  Location of the center in y pixels
		dir = fxpar(hdr_r,'R_SUN');

		;Load Solar Coordinates
		P0r = fxpar(hdr_r, 'P_ANGLE')
		RDr = fxpar(hdr_r, 'OBS_DIST')/0.0046491
		B0r = fxpar(hdr_r, 'B0')
		L0r = fxpar(hdr_r, 'L0')

		;Observer Coordinates
		X_sclr = fxpar(hdr_r, 'XSCALE')
		Y_sclr = fxpar(hdr_r, 'YSCALE')	
	
	endif

endelse

while ( (scc_sw eq 0) and ( !mouse.button ne 4 )) do begin
  
	cursor, u, v, 3, /device

	;Check Window  
	if (u le pls+3.0*pad) then begin
		pos  = [2.0*pad, pad+plxy]
		refm = -1
		rad  = dil
		xcen = hfxl
		ycen = hfyl	

		scl  = (X_scll+Y_scll)/2.0 
		
	endif else begin
		pos = [pls+4.0*pad, pad+plxy]
		refm = 1
		rad  = dir
		xcen = hfxr
		ycen = hfyr

	    scl  = (X_sclr+Y_sclr)/2.0 
		
	end

	;Find positions
	
	refx = ((u-pos[0]))/d_zoom
	refy = ((v-pos[1]))/d_zoom
	
	;Find distances
    disp = sqrt((refx-xcen)^2.0 + (refy-ycen)^2.0)

	;Making sure the point is inside the solar disk
	if disp gt rad then begin
		print, 'Please click inside the solar disk'
		print, 'Left click once on disk  ||  Rigth Click to cancel'
	endif else begin
	
        ; Differential rotation profile from Snodgrass (1983), gives rotation
        ; rate vs. lat. in microradians per sec:
        ;
        ;  omega = snod_A + snod_B*sin(latitude)^2 + snod_C*sin(latitude)^4
        ;====================================================================
        ;snod_A =  2.902  ; magnetic rot. coeffs, in microrad.    
        snod_A =  0.0367;2.902  ; magnetic rot. coeffs, in microrad.    
        ;Set to 0.0367 because the heliographic coordinates include the carrington rotation
        
        snod_B = -0.464
        snod_C = -0.328   

		Lath = 90-findgen(1,181)
		
        ;Accounting for differential rotation
        omegap = 1e-6*(snod_A + $ ; Omega at each pixel's lat., in radians
                      snod_B*sin(abs(Lath)*!dtor)^2 + $
                      snod_C*sin(abs(Lath)*!dtor)^4 )/!dtor	
		
	
		;Calculating longitude of clicked point and equivalent on the other magnetogram window
		
		;LEFT MAGNETOGRAM-------------------------------------------------------
		if refm eq -1 then begin

			helio = arcmin2hel((refx-xcen)*scl, (refy-ycen)*scl, date = datel, p = P0l, b0 = B0l, l0 = L0l, sphere = 1, rsun = rad*scl*60)
			dtim = anytim(datel,/utime)-anytim(dater,/utime)

			Lonhl = helio[1] + omegap*0.0			
			Lonhr = helio[1] - omegap*dtim
		 
		;RIGHT MAGNETOGRAM-------------------------------------------------------
		endif else begin

			helio = arcmin2hel((refx-xcen)*scl, (refy-ycen)*scl, date = dater, p = P0r, b0 = B0r, l0 = L0r, sphere = 1, rsun = rad*scl*60)
			dtim = anytim(dater,/utime)-anytim(datel,/utime)

			Lonhr = helio[1] + omegap*0.0			
			Lonhl = helio[1] - omegap*dtim		
		
		endelse
		
		latlr = helio[0]
		
		print, refx, refy, helio[0], helio[1], dtim
	
	
		scc_sw = 1	

	endelse
		
endwhile

if (!mouse.button eq 4) then inx = -2
                
return
END

;----------------------------------------------------------------------------------------------------------





PRO amj_pick_pnr, mdi_ir, PRs, NRs, inx, pn_sw, pls, pad, plxy, plxx, d_zoom, all = all

print, 'Left click once on PNR  ||  Rigth Click to cancel'
scc_sw = 0 ;Switch to indicate success

while ( (scc_sw eq 0) and ( !mouse.button ne 4 )) do begin
  
  cursor, u, v, 3, /device
  
  ;Check Window  
  if (u gt pls+3.0*pad) then begin
    pos = [pls+4.0*pad, pad+plxy]
    mdi_t = mdi_ir
    
    ;finding Available PRs
    if keyword_set(all) then begin
      pr_ix = where(PRs.mdi_i eq mdi_t, n_regp)      
    endif else begin
      pr_ix = where((PRs.mdi_i eq mdi_t) and (PRs.lnk_sw eq 0), n_regp)
    endelse
    
    ;Finding closest positive region    
    mdisp = 9999.0
    inp = -1
    if (n_regp gt 0) then begin

      ;Positions of regions
      xp = PRs[pr_ix].fcenxp*d_zoom + pos[0]
      yp = PRs[pr_ix].fcenyp*d_zoom + pos[1]      

      ;Find distances
      disp = sqrt((u-xp)^2.0 + (v-yp)^2.0)

      ;Find Minima
      minp = min(disp,inp)
      
      if (minp lt 2.0*PRs[pr_ix[inp]].dcenp*d_zoom) then begin
        scc_sw = 1
      endif
      
    endif
    
    ;finding Available PRs
    if keyword_set(all) then begin
      nr_ix = where(NRs.mdi_i eq mdi_t, n_regn)
    endif else begin
      nr_ix = where((NRs.mdi_i eq mdi_t) and (NRs.lnk_sw eq 0), n_regn)
    endelse

    ;Finding closest negative region 
    mdisn = 9999.0
    inn = -1
    if (n_regn gt 0) then begin

      ;Positions of regions
      xp = NRs[nr_ix].fcenxp*d_zoom + pos[0]
      yp = NRs[nr_ix].fcenyp*d_zoom + pos[1]      

      ;Find distances
      disn = sqrt((u-xp)^2.0 + (v-yp)^2.0)

      ;Find Minima
      minn = min(disn,inn)

      if (minn lt 2.0*NRs[nr_ix[inn]].dcenp*d_zoom) then begin
        scc_sw = 1
      endif
      
    endif

    if (scc_sw eq 1) then begin
      
      ;Choosing closest region
      if (minp lt minn) then begin
         inx = pr_ix[inp]
         pn_sw = 1
      endif else begin
         inx = nr_ix[inn]
         pn_sw = -1
      endelse
      
    endif else begin

      print, 'No available PNR selected'
      print, 'Left click once on PNR  ||  Rigth Click to cancel'
    
    endelse
    
  endif else begin
    print, 'PNRs must be selected on the right window'
    print, 'Left click once on PNR  ||  Rigth Click to cancel'    

  end
    
endwhile

if (!mouse.button eq 4) then inx = -2
                
return
END



PRO amj_pnr_merge, CRD_in, mdi_i, inar, innar, Rs, ARs, pl_sw, date

;Extracting indices of Growing region
inpg = long(strsplit(Rs[inar].indx,/extract))

;Extracting indices of Killed region
inpk = long(strsplit(Rs[innar].indx,/extract))

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
Rs[inar].indx = strjoin(string(ind_gn))                    

;Updating relevant quantities
tmp_r = {flux: flux, area: area, $
          fcn_lt: fcenlt, fcn_ln: fcenln, dcen: dcen, $ 
          fcenxp: fcenxp, fcenyp: fcenyp, dcenp: dcenp, dm:dm, qm:qm, fr_lbl:-999} 
               
prtmp = Rs[inar]                        
STRUCT_ASSIGN, tmp_r, prtmp, /nozero
 
Rs[inar] = prtmp       
                   

;Modifying AR properties
if (Rs[inar].lnk_sw eq 1) then begin
  
  tmp_in = where((ARs.mdi_i eq mdi_i) and (ARs.labl eq Rs[inar].ar_lbl))
  print, tmp_in
    
  ;For positive regions
  if (pl_sw eq 1) then begin
      
    ;Adding region indexes
    ARs[tmp_in].indxp = strjoin(string(ind_gn))                    

    ;Making Sure all visible longitudes are in the same hemisphere
    ;---------------------------------------------------------
    sun_data = get_sun(date,carr=carr,he_lon=he_lon)

    Lonp = fcenln
    Lonn = ARs[tmp_in].fcn_lnn
                            
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

        LatF  = ARs[tmp_in].fcn_ltn*!dtor
        LongF = Lonn*!dtor
        
        lp = 1;
    
    endif else begin

        LatF  = fcenlt*!dtor
        LongF = Lonp*!dtor

        LatL  = ARs[tmp_in].fcn_ltn*!dtor
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
                   
    ar = ARs[tmp_in]                        
    STRUCT_ASSIGN, tmp_ar, ar, /nozero 
    ARs[tmp_in] = ar  
  
  ;For negative regions  
  endif else begin
    
    ;Adding region indexes
    ARs[tmp_in].indxn = strjoin(string(ind_gn))                    

    ;Making Sure all visible longitudes are in the same hemisphere
    ;---------------------------------------------------------
    sun_data = get_sun(date,carr=carr,he_lon=he_lon)

    Lonp = ARs[tmp_in].fcn_lnp
    Lonn = fcenln
                            
    if ( (he_lon ge 90.0) and (he_lon le 270.0) ) then begin
            
      tmpin = where(Lonp lt 0, n_tmp)
      if (n_tmp gt 0) then Lonp[tmpin] = Lonp[tmpin]+360.0
            
      tmpin = where(Lonn lt 0, n_tmp)
      if (n_tmp gt 0) then Lonn[tmpin] = Lonn[tmpin]+360.0
            
    endif 
    
    
    ;Finding leading polarity
    if ( Lonp ge Lonn ) then begin
    
        LatL  = ARs[tmp_in].fcn_ltn*!dtor
        LongL = Lonp*!dtor

        LatF  = fcenlt*!dtor
        LongF = Lonn*!dtor
        
        lp = 1;
    
    endif else begin

        LatF  = ARs[tmp_in].fcn_ltn*!dtor
        LongF = Lonp*!dtor

        LatL  = fcenlt*!dtor
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
    tmp_ar = {fluxn: flux, arean: area, $
              fcn_ltn: fcenlt, fcn_lnn: fcenln, dcenn: dcen, $ 
              fcenxpn: fcenxp, fcenypn: fcenyp, dcenpn: dcenp, $
              dis: Dis, tilt: Tilt, lp: lp, dm:dm, qm:qm} 
                   
    ar = ARs[tmp_in]                        
    STRUCT_ASSIGN, tmp_ar, ar, /nozero 
    ARs[tmp_in] = ar      
      
  endelse 


endif

;Removing merged region
first = 0
last = N_Elements(Rs)-1
CASE innar OF
  first: Rs = Rs[1:*]
  last: Rs = Rs[first:last-1]
  ELSE: Rs = [ Rs[first:innar-1], Rs[innar+1:last] ]
ENDCASE
    


return
END


;-----------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------

;Control constant definitions


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


;EDITS
PRO main_ar_id, start_date=start_date, continue = continue, restoref = restoref, pnr_file = pnr_file, labl = labl, buff = buff, nobuff = nobuff, kerth1=kerth1, kerth2=kerth2, arth=arth, eros=eros, dila=dila, dislim=dislim, ovrlim=ovrlim, ardislim1=ardislim1, exp_f=exp_f, exp_d=exp_d, exp_s=exp_s, MxFlxim=MxFlxim, instr = instr, rtrn_chk = rtrn_chk

device,retain=2
  
;Creating backups------------------------

;Final State
n_bkps = 2
spawn, 'ls Ar_id_Sav*.sav', result 
s = size(result)
if s[0] ne 0 then begin

  if (s[1] lt n_bkps) then n_bkps = s[1]
  
  for i = n_bkps, 1, -1 do begin
    spawn, '\cp -f Ar_id_Save' + strtrim(string(i-1),2) + '.sav Ar_id_Save' + strtrim(string(i),2) + '.sav'    
  endfor

endif


;;ARs
;n_bkps = 2
;spawn, 'ls ARs*.sav', result 
;s = size(result)
;if s[0] ne 0 then begin
;
;  if (s[1] lt n_bkps) then n_bkps = s[1]
;  
;  for i = n_bkps, 1, -1 do begin
;    spawn, '\cp -f ARs' + strtrim(string(i-1),2) + '.sav ARs' + strtrim(string(i),2) + '.sav'    
;  endfor
;
;endif



;initializing run and variables
;-------------------------------------------------------------------------------------

;Return Check

rchk_sw = 0		;Switch that plots the operational active regions on the reference magnetogram
dy_skp  = 0	    ;Default amount of days to skip ahead to see the return

if keyword_set(rtrn_chk) then begin
	rchk_sw = 1
	dy_skp  = 24
endif

;Reference longitude and latitude
ref_sw = 0   			;Switch that turns on the reference lines.
Lonhrl = !values.f_nan	;Longitude of the reference line in the left magnetogram
Lonhrr = !values.f_nan  ;Longitude of the reference line in the right magnetogram
latlr  = !values.f_nan  ;Latitude of the reference line in the both magnetograms

;Buffer definition
buff_sw = 0; Switch indicating the run is using a buffer
bfr_ds = 60;  maximum number of buffer days
bfr_lm = 10; Number of days that trigger a buffer update
sv_bfr = 10; Number of days used to save
sv_sw = 1 ; Saving switch for the buffer run

iskip = 10; number of days to skip

ds_swl = 2; Switch AR overlay
ds_swr = 1; Switch AR overlay

mdi_ir_vis = 0 ; Variable that keeps track of the maximum date where AR detection has taken place


und_sw = 1; Save undo state switch

max_sat = 600;   Default image contrast
; constants to change max_sat
incrmnt = 100
max_max_sat = 2000
min_max_sat = 0

prepnr_sw = 0; Pre-detection of PNRs




if ( (not keyword_set( continue )) and (not keyword_set(restoref)) ) then begin

	IF( keyword_set( start_date ) ) THEN date0 = start_date $
	ELSE begin
		if instr eq 1 then date0 = '1981-06-08';	KPVT 512
		if instr eq 2 then date0 = '1992-04-21';	KPVT SPMG
		if instr eq 3 then date0 = '1997-01-01';	MDI
		if instr eq 4 then date0 = '2010-11-21';	HMI
	endelse

	;First reference day for keeping track time easily
	if instr eq 1 then DayOff = julday(1,1,1970);	KPVT 512
	if instr eq 2 then DayOff = julday(1,1,1970);	KPVT SPMG
	if instr eq 3 then DayOff = julday(1,1,1993);	MDI
	if instr eq 4 then DayOff = julday(1,1,2009);	HMI

      
    mdi_ir  = julday(strmid(date0,5,2),strmid(date0,8,2),strmid(date0,0,4))-DayOff
    mdi_il = mdi_ir - 1
    
    
    lbl = 0;  AR label
    detpn_sw  = 1; %PNR Detection switch
    detar_sw  = 1; %AR Detection switch
        
    ;Initializing AR array
    ARs = {ar, mdi_i: 0L, date: '', labl: 0, clr: 0, indxp: '', indxn: '', fluxp: !values.f_nan, fluxn: !values.f_nan, areap:!values.f_nan, arean:!values.f_nan, $
              fcn_ltp: !values.f_nan, fcn_lnp: !values.f_nan, dcenp: !values.f_nan, $ 
              fcn_ltn: !values.f_nan, fcn_lnn: !values.f_nan, dcenn: !values.f_nan, $
              fcenxpp: !values.f_nan, fcenypp: !values.f_nan, dcenpp: !values.f_nan, $
              fcenxpn: !values.f_nan, fcenypn: !values.f_nan, dcenpn: !values.f_nan, $               
              dis: !values.f_nan, tilt: !values.f_nan, lp: !values.f_nan, dm: !values.f_nan, qm: !values.f_nan}  
    
    ;Initializing PR array
    PRs = {rgn,  lnk_sw: 0,  mdi_i: 0L, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}
    
    ;Initializing NR array
    NRs = PRs
    
endif else begin
  
  detpn_sw  = 0; %PNR Detection switch
  detar_sw  = 0; %AR Detection switch  
  
  
  if keyword_set(continue) then restore, 'Ar_id_Save0.sav'
  if keyword_set(restoref) then begin 
    restore, restoref
    mdi_ir = min(ARs[where(ARs.mdi_i ne 0)].mdi_i)
    mdi_il = mdi_ir - 1    
  endif
  
	IF( keyword_set( start_date ) ) THEN date0 = start_date $
	ELSE begin
		if instr eq 1 then date0 = '1981-06-08';	KPVT 512
		if instr eq 2 then date0 = '1992-04-21';	KPVT SPMG
		if instr eq 3 then date0 = '1997-01-01';	MDI
		if instr eq 4 then date0 = '2010-11-21';	HMI
	endelse

	;First reference day for keeping track time easily
	if instr eq 1 then DayOff = julday(1,1,1970);	KPVT 512
	if instr eq 2 then DayOff = julday(1,1,1970);	KPVT SPMG
	if instr eq 3 then DayOff = julday(1,1,1993);	MDI
	if instr eq 4 then DayOff = julday(1,1,2009);	HMI

   
  ;-----------------------------------------------------------------
  ;CLEANUP
  
  pr_ix = where((PRs.lnk_sw eq 0) and (PRs.ar_lbl ne 0), n_regp)
  if (n_regp gt 0) then PRs[pr_ix].ar_lbl = 0

  nr_ix = where((NRs.lnk_sw eq 0) and (NRs.ar_lbl ne 0), n_regn)
  if (n_regn gt 0) then NRs[nr_ix].ar_lbl = 0
  
  
endelse


if keyword_set(pnr_file) then begin
  
  restore, pnr_file
  mdi_ir = min([min(PRs[where(PRs.mdi_i ne 0)].mdi_i),min(NRs[where(NRs.mdi_i ne 0)].mdi_i)])
  mdi_il = mdi_ir - 1    
  prepnr_sw = 1
  detpn_sw  = 0

endif

if (prepnr_sw eq 1) then detpn_sw  = 0

if keyword_set(continue) and keyword_set(buff) then buff_sw=1
if keyword_set(continue) and keyword_set(nobuff) then buff_sw=0


;Active region detection constants
ar_cnst={dis_lim1:5.0, dis_lim2:4.0, exp_f: 1.0,  exp_d: 4.0, exp_s: 1.0, mlth: 40.0, mxB: 180.0, MxFlxim:3.0, Imb_tol: 0.10, Imb_it: 10, lim_lon: -90.0, k_sig:15.0, npr: 5, nmgnt: 5 , vld_thr: 0.69, valid_range:[-20000.,20000.]}

;KPVT 512
if instr eq 1 then begin
	seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;KPVT SMPG
if instr eq 2 then begin
	seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;MDI
if instr eq 3 then begin
	seg_const={ker_th:500.0, ker_th2:275.0, ar_th:325.0, ar_th2:50.0, eros_size:9.0, dila_size:18.0, npssu: 3, dis_lim:3.0, ovr_lim: 0.3, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;HMI
if instr eq 4 then begin
	seg_const={ker_th:200.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:30.0, npssu:1, dis_lim:2.0, ovr_lim: 0.0, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;EDITS
if ( n_elements( kerth1 ) ne 0 ) then seg_const.ker_th = kerth1
if ( n_elements( kerth2 ) ne 0 ) then seg_const.ker_th2 = kerth2
if ( n_elements( arth ) ) ne 0 then seg_const.ar_th = arth
if ( n_elements( eros ) ) ne 0 then seg_const.eros_size = eros
if ( n_elements( dila ) ) ne 0 then seg_const.dila_size = dila
if ( n_elements( dislim ) ne 0 ) then seg_const.dis_lim = dislim
if ( n_elements( ovrlim ) ne 0 ) then seg_const.ovr_lim = ovrlim
if ( n_elements( ardislim1 ) ne 0 ) then ar_cnst.dis_lim1 = ardislim1
if ( n_elements( exp_f ) ne 0 ) then ar_cnst.exp_f = exp_f
if ( n_elements( exp_d ) ne 0 ) then ar_cnst.exp_d = exp_d
if ( n_elements( exp_s ) ne 0 ) then ar_cnst.exp_s = exp_s
if ( n_elements(MxFlxim) ne 0 ) then ar_cnst.MxFlxim = MxFlxim


;Creating buffer if run is continued
if (buff_sw eq 1) then begin

  print, 'Creating buffer' 
  
  PRs_Big = PRs
  NRs_Big = NRs
  ARs_Big = ARs

  if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin
      PRsFr_Big = PRsFr
      NRsFr_Big = NRsFr      
  end

  ;Updating Variables that keep track of the limits in which the buffer need to be extended
  bff_dn = mdi_ir - 2*bfr_lm
  bff_up = mdi_ir - 2*bfr_lm + bfr_ds

  mdi_bf_lm1 = bff_dn + bfr_lm
  mdi_bf_lm2 = bff_up - bfr_lm
  
  ;Updating buffer
  tmp_in = where((PRs_Big.mdi_i le bff_up) and (PRs_Big.mdi_i ge bff_dn), nprs)
  if (nprs ne 0) then begin
    PRs = PRs_Big[tmp_in]          
  endif else begin
    PRs = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
  endelse
        
  tmp_in = where((NRs_Big.mdi_i le bff_up) and (NRs_Big.mdi_i ge bff_dn), nprs)
  if (nprs ne 0) then begin
    NRs = NRs_Big[tmp_in]          
  endif else begin
    NRs = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
  endelse
 
  tmp_in = where((ARs_Big.mdi_i le bff_up) and (ARs_Big.mdi_i ge bff_dn), nprs)
  if (nprs ne 0) then begin
    ARs = ARs_Big[tmp_in]          
  endif else begin
    ARs = {ar, mdi_i: 0, date: '', labl: 0, clr: 0, indxp: '', indxn: '', fluxp: !values.f_nan, fluxn: !values.f_nan, areap:!values.f_nan, arean:!values.f_nan, $
              fcn_ltp: !values.f_nan, fcn_lnp: !values.f_nan, dcenp: !values.f_nan, $ 
              fcn_ltn: !values.f_nan, fcn_lnn: !values.f_nan, dcenn: !values.f_nan, $
              fcenxpp: !values.f_nan, fcenypp: !values.f_nan, dcenpp: !values.f_nan, $
              fcenxpn: !values.f_nan, fcenypn: !values.f_nan, dcenpn: !values.f_nan, $               
              dis: !values.f_nan, tilt: !values.f_nan, lp: !values.f_nan, dm: !values.f_nan, qm: !values.f_nan}  
  endelse
  

  if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin
  
      tmp_in = where((PRsFr_Big.mdi_i le bff_up) and (PRsFr_Big.mdi_i ge bff_dn), nprs)
      if (nprs ne 0) then begin
        PRsFr = PRsFr_Big[tmp_in]          
      endif else begin
        PRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
      endelse
            
      tmp_in = where((NRsFr_Big.mdi_i le bff_up) and (NRsFr_Big.mdi_i ge bff_dn), nprs)
      if (nprs ne 0) then begin
        NRsFr = NRsFr_Big[tmp_in]          
      endif else begin
        NRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
      endelse
      
   endif  
        
endif



;Using buffer for the first time
if keyword_set(buff) and (buff_sw eq 0) then begin
  
  buff_sw = 1
   
  print, 'Creating buffer' 
  mdi_imin = min([min(PRs[where(PRs.mdi_i ne 0)].mdi_i),min(NRs[where(NRs.mdi_i ne 0)].mdi_i)])
  mdi_imax = max([max(PRs[where(PRs.mdi_i ne 0)].mdi_i),max(NRs[where(NRs.mdi_i ne 0)].mdi_i)])  
  
  PRs_Big = PRs
  NRs_Big = NRs
  ARs_Big = ARs
  
  if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin
      PRsFr_Big = PRsFr
      NRsFr_Big = NRsFr      
  end
  
  bff_dn = mdi_imin
  bff_up = mdi_imin+bfr_ds
  
  PRs = PRs_Big[where(PRs_Big.mdi_i le bff_up)]
  NRs = NRs_Big[where(NRs_Big.mdi_i le bff_up)]
  ARs = ARs_Big[where(ARs_Big.mdi_i le bff_up)]
    
  if (n_elements(PRsFr_Big) gt 0) or (n_elements(NRsFr_Big) gt 0) then begin
  
      tmp_in = where((PRsFr_Big.mdi_i le bff_up), nprs)
      if (nprs ne 0) then begin
        PRsFr = PRsFr_Big[tmp_in]          
      endif else begin
        PRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
      endelse
            
      tmp_in = where((NRsFr_Big.mdi_i le bff_up), nprs)
      if (nprs ne 0) then begin
        NRsFr = NRsFr_Big[tmp_in]          
      endif else begin
        NRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
      endelse
      
   endif    
  
  ;Variables that keep track of the limits in which the buffer need to be extended
  mdi_bf_lm1 = mdi_imin
  mdi_bf_lm2 = mdi_imin + bfr_ds - bfr_lm
  
  
endif



;updating label
if (n_elements(labl) ne 0) then lbl = labl


undo = 1

;Initializing movement switches.  These mark the direction of the search for the next valid magnetograms
rw_msw =  1
lw_msw = -1


REPEAT BEGIN

  stat = 1   ;Variable that indicates the need for the program to end and save.

  ;Detection Window
  n_itr = 0
  if rw_msw ne 0 then begin   
    print, 'Looking for valid Working Window magnetogram'
    repeat begin
      
      ;Reading files
      caldat, mdi_ir+DayOff, Month, Day, Year
      dater = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
      print, dater
      mgr = amj_file_read( dater, hdrr, instr )
          
      sr = size(mgr)
    
      mdi_ir = mdi_ir + rw_msw
      n_itr = n_itr+1
      
    endrep until ( (sr[2] eq 8) or (n_itr gt 365) )
    
    mdi_ir = mdi_ir - rw_msw
    rw_msw = 0
    
    mdi_il = mdi_ir - 1 + dy_skp
    lw_msw = -1
  
  endif
  
  if (n_itr gt 365) then begin
     stat = 0
  endif
  
  
  ;Reference window
  n_itr = 0
  if lw_msw ne 0 then begin
    print, 'Looking for valid Reference Window magnetogram' 
    repeat begin
      ;Reading files 
      caldat, mdi_il+DayOff, Month, Day, Year
      datel = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
      print, datel
      mgl = amj_file_read( datel, hdrl, instr )
    
      sl = size(mgl)
    
      mdi_il = mdi_il + lw_msw
      n_itr = n_itr+1
      
    endrep until ( (sl[2] eq 8) or (n_itr gt 365) )
    
    mdi_il = mdi_il - lw_msw
    lw_msw = 0

  endif
      
  if (n_itr gt 365) then begin
     stat = 0
  endif

  if (stat eq 1) then begin

	  ;If doing a checkup run make sure that the detection of PNRs and ARs is deactivated
	  if keyword_set(rtrn_chk) then begin
		detpn_sw = 0
		detar_sw = 0
	  endif
  
      ;Performing detection of possitive and negative regions
      if ( (detpn_sw eq 1) ) then begin
        amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const;, /disp		
        amj_pnr_dt, CRD, mdi_ir, PRs, NRs, instr, seg_const=seg_const;, /disp,/not_merge, /detdisp

        detpn_sw = 0   
      endif
      
      ;Performing detection of active regions
      tmp_in = where((PRs.mdi_i eq mdi_ir) and (PRs.lnk_sw eq 0), n_prs)
      tmp_in = where((NRs.mdi_i eq mdi_ir) and (NRs.lnk_sw eq 0), n_nrs)
;      tmp_in = where((ARs.mdi_i eq mdi_ir), n_ars)     
;      detar_sw = 0;
    
;      if ( (detar_sw eq 1) and ( n_prs gt 0 ) and ( n_nrs gt 0 ) and ( n_ars eq 0 ) ) then begin
      if ( (detar_sw eq 1) and ( n_prs gt 0 ) and ( n_nrs gt 0 ) ) then begin
        amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const                       
        amj_ar_dt_track_dr, CRD, mdi_ir, lbl, PRs, NRs, ARs, mgr.date, ar_cnst=ar_cnst, seg_const=seg_const;,/display;, /bck_trck;, /display
        detar_sw = 0
        mdi_ir_vis = [mdi_ir_vis,mdi_ir]
      endif
      
      ;Definition of window parameters and window initialization
;      pls = 800;
      pls = 600;
      pad = 40;
      plxy = 100;
      plxx = 200;
      
      sz=size(mgr.img)
      d_zoom = pls/float(sz[1])
      
      window,0,xsize=2.0*pls+4.0*pad+plxx,ysize=pls+2.0*pad+plxy
      
      ;Exploration Window
      ;AR overlay
      if (ds_swl eq 1) then begin
        amj_mgplot, mgl.img, mdi_il, instr, hdr_i=hdrl, mdi_rf = mdi_ir, hdr_f=hdrr, ref_sw = ref_sw, lath = latlr, Lonh=Lonhrl, rchk_sw = rchk_sw, seg_const = seg_const, /tag, ARs = ARs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)]    
      endif
      ;Mixed overlay
      if (ds_swl eq 2) then begin
        amj_mgplot, mgl.img, mdi_il, instr, hdr_i=hdrl, mdi_rf = mdi_ir, hdr_f=hdrr, ref_sw = ref_sw, lath = latlr, Lonh=Lonhrl, rchk_sw = rchk_sw, seg_const = seg_const, /tag, PRs = PRs, NRs = NRs, ARs = ARs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)]    
      endif
      ;PNR overlay
      if (ds_swl eq 3) then begin
        amj_mgplot, mgl.img, mdi_il, instr, hdr_i=hdrl, mdi_rf = mdi_ir, hdr_f=hdrr, ref_sw = ref_sw, lath = latlr, Lonh=Lonhrl, rchk_sw = rchk_sw, seg_const = seg_const, /tag, PRs = PRs, NRs = NRs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)]        
      endif
      ;No overlay
      if (ds_swl eq 4) then begin
        amj_mgplot, mgl.img, mdi_il, instr, hdr_i=hdrl, mdi_rf = mdi_ir, hdr_f=hdrr, ref_sw = ref_sw, lath = latlr, Lonh=Lonhrl, seg_const = seg_const, /tag, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)]        
      endif
       
      ;Detection Window
      ;Mixed Overlay
      if (ds_swr eq 1) then begin
        amj_mgplot, mgr.img, mdi_ir, instr, hdr_i=hdrr, seg_const = seg_const, /tag, PRs = PRs, NRs = NRs, ARs = ARs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [pls+4.0*pad, pad+plxy, 2.0*pls+4.0*pad, pls+pad+plxy], title = dater, xrange = [min(mgr.x), max(mgr.x)] , yrange = [min(mgr.y), max(mgr.y)]    
      endif
    ;  ;AR overlay
    ;  if (ds_swr eq 2) then begin
    ;    amj_mgplot, mgr.img, mdi_ir, instr, ARs = ARs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [pls+4.0*pad, pad+plxy, 2.0*pls+4.0*pad, pls+pad+plxy], title = dater, xrange = [min(mgr.x), max(mgr.x)] , yrange = [min(mgr.y), max(mgr.y)]    
    ;  endif
      ;No overlay
      if (ds_swr eq 2) then begin
        amj_mgplot, mgr.img, mdi_ir, instr, hdr_i=hdrr, seg_const = seg_const, /tag, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [pls+4.0*pad, pad+plxy, 2.0*pls+4.0*pad, pls+pad+plxy], title = dater, xrange = [min(mgr.x), max(mgr.x)] , yrange = [min(mgr.y), max(mgr.y)]    
      endif
      ;Mixed Overlay with reference
      if (ds_swr eq 3) then begin
        amj_mgplot, mgr.img, mdi_ir, instr, hdr_i=hdrr, ref_sw = ref_sw, lath = latlr, Lonh=Lonhrr, seg_const = seg_const, /tag, PRs = PRs, NRs = NRs, ARs = ARs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [pls+4.0*pad, pad+plxy, 2.0*pls+4.0*pad, pls+pad+plxy], title = dater, xrange = [min(mgr.x), max(mgr.x)] , yrange = [min(mgr.y), max(mgr.y)]    
      endif
     
  endif

                       

  REPEAT BEGIN 
    redraw = 0
    print, 'Please make a choice'
     ;ltmpcs = 13
     ;CASE (tmpcs) of ;button_choice(pls, pad, plxy, plxx) OF
   CASE button_choice(pls, pad, plxy, plxx) OF
      0: BEGIN
            print, 'Please click inside one of the buttons'
          
         END
      1: BEGIN
            print, 'Quit'
            stat = 0
          
         END
      2: BEGIN
            print, 'Delete One'
            
            amj_pick_ar, mdi_ir, mdi_il, ARs, ar_in, pls, pad, plxy, plxx, d_zoom
            
            if (ar_in eq -2) then begin
              print, 'Cancelled. No ARs were erased'
            endif else begin
              
              ;finding corresponding PRs
              pr_ix = where(PRs.mdi_i eq ARs[ar_in].mdi_i)
              nr_ix = where(NRs.mdi_i eq ARs[ar_in].mdi_i)
              
              ;Find distances                            
              disp = sqrt((ARs[ar_in].fcn_ltp - PRs[pr_ix].fcn_lt)^2.0 + (ARs[ar_in].fcn_lnp - PRs[pr_ix].fcn_ln)^2.0)
              disn = sqrt((ARs[ar_in].fcn_ltn - NRs[nr_ix].fcn_lt)^2.0 + (ARs[ar_in].fcn_lnn - NRs[nr_ix].fcn_ln)^2.0)
              
              ;Find Minima
              minp = min(disp,inp)
              minn = min(disn,inn)
              
              ;Restoring indices for availability             
              PRs[pr_ix[inp]].lnk_sw = 0
              PRs[pr_ix[inp]].ar_lbl = 0
              NRs[nr_ix[inn]].lnk_sw = 0        
              NRs[nr_ix[inn]].ar_lbl = 0        

              ;Removing region              
              first = 0
              last = N_Elements(ARs)-1
              CASE ar_in OF
                first: ARs = ARs[1:*]
                last: ARs = ARs[first:last-1]
                ELSE: ARs = [ ARs[first:ar_in-1], ARs[ar_in+1:last] ]
              ENDCASE

              redraw = 1
              und_sw = 1               
            endelse         
         END         
      3: BEGIN
            print, 'Delete All'
            
            amj_pick_ar, mdi_ir, mdi_il, ARs, ar_in, pls, pad, plxy, plxx, d_zoom
            
            if (ar_in eq -2) then begin
              print, 'Cancelled. No ARs were erased'
            endif else begin
              
              tmp_lbl = ARs[ar_in].labl
              lblin = where(ARs.labl eq tmp_lbl, nars)
                            
              while (nars ne 0) do begin
              
                ;finding corresponding PRs
                pr_ix = where(PRs.mdi_i eq ARs[lblin[0]].mdi_i)
                nr_ix = where(NRs.mdi_i eq ARs[lblin[0]].mdi_i)
                
                ;Find distances                            
                disp = sqrt((ARs[lblin[0]].fcn_ltp - PRs[pr_ix].fcn_lt)^2.0 + (ARs[lblin[0]].fcn_lnp - PRs[pr_ix].fcn_ln)^2.0)
                disn = sqrt((ARs[lblin[0]].fcn_ltn - NRs[nr_ix].fcn_lt)^2.0 + (ARs[lblin[0]].fcn_lnn - NRs[nr_ix].fcn_ln)^2.0)
                
                ;Find Minima
                minp = min(disp,inp)
                minn = min(disn,inn)
                
                ;Restoring indices for availability             
                PRs[pr_ix[inp]].lnk_sw = 0
                PRs[pr_ix[inp]].ar_lbl = 0
                NRs[nr_ix[inn]].lnk_sw = 0        
                NRs[nr_ix[inn]].ar_lbl = 0        
  
                ;Removing region              
                first = 0
                last = N_Elements(ARs)-1
                CASE lblin[0] OF
                  first: ARs = ARs[1:*]
                  last: ARs = ARs[first:last-1]
                  ELSE: ARs = [ ARs[first:lblin[0]-1], ARs[lblin[0]+1:last] ]
                ENDCASE
              
                lblin = where(ARs.labl eq tmp_lbl, nars)
              
              endwhile

              redraw = 1
              und_sw = 1              
            endelse          
         END         
      4: BEGIN
            print, 'CREATE'
            
            ;finding current PRs
            pr_ix = where((PRs.mdi_i eq mdi_ir) and (PRs.lnk_sw eq 0), num_p)
            nr_ix = where((NRs.mdi_i eq mdi_ir) and (NRs.lnk_sw eq 0), num_n)            
            
            ;Checking if there are at least one region available of each kind
            if ( (num_p gt 0) and (num_n gt 0) ) then begin
            
              ;Detecting first region
              amj_pick_pnr, mdi_ir, PRs, NRs, r1_in, r1_sw, pls, pad, plxy, plxx, d_zoom
  
              if (r1_in eq -2) then begin
                print, 'Cancelled. No ARs were created'
              endif else begin
                
                ;Detecting second region
                amj_pick_pnr, mdi_ir, PRs, NRs, r2_in, r2_sw, pls, pad, plxy, plxx, d_zoom
                
                if (r2_in eq -2) then begin
                  print, 'Cancelled. No ARs were created'
                endif else begin
                  
                  ;Making sure they have opposite signs
                  if (r1_sw*r2_sw lt 0) then begin
                    
                    ;Choosing positive and negative parts
                    if (r1_sw gt 0) then begin
                      inp = r1_in
                      inn = r2_in
                    endif else begin
                      inp = r2_in
                      inn = r1_in
                    endelse
                    
                    ;Adding region to list, tracking, and updating PNR indices
					print, 'Creating'
                    amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const                                           
                    amj_ar_manual, CRD, mdi_ir, inp, inn, lbl, PRs, NRs, ARs, mgr.date, ar_cnst=ar_cnst, seg_const=seg_const
                    redraw = 1    
                    und_sw = 1                    
                  endif else begin
                    
                    print, 'Cancelled. Selected PNRs must be of opposite signs'
                    
                  endelse
                  
                  
                endelse
                
                
              endelse  

            endif else begin
              
              print, 'Cancelled. At least one free positive and negative regions are necessary'

            endelse                      
         END         
      5: BEGIN
            print, 'Create and Track'
            
            ;finding current PRs
            pr_ix = where((PRs.mdi_i eq mdi_ir) and (PRs.lnk_sw eq 0), num_p)
            nr_ix = where((NRs.mdi_i eq mdi_ir) and (NRs.lnk_sw eq 0), num_n)            
            
            ;Checking if there are at least one region available of each kind
            if ( (num_p gt 0) and (num_n gt 0) ) then begin
            
              ;Detecting first region
              amj_pick_pnr, mdi_ir, PRs, NRs, r1_in, r1_sw, pls, pad, plxy, plxx, d_zoom
  
              if (r1_in eq -2) then begin
                print, 'Cancelled. No ARs were created'
              endif else begin
                
                ;Detecting second region
                amj_pick_pnr, mdi_ir, PRs, NRs, r2_in, r2_sw, pls, pad, plxy, plxx, d_zoom
                
                if (r2_in eq -2) then begin
                  print, 'Cancelled. No ARs were created'
                endif else begin
                  
                  ;Making sure they have opposite signs
                  if (r1_sw*r2_sw lt 0) then begin
                    
                    ;Choosing positive and negative parts
                    if (r1_sw gt 0) then begin
                      inp = r1_in
                      inn = r2_in
                    endif else begin
                      inp = r2_in
                      inn = r1_in
                    endelse
                    
                    ;Adding region to list, tracking, and updating PNR indices
					print, 'Creating and Tracking'
                    amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const                                           
                    amj_ar_manual, CRD, mdi_ir, inp, inn, lbl, PRs, NRs, ARs, mgr.date, ar_cnst=ar_cnst, seg_const=seg_const, /track
                    redraw = 1    
                    und_sw = 1                    
                  endif else begin
                    
                    print, 'Cancelled. Selected PNRs must be of opposite signs'
                    
                  endelse
                  
                  
                endelse
                
                
              endelse  

            endif else begin
              
              print, 'Cancelled. At least one free positive and negative regions are necessary'

            endelse                      
         END   
      6: BEGIN
            print, 'Re-label'
            
            print, 'Select AR to be labeled'
            amj_pick_ar, mdi_ir, mdi_il, ARs, ar_in1, pls, pad, plxy, plxx, d_zoom
            
            if (ar_in1 eq -2) then begin
              print, 'Cancelled. No ARs were re-labeled'
            endif else begin

              print, 'Select reference AR'
              amj_pick_ar, mdi_ir, mdi_il, ARs, ar_in2, pls, pad, plxy, plxx, d_zoom
              
              if (ar_in2 eq -2) then begin
                print, 'Cancelled. No ARs were re-labeled'
              endif else begin
                

                lblin = where(PRs.ar_lbl eq ARs[ar_in1].labl)
                PRs[lblin].ar_lbl = ARs[ar_in2].labl
                
                lblin = where(NRs.ar_lbl eq ARs[ar_in1].labl)
                NRs[lblin].ar_lbl = ARs[ar_in2].labl
                
                lblin = where(ARs.labl eq ARs[ar_in1].labl) 
                ARs[lblin].labl = ARs[ar_in2].labl
                ARs[lblin].clr = ARs[ar_in2].clr
               


                redraw = 1
                und_sw = 1                 
              endelse            
              
            endelse                     
         END         
      7: BEGIN
            print, 'Merge'

            ;finding current PRs
            pr_ix = where((PRs.mdi_i eq mdi_ir) and (PRs.lnk_sw eq 0), num_p)
            nr_ix = where((NRs.mdi_i eq mdi_ir) and (NRs.lnk_sw eq 0), num_n)            
            
            ;Detecting first region
            amj_pick_pnr, mdi_ir, PRs, NRs, r1_in, r1_sw, pls, pad, plxy, plxx, d_zoom, /all

            if (r1_in eq -2) then begin
              print, 'Cancelled. No regions were merged'
            endif else begin
              
              ;Detecting second region
              amj_pick_pnr, mdi_ir, PRs, NRs, r2_in, r2_sw, pls, pad, plxy, plxx, d_zoom, /all
              
              if (r2_in eq -2) then begin
                print, 'Cancelled. No regions were merged'
              endif else begin
                
                ;Making sure they have the same sign
                if (r1_sw*r2_sw gt 0) then begin
                  
                  ;Positive regions
                  if (r1_sw gt 0) then begin
                    
                    pl_sw = 1
                    ;Making sure at most one of them belongs to an active region
                    if( (PRs[r1_in].lnk_sw eq 1) and (PRs[r2_in].lnk_sw eq 1) ) then begin
                      
                      print, 'Cancelled. Selected PNRs cannot be both part of active regions'
                      
                    endif else begin
                    
                      ;Choosing the region that will remain based on which of them belongs to an active region
                      if (PRs[r1_in].lnk_sw eq 1) then begin
                        inar  = r1_in
                        innar = r2_in
                      endif else begin
                        inar = r2_in
                        innar = r1_in
                      endelse
                      
                      ;Adding region to list and updating PNR indices
					  print, 'Merging'
                      amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const                                           
                      amj_pnr_merge, CRD, mdi_ir, inar, innar, PRs, ARs, pl_sw, date


                      redraw = 1    
                      und_sw = 1
                      
                    endelse
                  
                  ;Negative regions
                  endif else begin
                    
                    pl_sw = -1
                    ;Making sure at most one of them belongs to an active region
                    if( (NRs[r1_in].lnk_sw eq 1) and (NRs[r2_in].lnk_sw eq 1) ) then begin
                      
                      print, 'Cancelled. Selected PNRs cannot be both part of active regions'
                      
                    endif else begin
                    
                      ;Choosing the region that will remain based on which of them belongs to an active region
                      if (NRs[r1_in].lnk_sw eq 1) then begin
                        inar  = r1_in
                        innar = r2_in
                      endif else begin
                        inar = r2_in
                        innar = r1_in
                      endelse
                      
                      ;Adding region to list and updating PNR indices
                      amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const                                           
                      amj_pnr_merge, CRD, mdi_ir, inar, innar, NRs, ARs, pl_sw, date


                      redraw = 1    
                      und_sw = 1
                      
                    endelse
                    
                    
                  endelse
                  
                                      
                endif else begin
                  
                  print, 'Cancelled. Selected PNRs must be of the same sign'
                  
                endelse
                
                
              endelse
              
              
            endelse  
          
         END
      8: BEGIN
            print, 'Fragment'

            ;finding current PRs
            pr_ix = where((PRs.mdi_i eq mdi_ir) and (PRs.lnk_sw eq 0), num_p)
            nr_ix = where((NRs.mdi_i eq mdi_ir) and (NRs.lnk_sw eq 0), num_n)            
            
            ;Detecting region
            amj_pick_pnr, mdi_ir, PRs, NRs, r1_in, r1_sw, pls, pad, plxy, plxx, d_zoom

            if (r1_in eq -2) then begin
              print, 'Cancelled. No PNR was fragmented'
            endif else begin
			
			  print, 'Fragmenting'
              prfrg_sw = 0 ;Swiths that indicates that a prefragmented set of regions has been found

              ;Creating temporary PNR arrays
              tmp_PRs = {rgn,  lnk_sw: 0,  mdi_i: 0L, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}    
              tmp_NRs = tmp_PRs

                    
              ;Positive regions                    
              if (r1_sw gt 0) and (n_elements(PRsFr) gt 0) then begin
                  
                  ;Finding if there are associated pre-fragmented regions that day
                  frin = where((PRsFr.mdi_i eq PRs[r1_in].mdi_i) and (PRsFr.fr_lbl eq PRs[r1_in].fr_lbl), n_fr)
                  
                  if n_fr gt 0 then begin
                    
                      tmp_PRs = [tmp_PRs, PRsFr[frin]]
                      prfrg_sw = 1
                    
                  endif                
                
              endif

              ;Negative regions                    
              if (r1_sw lt 0) and (n_elements(NRsFr) gt 0) then begin
                  
                  ;Finding if there are associated pre-fragmented regions that day
                  frin = where((NRsFr.mdi_i eq NRs[r1_in].mdi_i) and (NRsFr.fr_lbl eq NRs[r1_in].fr_lbl), n_fr)
                  
                  if n_fr gt 0 then begin
                    
                       tmp_NRs = [tmp_NRs, NRsFr[frin]]
                       prfrg_sw = 1
                    
                  endif                
                
              endif              
                                                    
              if (prfrg_sw eq 0) then begin   
                          
                  ;Extracting necessary pixels
                  if (r1_sw gt 0) then begin
                    ex_in = long(strsplit(PRs[r1_in].indx,/extract))
                  endif else begin
                    ex_in = long(strsplit(NRs[r1_in].indx,/extract))
                  endelse
                  
                  ;Creating temporary image
                  tmp_im = mgr.img*0.0
                  tmp_im[ex_in] = mgr.img[ex_in]
                  
                  
                  ;Creating temporary seg_const
                  tmp_seg_const = seg_const
                  tmp_seg_const.npssu = 2
          			  tmp_seg_const.ker_th=200
          			  tmp_seg_const.ker_th2=150
          			  tmp_seg_const.dis_lim = 2
                  tmp_seg_const.ovr_lim = 0.7
          			  ;tmp_seg_const.dila_size=5
          			  
                  amj_coord, tmp_im, hdrr, CRD              
                  amj_pnr_dt, CRD, mdi_ir, tmp_PRs, tmp_NRs, seg_const=tmp_seg_const, pnr_lbl = -999;,/not_merge;, /detdisp

              endif
              
              
              if ( ((n_elements(tmp_PRs) gt 1) and (r1_sw gt 0)) or ((n_elements(tmp_NRs) gt 1) and (r1_sw lt 0)) ) then begin

                  ;Positive regions
                  if (r1_sw gt 0) then begin
                    
                    ;Removing fragmented region and adding the new ones in its place
                    first = 0
                    last = N_Elements(PRs)-1
                    CASE r1_in OF
                      first: PRs = [tmp_PRs, PRs[1:*]]
                      last:  PRs = [PRs[first:last-1], tmp_PRs[1:*]]
                      ELSE:  PRs = [PRs[first:r1_in-1], tmp_PRs[1:*], PRs[r1_in+1:last] ]
                    ENDCASE
                    
                  endif else begin
                    
                    ;Removing fragmented region and adding the new ones in its place
                    first = 0
                    last = N_Elements(NRs)-1
                    CASE r1_in OF
                      first: NRs = [tmp_NRs, NRs[1:*]]
                      last:  NRs = [NRs[first:last-1], tmp_NRs[1:*]]
                      ELSE:  NRs = [NRs[first:r1_in-1], tmp_NRs[1:*], NRs[r1_in+1:last] ]
                    ENDCASE
                    
                  endelse
                  
                  und_sw = 1  
                  redraw = 1
                               
              endif else begin
                
                  print, 'Cancelled. Region could not be fragmented'
                
              endelse
                                  
    

                                   
            endelse  
          
         END
      9: BEGIN
            print, 'Re-detect'
                          
            ;Finding Regions not present in the detection window
            pr_ix = where(PRs.mdi_i ne mdi_ir)
            PRs = PRs[pr_ix]
            
            nr_ix = where(NRs.mdi_i ne mdi_ir)                        
            NRs = NRs[nr_ix]              

            ar_ix = where(ARs.mdi_i ne mdi_ir)
            ARs = ARs[ar_ix]              
              
            detpn_sw  = 1; %PNR Detection switch
            detar_sw  = 1; %AR Detection switch
            redraw = 1
            und_sw = 1
         END         
     10: BEGIN
            print, 'Synchronize'
            mdi_il = mdi_ir - 1
            redraw = 1
			lw_msw = -1			
         END		 
     11: BEGIN
            print, 'Contrast Up'
            max_sat = max_sat - incrmnt
            if (max_sat le min_max_sat) then begin
              max_sat = min_max_sat
              print, 'Maximum contrast reached'
            endif
            print, 'max_sat = ', max_sat
            redraw = 1
  		   END
     12: BEGIN
            print, 'Contrast Down'
            max_sat = max_sat + incrmnt
            if (max_sat ge max_max_sat) then begin
              max_sat = max_max_sat 
             	print, 'Minimum contrast reached'
            endif
            print, 'max_sat = ', max_sat
            redraw = 1
         END
     13: BEGIN
            print, 'Print Screen'
;            tmp_fl = 'prnt_' + strtrim(string(mdi_ir,'(I)'),2) + '_' + strtrim(string(seg_const.ker_th,'(I)'),2) + '_' + strtrim(string(seg_const.ker_th2,'(I)'),2) + '_' + strtrim(string(seg_const.ar_th,'(I)'),2) + '_' + strtrim(string(seg_const.eros_size,'(I)'),2) + '_' + strtrim(string(seg_const.dila_size,'(I)'),2) + '_' + strtrim(string(seg_const.dis_lim,'(I)'),2) + '_' + strtrim(string(seg_const.ovr_lim,'(F3.1)'),2) + '.png'
            tmp_fl = 'prnt_' + strtrim(string(mdi_ir,'(I)'),2) + '_' + strtrim(string(seg_const.ker_th,'(I)'),2) + '_' + strtrim(string(seg_const.eros_size,'(I)'),2) + '_' + strtrim(string(seg_const.dila_size,'(I)'),2) + '_'+ strtrim(string(seg_const.ker_th2,'(I)'),2) + '_' + strtrim(string(seg_const.ar_th,'(I)'),2) + '_'  + strtrim(string(seg_const.dis_lim,'(I)'),2) + '_' + strtrim(string(seg_const.ovr_lim,'(F3.1)'),2) + '.png'
;            tmp_fl = 'prnt_AR_not_merge' + strtrim(string(mdi_ir,'(I)'),2) + '_' + strtrim(string(ar_cnst.dis_lim1,'(F4.1)'),2) + '_' + strtrim(string(ar_cnst.MxFlxim,'(F4.1)'),2) + '_' + strtrim(string(ar_cnst.exp_f,'(F4.1)'),2) + '_' + strtrim(string(ar_cnst.exp_d,'(F4.1)'),2) + '_'+ strtrim(string(ar_cnst.exp_s,'(F4.1)'),2) + '_' + strtrim(string(seg_const.ar_th,'(I)'),2) + '.png'
            write_png, tmp_fl, TVRD(/true)
            stat = 0
            sv_sw = 0
         END	  
     14: BEGIN
            print, 'Longitudinal Guide'
			amj_pick_long, hdrr, hdrl, ARs, latlr, Lonhrl, Lonhrr, pls, pad, plxy, plxx, d_zoom, instr
			
			;print, Lonhrl
			;print, Lonhrr
			
			ref_sw = 1
            ds_swr = 3
            redraw = 1			
			
         END
     15: BEGIN
            print, '<< Reference'
            mdi_il = mdi_il - iskip
            redraw =  1
            lw_msw = -1
			
			;Hiding reference lines
			if (ds_swr gt 2) then begin
				ds_swr = 1
				ref_sw = 0
			endif			
			
         END
     16: BEGIN
            print, '< Reference'
            mdi_il = mdi_il - 1
            redraw =  1
            lw_msw = -1
			
			;Hiding reference lines
			if (ds_swr gt 2) then begin
				ds_swr = 1
				ref_sw = 0
			endif			
			
         END
  	 17: BEGIN
  	        print, 'Jump to'
      			date = ''
      			read, 'Enter date (e.g. 2001-08-31): ', date
      			
      			mg_tmp = amj_file_read( date, hdrr, instr)
            mdi_tmp  = julday(strmid(date,5,2),strmid(date,8,2),strmid(date,0,4))-DayOff
      			
      			tmpsize = size(mg_tmp)
      			if (tmpsize[2] ne 8) then begin
      				print, 'No available MDI magnetogram.'
      			endif else begin
      				mdi_il = mdi_tmp
      				redraw = 1
      				lw_msw = 1
					
					;Hiding reference lines
					if (ds_swr gt 2) then begin
						ds_swr = 1
						ref_sw = 0
					endif					
					
      			endelse  			
  		   END
     18: BEGIN
            print, '> Reference'
            mdi_il = mdi_il + 1
            redraw = 1
            lw_msw = 1
			
			;Hiding reference lines
			if (ds_swr gt 2) then begin
				ds_swr = 1
				ref_sw = 0
			endif			
			
         END
     19: BEGIN
            print, '>> Reference'
            mdi_il = mdi_il + iskip
            redraw = 1 
            lw_msw = 1

			;Hiding reference lines
			if (ds_swr gt 2) then begin
				ds_swr = 1
				ref_sw = 0
			endif			
			
         END
     20: BEGIN
            print, '< Control'
            mdi_ir = mdi_ir - 1
            mdi_il = mdi_il - 1 + dy_skp
            redraw =  1
            rw_msw = -1 

			;Hiding reference lines
			if (ds_swr gt 2) then begin
				ds_swr = 1
				ref_sw = 0
			endif
				
         END
     21: BEGIN
            print, 'Jump to'
            date = ''
            read, 'Enter date (e.g. 2001-08-31): ', date

            mg_tmp = amj_file_read( date, hdrr, instr)
            mdi_tmp  = julday(strmid(date,5,2),strmid(date,8,2),strmid(date,0,4))-DayOff
            
            tmpsize = size(mg_tmp)
            if (tmpsize[2] ne 8) then begin
              print, 'No available MDI magnetogram.'
            endif else begin
				mdi_ir = mdi_tmp
				mdi_il = mdi_ir - 1 + dy_skp
				redraw =  1
				rw_msw =  1
				lw_msw = -1
			  
				;Hiding reference lines
				if (ds_swr gt 2) then begin
					ds_swr = 1
					ref_sw = 0
				endif			  
			  
            endelse       
         END
     22: BEGIN
            print, '> Control'
            mdi_ir = mdi_ir + 1			
            mdi_il = mdi_il - 1 + dy_skp
            sv_sw = 1
            
            if (total(mdi_ir_vis eq mdi_ir) eq 0) then begin
              detar_sw = 1
              if (prepnr_sw eq 0) then detpn_sw = 1
              und_sw = 1
              
;              if not keyword_set(pnr_file) then begin
;                ;Ensuring there are not repeated PNR detections
;                pr_ix = where(PRs.mdi_i ne mdi_ir)
;                PRs = PRs[pr_ix]
;                
;                nr_ix = where(NRs.mdi_i ne mdi_ir)                        
;                NRs = NRs[nr_ix]                            
;              endif
              
            endif
			
			;Hiding reference lines
            if (ds_swr gt 2) then begin
				ds_swr = 1
				ref_sw = 0
			endif
			
            
            redraw = 1
            rw_msw = 1          
          
         END
     23: BEGIN
            print, 'Left Overlay'
            ds_swl = ds_swl + 1
            if (ds_swl gt 4) then ds_swl = 1
            redraw = 1                    
         END
     24: BEGIN
            print, 'Right Overlay'
            ds_swr = ds_swr + 1			
			;Hiding reference lines
            if (ds_swr gt 2) then begin
				ds_swr = 1
				ref_sw = 0
			endif
            redraw = 1
			
         END
    ENDCASE
    
;    if und_sw eq 1 then begin  
;      	spawn, '\cp -f Undo_State.sav Undo_StateOld.sav'    	
;        SAVE, ARs, PRs, NRs, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, FILENAME = 'Undo_State.sav'
;        und_sw = 0
;  	endif
  	
    if (buff_sw eq 0) then begin 
      
;        if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin         
;            SAVE, ARs, PRs, NRs, PRsFr, NRsFr, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, FILENAME = 'Ar_id_Save0.sav'
;        endif else begin
;            SAVE, ARs, PRs, NRs, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, FILENAME = 'Ar_id_Save0.sav'                    
;        endelse 
      
        if ( ((mdi_ir mod sv_bfr) eq 0) and (sv_sw eq 1) ) then begin
            spawn, '\cp -f Ar_id_Save0.sav Ar_id_Save_Bckp.sav'
                  
            if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin         
                SAVE, ARs, PRs, NRs, PRsFr, NRsFr, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, instr, FILENAME = 'Ar_id_Save0.sav'        
            endif else begin
                SAVE, ARs, PRs, NRs, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, instr, FILENAME = 'Ar_id_Save0.sav'                    
            endelse
        endif

    ;using buffer
    endif else begin
      ;If either of the limits is overtaken
      if ( (mdi_ir lt mdi_bf_lm1) or (mdi_ir gt mdi_bf_lm2) ) then begin
        
        
        print, 'Updating Buffer'    
            
        ;Updating the big files
        ;----------------------
              
        ;Positive regions
        tmp_in = where(PRs_Big.mdi_i lt bff_dn, nprs1 ) 
        tmp_in2 = where(PRs_Big.mdi_i gt bff_up, nprs2 ) 
        
        if ( (nprs1 eq 0) and (nprs2 eq 0) ) then begin
          PRs_Big = PRs
        endif else if (nprs1 eq 0) then begin
          PRs_Big = [PRs, PRs_Big[tmp_in2]]
        endif else if (nprs2 eq 0) then begin
          PRs_Big = [PRs_Big[tmp_in], PRs]
        endif else begin
          PRs_Big = [PRs_Big[tmp_in], PRs, PRs_Big[tmp_in2]]
        endelse
              
  
        ;Negative regions
        tmp_in = where(NRs_Big.mdi_i lt bff_dn, nprs1 ) 
        tmp_in2 = where(NRs_Big.mdi_i gt bff_up, nprs2 ) 
        
        if ( (nprs1 eq 0) and (nprs2 eq 0) ) then begin
          NRs_Big = NRs
        endif else if (nprs1 eq 0) then begin
          NRs_Big = [NRs, NRs_Big[tmp_in2]]
        endif else if (nprs2 eq 0) then begin
          NRs_Big = [NRs_Big[tmp_in], NRs]
        endif else begin
          NRs_Big = [NRs_Big[tmp_in], NRs, NRs_Big[tmp_in2]]
        endelse  
  
        ;Active regions
        tmp_in = where(ARs_Big.mdi_i lt bff_dn, nprs1 ) 
        tmp_in2 = where(ARs_Big.mdi_i gt bff_up, nprs2 ) 
        
        if ( (nprs1 eq 0) and (nprs2 eq 0) ) then begin
          ARs_Big = ARs
        endif else if (nprs1 eq 0) then begin
          ARs_Big = [ARs, ARs_Big[tmp_in2]]
        endif else if (nprs2 eq 0) then begin
          ARs_Big = [ARs_Big[tmp_in], ARs]
        endif else begin
          ARs_Big = [ARs_Big[tmp_in], ARs, ARs_Big[tmp_in2]]
        endelse  
            
        ;Updating Variables that keep track of the limits in which the buffer need to be extended
        bff_dn = mdi_ir - 2*bfr_lm
        bff_up = mdi_ir - 2*bfr_lm + bfr_ds
    
        mdi_bf_lm1 = bff_dn + bfr_lm
        mdi_bf_lm2 = bff_up - bfr_lm        
  
        ;Updating buffer        
        tmp_in = where((PRs_Big.mdi_i le bff_up) and (PRs_Big.mdi_i ge bff_dn), nprs)
        if (nprs ne 0) then begin
          PRs = PRs_Big[tmp_in]          
        endif else begin
          PRs = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
        endelse
        
        tmp_in = where((NRs_Big.mdi_i le bff_up) and (NRs_Big.mdi_i ge bff_dn), nprs)
        if (nprs ne 0) then begin
          NRs = NRs_Big[tmp_in]          
        endif else begin
          NRs = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
        endelse
 
        tmp_in = where((ARs_Big.mdi_i le bff_up) and (ARs_Big.mdi_i ge bff_dn), nprs)
        if (nprs ne 0) then begin
          ARs = ARs_Big[tmp_in]          
        endif else begin
          ARs = {ar, mdi_i: 0, date: '', labl: 0, clr: 0, indxp: '', indxn: '', fluxp: !values.f_nan, fluxn: !values.f_nan, areap:!values.f_nan, arean:!values.f_nan, $
              fcn_ltp: !values.f_nan, fcn_lnp: !values.f_nan, dcenp: !values.f_nan, $ 
              fcn_ltn: !values.f_nan, fcn_lnn: !values.f_nan, dcenn: !values.f_nan, $
              fcenxpp: !values.f_nan, fcenypp: !values.f_nan, dcenpp: !values.f_nan, $
              fcenxpn: !values.f_nan, fcenypn: !values.f_nan, dcenpn: !values.f_nan, $               
              dis: !values.f_nan, tilt: !values.f_nan, lp: !values.f_nan, dm: !values.f_nan, qm: !values.f_nan}  
        endelse

        if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin
            ;Updating buffer        
            tmp_in = where((PRsFr_Big.mdi_i le bff_up) and (PRsFr_Big.mdi_i ge bff_dn), nprs)
            if (nprs ne 0) then begin
              PRsFr = PRsFr_Big[tmp_in]          
            endif else begin
              PRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
            endelse
            
            tmp_in = where((NRsFr_Big.mdi_i le bff_up) and (NRsFr_Big.mdi_i ge bff_dn), nprs)
            if (nprs ne 0) then begin
              NRsFr = NRsFr_Big[tmp_in]          
            endif else begin
              NRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}          
            endelse

        endif


      endif
      
      
      ;saving periodically
      if ( ((mdi_ir mod sv_bfr) eq 0) and (sv_sw eq 1) ) then begin
      
      
        print, 'Saving run'
        
        
        
        ;Creating buffer backup------------------------

        spawn, '\cp -f Ar_id_Save0.sav Ar_id_Save_Buff_Bckp.sav'    
        
        
        
        tmpPRsb = PRs
        tmpNRsb = NRs
        tmpARsb = ARs
        
        PRs = PRs_Big
        NRs = NRs_Big
        ARs = ARs_Big
 
        if (n_elements(PRsFr_Big) gt 0) or (n_elements(NRsFr_Big) gt 0) then begin            
            tmpPRsFrb = PRsFr
            tmpNRsFrb = NRsFr
    
            PRsFr = PRsFr_Big
            NRsFr = NRsFr_Big
        endif
        
        ;Saving file
        if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin
            SAVE, ARs, PRs, NRs, PRsFr, NRsFr, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, instr, FILENAME = 'Ar_id_Save0.sav'          
        endif else begin  
            SAVE, ARs, PRs, NRs, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, instr, FILENAME = 'Ar_id_Save0.sav'
        endelse

        PRs = tmpPRsb
        NRs = tmpNRsb
        ARs = tmpARsb

        if (n_elements(PRsFr_Big) gt 0) or (n_elements(NRsFr_Big) gt 0) then begin
            PRsFr = tmpPRsFrb
            NRsFr = tmpNRsFrb
        endif
        
        sv_sw = 0

        
      endif
      
    endelse 
	
  ENDREP UNTIL( redraw or ( stat eq 0 ) )
ENDREP UNTIL ( stat EQ 0 )

;Final Save in case Buffer is being used
if (buff_sw eq 1) then begin
  
  ;Updating the big files
  ;----------------------
        
  ;Positive regions
  tmp_in = where(PRs_Big.mdi_i lt bff_dn, nprs1 ) 
  tmp_in2 = where(PRs_Big.mdi_i gt bff_up, nprs2 ) 
  
  if ( (nprs1 eq 0) and (nprs2 eq 0) ) then begin
    PRs_Big = PRs
  endif else if (nprs1 eq 0) then begin
    PRs_Big = [PRs, PRs_Big[tmp_in2]]
  endif else if (nprs2 eq 0) then begin
    PRs_Big = [PRs_Big[tmp_in], PRs]
  endif else begin
    PRs_Big = [PRs_Big[tmp_in], PRs, PRs_Big[tmp_in2]]
  endelse
        

  ;Negative regions
  tmp_in = where(NRs_Big.mdi_i lt bff_dn, nprs1 ) 
  tmp_in2 = where(NRs_Big.mdi_i gt bff_up, nprs2 ) 
  
  if ( (nprs1 eq 0) and (nprs2 eq 0) ) then begin
    NRs_Big = NRs
  endif else if (nprs1 eq 0) then begin
    NRs_Big = [NRs, NRs_Big[tmp_in2]]
  endif else if (nprs2 eq 0) then begin
    NRs_Big = [NRs_Big[tmp_in], NRs]
  endif else begin
    NRs_Big = [NRs_Big[tmp_in], NRs, NRs_Big[tmp_in2]]
  endelse  

  ;Active regions
  tmp_in = where(ARs_Big.mdi_i lt bff_dn, nprs1 ) 
  tmp_in2 = where(ARs_Big.mdi_i gt bff_up, nprs2 ) 
  
  if ( (nprs1 eq 0) and (nprs2 eq 0) ) then begin
    ARs_Big = ARs
  endif else if (nprs1 eq 0) then begin
    ARs_Big = [ARs, ARs_Big[tmp_in2]]
  endif else if (nprs2 eq 0) then begin
    ARs_Big = [ARs_Big[tmp_in], ARs]
  endif else begin
    ARs_Big = [ARs_Big[tmp_in], ARs, ARs_Big[tmp_in2]]
  endelse  

      
  ;Updating Variables that keep track of the limits in which the buffer need to be extended
  bff_dn = mdi_ir - 2*bfr_lm
  bff_up = mdi_ir - 2*bfr_lm + bfr_ds

  mdi_bf_lm1 = bff_dn + bfr_lm
  mdi_bf_lm2 = bff_up - bfr_lm
  
  PRs = PRs_Big
  NRs = NRs_Big
  ARs = ARs_Big
  
  if (n_elements(PRsFr_Big) gt 0) or (n_elements(NRsFr_Big) gt 0) then begin
      PRsFr = PRsFr_Big
      NRsFr = NRsFr_Big
  endif  
    
endif

;Saving file
if (n_elements(PRsFr) gt 0) or (n_elements(NRsFr) gt 0) then begin
    SAVE, ARs, PRs, NRs, PRsFr, NRsFr, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, instr, FILENAME = 'Ar_id_Save0.sav'
endif else begin
    SAVE, ARs, PRs, NRs, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, instr, FILENAME = 'Ar_id_Save0.sav'
endelse



print, 'Next Active Region Label is ', lbl

;SAVE, ARs, FILENAME = 'ARs0.sav'

return
END


;----------------------------------------------------------------------------------------------------------

 
