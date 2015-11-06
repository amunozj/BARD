PRO amj_mgplot, im, mdi_i, instr, hdr_i=hdr_i, seg_const=seg_const, PRs = PRs, NRs = NRs, ARs = ARs, tag = tag, d_xsize = d_xsize, d_ysize = d_ysize, max = max, pos = pos, title = title, xrange = xrange, yrange = yrange, prt=prt, sqs_nm = sqs_nm, shw_lbl = shw_lbl, lath = lath, Lonh=Lonh, ref_sw = ref_sw


    ;Constant with parameters for plotting, make sure the values are the same as in kpvt_pnr_dt.rpo
    seg_const_def={k_sig:15.0, valid_range:[-20000.,20000.], deg_lim:70.0}
    
    if not keyword_set(seg_const) then begin
        seg_const=seg_const_def
    endif


;    plot,[1,1],/nodata,xstyle=5,ystyle=5
    loadct, 0, /silent
    tv,bytscl(congrid(im,d_xsize,d_ysize),min=-max,max=max), pos[0], pos[1] 
    loadct, 0, /silent

    
    sz=size(im)
    d_zoom = d_xsize/float(sz[1])
    
    ;plotting positive regions
    if keyword_set(PRs) then begin
      
        if keyword_set(ARs) then begin
            ;finding Available PRs
            pr_ix = where((PRs.mdi_i eq mdi_i) and (PRs.lnk_sw eq 0), n_regp)
        endif else begin
            ;finding All PRs
            pr_ix = where(PRs.mdi_i eq mdi_i, n_regp)          
        endelse
        
        if (n_regp gt 0) then begin
            pregions = PRs[pr_ix]

            for i=0,n_regp-1 do begin
                if keyword_set(shw_lbl) then begin
                    str='P'+strtrim(string((i+1) mod 100,'(I2)'),2)
                    xyouts,pregions[i].fcenxp*d_zoom,pregions[i].fcenyp*d_zoom,str,charsize=1.5, color=255,/device
                endif
                tvcircle,2.0*pregions[i].dcenp*d_zoom, pregions[i].fcenxp*d_zoom + pos[0], pregions[i].fcenyp*d_zoom + pos[1], color=255, /device, thick = 1              
            endfor
        endif
        
    endif

    ;plotting available negative regions
    if keyword_set(NRs) then begin
      
        if keyword_set(ARs) then begin
            ;finding Available NRs
            nr_ix = where((NRs.mdi_i eq mdi_i) and (NRs.lnk_sw eq 0), n_regn)
        endif else begin
            ;finding All NRs
            nr_ix = where(NRs.mdi_i eq mdi_i, n_regn)          
        endelse
        
        if (n_regn gt 0) then begin
            nregions = NRs[nr_ix]

            for i=0,n_regn-1 do begin
				if keyword_set(shw_lbl) then begin
               	 str='N'+strtrim(string((i+1) mod 100,'(I2)'),2)
               	 xyouts,nregions[i].fcenxp*d_zoom,nregions[i].fcenyp*d_zoom,str,charsize=1.5, color=0,/device
				endif
                tvcircle,2.0*nregions[i].dcenp*d_zoom, nregions[i].fcenxp*d_zoom + pos[0], nregions[i].fcenyp*d_zoom + pos[1], color=0,/device, thick = 1
            endfor
        endif
        
    endif
    
    
    if keyword_set(ARs) then begin
      
        arin = where(ARs.mdi_i eq mdi_i, nars)        
        if (nars gt 0) then begin

			set_plot,'X'
            tmpArs = ARs[arin]
            
;            print, 'ARs.labl',tmpArs.labl
;            print, 'ARs.fcn_ltp',tmpArs.fcn_ltp
;            print, 'ARs.fcn_lnp',tmpArs.fcn_lnp
;            print, 'ARs.fcn_ltn',tmpArs.fcn_ltn
;            print, 'ARs.fcn_lnn',tmpArs.fcn_lnn
;            
;            print, ' '
;            
;            pr_ix = where(PRs.mdi_i eq mdi_i, n_regp)            
;            print, 'PRs.ar_lbl',PRs[pr_ix].ar_lbl
;            print, 'PRs.lnk_sw',PRs[pr_ix].lnk_sw
;            print, 'PRs.fcn_lt',PRs[pr_ix].fcn_lt
;            print, 'PRs.fcn_ln',PRs[pr_ix].fcn_ln
;
;            print, ' '
;            
;            nr_ix = where(NRs.mdi_i eq mdi_i, n_regn)
;            print, 'NRs.ar_lbl',NRs[nr_ix].ar_lbl
;            print, 'NRs.lnk_sw',NRs[nr_ix].lnk_sw
;            print, 'NRs.fcn_lt',NRs[nr_ix].fcn_lt
;            print, 'NRs.fcn_ln',NRs[nr_ix].fcn_ln
                        
            for i=0,nars-1 do begin
                 
                loadct, 13, /silent
                
                if keyword_set(tag) then begin
                  str='P'+strtrim(string((tmpArs[i].labl)mod 100,'(I2)'),2)
                  xyouts,tmpArs[i].fcenxpp*d_zoom + pos[0], tmpArs[i].fcenypp*d_zoom + pos[1],str,charsize=1.5,color= long(tmpArs[i].clr),/device
                endif
                tvcircle,2.0*tmpArs[i].dcenpp*d_zoom, tmpArs[i].fcenxpp*d_zoom + pos[0], tmpArs[i].fcenypp*d_zoom + pos[1], color = strtrim(string(tmpArs[i].clr,'(I3)'),2),/device, thick = 2
                
                if keyword_set(tag) then begin
                  str='N'+strtrim(string((tmpArs[i].labl)mod 100,'(I2)'),2)
                  xyouts,tmpArs[i].fcenxpn*d_zoom + pos[0], tmpArs[i].fcenypn*d_zoom + pos[1],str,charsize=1.5,color= long(tmpArs[i].clr),/device
                endif
                tvcircle,2.0*tmpArs[i].dcenpn*d_zoom, tmpArs[i].fcenxpn*d_zoom + pos[0], tmpArs[i].fcenypn*d_zoom + pos[1], color = strtrim(string(tmpArs[i].clr,'(I3)'),2),/device, thick = 2
        
                loadct, 0, /silent
            endfor    
        endif 
    
    endif
    
    if keyword_set(Lath) and keyword_set(Lonh) and keyword_set(hdr_i) and keyword_set(ref_sw) then begin
     
        ;
        ;Calculation of heliospheric coordinates
		
		;Magnetogram being draw-------------------------------------------------------
		
		;HMI uses structures for header values
		if instr eq 4 then begin
			datei = hdr_i.DATE_OBS

			;Define center and radius
			hfxi = hdr_i.CRPIX1 ;  Location of the center in x pixels 
			hfyi = hdr_i.CRPIX2 ;    Location of the center in y pixels
			dii = hdr_i.RSUN_OBS/hdr_i.CDELT1;

			;Load Solar Coordinates
			P0i = 0.0
			RDi = hdr_i.DSUN_OBS/hdr_i.RSUN_REF
			B0i = hdr_i.CRLT_OBS
			L0i = hdr_i.CRLN_OBS

			;Observer Coordinates
			X_scli = hdr_i.CDELT1/60.0
			Y_scli = hdr_i.CDELT2/60.0
			
		endif else begin


			datei = fxpar(hdr_i, 'DATE_OBS')

			;KPVT-512
			if instr eq 1 then begin
			
				;Define center and radius
				hfxi = fxpar(hdr_i, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
				hfyi = fxpar(hdr_i, 'CRPIX2A');+1.0;    Location of the center in y pixels
				dii = fxpar(hdr_i,'EPH_R0');

				;Load Solar Coordinates
				P0i = 0.0
				RDi = !values.f_nan
				B0i = fxpar(hdr_i, 'EPH_B0')
				L0i = fxpar(hdr_i, 'EPH_L0')

				;Observer Coordinates
				X_scli = fxpar(hdr_i, 'CDELT1')*fxpar(hdr_i, 'CRR_SCLX')/60.0
				Y_scli = fxpar(hdr_i, 'CDELT2')*fxpar(hdr_i, 'CRR_SCLY')/60.0

			endif

			;MDI
			if instr eq 3 then begin
			
				;Define center and radius
				hfxi = fxpar(hdr_i, 'X0');  Location of the center in x pixels 
				hfyi = fxpar(hdr_i, 'Y0');  Location of the center in y pixels
				dii = fxpar(hdr_i,'R_SUN');

				;Load Solar Coordinates
				P0i = fxpar(hdr_i, 'P_ANGLE')
				RDi = fxpar(hdr_i, 'OBS_DIST')/0.0046491
				B0i = fxpar(hdr_i, 'B0')
				L0i = fxpar(hdr_i, 'L0')

				;Observer Coordinates
				X_scli = fxpar(hdr_i, 'XSCALE')
				Y_scli = fxpar(hdr_i, 'YSCALE')	
			
			endif

		endelse
		
		Latht = 90-findgen(1,181)
            
        Vsbl = hel2arcmin(Latht, Lonh, vsblN, p = P0i, b0 = B0i, l0 = L0i , rsun = dii*X_scli*60.0)/X_scli
		
		;print, Vsbl
		
        Vsbl_ind = where(vsblN eq 0)
        Vsbl[0,Vsbl_ind] = !values.f_nan
        Vsbl[1,Vsbl_ind] = !values.f_nan
		
		
		Lonht = 180-findgen(1,361)
            
        Vsbl2 = hel2arcmin(Lath, Lonht, vsblN, p = P0i, b0 = B0i, l0 = L0i , rsun = dii*X_scli*60.0)/X_scli
		
		;print, Vsbl
		
        Vsbl_ind = where(vsblN eq 0)
        Vsbl2[0,Vsbl_ind] = !values.f_nan
        Vsbl2[1,Vsbl_ind] = !values.f_nan		
		
		
                
        loadct, 13, /silent
        
        plots, (Vsbl[0,*] + hfxi)*d_zoom + pos[0], (Vsbl[1,*] + hfyi)*d_zoom + pos[1] ,color='00FFFF'x,/device,thick=2,LINESTYLE=2
        plots, (Vsbl2[0,*] + hfxi)*d_zoom + pos[0], (Vsbl2[1,*] + hfyi)*d_zoom + pos[1] ,color='00FFFF'x,/device,thick=2,LINESTYLE=2
      
        loadct, 0, /silent
   
    endif

    Plot, [0,1], /NoErase, /NoData, XStyle=1, YStyle=1, $
            /Device, Position=pos, $
            XRange=xrange, YRange=yrange, $
            Title=title, charsize=2

;    stop

			
	;----------------------------------		

; code to save each frame of the movie
; parameters to add
sz = size(im)
print_zoom = 1.0
print_xsize=sz[1]*print_zoom
print_ysize=sz[2]*print_zoom 


if keyword_set(sqs_nm) then begin
    flnm = '/nfs/hl0/data/mdeluca/work/movie/AR' + strtrim(string(sqs_nm/1000.,'(F5.3)'),2) + '.eps'
endif


if keyword_set(prt) and keyword_set(ARs) then begin

    arin = where(ARs.mdi_i eq mdi_i, n_ars)        

    set_plot,'PS'
    device, filename=flnm,/color, encapsulated = 1, decomposed = 0, bits_per_pixel=8,/portr,/inches,xsize=6.0*print_zoom,ysize=6.0*print_zoom,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    px = !d.x_vsize ;Get size of window in device units
    py = !d.y_vsize

    plot,[1,1],/nodata,xstyle=5,ystyle=5
    ;tv,bytscl(congrid(im_s_gp+im_s_gn,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu) 
    loadct, 0, /silent
    tv,bytscl(congrid(im,print_xsize,print_ysize),min=-max,max=max)
    
    if (n_ars ge 1) then begin
	
		tmpArs = ARs[arin]
	
        for i=0,n_ars-1 do begin
            str='P'+strtrim(string(tmpArs[i].labl,'(I5)'),2)
            
            loadct, 13, /silent
            
            !P.COLOR=long(tmpArs[i].clr)
            if keyword_set(tag) then xyouts,tmpArs[i].fcenxpp*px/sz[1],tmpArs[i].fcenypp*py/sz[2],str,charsize=1.5,charthick = 4,/device
;            xyouts,tmpArs[i].fcenxpp*px/sz[1],tmpArs[i].fcenypp*py/sz[2],str,charsize=1.5,charthick = 2,color = long(tmpArs[i].clr),/device
            tvcircle,2.0*tmpArs[i].dcenpp*px/sz[1], tmpArs[i].fcenxpp*px/sz[1],tmpArs[i].fcenypp*py/sz[2], color = strtrim(string(tmpArs[i].clr,'(I3)'),2), thick=4,/device
            
            str='N'+strtrim(string(tmpArs[i].labl,'(I5)'),2)
            !P.COLOR=long(tmpArs[i].clr)
            if keyword_set(tag) then xyouts,tmpArs[i].fcenxpn*px/sz[1],tmpArs[i].fcenypn*py/sz[2],str,charsize=1.5,charthick = 4,/device
;            xyouts,tmpArs[i].fcenxpn*px/sz[1],tmpArs[i].fcenypn*py/sz[2],str,charsize=1.5,charthick = 2,color= long(tmpArs[i].clr),/device
            tvcircle,2.0*tmpArs[i].dcenpn*px/sz[1], tmpArs[i].fcenxpn*px/sz[1],tmpArs[i].fcenypn*py/sz[2], color = strtrim(string(tmpArs[i].clr,'(I3)'),2), thick=4,/device
    
            
        endfor    
    endif
    
    loadct, 0, /silent
    
    xyouts,0.0*px,0.0*py,title,charsize=1.75, charthick = 6,color= 254,/device
    
    device, /close_file
    set_plot,'X'
    loadct, 0, /silent           
    
    
endif
			

return
end
