PRO kpvt_pnr_pre_dt, start_date=start_date, end_date=end_date, ars_dt = ars_dt, pre_frag = pre_frag

SET_PLOT, 'Z'
thisDevice = !D.Name

;initializing run and variables
;-------------------------------------------------------------------------------------

;Active region detection constants

;MDI
;ar_cnst={dis_lim1:4.0, dis_lim2:4.0, exp_f: 1.0,  exp_d: 4.0, exp_s: 0.5, mlth: 40.0, mxB: 180.0, MxFlxim:3.0, Imb_tol: 0.10, Imb_it: 10, lim_lon: 0.0, k_sig:15.0, npr: 5, nmgnt: 5 , vld_thr: 0.69, valid_range:[-20000.,20000.]}
;seg_const={ker_th:500.0, ker_th2:275.0, ar_th:325.0, ar_th2:50.0, eros_size:9.0, dila_size:18.0, npssu: 4, dis_lim:3.0, ovr_lim: 0.3, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}

;KPVT
ar_cnst={dis_lim1:4.0, dis_lim2:4.0, exp_f: 1.0,  exp_d: 4.0, exp_s: 1.0, mlth: 40.0, mxB: 180.0, MxFlxim:3.0, Imb_tol: 0.10, Imb_it: 10, lim_lon: -90.0, k_sig:15.0, npr: 5, nmgnt: 5 , vld_thr: 0.69, valid_range:[-20000.,20000.]}
seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}



ar_cnst_lg = ar_cnst
seg_const_lg = seg_const
;EDITS

DayOff = julday(1,1,1970);  First reference day for keeping track time easily


mdi_ir  = julday(strmid(start_date,5,2),strmid(start_date,8,2),strmid(start_date,0,4))-DayOff
mdi_if  = julday(strmid(end_date,5,2),strmid(end_date,8,2),strmid(end_date,0,4))-DayOff

 
;Initializing PR array
PRs = {rgn,  lnk_sw: 0,  mdi_i: 0L, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm:!values.f_nan, qm: !values.f_nan}

;Initializing NR array
NRs = PRs


if keyword_set(pre_frag) then begin
   
    ;Initializing PR frag array
    PRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm:!values.f_nan, qm: !values.f_nan}
    ;Initializing NR frag array
    NRsFr = PRsFr
  
endif  


if keyword_set(ars_dt) then begin
  
    mdi_ir_vis = 0 ; Variable that keeps track of the maximum date where AR detection has taken place
    lbl = 0;  AR label
      
    ;Initializing AR array
    ARs = {ar, mdi_i: 0, date: '', labl: 0, clr: 0, indxp: '', indxn: '', fluxp: !values.f_nan, fluxn: !values.f_nan, areap:!values.f_nan, arean:!values.f_nan, $
              fcn_ltp: !values.f_nan, fcn_lnp: !values.f_nan, dcenp: !values.f_nan, $ 
              fcn_ltn: !values.f_nan, fcn_lnn: !values.f_nan, dcenn: !values.f_nan, $
              fcenxpp: !values.f_nan, fcenypp: !values.f_nan, dcenpp: !values.f_nan, $
              fcenxpn: !values.f_nan, fcenypn: !values.f_nan, dcenpn: !values.f_nan, $               
              dis: !values.f_nan, tilt: !values.f_nan, lp: !values.f_nan, dm:!values.f_nan, qm: !values.f_nan}   
  
endif
    

while (mdi_ir le mdi_if) do begin

  repeat begin
    ;Reading files
    caldat, mdi_ir+DayOff, Month, Day, Year
    dater = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
    mgr = amj_kpvt_archive( dater, hdrr )    
            
    sr = size(mgr)

    mdi_ir = mdi_ir + 1
    print, mdi_ir
    
  endrep until ( (sr[2] eq 8) or (mdi_ir gt mdi_if) )  
  mdi_ir = mdi_ir - 1


  if (sr[2] eq 8) then begin
      print, dater
      ;Performing detection of possitive and negative regions
      amj_kpvt_coord, mgr.img, hdrr, CRD, seg_const=seg_const                       
      amj_pnr_dt, CRD, mdi_ir, PRs, NRs, seg_const=seg_const;, /detdisp
      
      if keyword_set(ars_dt) then begin
        amj_ar_dt_track_dr, CRD, mdi_ir, lbl, PRs, NRs, ARs, ar_cnst=ar_cnst, seg_const=seg_const, /bck_trck;, /display
        detar_sw = 0
        mdi_ir_vis = [mdi_ir_vis,mdi_ir]    
      endif
      
      sz_fr = 30  ; Minimum size for pre-fragmentation
      if keyword_set(pre_frag) then begin
        
          ;Creating temporary seg_const
          tmp_seg_const = seg_const
          tmp_seg_const.npssu = 2
          tmp_seg_const.ker_th=200
          tmp_seg_const.ker_th2=150
          tmp_seg_const.dis_lim = 2
          tmp_seg_const.ovr_lim = 0.7
          ;tmp_seg_const.dila_size=5
          
          
          ; Positive regions------------------------------------
          frin = where((PRs.mdi_i eq mdi_ir) and (PRs.dcenp gt sz_fr), n_fr)      
          if n_fr gt 0 then begin
    
              for i=0,n_fr-1 do begin
    
                  print, PRs[frin[i]].fr_lbl
                  
                  ;Extracting necessary pixels
                  ex_in = long(strsplit(PRs[frin[i]].indx,/extract))
                  
                  ;Creating temporary image
                  tmp_im = mgr.img*0.0
                  tmp_im[ex_in] = mgr.img[ex_in]
                  
                  ;Creating temporary PNR arrays
                  tmp_PRs = {rgn,  lnk_sw: 0,  mdi_i: 0L, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}    
                  tmp_NRs = tmp_PRs
                  
                  amj_kpvt_coord, tmp_im, hdrr, CRD              
                  amj_pnr_dt, CRD, mdi_ir, tmp_PRs, tmp_NRs, seg_const=tmp_seg_const, pnr_lbl = PRs[frin[i]].fr_lbl;,/not_merge;, /detdisp
                  
                  if n_elements(tmp_PRs) gt 1 then begin              
                      PRsFr = [PRsFr, tmp_PRs[1:n_elements(tmp_PRs)-1]]                  
                  endif
      
              endfor
      
          endif
            
          ; Negative regions------------------------------------
          frin = where((NRs.mdi_i eq mdi_ir) and (NRs.dcenp gt sz_fr), n_fr)      
          if n_fr gt 0 then begin
    
              for i=0,n_fr-1 do begin
    
                  print, NRs[frin[i]].fr_lbl
                  
                  ;Extracting necessary pixels
                  ex_in = long(strsplit(NRs[frin[i]].indx,/extract))
                  
                  ;Creating temporary image
                  tmp_im = mgr.img*0.0
                  tmp_im[ex_in] = mgr.img[ex_in]
                  
                  ;Creating temporary PNR arrays
                  tmp_PRs = {rgn,  lnk_sw: 0,  mdi_i: 0L, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm: !values.f_nan, qm: !values.f_nan}    
                  tmp_NRs = tmp_PRs
                  
                  amj_kpvt_coord, tmp_im, hdrr, CRD              
                  amj_pnr_dt, CRD, mdi_ir, tmp_PRs, tmp_NRs, seg_const=tmp_seg_const, pnr_lbl = NRs[frin[i]].fr_lbl;,/not_merge;, /detdisp
                  
                  if n_elements(tmp_NRs) gt 1 then begin              
                      NRsFr = [NRsFr, tmp_NRs[1:n_elements(tmp_NRs)-1]]                  
                  endif
      
              endfor
      
          endif
        
      endif  
                  
  endif
  
  mdi_ir = mdi_ir + 1      

  if keyword_set(ars_dt) then begin
  SAVE, PRs, NRs, ARs, detar_sw, mdi_ir_vis, lbl, seg_const_lg, ar_cnst_lg, FILENAME = 'PNRs_ARs_' + start_date + '_to_' + end_date + '.sav'  

  endif else begin
    
      if keyword_set(pre_frag) then begin
        
          SAVE, PRs, NRs, PRsFr, NRsFr, seg_const_lg, FILENAME = 'PNRs_PrFr' + start_date + '_to_' + end_date + '.sav'  
        
      endif else begin
    
          SAVE, PRs, NRs, seg_const_lg, FILENAME = 'PNRs_' + start_date + '_to_' + end_date + '.sav'  
  
      endelse
  
  endelse
  
  
endwhile


return
END


;----------------------------------------------------------------------------------------------------------

 
