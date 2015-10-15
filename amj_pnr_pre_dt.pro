PRO amj_pnr_pre_dt, start_date=start_date, end_date=end_date, ars_dt = ars_dt, pre_frag = pre_frag


;initializing run and variables
;-------------------------------------------------------------------------------------

;Active region detection constants
ar_cnst={dis_lim1:4.0, dis_lim2:4.0, exp_f: 1.0,  exp_d: 4.0, exp_s: 0.5, mlth: 40.0, mxB: 180.0, MxFlxim:3.0, Imb_tol: 0.10, Imb_it: 10, lim_lon: 0.0, k_sig:15.0, npr: 5, nmgnt: 5 , vld_thr: 0.69, valid_range:[-20000.,20000.]}
seg_const={ker_th:500.0, ker_th2:275.0, ar_th:325.0, ar_th2:50.0, eros_size:9.0, dila_size:18.0, npssu: 4, dis_lim:3.0, ovr_lim: 0.3, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}


ar_cnst_lg = ar_cnst
seg_const_lg = seg_const
;EDITS




mdi_ir  = long( mdi_datestr( start_date ) )
mdi_if  = long( mdi_datestr( end_date ) )
 
;Initializing PR array
PRs = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm:!values.f_nan, qm: !values.f_nan}
;Initializing NR array
NRs = PRs

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
    dater = mdi_datestr( string( mdi_ir ), /inverse ) 
    mgr = amj_mdi_archive( dater, i=1, hdrr )
        
    sr = size(mgr)

    mdi_ir = mdi_ir + 1
    
  endrep until (sr[2] eq 8)  
  mdi_ir = mdi_ir - 1

  print, dater
  ;Performing detection of possitive and negative regions
  amj_coord, mgr.img, hdrr, CRD                       
  amj_pnr_dt, CRD, mdi_ir, PRs, NRs, seg_const=seg_const;, /detdisp
  
  if keyword_set(ars_dt) then begin
    amj_ar_dt_track_dr, CRD, mdi_ir, lbl, PRs, NRs, ARs, ar_cnst=ar_cnst, seg_const=seg_const, /bck_trck;, /display
    detar_sw = 0
    mdi_ir_vis = [mdi_ir_vis,mdi_ir]    
  endif
  
  if keyword_set(pre_frag) then begin

        ;Creating temporary seg_const
        tmp_seg_const = seg_const
        tmp_seg_const.npssu   = 2
        tmp_seg_const.ker_th  = 200
        tmp_seg_const.ker_th2 = 150
        tmp_seg_const.dis_lim = 2
        tmp_seg_const.ovr_lim = 0.7
        ;tmp_seg_const.dila_size=5
        
        amj_pnr_dt, CRD, mdi_ir, PRsFr, NRsFr, seg_const=seg_const            
    
  endif  
  
  
  mdi_ir = mdi_ir + 1

  if keyword_set(ars_dt) then begin
  SAVE, PRs, NRs, ARs, detar_sw, mdi_ir_vis, lbl, seg_const_lg, ar_cnst_lg, FILENAME = 'PNRs_ARs_' + start_date + '_to_' + end_date + '.sav'  

  endif else begin
  SAVE, PRs, NRs, seg_const_lg, FILENAME = 'PNRs_' + start_date + '_to_' + end_date + '.sav'  
  endelse
  
  
endwhile


return
END


;----------------------------------------------------------------------------------------------------------

 
