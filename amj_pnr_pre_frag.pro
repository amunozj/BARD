PRO amj_pnr_pre_frag, pnr_file, O_fldr = O_fldr

SET_PLOT, 'Z'
thisDevice = !D.Name


;initializing run and variables
;-------------------------------------------------------------------------------------

;Region detection constants

;MDI
;seg_const={ker_th:500.0, ker_th2:275.0, ar_th:325.0, ar_th2:50.0, eros_size:9.0, dila_size:18.0, npssu: 4, dis_lim:3.0, ovr_lim: 0.3, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}

;KPVT
seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}


restore, pnr_file

mdi_ir  = min([min(PRs[where(PRs.mdi_i ne 0)].mdi_i),min(NRs[where(NRs.mdi_i ne 0)].mdi_i)])
mdi_if  = max([max(PRs[where(PRs.mdi_i ne 0)].mdi_i),max(NRs[where(NRs.mdi_i ne 0)].mdi_i)])

DayOff = julday(1,1,1970);  First reference day for keeping track time easily

caldat, mdi_ir+DayOff, Month, Day, Year
start_date = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)

caldat, mdi_if+DayOff, Month, Day, Year
end_date = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)


;Initializing PR frag array
PRsFr = {rgn,  lnk_sw: 0,  mdi_i: 0, ar_lbl: 0, fr_lbl: 0L, date: '', indx:'', flux:!values.f_nan, area:!values.f_nan, fcn_lt:!values.f_nan, fcn_ln:!values.f_nan, dcen:!values.f_nan, fcenxp: !values.f_nan, fcenyp: !values.f_nan, dcenp:!values.f_nan, dm:!values.f_nan, qm: !values.f_nan}
;Initializing NR frag array
NRsFr = PRsFr

    

while (mdi_ir le mdi_if) do begin

  repeat begin
    ;Reading files
    caldat, mdi_ir+DayOff, Month, Day, Year
    dater = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
    print, dater
    mgr = amj_kpvt_archive( dater, hdrr )    
            
    sr = size(mgr)

    mdi_ir = mdi_ir + 1
;    print, mdi_ir
    
  endrep until ( (sr[2] eq 8) or (mdi_ir gt mdi_if) )  
  mdi_ir = mdi_ir - 1


  if (sr[2] eq 8) then begin
      
      sz_fr = 30  ; Minimum size for pre-fragmentation
        
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

;              print, PRs[frin[i]].fr_lbl
              
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

;              print, NRs[frin[i]].fr_lbl
              
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
  
  mdi_ir = mdi_ir + 1      
  
  if keyword_set(O_fldr) then begin
      fl_nm = O_fldr + 'PNRs_PrFr' + start_date + '_to_' + end_date + '.sav'
  endif else begin
      fl_nm = 'PNRs_PrFr' + start_date + '_to_' + end_date + '.sav'    
  endelse


  SAVE, PRs, NRs, PRsFr, NRsFr, seg_const_lg, FILENAME = fl_nm    
  
endwhile


return
END


;----------------------------------------------------------------------------------------------------------
