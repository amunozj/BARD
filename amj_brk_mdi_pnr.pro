
;EDITS
PRO amj_brk_mdi_pnr, pnr_file
 
restore, pnr_file


mdi_imin = min([min(PRs[where(PRs.mdi_i ne 0)].mdi_i),min(NRs[where(NRs.mdi_i ne 0)].mdi_i)])
mdi_imax = max([max(PRs[where(PRs.mdi_i ne 0)].mdi_i),max(NRs[where(NRs.mdi_i ne 0)].mdi_i)])

prsmdi = PRs.mdi_i
nrsmdi = NRs.mdi_i

nopnrs = 0;
yespnrs = 0;

bndrsa = 0;
bndrsq = 0;
whl_sw = 0

actv_sw = 0;

n_qd = 0;
ndsqtlim = 6


for i = mdi_imin, mdi_imax do begin
  
  
  tmpin = where(prsmdi eq i, nprs)
  tmpin = where(nrsmdi eq i, nnrs)
  
  ;if there are regions and we are in quiet period
  if ( (nprs gt 0) and (nnrs gt 0) and (actv_sw eq 0) ) then begin
    bndrsa = [bndrsa, i]
    actv_sw = 1
    n_qd = 0
  endif
  
  ;if we are in active period and hit a day without both positive and negative regions
  if ( (nprs eq 0) or (nnrs eq 0) and (actv_sw eq 1) ) then begin
    
    ;If we run into the first day without regions
    if (n_qd eq 0) then begin
      tmp_i = i
      n_qd = n_qd + 1
    endif else begin
      ;Otherwise check if day without regions is consecutive and increase the counter
      if (i-tmp_i eq 1) then begin
        tmp_i = i
        n_qd = n_qd + 1
          
        ;If a quiet period of certain length is reached then end the active period  
        if (n_qd ge ndsqtlim) then begin
          bndrsq = [bndrsq, i-ndsqtlim+1]
          actv_sw = 0
        endif
      
      ;If at least a day with regions has passer reset the counter  
      endif else begin
        n_qd = 0
      endelse
      
    endelse

  endif  
  
  
  
;  if ( (nprs gt 0) and (nprs gt 0) ) then begin
;    yespnrs = [yespnrs, i]
;  endif
;  if ( (nprs eq 0) or (nprs eq 0) ) then begin
;    nopnrs = [nopnrs, i]
;  endif
  
endfor


bndrsa = bndrsa[1:n_elements(bndrsa)-1]
bndrsa = bndrsa[0:n_elements(bndrsa)-2]

bndrsq = bndrsq[1:n_elements(bndrsq)-1]

print, 'Starting of active phase'
print, bndrsa

print, 'End of active phase'
print, bndrsq

print, 'Length of active phase'
print, bndrsq-bndrsa

print, 'Length of quiet phase'
print, bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]


lthactv = bndrsq-bndrsa
lthaqt = bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]

mdi_icr = mdi_imin  ; variable storing the current day of the mission
max_dys = 350

mrg_sw = 1
dyssep = ndsqtlim
dyssepmx = 10

while (mrg_sw eq 1) do begin
  
  tmp_in = where(lthaqt eq dyssep, nds)
  while (nds gt 0) do begin
      
    if (tmp_in[0] lt n_elements(bndrsa)-2) then begin
      
      if ( (lthactv[tmp_in[0]]+lthaqt[tmp_in[0]]+lthactv[tmp_in[0]+1]) le max_dys) then begin 
        bndrsa = [bndrsa[0:tmp_in[0]], bndrsa[tmp_in[0]+2:(n_elements(bndrsa)-1)]]
        bndrsq = [bndrsq[0:(tmp_in[0]-1)], bndrsq[tmp_in[0]+1:(n_elements(bndrsq)-1)]]
      
        lthactv = bndrsq-bndrsa
        lthaqt = bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]
        
      endif else begin
        lthaqt[tmp_in[0]] = 0      
      endelse
    
    endif else begin
      lthaqt[tmp_in[0]] = 0
    endelse

    tmp_in = where(lthaqt eq dyssep, nds)    
    
  endwhile

  dyssep = dyssep + 1
  if (dyssep gt dyssepmx) then mrg_sw = 0
   
  
endwhile

lthactv = bndrsq-bndrsa
lthaqt = bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]


print, '-------------------------------------------------------------'

print, 'Starting of active phase'
print, bndrsa

print, 'End of active phase'
print, bndrsq

print, 'Length of active phase'
print, bndrsq-bndrsa

print, 'Length of quiet phase'
print, bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]


stop


tmpPRs = PRs
tmpNRs = NRs


tmp_in = where(lthaqt eq 7);ndsqtlim)

;for i = 0, tmp_in[0] do begin

for i = 0, n_elements(bndrsa)-1 do begin
  
  mdi_i1 = bndrsa[i]
  mdi_i2 = bndrsq[i]
  
  start_date = mdi_datestr( string( mdi_i1 ), /inverse )
  end_date = mdi_datestr( string( mdi_i2 ), /inverse )

  print, start_date + ' to ' + end_date

  PRs = tmpPRs[where( (prsmdi ge mdi_i1) and (prsmdi le mdi_i2) )]
  NRs = tmpNRs[where( (nrsmdi ge mdi_i1) and (nrsmdi le mdi_i2) )]

  SAVE, PRs, NRs, seg_const_lg, FILENAME = 'PNRs_' + start_date + '_to_' + end_date + '.sav'  
  
  
endfor

return
END


;----------------------------------------------------------------------------------------------------------

 
