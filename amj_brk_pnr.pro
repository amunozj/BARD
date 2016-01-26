
;EDITS
PRO amj_brk_pnr, pnr_file, instr, max_intr = max_intr

SET_PLOT, 'Z'
thisDevice = !D.Name
 
restore, pnr_file


mdi_imin = min([min(PRs[where(PRs.mdi_i ne 0)].mdi_i),min(NRs[where(NRs.mdi_i ne 0)].mdi_i)])
mdi_imax = max([max(PRs[where(PRs.mdi_i ne 0)].mdi_i),max(NRs[where(NRs.mdi_i ne 0)].mdi_i)])

prsmdi = PRs.mdi_i
nrsmdi = NRs.mdi_i

;nopnrs = 0;
;yespnrs = 0;

bndrsa = 0    ;Boundaries of active periods
bndrsq = 0    ;Boundaries of quiet periods

actv_sw = 0;  ;Switch that keeps track of whether the current period is active or quiet

n_qd = 0;     ;Number of days of quiet period
ndsqtlim = 7  ;Minimum number of quiet days necessary to end a quiet period


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

;If we reach the end of the interval without finding a quiet period, end the current active period
if (actv_sw eq 1) then begin
  bndrsq = [bndrsq, i]
  actv_sw = 0
endif

;stop

bndrsa = bndrsa[1:n_elements(bndrsa)-1]
;bndrsa = bndrsa[0:n_elements(bndrsa)-2]

bndrsq = bndrsq[1:n_elements(bndrsq)-1]

print, 'Starting of active phase'
print, bndrsa

print, 'End of active phase'
print, bndrsq

print, 'Length of active phase'
print, bndrsq-bndrsa


if n_elements(bndrsa) gt 1 then begin
	print, 'Length of quiet phase'
	print, bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]


	lthactv = bndrsq-bndrsa
	lthaqt = bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]

	mdi_icr = mdi_imin  ; variable storing the current day of the mission
	max_dys = 5000       ; Maximum length of each active period

	mrg_sw = 1          ;Switch that determines when the process of mergin is over
	dyssep = ndsqtlim   ;Minimum number of quiet days necessary to end a quiet period
	dyssepmx = 10

	while (mrg_sw eq 1) do begin
	  
	  tmp_in = where(lthaqt eq dyssep, nds)
	  while (nds gt 0) do begin
		  
		if (tmp_in[0] lt n_elements(bndrsa)-2) then begin
		  
		  if ( (lthactv[tmp_in[0]]+lthaqt[tmp_in[0]]+lthactv[tmp_in[0]+1]) le max_dys) then begin 
			bndrsa = [bndrsa[0:tmp_in[0]], bndrsa[tmp_in[0]+2:(n_elements(bndrsa)-1)]]
			
			if tmp_in[0] eq 0 then begin
				bndrsq = [bndrsq[tmp_in[0]+1:(n_elements(bndrsq)-1)]]
			endif else begin  
				bndrsq = [bndrsq[0:(tmp_in[0]-1)], bndrsq[tmp_in[0]+1:(n_elements(bndrsq)-1)]]
			endelse
		  
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

endif


print, '-------------------------------------------------------------'

print, 'Starting of active phase'
print, bndrsa

print, 'End of active phase'
print, bndrsq

print, 'Length of active phase'
print, bndrsq-bndrsa

if n_elements(bndrsa) gt 1 then begin
	print, 'Length of quiet phase'
	print, bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]
endif



if keyword_set(max_intr) then begin

	for i = 0, n_elements(bndrsa)-1 do begin
	
		if (bndrsq[i]-bndrsa[i]) gt 2*max_intr then begin
			
			bndrsatmp = indgen(uint((bndrsq[i]-bndrsa[i])/max_intr))*max_intr+bndrsa[i]
			bndrsqtmp = (indgen(uint((bndrsq[i]-bndrsa[i])/max_intr)-1)+1)*max_intr+bndrsa[i]
			
			if n_elements(bndrsa) gt 1 then begin
			
				;Adding additional intervals
				first = 0
				last = n_elements(bndrsa)-1
				
				CASE i OF
					first: begin
						bndrsa = [bndrsatmp, bndrsa[1:*]]
						bndrsq = [bndrsqtmp, bndrsq]
					end
				  
					last: begin  
						bndrsa = [bndrsa[first:last-1], bndrsatmp]
						bndrsq = [bndrsq[first:last-1], bndrsqtmp, bndrsq[last]]
					end	
						
					ELSE: begin 
						bndrsa = [bndrsa[first:i-1], bndrsatmp, bndrsa[i+1:last] ]
						bndrsq = [bndrsq[first:i-1], bndrsqtmp, bndrsq[i:last] ]
					end	
						
				ENDCASE
				
			endif else begin

				bndrsa = bndrsatmp
				bndrsq = [bndrsqtmp, bndrsq]
			
			endelse
		
		endif
	
	
	endfor

endif

print, '-------------------------------------------------------------'

print, 'Starting of active phase'
print, bndrsa

print, 'End of active phase'
print, bndrsq

print, 'Length of active phase'
print, bndrsq-bndrsa

if n_elements(bndrsa) gt 1 then begin
	print, 'Length of quiet phase'
	print, bndrsa[1:n_elements(bndrsa)-1] - bndrsq[0:n_elements(bndrsq)-2]
endif

;stop

tmpPRs = PRs
tmpNRs = NRs

szPfr = size(PRsFr)
szNfr = size(NRsFr)

if (szPfr[0] gt 0) or (szNfr[0] gt 0) then begin
	tmpPRsFr = PRsFr
	tmpNRsFr = NRsFr
endif

;First reference day for keeping track time easily
if instr eq 1 then DayOff = julday(1,1,1970);	KPVT 512
if instr eq 2 then DayOff = julday(1,1,1970);	KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993);	MDI
if instr eq 4 then DayOff = julday(1,1,2009);	HMI


for i = 0, n_elements(bndrsa)-1 do begin
  
  mdi_i1 = bndrsa[i]
  mdi_i2 = bndrsq[i]
  
  caldat, mdi_i1+DayOff, Month, Day, Year
  start_date = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)

  caldat, mdi_i2+DayOff, Month, Day, Year
  end_date = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)

  print, start_date + ' to ' + end_date

  PRs = tmpPRs[where( (prsmdi ge mdi_i1) and (prsmdi le mdi_i2) )]
  NRs = tmpNRs[where( (nrsmdi ge mdi_i1) and (nrsmdi le mdi_i2) )]

  if (szPfr[0] gt 0) or (szNfr[0] gt 0) then begin

	PRsFr = tmpPRsFr[where( ((tmpPRsFr.mdi_i ge mdi_i1) and (tmpPRsFr.mdi_i le mdi_i2) ) or tmpPRsFr.mdi_i eq 0)]
	NRsFr = tmpNRsFr[where( ((tmpNRsFr.mdi_i ge mdi_i1) and (tmpNRsFr.mdi_i le mdi_i2) ) or tmpNRsFr.mdi_i eq 0)]
  
	SAVE, PRs, NRs, PRsFr, NRsFr, seg_const_lg, FILENAME = 'PNRs_PrFr_' + start_date + '_to_' + end_date + '.sav' 
  
  endif else begin
  
	SAVE, PRs, NRs, seg_const_lg, FILENAME = 'PNRs_' + start_date + '_to_' + end_date + '.sav'  
  
  endelse
  
endfor

return
END


;----------------------------------------------------------------------------------------------------------

 
