PRO zaw_runall_ar

dates = ['1989-09-05','1988-12-28','1989-09-14','1989-06-17','1989-11-14']

;expf_arr=[0,1,2]
;expd_arr=[0,4,8]
;exps_arr=[0,0.5,1]
;ardislim_arr=[5,7,9]
MxFlxim_arr=[2.3,2.4,2.6,2.7]

for i = 0,4 do begin

	date = dates[i]
	print, 'RUNNING FOR ' + date	
	
;	for k = 0,2 do begin
;		for j = 0,2 do begin
			for h = 0,3 do begin
	kpvt_ar_id, kerth1=400,kerth2=200, eros=10,dila=20, start_date = date, exp_f=1, exp_d=4, exp_s=1, ardislim1=5, MxFlxim=MxFlxim_arr[h]
			endfor
;		endfor
;	endfor

endfor 

END
