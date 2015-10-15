PRO zaw_runall

dates = ['1984-02-09','1984-02-25','1984-04-29','1984-05-13','1984-05-24','1984-11-11','1985-06-25','1985-10-23','1985-12-14','1986-07-10','1986-02-23','1988-07-17','1988-12-28','1989-03-13','1989-06-09','1989-06-17','1989-09-05','1989-09-13','1989-11-14']

;kerth1_arr = [200.0,400.0,600.0]
;kerth2_arr=[350,360,370]
;eros_arr = [5,10,15,20]
;dila_arr = [10,20,30,40]

for i = 0,18 do begin

	date = dates[i]
	print, 'RUNNING FOR ' + date	
	
;	for k = 0,2 do begin
;		for j = 0,3 do begin
;			for h = 0,3 do begin
;	kpvt_ar_id, kerth1=kerth1_arr, eros=eros_arr[j],dila=dila_arr[h], start_date = date
	kpvt_ar_id, kerth1=400,kerth2=200, eros=10,dila=20, start_date = date
;			endfor
;		endfor
;	endfor

endfor 

END
