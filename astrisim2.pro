
pro astrisim2

;;;;;;;;;  Run parameters

@astrisim.par

tsec=t_hours*3600.

ins_pile=ins_path+ins_file

if size(src_list,/n_dim) eq 0 then readcol,src_path+src_list,src_list,for='a'

src_pist=src_path+src_list

out_pile=out_path+out_file

print
print,'  ______________   ASTRI SCIENTIFIC SIMULATOR (v2.0) _____________ '
print
print,'Instrument Response File : ',ins_pile
print,'Simulated Sources        : ',transpose(src_pist)
print
print,'Livetime (hours)         : ',t_hours
print,'Pointing Direction (l,b) : ',l_asse,b_asse
print



;;;;;;;;;; Instrument response

readcol,ins_path+ins_file,lene_c,aeff_c,r68_c,r80_c,eres_c,bgrate_c,skip=1,/sile


;;; HardCoded

fov=3.2
bgminene=0.8 

;;;;;;;  Background

nch_p=n_elements(bgrate_c)
cts_p=fltarr(nch_p)
lelow =lene_c(0:-2)
lehigh=lene_c(1:*)

lene_p=0

for i=0,nch_p-2 do begin
  
  if 10.^lehigh(i) lt bgminene then bgrate_c(i)=0       ; Soglia 

  bgrate_i=sqrt(bgrate_c(i)*(fov/r80_c(i))^2. *  bgrate_c(i+1)*(fov/r80_c(i+1))^2.) ;  prot / sec TeV (on the whole FOV)

  if bgrate_i gt 0 then cts_p(i)=randomu(see,poiss=bgrate_i * tsec)

  if cts_p(i) gt 0 then lene_p=[lene_p,randomu(see2,cts_p(i))*(lehigh(i)-lelow(i))+lelow(i)]

endfor

ene_p=10^lene_p(1:*)
np= n_elements(ene_p)

theta_p=sqrt(abs(randomu(see4,np,/normal)*fov^2.))
phi_p=randomu(see3,np)*360.

l_p= ( l_asse + theta_p*cos(phi_p/!radeg) / cos(b_asse/!radeg) ) 
b_p=   b_asse + theta_p*sin(phi_p/!radeg)

print,'Total Background events          :',n_elements(ene_p), '    at ['+string(min(ene_p))+';'+string(max(ene_p))+'] TeV'




;;;;;;;;;  Sources

lene=0	; log(energy) of the event
l=0	; l of the event
b=0	; b of the event


for is=0,n_elements(src_pist)-1 do begin
	srcfile=src_pist[is]
	print, '>>> Calculating: ',srcfile
	
    	ebounds=mrdfits(srcfile,2,/fsc,/silen)

    	lemin=ebounds.l_emin
    	lemax=ebounds.l_emax

    	nch=n_elements(lemin)
    	cts=fltarr(nch)

	;;; Parameters of the fits-file
	
	tot_flux=mrdfits(srcfile,1,h1,/silen,/fsc) ; input total array of the flux
	f_size = size(tot_flux)

	if ((fxpar(h1,'CRPIX1') eq 0) && (fxpar(h1,'CRPIX2') eq 0)) then begin
		delt_fix=-1.0
	endif else delt_fix=0.0

	;;; Two for-loops to go through all pixels of the map

	for ind0=0,f_size[1]-1 do begin
    		for ind1=0,f_size[2]-1 do begin			
			w=where(tot_flux[ind0,ind1,*] ne 0,w_num)
			;print,w
			if (w_num ne 0) then begin
				xyad, h1, ind0, ind1, l_src, b_src, /GALACTIC	
				l_src=l_src+delt_fix*fxpar(h1,'CDELT1')
				b_src=b_src+delt_fix*fxpar(h1,'CDELT2')
				flux=tot_flux[ind0,ind1,*] 	; flux of non-zero (ind0,ind1)-pixel

				gcirc,1,l_asse/15.,b_asse,l_src/15.,b_src,dsec	;  Off-axis angle
				theta_src=dsec/3600.

				;;;;;;;;;;;; Photons Generation

				for ie=nch-1,0,-1 do begin
				 ene_ie=10.^(0.5*(lemax(ie)+lemin(ie)))

				 ;aeff_ie=interpol(aeff_c, lene_c, alog10(ene_ie)) * interpol([1,1,0.5,0,0],[0,8,9,10,20.]/2.,theta_src)
				 aeff_ie=interpol(aeff_c, lene_c, alog10(ene_ie)) * exp(-((theta_src/fov)^4.)/2.)

				 cts_mu=flux(ie)*tsec  * aeff_ie*1d4
				 ;print,flux(ie),10.^lemax(ie),10.^lemin(ie) , aeff_ie*1d4
		    
				 if cts_mu gt 0 then cts(ie)=randomu(see,/double,poiss=cts_mu)    
				 if cts(ie) gt 0 then begin   
				  lene=[lene, randomu(see2,cts(ie)) * (lemax(ie)-lemin(ie)) + lemin(ie)]
				  l= [l, fltarr(cts(ie))+l_src]
				  b= [b, fltarr(cts(ie))+b_src] 
				 endif
				;print, cts(ie)
				endfor
			endif
		endfor
	endfor
endfor

ene=10.^lene	; energy of the event

print,'Total Sources photons :',n_elements(ene)-1, '    at ['+string(min(ene))+';'+string(max(ene))+'] TeV'


;;;   Instrument Signature


;gcirc,1,l_asse/15.,b_asse,l/15.,b,dsec 
;theta=dsec/3600.

;ciccia=psf(l,b,lr,br,ins_file=ins_pile,ene=ene,theta=theta)
ener=sfoca(l,b,lr,br,ins_file=ins_pile,ene=ene,theta=theta)


;;;; Merge

ltot=[lr,l_p]
btot=[br,b_p]
enetot=[ener,ene_p]

t=randomu(ss,n_elements(ltot),/double)*tsec
t(0:n_elements(ene)-1)=long(t(0:n_elements(ene)-1))+0.0       ;  Photon signature

enetot=enetot(1:*)
ltot=ltot(1:*) 
btot=btot(1:*)


; S/N 
;
;
;for is=0,n_elements(src_pist)-1 do begin
;for is=0,0 do begin
;
;    srcfile=src_pist[is]
;    !p.title=src_list(is)+'   '+strtrim(fix(tsec/3600.),1)+' hours'
;
;    l_src=fxpar(h1,'CRVAL1')
;    b_src=fxpar(h1,'CRVAL2')
;
;    gcirc,1,l_src/15.,b_src,l_p/15.,b_p,dsec_p           ;  Background around the source
;    gcirc,1,l_src/15.,b_src,ltot/15.,btot,dsectot        ;  Events around the source
;
;    sig=sig(dsectot,enetot,dsec_p,ene_p,tsec,ins_file=ins_pile)
;
;endfor



;;;; Output

plot,l_p,b_p,psym=3,xrange=[10,-10]+l_asse,yrange=[-10,10]+b_asse
oplot,lr,br,psym=3,col=255

print
print,' Saving files...   ('+out_pile+')'

write_ctalike,l_asse,b_asse, ltot mod  360.,btot,t,enetot,tsec,outfile=out_pile+'.evt'
spawn,'cp astrisim.par '+out_pile+'.par'


; aitoff_grid,  /new , lines=2,col=255
; aitoff, l_p,b_p , xait, yait
; plots,xait,yait,psym=3 
; aitoff, l,b , xait, yait
; plots,xait,yait,psym=3 ,col=verde


;sig2,out_pile+'.evt',ins_path+ins_file,l_src,b_src,tsec

print,' Done!'

end
