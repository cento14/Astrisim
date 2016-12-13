pro  write_ctalike,l_pnt,b_pnt,l,b,t,ene,tsec,outfile=outfile

infile2='$ASTRISIM/template.fits'

na=n_elements(l) *1L

m=mrdfits(infile2,0,h0,/silen)
m1=mrdfits(infile2,1,h1,/silen)
m2=mrdfits(infile2,2,h2,/silen)

fxwrite,outfile,h0       ; creo il file fits con l'header primario 
fxaddpar,h1,'NAXIS2',na,'number of rows in table'

fxaddpar,h1,'ONTIME',tsec,'[s] Total good time including deadtime'
fxaddpar,h1,'LIVETIME',tsec,'[s] Total livetime'

;;; Adjust outfile for current observations
euler,l_pnt,b_pnt,ra_pnt,dec_pnt,2						;ab; Converting (l,b) to (ra,dec)
fxaddpar,h1,'DSVAL2','CIRCLE('+string(ra_pnt)+','+string(dec_pnt)+',10)'	;ab; Data selection value: coordinates
fxaddpar,h1,'DSVAL3',string(min(ene))+':'+string(max(ene))           		;ab; Data selection value: energy
fxaddpar,h1,'RA_PNT',ra_pnt							;ab; [deg] Pointing Right Ascension
fxaddpar,h1,'DEC_PNT',dec_pnt							;ab; [deg] Pointing Declination

 
fxbcreate,n_file,outfile,h1  ; creo nel file fits l'header della 1a tabella 

euler,l,b,ra,dec,2
 

  for i=0,na-1 do begin
  fxbwrite,n_file,i+1L,1,i+1
  fxbwrite,n_file,t(i)+0.d0,3,i+1
  fxbwrite,n_file,float(ra(i)),7,i+1
  fxbwrite,n_file,float(dec(i)),8,i+1
  fxbwrite,n_file,float(ene(i)),21,i+1

  endfor


fxbcreate,unit,outfile,h2


fxbwrite,unit,0d0,1,1
fxbwrite,unit,double(tsec),2,1


fxbfinish,unit
close,/all


end
