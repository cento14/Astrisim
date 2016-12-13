

;;; VelaX: (ra=128.75,dec=-45.6) <-> (l=263.86,b=-3.09)

l=263.86
b=-3.09


;;;  1)  Old mksrc, only for point-like sources. As input it takes 3 vectors describing the spectrum.


lemin = findgen(30.)/10.-1                      ; Lower boundary for energy channels   [ Log(TeV) ]
lemax = lemin+.1				; Higher boundary for energy channels  [ Log(TeV) ]
flux = (10.^(lemin+0.05))^(-1.36) * 11.6e-12    ; Flux                                 [ph / cm sec]

mksrc, 'VelaX_1.fits', l, b, lemin, lemax, flux



;;; 2)  For both point or diffuse sources. The source morphology is passed through fits files.
        Several analytical models are availables for spectra.
        EBL absorption evaluated if z_src is given (from Dominguez+11).


e_min= 0.04    ; minimum limit of the energy range [TeV]
e_max= 160.0   ; maximum limit of the energy range [TeV]


mksrc2, 'VelaX_2.fits', l, b, e_min, e_max, NEbins=30, $
		spec_type='ExpCutoff',k=11.6e-12,index1=-1.36,E0=1.0,Ecut=13.6, $
		DiffMap='VelaX_morphology.fits',exten_no=0, z_src=0.0


		
;;; 3)  Takes as input a Gsed file  ( Gsed files are availables at http://www.iasf-milano.inaf.it/~giuliani/astrisim/gsed )


mksrc3,'HESSJ1418-609.fits', l, b, gsed_file='1418-609.fits', map_ext='Map_TeVCat', sed_ext='MODEL_HESS2006'




end


