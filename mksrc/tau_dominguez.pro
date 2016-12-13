function tau_dominguez,e,z

; z
red=[0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 0.10157895, 0.11684211, 0.13210526,  0.14736842,    0.16263158, 0.17789474,  0.19315789,  0.20842105,  0.22368421, 0.23894737, 0.25421053,  0.26947368,  0.28473684, 0.3 ,  0.35,  0.4, 0.45,  0.5,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,  0.9 ,   0.95,  1.,1.2,1.4,1.6,1.8,2.]
dimred=n_elements(red)

;energy is read in GeV
readcol,'energy_dominguez.dat',energy,/silent
dimene=n_elements(energy)

etev=e*1.d-3

openr,1,'tau_dominguez11.dat'
matrixa=dblarr(dimred,dimene)
readf,1,matrixa
close,1

yntmp=dblarr(dimene)
ymtmp=dblarr(dimred)

for jj=0, dimred-1 do begin
    for kk=0, dimene-1 do begin
        yntmp(kk)=matrixa(jj,kk)
    endfor
    ymtmp(jj)=interpol(yntmp,energy,etev,/spline)
endfor
y=interpol(ymtmp,red,z,/spline)

pippo=[y,0]

ans=max(pippo)
return,ans
end
