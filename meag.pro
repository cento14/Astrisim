function meag,v

; Da la media geometrica di v

n=n_elements(v)
m=1.0d0

for i=0,n-1 do begin
m=m*(v(i))^(1.d0/n)
endfor

return,m

end