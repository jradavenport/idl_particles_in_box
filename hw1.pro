pro hw1, N, cs, v_init, display=display
; James Davenport
; Astronomy 507
; Autumn 2009

;------------------------------------------------------
; 2D MONATOMIC IDEAL GAS IN A BOX
; All gas particles have the same mass, cross section, & KE_initial
; The box has size 1x1
;

  
if not keyword_set(cs) then $
   cs = 0.01                 ; the radial cross section of the gas particles, in box units
A = FINDGEN(17) * (!PI*2/16.)  
USERSYM, COS(A), SIN(A), /FILL 
!p.multi=[0,1,1] 

mass = fltarr(n) * 0d0 + 1d0


; The box is bound at x=[0,1] and y=[0,1]
; the N particles have random initial positions in X and Y
x_init = randomu(s,n)
y_init = randomu(s,n)
; the N particles have the same velocity (i.e. same KE)
;    but in random directions
;                              *recall* KE = 1/2 m*v^2

if not keyword_set(v_init) then $
   v_init = 0.5d-2
;; thought problem: what happens if random initial velocities are given?
;; v_init = randomu(s,n)*0.02 
v_theta = randomu(s,n)*2d0*!dpi-!dpi ; uniformly random angle of velocity
vx_init = v_init * cos(v_theta)
vy_init = v_init * sin(v_theta)


; now let the system evolve with time
ntime = 2000
x = dblarr(n,ntime+1)
y = dblarr(n,ntime+1)
vx = dblarr(n,ntime+1)
vy = dblarr(n,ntime+1)
;fill these arrays with the initial values
x[*,0] = x_init
y[*,0] = y_init
vx[*,0] = vx_init
vy[*,0] = vy_init

;make array for pressure @ each timestep to fill
Pres = dblarr(ntime+1) * 0d0

FOR j=0L,ntime-1L DO BEGIN
   if keyword_set(display) then begin
      ;; if stddev(sqrt(vx[*,j]^2+vy[*,j]^2)) ne 0 then begin
         plothist,sqrt(vx[*,j]^2+vy[*,j]^2),/autobin,$
                  position=[.6,.1,.9,.9],xtitle='velocity',/xstyle,/fill,$
                  ytickname=[replicate(' ',7)],xrange=[0,0.025]
         plot,x[*,j],y[*,j],psym=4,xrange=[0,1],yrange=[0,1],/xstyle,/ystyle $
              ,title=string(j),position=[.1,.1,.6,.9],/noerase
      ;; endif
      plot,x[*,j],y[*,j],psym=4,xrange=[0,1],yrange=[0,1],/xstyle,/ystyle,title=string(j),position=[.1,.1,.6,.9],/noerase
      wait,0.025                 ; so i can watch the movie
   endif

;; thought problem: at some time excite the hell out of a single
;; particle (add a lot of velocity/energy)
;;   if j eq 20 then vx[0,j] = 1
 
; now deal with collisions
;Particle - Wall: Perfectly Reflective
   ;X=0
   IF TOTAL(x[*,j] LE CS) GE 1 THEN BEGIN
      col = where(x[*,j] le CS) ; the particles colliding
      FOR k=0L,TOTAL(x[*,j] LE CS)-1 DO BEGIN
         if vx[col[k],j] le 0 then $
            vx[col[k],j] = vx[col[k],j]*(-1.)
      ENDFOR
   ENDIF
   ;X=1
   IF TOTAL(x[*,j] GE 1.-CS) GE 1 THEN BEGIN
      col = where(x[*,j] ge 1.-CS)
      FOR k=0L,TOTAL(x[*,j] ge 1.-CS)-1 DO BEGIN
         if vx[col[k],j] ge 0 then $
            vx[col[k],j] = vx[col[k],j]*(-1.)
      ENDFOR
   ENDIF
   ;Y=0
   IF TOTAL(y[*,j] LE CS) GE 1 THEN BEGIN
      col = where(y[*,j] le CS) ; the particles colliding
      FOR k=0L,TOTAL(y[*,j] LE CS)-1 DO BEGIN
         if vy[col[k],j] le 0 then $
            vy[col[k],j] = vy[col[k],j]*(-1.)
      ENDFOR
   ENDIF
   ;Y=1
   IF TOTAL(y[*,j] GE 1.-CS) GE 1 THEN BEGIN
      col = where(y[*,j] ge 1.-CS)
      FOR k=0L,TOTAL(y[*,j] ge 1.-CS)-1 DO BEGIN
         if vy[col[k],j] ge 0 then begin
            vy[col[k],j] = vy[col[k],j]*(-1.)
            ;measure pressure on this wall
            Pres[j] = Pres[j] - (mass[col[k]]*vy[col[k],j])
         endif
      ENDFOR
   ENDIF

;Particle - Particle: Elastic
   ;run thru each particle once to see if its colliding w/ any other
   FOR k=0L,N-1L DO BEGIN
      ;find if any particles are colliding with particle k
      if total(sqrt((x[k,j]-x[*,j])^2.+(y[k,j]-y[*,j])^2.) le CS) $
         gt 1 then begin
         ;find which paricle(s) its colliding with
         col = where(sqrt((x[k,j]-x[*,j])^2.+(y[k,j]-y[*,j])^2.) le CS and $
                     sqrt((x[k,j]-x[*,j])^2.+(y[k,j]-y[*,j])^2.) gt 0d0)
         for i=0L,n_elements(col)-1 do begin
            x1  = x[k,j] 	& y1  = y[k,j]
            vx1 = vx[k,j]	& vy1 = vy[k,j]
            x2  = x[col[i],j] 	& y2  = y[col[i],j]
            vx2 = vx[col[i],j]	& vy2 = vy[col[i],j]
;            vx2 = vx2 - vx2	& vy2 = vy2 - vy2 ; obj2 @ rest
;            vx1 = vx1 - vx2	& vy1 = vy1 - vy2 ; obj1 relative to obj2
            m1 = mass[k] 	& m2 = mass[col[i]]
            
            theta_v = atan(vy1,vx1)  	  ; velocity impact anlge
            theta = atan((y2-y1),(x2-x1)) ; impact angle
            a = tan(theta + theta_v)

            dvx2 = 2d0*(vx1 - vx2 +a*(vy1-vy2)) / ((1d0+a^2.)*(1d0+m2/m1))
            ;make sure the particles are heading towards eachother!
            if sqrt((x1+vx1-x2)^2+(y1+vy1-y2)^2) lt $
               sqrt((x1-x2)^2+(y1-y2)^2) then begin
               vx2f = vx2 + dvx2
               vy2f = vy2 + a*dvx2
               vx1f = vx1 - m2/m1*dvx2
               vy1f = vy1 - a*m2/m1*dvx2
               
               vx[k,j] = vx1f 	& vx[col[i],j] = vx2f
               vy[k,j] = vy1f 	& vy[col[i],j] = vy2f     
         endif
         endfor
      endif
   ENDFOR

   x[*,j+1] = x[*,j]+vx[*,j]
   y[*,j+1] = y[*,j]+vy[*,j]
   
   vx[*,j+1] = vx[*,j]
   vy[*,j+1] = vy[*,j]

ENDFOR ; END OF TIME LOOP

;=================================================================
;MAKE SOME HANDY PLOTS
!p.multi=[0,2,2]
set_plot,'PS'
device,filename='hw1p1b.eps',/landscape
plot,x[*,0],y[*,0],psym=8,xrange=[0,1],yrange=[0,1],/xstyle,/ystyle,title='t!d0!n '+strtrim(string(N),2)+' particles',xtitle='X',ytitle='Y'
plot,x[*,1999],y[*,1999],psym=8,xrange=[0,1],yrange=[0,1],/xstyle,/ystyle,title='t!df!n',xtitle='X',ytitle='Y'
plothist,sqrt(vx_init^2.+vy_init^2.),/fill,bin=.002,xrange=[0,0.05],/xstyle,xtitle='Velocity',ytitle='Num',title='v!dinit!n=0.01'
plothist,sqrt(vx[*,1999]^2.+vy[*,1999]^2.),/fill,bin=0.002,xrange=[0,0.05],/xstyle,xtitle='Velocity',ytitle='Num'
device,/close_File


vel = sqrt(vx^2.+vy^2.)


;=================================================================
!p.multi=[0,1,1]
tmp=0d0
for j=0L,499 do $
   tmp = tmp + histogram(vel[*,j+1000],bin=0.001,min=0,max=0.1)

mhist = tmp/500d0
mhistb = findgen(n_elements(mhist))*0.001
device,filename='hw1p2b.eps',/portrait
plot,mhistb+0.0005,mhist,psym=10,xrange=[0,0.025],title=strtrim(string(N),2)+' particles',xtitle='velocity',ytitle='Num'
legend,['Averaged over 500 timesteps','!3Maxwell-Boltzmann'],linestyle=[0,2],/right

vv = findgen(4000)*0.00001
Temp = v_init^2. / 3.   ; m=1, k cancels out
max_boltz = 4d0 * !dpi * N * (1d0 / (2d0 * !dpi * Temp))^(3./2.) * vv^2. * exp(-1. * vv^2./(2. * temp)) * 0.00001
oplot,vv,max_boltz*n,linestyle=2
device,/close_File




;=================================================================

mvel = dblarr(ntime+1)
svel = dblarr(ntime+1)
for j=0L,ntime do begin
   mvel[j] = mean(vel[*,j])
   svel[j] = stddev(vel[*,j])
endfor

!p.multi=[0,1,2]

device,filename='hw1p3b.eps',/portrait
plot,mvel,yrange=[0.007,0.011],xtitle='time',ytitle='<velocity>',title=strtrim(string(N),2)+' particles'
plot,svel,xtitle='time',ytitle='!7r!3!dvel!n'
device,/close_file

print,'v_init',v_init, ', <v_final>',mean(mvel[1000:*]),', frac',mean(mvel[1000:*])/v_init
print,'stddev(v_final)',stddev(mvel[1000:*])


print,''
print,'Ave Pressure',mean(pres),' = mean(total(mass*vel))' 
!p.multi=[0,1,1]
set_plot,'X'
stop
return
end
