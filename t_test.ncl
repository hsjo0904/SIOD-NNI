;;=========================================================
;; significant_test.ncl
;;=========================================================
df = new((/dimz(1),dimz(2)/),"integer",-999)

do ix = 0, dimz(2)-1
do iy = 0, dimz(1)-1

xi = rm_single_dims(x)
yi = rm_single_dims(y(:,iy,ix))

if(any(ismissing(xi)).or.any(ismissing(yi))) then
df(iy,ix) = -999
else

; mean
mx = dim_avg_n_Wrap(xi,0)
my = dim_avg_n_Wrap(yi,0)

; variance
vx = dim_variance_n_Wrap(xi,0)
vy = dim_variance_n_Wrap(yi,0)

nx = new((/dimz(0)/),"float")
ny = new((/dimz(0)/),"float")

; autocorrelation
do lag = 0, dimz(0)-2
  nlag = dimz(0)-2 - lag
  xacr = new((/nlag+1/),"float")
  yacr = new((/nlag+1/),"float")

  do i = 0, nlag
  xacr(i) = ((xi(i)-mx)*(xi(i+lag+1)-mx))/dimz(0)
  yacr(i) = ((yi(i)-my)*(yi(i+lag+1)-my))/dimz(0)
  end do

  if(vx.ne.0) then
  nx(lag) = dim_sum_n(xacr,0)/vx
  else
  nx(lag) = dim_sum_n(xacr,0)
  end if

  if(vy.ne.0) then
  ny(lag) = dim_sum_n(yacr,0)/vy
  else
  ny(lag) = dim_sum_n(yacr,0)
  end if

  delete(xacr)
  delete(yacr)
end do

; effective degree of freedom
en = floattointeger(dimz(0)*0.8)
c  = dim_sum_n(nx(0:en)*ny(0:en),0)
d  = 1.+2.*c
df(iy,ix) = floattointeger(dimz(0)/d)

end if
end do
end do

tval = onedtond(rc@tval,dimsizes(rc))
tval = r*sqrt(df/(1-r^2+0.000001))
b    = tval
b    = 0.5
prob = betainc(df/(df+tval^2),df/2.0,b)

copy_VarCoords(rc,df)
copy_VarCoords(rc,prob)
