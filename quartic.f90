function vx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, vp4, vx

if (nd.eq.1) then
  vp4=0.046145894861193144
  vp2=-0.5
  vx=vp4*x**4+vp2*x**2
else
  vp2=0.5
  vx=vp2*x**2
end if

end function

function dvx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, vp4, dvx

if (nd.eq.1) then
  vp4=0.046145894861193144
  vp2=-0.5
  dvx=4d0*vp4*x**3+2d0*vp2*x
else
  vp2=0.5
  dvx=2d0*vp2*x
end if

end function

function d2vx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, vp4, d2vx

if (nd.eq.1) then
  vp4=0.046145894861193144
  vp2=-0.5
  d2vx=12d0*vp4*x**2+2d0*vp2
else
  vp2=0.5
  d2vx=2d0*vp2
end if

end function

function vcpl(x1,x2)

real*8, intent(in) :: x1, x2
real*8 :: vcpl

!vpc=0.2d0
vcpl=vcp*x1*x2

end function

function dvcpl(x1,x2,nd)
integer, intent(in) :: nd
real*8, intent(in) :: x1, x2
real*8 :: dvcpl

!vpc=0.2d0
if (nd.eq.1) then
  dvcpl=vcp*x2
else if (nd.eq.2) then
  dvcpl=vcp*x1
end if

end function

function d2vcpl(x1,x2)

real*8, intent(in) :: x1, x2
real*8 :: d2vcpl

!vpc=0.2d0
d2vcpl=vcp

end function
