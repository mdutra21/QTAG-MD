function vx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, vx

vp2=0.5d0
vx=vp2*x**2

end function

function dvx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, dvx

vp2=0.5d0
dvx=2d0*vp2*x

end function

function d2vx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, d2vx

vp2=0.5d0
d2vx=2d0*vp2

end function

function vcpl(x1,x2)

real*8, intent(in) :: x1, x2
real*8 :: vpc, vcpl

vpc=0.2d0
vcpl=vpc*x1*x2

end function

function dvcpl(x1,x2,nd)
integer, intent(in) :: nd
real*8, intent(in) :: x1, x2
real*8 :: vpc, dvcpl

vpc=0.2d0
if (nd.eq.1) then
  dvcpl=vpc*x2
else if (nd.eq.2) then
  dvcpl=vpc*x1
end if

end function

function d2vcpl(x1,x2)

real*8, intent(in) :: x1, x2
real*8 :: d2vcpl

vpc=0.2d0
d2vcpl=vpc

end function
