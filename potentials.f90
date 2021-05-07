module potentials

use constants
implicit none

contains

subroutine potcalc(qpasi,qpasj,V)

real*8, dimension(1,1:ndof,1:4), intent(in) :: qpasi, qpasj
complex*16, intent(out) :: V

select case (potvar)
case ('BM1')
  call bm1(qpasi,qpasj,V)
case ('BM2')
  call bm2(qpasi,qpasj,V)
case ('HARMONIC')
  call harmonic(qpasi,qpasj,V)
case ('LHA')
  call lha(qpasi,qpasj,V)
case default
  write(*,*) "Error: Unrecognized potential type!"
  stop
end select

end subroutine

subroutine bm1(qpasi,qpasj,V)

real*8, dimension(1,1:ndof,1:4), intent(in) :: qpasi, qpasj
complex*16, intent(out) :: V

integer :: n
real*8 :: qi, qj, pbi, pbj, ai, aj, vc1, vc2
real*8 :: qi2, qj2, pbi2, pbj2, ai2, aj2, qq, pp
complex*16 :: z1, z2, zz

!Exact integration for quartic potential (Batista Models) with QTAG bases.
V=(0d0,0d0)
do n=1,ndof
  qi=qpasi(1,n,1)
  qj=qpasj(1,n,1)
  pbi=qpasi(1,n,2)
  pbj=qpasj(1,n,2)
  ai=qpasi(1,n,3)
  aj=qpasj(1,n,3)

  if (n.eq.1) then
    vc1=0.046145894861193144
    vc2=-0.5
  else
    vc1=0d0
    vc2=0.5d0
  end if

  qq=(ai*qi+aj*qj)/(ai+aj)
  pp=(pbi-pbj)/(ai+aj)

  V=V+vc1*((qq-im*pp)**4-pp**4)+6d0*vc1/(ai+aj)*(qq-im*pp)**2+3d0*vc1/(ai+aj)**2+vc2*(qq-im*pp)**2+vc2/(ai+aj)

!Coupling...
  if (n.lt.ndof) then
    qi2=qpasi(1,n+1,1)
    qj2=qpasj(1,n+1,1)
    pbi2=qpasi(1,n+1,2)
    pbj2=qpasj(1,n+1,2)
    ai2=qpasi(1,n+1,3)
    aj2=qpasj(1,n+1,3)

    z1=(ai*qi+aj*qj-im*(pbi-pbj))/(ai+aj)
    z2=(ai2*qi2+aj2*qj2-im*(pbi2-pbj2))/(ai2+aj2)
    zz=z1*z2

    V=V+vcp*zz
  end if
end do

end subroutine

subroutine bm2(qpasi,qpasj,V)

real*8, dimension(1,1:ndof,1:4), intent(in) :: qpasi, qpasj
complex*16, intent(out) :: V

integer :: n
real*8 :: qi, qj, pbi, pbj, ai, aj, vc1, vc2
real*8 :: qi2, qj2, pbi2, pbj2, ai2, aj2, qq, pp
complex*16 :: z1, z2, zz

qi=qpasi(1,1,1)
qj=qpasj(1,1,1)
pbi=qpasi(1,1,2)
pbj=qpasj(1,1,2)
ai=qpasi(1,1,3)
aj=qpasj(1,1,3)
z1=(ai*qi+aj*qj-im*(pbi-pbj))/(ai+aj)

!Exact integration for quartic potential (Batista Models) with QTAG bases.
V=(0d0,0d0)
do n=1,ndof
  qi=qpasi(1,n,1)
  qj=qpasj(1,n,1)
  pbi=qpasi(1,n,2)
  pbj=qpasj(1,n,2)
  ai=qpasi(1,n,3)
  aj=qpasj(1,n,3)

  if (n.eq.1) then
    vc1=0.046145894861193144
    vc2=-0.5
  else
    vc1=0d0
    vc2=0.5d0
  end if

  qq=(ai*qi+aj*qj)/(ai+aj)
  pp=(pbi-pbj)/(ai+aj)

  V=V+vc1*((qq-im*pp)**4-pp**4)+6d0*vc1/(ai+aj)*(qq-im*pp)**2+3d0*vc1/(ai+aj)**2+vc2*(qq-im*pp)**2+vc2/(ai+aj)

!Coupling...
  if (n.gt.1) then
    z2=((ai*qi+aj*qj-im*(pbi-pbj))/(ai+aj))**2+1d0/(ai+aj)
    zz=z1*z2
    V=V+vcp*zz
  end if
end do

end subroutine

subroutine harmonic(qpasi,qpasj,V)

real*8, dimension(1,1:ndof,1:4), intent(in) :: qpasi, qpasj
complex*16, intent(out) :: V

integer :: n
real*8 :: qi, qj, pbi, pbj, ai, aj, vc1, vc2
real*8 :: qi2, qj2, pbi2, pbj2, ai2, aj2, qq, pp
complex*16 :: z1, z2, zz

!Exact integration for quartic potential (Batista Models) with QTAG bases.
V=(0d0,0d0)
do n=1,ndof
  qi=qpasi(1,n,1)
  qj=qpasj(1,n,1)
  pbi=qpasi(1,n,2)
  pbj=qpasj(1,n,2)
  ai=qpasi(1,n,3)
  aj=qpasj(1,n,3)

  vc1=0d0
  vc2=0.5d0

  qq=(ai*qi+aj*qj)/(ai+aj)
  pp=(pbi-pbj)/(ai+aj)

  V=V+vc1*((qq-im*pp)**4-pp**4)+6d0*vc1/(ai+aj)*(qq-im*pp)**2+3d0*vc1/(ai+aj)**2+vc2*(qq-im*pp)**2+vc2/(ai+aj)

  if (n.lt.ndof) then
    qi2=qpasi(1,n+1,1)
    qj2=qpasj(1,n+1,1)
    pbi2=qpasi(1,n+1,2)
    pbj2=qpasj(1,n+1,2)
    ai2=qpasi(1,n+1,3)
    aj2=qpasj(1,n+1,3)

    z1=(ai*qi+aj*qj-im*(pbi-pbj))/(ai+aj)
    z2=(ai2*qi2+aj2*qj2-im*(pbi2-pbj2))/(ai2+aj2)
    zz=z1*z2

    V=V+vcp*zz
  end if
end do

end subroutine

subroutine lha(qpasi,qpasj,V)

real*8, dimension(1,1:ndof,1:4), intent(in) :: qpasi, qpasj
complex*16, intent(out) :: V

integer :: n
real*8 :: vv0, vv1, vv2
real*8 :: qi, qj, pbi, pbj, ai, aj
real*8 :: qi2, qj2, pbi2, pbj2, ai2, aj2
complex*16 :: vi, vj, z, z1, z2, zz

vi=(0d0,0d0)
vj=(0d0,0d0)
!LHA for any potential, needs set of functions vx..d2vx, etc.
do n=1,ndof
  qi=qpasi(1,n,1)
  qj=qpasj(1,n,1)
  pbi=qpasi(1,n,2)
  pbj=qpasj(1,n,2)
  ai=qpasi(1,n,3)
  aj=qpasj(1,n,3)

  z=(ai*qi+aj*qj+im*(pbj-pbi))/(ai+aj)

  vv0=vx(qj,n)-dvx(qj,n)*qj+d2vx(qj,n)/2d0*qj**2
  vv1=-d2vx(qj,n)*qj+dvx(qj,n)
  vv2=d2vx(qj,n)/2d0
  vj=vj+vv0+vv1*z+vv2*(z**2+1d0/(ai+aj))

  vv0=vx(qi,n)-dvx(qi,n)*qi+d2vx(qi,n)/2d0*qi**2
  vv1=-d2vx(qi,n)*qi+dvx(qi,n)
  vv2=d2vx(qi,n)/2d0
  vi=vi+vv0+vv1*z+vv2*(z**2+1d0/(ai+aj))

!Coupling...
  if (n.lt.ndof) then
    qi2=qpasi(1,n+1,1)
    qj2=qpasj(1,n+1,1)
    pbi2=qpasi(1,n+1,2)
    pbj2=qpasj(1,n+1,2)
    ai2=qpasi(1,n+1,3)
    aj2=qpasj(1,n+1,3)

    z1=z
    z2=(ai2*qi2+aj2*qj2-im*(pbi2-pbj2))/(ai2+aj2)
    zz=z1*z2

    vv0=vcpl(qi,qi2)
    vv1=-d2vcpl(qi,qi2)*z1*qi2-d2vcpl(qi,qi2)*z2*qi
    vv2=d2vcpl(qi,qi2)*zz

    vi=vi+2d0*vv0+vv1+vv2

    vv0=vcpl(qj,qj2)
    vv1=-d2vcpl(qi,qi2)*z1*qj2-d2vcpl(qi,qi2)*z2*qj
!not sure about the above line... vv1j as a function of i?
    vv2=d2vcpl(qj,qj2)*zz

    vj=vj+2d0*vv0+vv1+vv2
  end if
end do

V=(vi+vj)/2d0

end subroutine

function vx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, vp4, vx

!if (nd.eq.1) then
!!  vp4=0.046145894861193144
!!  vp2=-0.5
!  vp4=0d0
!  vp2=0.5
!  vx=vp4*x**4+vp2*x**2
!else
!  vp2=0.5
!  vx=vp2*x**2
!end if

vx=0d0

end function

function dvx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, vp4, dvx

!if (nd.eq.1) then
!  vp4=0.046145894861193144
!  vp2=-0.5
!  dvx=4d0*vp4*x**3+2d0*vp2*x
!else
!  vp2=0.5d0
!  dvx=2d0*vp2*x
!end if

dvx=0d0

end function

function d2vx(x,nd)

integer, intent(in) :: nd
real*8, intent(in) :: x
real*8 :: vp2, vp4, d2vx

!if (nd.eq.1) then
!  vp4=0.046145894861193144
!  vp2=-0.5
!  d2vx=12d0*vp4*x**2+2d0*vp2
!else
!  vp2=0.5
!  d2vx=2d0*vp2
!end if

d2vx=0d0

end function

function vcpl(x1,x2)

real*8, intent(in) :: x1, x2
real*8 :: vcpl

vcpl=vcp*x1*x2

end function

function dvcpl(x1,x2,nd)
integer, intent(in) :: nd
real*8, intent(in) :: x1, x2
real*8 :: dvcpl

if (nd.eq.1) then
  dvcpl=vcp*x2
else if (nd.eq.2) then
  dvcpl=vcp*x1
end if

end function

function d2vcpl(x1,x2)

real*8, intent(in) :: x1, x2
real*8 :: d2vcpl

d2vcpl=vcp

end function

!Exact integration for quartic potential assuming fully real basis.
!  if (n.eq.1) then
!    vc1=0.046145894861193144
!    vc2=-0.5
!  else
!    vc1=0d0
!    vc2=0.5d0
!  end if
!  V=V+(vc1*aqp**2+vc2*as**2)*aqp**2/as**4+3d0*vc1/as**2+(6d0*vc1*aqp**2+vc2*as**2)/as**3
!End quartic

end module
