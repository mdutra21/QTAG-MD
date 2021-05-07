!------------------------------------------------------------------------------
!Module containing subroutines needed during trajectory calculations and
!updates. Included subroutines are as follows:
!prop_fb
!------------------------------------------------------------------------------

module propagators

use constants
use tier2
use potentials
implicit none

contains

!------------------------------------------------------------------------------
subroutine prop_diag(nbt, qpas, p, c, ov_no)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(inout) :: qpas
real*8, dimension(1:nbt,1:ndof), intent(inout) :: p
complex*16, dimension(1:nbt), intent(in) :: c
complex*16, dimension(1:nbt,1:nbt), intent(out) :: ov_no

integer :: i, j, k, n, n1
real*8 :: dq, dpb, ds, fac, V, qmax, qmin
complex*16 :: ovtot, zn

real*8, dimension(1:nbt,1:ndof,1:4) :: qpasn
real*8, dimension(1:nbt,1:ndof) :: dp, r, dr

!complex*16, dimension(1:nbt,1:nbt) :: mat

!Calculate new p and r according to Im(dpsi/psi) and Re(dpsi/psi).
call momentum(nbt, qpas, c, p, r, dp, dr)

!Momentum output...
!do i=1,nbt
!write(208,*) (p(i,n), n=1,ndof)
!end do

!All of p and r should be scaled if not using classical DoFs and tunneling
!system.
if (cls_chk.eqv..FALSE.) then
  p=p/2d0
  dp=dp/2d0
  r=r/2d0
  dr=r/2d0
else if (cls_chk.eqv..TRUE.) then
  p(:,1)=p(:,1)/2d0
  dp(:,1)=dp(:,1)/2d0
  r(:,1)=r(:,1)/2d0
  dr(:,1)=dr(:,1)/2d0
end if

!open(222,file='mom.txt',access='append')
!do k=1,nbt
!  write(222,form1) (iter*dt, qpas(k,n,1), p(k,n), n=1,ndof)
!end do
!close(222)

select case(avar)
case ("ADPT")
  do k=1, nbt
    do n=1, ndof
      qpasn(k,n,3) = qpas(k,n,3)-2d0*qpas(k,n,3)*dp(k,n)*dt/ms(n)
    end do
  end do
case ("FRZN")
  do k=1, nbt
    do n=1, ndof
      qpasn(k,n,3)=a0(n)
    end do
  end do
end select

select case(svar)
case ("ADPT")
  do k=1, nbt
    V=0d0
    do n=1, ndof-1
      V=V+vx(qpas(k,n,1),n)+vcpl(qpas(k,n,1),qpas(k,n+1,1))
    end do
    V=V+vx(qpas(k,ndof,1),ndof)
    do n=1,ndof
      qpasn(k,n,4) = qpas(k,n,4)+dt*((p(k,n)**2+r(k,n)**2+dr(k,n))/(2d0*ms(n))-V)
    end do
  end do
case ("FRZN")
  do k=1, nbt
    do n=1, ndof
      qpasn(k,n,4)=qpas(k,n,4)
    end do
  end do
end select

select case(pvar)
case ("ADPT")
  do k=1,nbt
    do n=1, ndof
!        qpasn(k,n,2)=p(k,n)/2d0
        qpasn(k,n,2)=p(k,n)
    end do
  end do
case ("FRZN")
  do k=1,nbt
    do n=1, ndof
      qpasn(k,n,2)=qpas(k,n,2)
    end do
  end do
end select

!Update positions according to q_j=p_j*dt/m.
select case(qvar)
case ("ADPT")
do k=1, nbt
  do n=1, ndof
    qpasn(k,n,1) = qpas(k,n,1)+p(k,n)*dt/(ms(n))
  end do
end do
!qmax=maxval(qpasn(:,1,1))
!qmin=minval(qpasn(:,1,1))
!if (abs(abs(qmax)-abs(qmin)).lt.0.5)  qvar="FRZN"

case ("FRZN")
do k=1,nbt
  qpasn(k,1,1)=qpas(k,1,1)
  do n=2,ndof
    qpasn(k,n,1) = qpas(k,n,1)+p(k,n)*dt/ms(n)
  end do
end do
end select

!Compute new-old overlap between iterations...
do j=1,nbt
  do i=1,nbt
    ovtot = (1d0,0d0)
    do k=1,ndof
      fac=2d0*(qpasn(i,k,3)*qpas(j,k,3))**0.25/dsqrt(2d0*qpasn(i,k,3)+2d0*(qpas(j,k,3)))
      dq=qpasn(i,k,1)-qpas(j,k,1)
      dpb=qpasn(i,k,2)-qpas(j,k,2)
      ds=qpasn(i,k,4)-qpas(j,k,4)
      zn=dpb**2+qpasn(i,k,3)*qpas(j,k,3)*dq**2+2d0*im*(qpasn(i,k,3)*(ds-qpas(j,k,2)*dq)+qpas(j,k,3)*(ds-qpasn(i,k,2)*dq))
      ovtot = fac*ovtot*exp(-zn/(2d0*(qpasn(i,k,3)+qpas(j,k,3))))
    end do
    ov_no(i,j)=ovtot
  end do
end do

do k=1,nbt
  do n=1,ndof
    do n1=1,4
      qpas(k,n,n1)=qpasn(k,n,n1)
    end do
  end do
end do

end subroutine

end module
