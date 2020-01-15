!------------------------------------------------------------------------------
!Module containing subroutines needed during trajectory calculations and
!updates. Included subroutines are as follows:
!overlap
!matrices
!momentum
!normalize
!autocorr
!reexpand
!------------------------------------------------------------------------------

module initialization

use constants
use tier1
use tier2

implicit none

contains

!------------------------------------------------------------------------------

subroutine initialize(nbt, qpas, c, co, sobolfile)

integer, intent(out) :: nbt
real*8, dimension(1:nbt_max,1:ndof,1:4), intent(out) :: qpas
complex*16, dimension(1:nbt_max), intent(out) :: c, co
character(len=255) :: sobolfile
logical :: cls_chk

if (init.eq.'GRID') then
  call classical_chk(cls_chk)
  if (cls_chk.eq..TRUE.) then
    call initialize_cls(nbt, qpas, c, co)
    write(*,*) "SINGLE BASIS DETECTED - CLASSICAL BATH ASSUMED"
  else
    call initialize_grid(nbt, qpas, c, co)
  end if
else if (init.eq.'SOBL') then
  call initialize_sobl(nbt, qpas, c, co, sobolfile)
else
  write(*,*) "UNRECOGNIZED INITIALIZATION KEYWORD!"
  stop
end if

return
end subroutine
!--------------------------------------------------------------------------

subroutine initialize_grid(nbt, qpas, c, co)

integer, intent(out) :: nbt
real*8, dimension(1:nbt_max,1:ndof,1:4), intent(out) :: qpas
complex*16, dimension(1:nbt_max), intent(out) :: c, co

integer :: i, i1, i2, i3, i4, i5, i6, i7, i8, j, k, d, l
real*8 :: x11, x12, x13, x2, x3, dq1, dq2, dq3
complex*16 :: normchk

!Initialize the values of the basis coordinates, momenta, and alphas at the
!established gridpoints.
do d=1,ndof
  qpas(:,d,2)=pwf(d)
  qpas(:,d,3)=a0(d)
  qpas(:,d,4)=0d0
end do

k = 0

if (ndof.eq.1) then
!FOR 1D CASE:
  if (nb(1).gt.1) then
  x11=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
  x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
  dq1=(x2-x11)/((nb(1)-1))
  do i = 1, nb(1)
    k = k+1
    qpas(k,1,1) = x11+(k-1)*dq1
  end do

  else if (nb(1).eq.1) then
    k=k+1
    qpas(1,1,1) = qwf(1)
  else
    write(*,*) "MUST HAVE POSITIVE NUMBER OF BASIS FUNCTIONS IN ALL DIMENSIONS!"
    stop
  end if  
else if (ndof.eq.2) then
!FOR 2D CASE:
  if (nb(1).gt.1.and.nb(2).gt.1) then
  x11=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
  x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
  dq1=(x2-x11)/((nb(1)-1))
  x12=qwf(2)-dsqrt(-0.5d0/awf(2)*log(tol))
  x2=qwf(2)+dsqrt(-0.5d0/awf(2)*log(tol))
!  dq2=(x2-x12)/((nb(2)-1))

!  x12=qwf(2)-dsqrt(-0.5d0/awf(2)*log(tol))
!  x2=qwf(2)+dsqrt(-0.5d0/awf(2)*log(tol))
  dq2=(x2-x12)/20d0
  x13=qwf(2)-dq2
  do i = 1, nb(1)
    do j = 1, nb(2)
      k = k+1
!      qpas(k,1,1) = x11+(i-1)*dq1
!      qpas(k,2,1) = x12+(j-1)*dq2

      qpas(k,1,1)=x11+(i-1)*dq1
      qpas(k,2,1)=x13+(j-1)*2d0*dq2
    end do
  end do
  else if (nb(2).eq.1) then
    x11=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
    x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
    dq1=(x2-x11)/((nb(1)-1))

    do i = 1, nb(1)
      k = k+1
      qpas(k,1,1)=x11+(i-1)*dq1
!      qpas(k,2,1)=-vcp*qpas(k,1,1)
      qpas(k,2,1)=qwf(2)
    end do
  else
    write(*,*) "MUST HAVE POSITIVE NUMBER OF BASIS FUNCTIONS IN ALL DIMENSIONS!"
    stop
  end if
else if (ndof.eq.3) then
!FOR 3D CASE:
  x11=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
  x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
  dq1=(x2-x11)/((nb(1)-1))
  x12=qwf(2)-dsqrt(-0.5d0/awf(2)*log(tol))
  x2=qwf(2)+dsqrt(-0.5d0/awf(2)*log(tol))
  dq2=(x2-x12)/((nb(2)-1))
  x13=qwf(3)-dsqrt(-0.5d0/awf(3)*log(tol))
  x3=qwf(3)+dsqrt(-0.5d0/awf(3)*log(tol))
  dq3=(x2-x13)/((nb(3)-1))
  do i = 1, nb(1)
    do j = 1, nb(2)
      do l=1, nb(3)
        k=k+1
        qpas(k,1,1) = x11+(i-1)*dq1
        qpas(k,2,1) = x12+(j-1)*dq2
        qpas(k,3,1) = x13+(l-1)*dq3
      end do
    end do
  end do
else if (ndof.eq.4) then
  x11=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
  x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
  dq1=(x2-x11)/((nb(1)-1))
!  dq1=(x2-x11)/20d0

  x12=qwf(2)-dsqrt(-0.5d0/awf(2)*log(tol))
  x2=qwf(2)+dsqrt(-0.5d0/awf(2)*log(tol))
  dq2=(x2-x12)/20d0

!  x13=qwf(2)-dq2

  do i1=1,nb(1)
  do i2=1,nb(2)
  do i3=1,nb(3)
  do i4=1,nb(4)
    k=k+1
    qpas(k,1,1)=x11+(i1-1)*dq1
    qpas(k,2,1)=qwf(2)-dq2+(i2-1)*2d0*dq2
    qpas(k,3,1)=qwf(3)-dq2+(i3-1)*2d0*dq2
    qpas(k,4,1)=qwf(4)-dq2+(i4-1)*2d0*dq2

!    qpas(k,1,1)=x11+(i1-1)*dq1
!    qpas(k,2,1)=qwf(2)-dq2-0.2*qpas(k,1,1)+(i2-1)*2d0*dq2
!    qpas(k,3,1)=qwf(3)-dq2-0.8*qpas(k,2,1)+(i3-1)*2d0*dq2
!    qpas(k,4,1)=qwf(4)-dq2-0.8*qpas(k,2,1)+(i4-1)*2d0*dq2
  end do
  end do
  end do
  end do
else if (ndof.eq.6) then
  x11=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
  x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
  dq1=(x2-x11)/((nb(1)-1))
!  dq1=(x2-x11)/20d0

  x12=qwf(2)-dsqrt(-0.5d0/awf(2)*log(tol))
  x2=qwf(2)+dsqrt(-0.5d0/awf(2)*log(tol))
  dq2=(x2-x12)/20d0

  x13=qwf(2)-dq2

  do i1=1,nb(1)
  do i2=1,nb(2)
  do i3=1,nb(3)
  do i4=1,nb(4)
  do i5=1,nb(5)
  do i6=1,nb(6)
    k=k+1
!    qpas(k,1,1)=qwf(1)-dq1+(i1-1)*2d0*dq1
    qpas(k,1,1)=x11+(i1-1)*dq1
    qpas(k,2,1)=x13+(i2-1)*2d0*dq2
    qpas(k,3,1)=x13+(i3-1)*2d0*dq2
    qpas(k,4,1)=x13+(i4-1)*2d0*dq2
    qpas(k,5,1)=x13+(i5-1)*2d0*dq2
    qpas(k,6,1)=x13+(i6-1)*2d0*dq2
  end do
  end do
  end do
  end do
  end do
  end do
end if

nbt = k

open(111,file='coords.txt')
do k=1,nbt
  write(111,form1) (qpas(k,j,1),j=1,ndof)
end do
close(111)

!Create the initial array of c(t=0) values based on the equation <g_i|c_k|g_k> =
!<g_k|psi_0>.
call initialize_coeffs(nbt, qpas(1:nbt,:,:), co(1:nbt), c(1:nbt))

return
end subroutine
!--------------------------------------------------------------------------
subroutine initialize_cls(nbt, qpas, c, co)

integer, intent(out) :: nbt
real*8, dimension(1:nbt_max,1:ndof,1:4), intent(out) :: qpas
complex*16, dimension(1:nbt_max), intent(out) :: c, co

integer :: i, k, n
real*8 :: x11, x2, dq1

do n=1,ndof
  qpas(:,n,2)=pwf(n)
  qpas(:,n,3)=a0(n)
  qpas(:,n,4)=0d0
end do

x11=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
dq1=(x2-x11)/((nb(1)-1))

k=0
do i = 1, nb(1)
  k = k+1
  qpas(k,1,1)=x11+(i-1)*dq1
  do n=2,ndof
    qpas(k,n,1)=-vcp*qpas(k,n-1,1)
!    qpas(k,n,1)=qwf(n)
  end do
end do

nbt = k
call initialize_coeffs(nbt, qpas(1:nbt,:,:), co(1:nbt), c(1:nbt))

return
end subroutine
!--------------------------------------------------------------------------
subroutine initialize_sobl(nbt, qpas, c, co, sobolfile)

integer, parameter :: lwork = 1000

real*8, dimension(1:nbt_max,1:ndof,1:4), intent(out) :: qpas
complex*16, dimension(1:nbt_max), intent(out) :: c, co
integer, intent(inout) :: nbt

integer :: n, nmax, d, nbt0, j, k, INFO
real*8 :: n0, dxv, nrm, a_norm, fac, d0, d1, d2, d3, d4
real*8 :: dq, dp, ds
complex*16 :: wf0, zn, ztot

integer, dimension(1:nbt_max) :: ipiv
real*8, dimension(1:ndof) :: xi, awf1, dpx, dx, a_sum
real*8, dimension(1:ndof,2) :: bds
complex*16, dimension(1:nbt_max) :: ctemp
complex*16, dimension(1:nbt_max,1:nbt_max) :: ov, ovtemp
complex*16, dimension(1:lwork) :: work
character(len=255) :: sobolfile

write(*,*) 'Random initialization selected: Sobol dataset generated'

do d=1,ndof
  qpas(:,d,2)=pwf(d)
  qpas(:,d,3)=a0(d)
  qpas(:,d,4)=0d0
end do

do k=1,ndof
  bds(k,1)=qwf(k)-dsqrt(-0.5/awf(k)*log(tol))
  bds(k,2)=qwf(k)+dsqrt(-0.5/awf(k)*log(tol))
end do

do k=1,ndof
  awf1(k) = 2.0d0*awf(k)
end do

nbt=0
nbt0=sum(nb)

c=(0d0,0d0)
co=(0d0,0d0)
ov=(0d0,0d0)

k=0
open(105,file=sobolfile)
nmax=100000
do n=1,nmax
!Read the Sobol coordinates as q for the new basis function.
  read(105,*) xi
  do d=1,ndof
    xi(d)=xi(d)*(bds(d,2)-bds(d,1))+bds(d,1)
  end do
!End Sobol read.
!Check if the basis function is being placed in an area of appropriate wf
!density. Continue if check passes.
  wf0=(1d0,0d0)
  do d=1,ndof
    n0 = dsqrt(dsqrt(2d0*awf(d)/pi))
    dxv = xi(d)-qwf(d)
    wf0 = wf0*n0*exp(-awf(d)*dxv**2+im*pwf(d)*dxv)
  end do
  if (abs(wf0)**2.gt.tol) then
!If k<nbt_max, place the new basis function.
    k=k+1
    if (k.gt.nbt_max) then
      write(*,*) 'INITIAL BASIS EXCEEDS NBT_MAX'
      stop
    end if
    do d=1,ndof
      qpas(k,d,1) = xi(d)
    end do
!End placement.
!Calculate co = <g_k|psi0> for the k-th basis function.
    d0=0d0
    d1=0d0
    d2=0d0
    d3=0d0
    d4=0d0
    a_norm=1d0
    do d=1,ndof
      a_sum(d)=qpas(k,d,3)+awf1(d)
      a_norm=a_norm*((2d0**3*qpas(k,d,3)*awf(d))**0.25/dsqrt(a_sum(d)))
      d2=d2-awf1(d)*qpas(k,d,3)*(qpas(k,d,1)-qwf(d))**2/(2d0*a_sum(d))
      d3=d3-(qpas(k,d,2)-pwf(d))**2/(2d0*a_sum(d))
      d0=d0+qpas(k,d,3)*pwf(d)*(qpas(k,d,1)-qwf(d))/a_sum(d)
      d1=d1+awf1(d)*qpas(k,d,2)*(qpas(k,d,1)-qwf(d))/a_sum(d)
      d4=d4-qpas(k,d,4)
    end do
    co(k)=a_norm*exp(d2+d3+im*(d0+d1+d4))

!Calculate  ovn = <g_k|g_j> for the k-th new basis function.
    do j=1,k
      ztot=(1d0,0d0)
      do d=1,ndof
        fac=2d0*(qpas(k,d,3)*qpas(j,d,3))**0.25/dsqrt(2d0*qpas(k,d,3)+2d0*(qpas(k,d,3)))
        dq=(qpas(k,d,1)-qpas(j,d,1))
        dp=(qpas(k,d,2)-qpas(j,d,2))
        ds=(qpas(k,d,4)-qpas(j,d,4))
        zn=dp**2+qpas(k,d,3)*qpas(j,d,3)*dq**2+2d0*im*(qpas(k,d,3)*(ds-qpas(j,d,2)*dq)+qpas(j,d,3)*(ds-qpas(k,d,2)*dq))
        ztot = fac*ztot*exp(-zn/(2d0*(qpas(k,d,3)+qpas(j,d,3))))
      end do
      ov(k,j) = ztot
    end do
    if (k.gt.1) then
      do j=1,k-1
        ov(j,k) = conjg(ov(k,j))
      end do
    end if

!If the number of new basis functions exceeds some minimum (nbt0), check to see if we
!can stop placement.
    if (mod(k,nbt0).eq.0) then
!...but first compute the new coefficients c and store them in ctemp.
    ctemp = co
    ovtemp = ov
    call zhesv('U',k,1,ovtemp(1:k,1:k),k,ipiv(1:k),ctemp(1:k),k,work,lwork,INFO)
    if(INFO.ne.0) then
      write(*,*) 'Reex: Matrix fails, INFO =',INFO
      stop
    end if
!Check the normalization...
    call normalize(k, qpas(1:k,:,:), ctemp(1:k), nrm)
    if (0.9995.lt.nrm.and.nrm.lt.1.005) exit

    end if
  end if
end do
close(105)
rewind(105)

nbt=k
c=ctemp

do n=nbt+1,nbt_max
  qpas(n,:,:)=0d0
end do

return
end subroutine

!--------------------------------------------------------------------------
!Evaluate b_k = <g_k|psi_0>, then compute the coefficient vector c through the
!matrix equation ov_mat*c = b. The coefficients contained in c are the values
!c_k for the Gaussian basis functions. Note that this will not work for re-
!expanded c_k, because the integral is solved analytically for psi(t=0), not
!psi(t>0).
subroutine initialize_coeffs(nbt, qpas, co, c)

integer, parameter :: lwork = 1000

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(out) :: c, co

integer :: i, k, INFO
real*8, dimension(1:ndof) :: awf1, dpx, dx, a_sum, a_avg
real*8 :: a_norm, d0, d1, d2, d3, d4

complex*16, dimension(1:nbt) :: b
complex*16, dimension(1:nbt,1:nbt) :: mat

integer, dimension(1:nbt) :: ipiv
complex*16, dimension(1:lwork) :: work

c = (0d0,0d0)
do i=1, ndof
  awf1(i) = 2.0d0*awf(i)
end do

!For QTAG basis functions
do k=1,nbt
  d0=0d0
  d1=0d0
  d2=0d0
  d3=0d0
  d4=0d0
  a_norm=1d0
  do i=1,ndof
    a_sum(i)=qpas(k,i,3)+awf1(i)
    a_norm=a_norm*((2d0**3*qpas(k,i,3)*awf(i))**0.25/dsqrt(a_sum(i)))
    d2=d2-awf1(i)*qpas(k,i,3)*(qpas(k,i,1)-qwf(i))**2/(2d0*a_sum(i))
    d3=d3-(qpas(k,i,2)-pwf(i))**2/(2d0*a_sum(i))
    d0=d0+qpas(k,i,3)*pwf(i)*(qpas(k,i,1)-qwf(i))/a_sum(i)
    d1=d1+awf1(i)*qpas(k,i,2)*(qpas(k,i,1)-qwf(i))/a_sum(i)
    d4=d4-qpas(k,i,4)
  end do
  b(k)=a_norm*exp(d2+d3+im*(d0+d1+d4))
  co(k)=b(k)
end do

!Calculate ov_mat.
call overlap(nbt,qpas(1:nbt,:,:),mat)

!Solve the matrix equation ov_mat*c = b, and store the results for c.
call zhesv('U', nbt, 1, mat, nbt, ipiv, b, nbt, work, lwork, INFO)
do k=1, nbt
  c(k) = b(k)
end do

return
end subroutine

subroutine classical_chk(cls_chk)

logical, intent(out) :: cls_chk

integer :: n

cls_chk=.FALSE.
do n=1,ndof
  if (nb(n).eq.1) cls_chk=.TRUE.
end do

return
end subroutine

end module
