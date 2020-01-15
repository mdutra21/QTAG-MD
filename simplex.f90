module simplex

use constants
use tier1

implicit none

contains

subroutine reex_traj(nbt, gbc, c, ntrk)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof), intent(in) :: gbc
complex*16, dimension(1:nbt), intent(in) :: c
integer, dimension(1:ndof,2), intent(out) :: ntrk

integer, dimension(1:nbt) :: nvt
integer, allocatable :: nvals(:)
real*8, dimension(1:ndof+1) :: y
real*8, dimension(1:ndof+1,1:ndof) :: q
real*8, dimension(1:nbt) :: nlt
real*8, allocatable :: nlocs(:)
integer :: n, n2, k, k2
real*8 :: ychk, dist

!call uphill_restart(nbt,c(1:nbt),gbc(1:nbt,1:ndof),q,y,.false.)
!call uphill_simplex(q,y,ndof+1,ndof,ndof,1d-10,nbt,gbc(1:nbt,1:ndof),c(1:nbt))
!ychk = y(1)
!3 call uphill_restart(nbt,c(1:nbt),gbc(1:nbt,1:ndof),q,y,.true.)
!call uphill_simplex(q,y,ndof+1,ndof,ndof,1d-10,nbt,gbc(1:nbt,1:ndof),c(1:nbt))
!if ((y(1)-ychk).ge.1d-6) then
!  write(*,*) "MULTIPLE MAXIMA FOUND IN SIMPLEX RESTART:", ychk, y(1)
!  ychk = y(1)
!  goto 3
!end if

!do n=1,ndof+1
!  write(*,*) q(n,:)
!end do
if (ndof.eq.1) then
  ntrk(1,1) = 1
  ntrk(1,2) = nbt
else
do k=1,ndof
  n2=0
  nvt=0
  nlt=0d0
  ychk = 2.3548/(2d0*dsqrt(a(k)))
  do n=1,nbt
    dist=0d0
    do k2=1,ndof
!q(1,k2)
      if (k.ne.k2) dist=dist+(gbc(n,k2)-qwf(k2))**2
    end do
    dist = dsqrt(dist)
    if (dist.lt.ychk) then
      n2=n2+1
      nvt(n2)=n
      nlt(n2)=gbc(n,k)
    end if
  end do
  if (n2.gt.0) then
  allocate(nvals(n2))
  allocate(nlocs(n2))
  do n=1,n2
    nvals(n)=nvt(n)
    nlocs(n)=nlt(n)
  end do
  ntrk(k,1)=nvals(minloc(nlocs,DIM=1))
  ntrk(k,2)=nvals(maxloc(nlocs,DIM=1))
  deallocate(nvals)
  deallocate(nlocs)
  end if
end do
end if

end subroutine

subroutine uphill_restart(nbt,c,gbc,q,y,restart)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof), intent(in) :: gbc
complex*16, dimension(1:nbt), intent(in) :: c
logical, intent(in) :: restart

real*8, dimension(1:ndof+1), intent(inout) :: y
real*8, dimension(1:ndof+1,1:ndof), intent(inout) :: q

integer :: nn, n, k, nndof
real*8, dimension(1:ndof) :: temp
complex*16 :: wf

nndof=ndof+1

if (.not.restart) then
  do nn=1,nndof
    do n=1,ndof
      q(nn,n) = 0d0
    end do
  end do
  k=1
  do nn=2,nndof
    q(nn,k)=0.7d0
    k=k+1
  end do
  do nn=1,nndof
    do n=1,ndof
      temp(n)=q(nn,n)
    end do
    call psi(temp, nbt, c, gbc, wf)
    y(nn)=real(conjg(wf)*wf)
  end do

else if (restart) then
  do nn=2,nndof
    do n=1,ndof
      q(nn,n)=q(1,n)
      if (nn.eq.n+1) then
        q(nn,n)=q(nn,n)+1.0
      end if
    end do
  end do

  do nn=1,nndof
    do n=1,ndof
      temp(n)=q(nn,n)
    end do
    call psi(temp, nbt, c, gbc, wf)
    y(nn)=real(conjg(wf)*wf)
  end do
end if

end subroutine

subroutine uphill_simplex(q,y,mq,nq,ndim,ftol,nbt,gbc,c)

integer, intent(in) :: mq, ndim, nq, nbt
real*8, intent(in) :: ftol
real*8, dimension(1:mq), intent(inout) :: y
real*8, dimension(1:mq,1:nq), intent(inout) :: q
real*8, dimension(1:nbt,1:ndof), intent(in) :: gbc
complex*16, dimension(1:nbt), intent(in) :: c

integer :: i, ihi, ilo, inhi, j, m, n, iter, itmax
real*8 :: a, rtol, swap, ysave, ytry
real*8, dimension(1:ndim) :: qsum
complex*16 :: wf

itmax = 10000
iter = 0
1 do n = 1, ndim
  a = 0.0
  do m = 1, ndim+1
    a = a + q(m,n)
  end do
  qsum(n) = a
end do

!Determine the highest (ihi, best), next-highest (inhi, second-best), and lowest
!(ilo, worst) points of the simplex by scanning through them.
2 ilo = 1
if (y(1).gt.y(2)) then
  ihi = 1
  inhi = 2
else
  ihi = 2
  inhi = 1
end if

do i = 1, ndim+1
  if (y(i).le.y(ilo)) ilo = i
  if (y(i).gt.y(ihi)) then
    inhi = ihi
    ihi = i
  else if (y(i).gt.y(inhi)) then
    if (i.ne.ihi) inhi = i
  end if
end do

!Check to see if the simplex has converged. If so, put the best (highest) point
!in the first index.
rtol = 2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
if (rtol.lt.ftol) then
  swap = y(1)
  y(1) = y(ihi)
  y(ihi) = swap
  do n = 1, ndim
    swap = q(1,n)
    q(1,n) = q(ihi,n)
    q(ihi,n) = swap
  end do
  return
end if

!Begin a new iteration of moving the simplex. First, extrapolate through the
!face across from the lowest point.
iter = iter + 2
call simtry(q,y,qsum,mq,nq,ndim,ilo,-1.0d0,nbt,gbc,c,ytry)
if (ytry.ge.y(ihi)) then
!If the extrapolation is successful, try again with an additional factor for efficiency.
!  ytry = simtry(q,y,qsum,mq,nq,ndim,ilo,2.0d0,nbt,gbc,c)
  continue
!If the point is worse than the second highest, check for a better intermediate
!point.
else if (ytry.le.y(inhi)) then
  ysave = y(ilo)
  call simtry(q,y,qsum,mq,nq,ndim,ilo,0.5d0,nbt,gbc,c,ytry)
!If the lowest point is stubborn (not replaced), contract around best (highest)
!point.
  if (ytry.le.ysave) then
    do i = 1, ndim+1
      if (i.ne.ihi) then
        do j = 1, ndim
          qsum(j) = 0.5*(q(i,j)+q(ihi,j))
          q(i,j) = qsum(j)
        end do
        call psi(qsum, nbt, c, gbc, wf)
        y(i) = real(conjg(wf)*wf)
      end if
    end do
    iter = iter + ndim
    goto 1
  end if
else
  iter = iter - 1
end if
!if (iter.gt.itmax) goto 999
goto 2

!999 continue

end subroutine

subroutine simtry(q,y,qsum,mq,nq,ndim,ilo,fac,nbt,gbc,c, ytry)

integer :: ilo, mq, ndim, nq, nbt
real*8 :: fac
real*8, dimension(1:nq) :: qsum
real*8, dimension(1:mq) :: y
real*8, dimension(1:mq,1:nq) :: q

real*8, dimension(1:nbt,1:ndof), intent(in) :: gbc
complex*16, dimension(1:nbt), intent(in) :: c

integer :: j
real*8 :: fac1, fac2, ytry
real*8, dimension(1:ndim) :: qtry
complex*16 :: wf

fac1 = (1.0d0 - fac)/ndim
fac2 = fac1 - fac
do j = 1, ndim
  qtry(j) = qsum(j)*fac1 - q(ilo,j)*fac2
end do
call psi(qtry, nbt, c, gbc, wf)
ytry = real(conjg(wf)*wf)
if (ytry.gt.y(ilo)) then
  y(ilo) = ytry
  do j = 1, ndim
    qsum(j) = qsum(j) - q(ilo,j) + qtry(j)
    q(ilo,j) = qtry(j)
  end do
end if

return
end subroutine

end module
