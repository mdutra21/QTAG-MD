module finalization

use constants
use tier1

implicit none

contains

subroutine finalize(nbt, gbc, p, c, wmod)

integer, intent(in) :: nbt, wmod
real*8, dimension(1:nbt,1:ndof), intent(in) :: gbc, p
complex*16, dimension(1:nbt), intent(in) :: c


call finalize_psi(nbt, gbc, c)
call finalize_mom(nbt, p)
call fourier(wmod)

end subroutine
!--------------------------------------------------------
subroutine finalize_psi(nbt, gbc, c)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof), intent(in) :: gbc
complex*16, dimension(1:nbt), intent(in) :: c

integer :: i, j
real*8 :: dq1, dq2, x1, x2
real*8, dimension(1:ndof) :: xi
complex*16 :: wf, normchk

normchk = (0d0,0d0)

if (ndof.eq.1) then
!FOR 1D CASE:
dq1=(xdim(1,2)-xdim(1,1))/(npx-1)
open(201,file = 'wf_final.txt', status='replace', action='write')
do i=1, npx
  x1 = xdim(1,1)+dq1*(i-1)
  xi(1)=x1
  call psi(xi(1:ndof),nbt,c,gbc,wf)
  normchk = normchk + conjg(wf)*wf*dq1
  write(201,*) x1, abs(wf)**2
end do
close(201)

else if (ndof.eq.2) then
!FOR 2D CASE:
dq1=(xdim(1,2)-xdim(1,1))/(npx-1)
dq2=(xdim(2,2)-xdim(2,1))/(npx-1)
open(201,file = 'wf_final.txt', status='replace', action='write')
do i=1, npx
  x1 = xdim(1,1)+dq1*(i-1)
  xi(1)=x1
  do j=1,npx
    x2 = xdim(2,1)+dq2*(j-1)
    xi(2)=x2
    call psi(xi(1:ndof),nbt,c,gbc,wf)
    normchk = normchk + conjg(wf)*wf*dq1*dq2
    write(201,*) x1, x2, abs(wf)**2
  end do
  write(201,*) ' '
end do
close(201)
else
  stop
end if

write(*,*) 'FINAL NORMALIZATION OF WF =' , real(normchk)

end subroutine
!-------------------------------------------------------
subroutine finalize_mom(nbt, p)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof), intent(in) :: p

integer :: k, n

open(203,file='mom_final.txt')
do k=1,nbt
  write(203,*) (p(k,n), n=1,ndof)
end do
close (203)

end subroutine
!-------------------------------------------------------
subroutine fourier(wmod)
integer, parameter :: epts=50000
integer :: i, j, n, nmax, wmod
real*8 :: a, b, e, emin, emax, de, dt, alfa
real*8, dimension(1:iter_max) :: t
complex*16, dimension(1:iter_max) :: corr
complex*16, dimension(1:epts) :: ft

nmax = iter_max/wmod
open(400,file='autocorr.txt')
do n=1, nmax
  read(400,*) t(n), a, b
  corr(n) = a+im*b
end do
close (400)
dt = t(2)-t(1)
alfa = -log(1d-4)/t(nmax)**2

emin=0d0
emax=0.005d0
de=(emax-emin)/(epts-1)

ft = (0d0,0d0)
do i=1,epts
  e = emin+(i-1)*de
  do j=1,nmax
    ft(i) = ft(i)+corr(j)*exp(2d0*pi*im*e*t(j))*exp(-alfa*t(j)**2)
  end do
end do
ft = ft/(2.d0*pi)

open(404,file='ac_ft.txt')
do i=1, epts
  e = emin+(i-1)*de
  write(404,1000) 2*pi*e, ft(i)
end do
close(404)

1000  format(20f14.6,1x)

end subroutine

end module

