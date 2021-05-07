!------------------------------------------------------------------------------
!Module containing subroutines needed during trajectory calculations and
!updates. Included subroutines are as follows:
!matrices
!momentum
!------------------------------------------------------------------------------

module tier2

use constants
use tier1

implicit none

contains

subroutine basis_diag(nbt, qpas, ov, H, co, cn, evals)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
complex*16, dimension(1:nbt,1:nbt), intent(in) :: ov, H
complex*16, dimension(1:nbt), intent(in) :: co
complex*16, dimension(1:nbt), intent(out) :: cn
real*8, dimension(1:nbt), intent(out) :: evals

integer :: i, j, n, nvals, lwork, INFO
real*8 :: abstol, dlamch

integer, dimension(1:nbt) :: IFAIL
integer, dimension(1:5*nbt) :: iwork
real*8, dimension(1:7*nbt) :: rwork
complex*16, dimension(1:nbt) :: c1
complex*16, dimension(1:nbt,1:nbt) :: err
complex*16, allocatable :: work(:)
complex*16, dimension(1:nbt,1:nbt) :: evecs

lwork=4500
allocate(work(lwork))

err=ov
abstol = 2d0*dlamch('S')
call zhegvx(1, 'V', 'A', 'L', nbt, H, nbt, ov, nbt, 1.0, 2.0, 1, 2, abstol, &
            nvals, evals, evecs, nbt, work, lwork, rwork, iwork, IFAIL, INFO)

if (INFO.ne.0) then
  write(*,*) "ERROR NUMBER:", INFO
  write(*,*) "basis_diag: Error in basis diagonalization!"
  write(*,*) 'POSITION INFO'
  write(*,*) (qpas(INFO-nbt,n,1), n=1,ndof)
  write(*,*) 'PHASE INFO'
  write(*,*) (qpas(INFO-nbt,n,2), n=1,ndof)
  write(*,*) 'WIDTH INFO'
  write(*,*) (qpas(INFO-nbt,n,3), n=1,ndof)
  write(*,*) 'OVERLAP INFO'
  do n=1,nbt
  write(*,*) err(INFO-nbt,n)
  end do
  stop
end if

!write(*,*) 'EVALS = ', evals(1), evals(2)
!write(*,*) 'EVALS2 = ', evals(3), evals(4)
write(402,form1) evals(1), evals(2), evals(3), evals(4)

c1=matmul(conjg(transpose(evecs)),co)
do n=1,nbt
  c1(n)=c1(n)*exp(-im*dt*evals(n))
end do
cn=matmul(evecs,c1)

deallocate(work)

return
end subroutine
!------------------------------------------------------------------------------
subroutine momentum(nbt, qpas, c, p, r, dp, dr)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
real*8, dimension(1:nbt,1:ndof), intent(out) :: p, r, dp, dr

if (pcalc.eq.'UMOD') then
  call mom_conv(nbt, qpas, c, p, r, dp, dr,'F')
else if (pcalc.eq.'PFIT') then
  call mom_fit(nbt, qpas, c, p, r, dp, dr)
else if (pcalc.eq.'PCNV') then
  call mom_conv(nbt, qpas, c, p, r, dp, dr,'T')
else if (pcalc.eq.'WCNV') then
  call psi_conv(nbt, qpas, c, p, r, dp, dr)
else
  write(*,*) "UNRECOGNIZED MOMENTUM KEYWORD!"
  stop
end if

!call basis_repulsion(nbt, qpas(1:nbt,:,:), p(1:nbt,:))

return
end subroutine

subroutine mom_fit(nbt, qpas, c, p, r, dp, dr)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
real*8, dimension(1:nbt,1:ndof), intent(inout) :: p, r, dp, dr

integer :: i, j, n
real*8 :: dq, fac
complex*16 :: z, prod1

real*8, dimension(1:nbt,1:ndof) :: zi, ri
complex*16, dimension(1:ndof) :: dz

do i=1,nbt
  z=(0d0,0d0)
  dz=(0d0,0d0)
  do j=1,nbt
    fac=1d0
    prod1=(1d0,0d0)
    do n=1,ndof
      dq=qpas(i,n,1)-qpas(j,n,1)
      fac=fac*(qpas(j,n,3)/pi)**0.25
      prod1=prod1*exp(-0.5*qpas(j,n,3)*dq**2+im*qpas(j,n,2)*dq+im*qpas(j,n,4))
    end do
    z=z+c(j)*fac*prod1
    do n=1,ndof
      dq=qpas(i,n,1)-qpas(j,n,1)
      dz(n)=dz(n)-(qpas(j,n,3)*dq-im*qpas(j,n,2))*c(j)*fac*prod1
    end do
  end do
  do n=1,ndof
    zi(i,n)=aimag(dz(n)/z)
    ri(i,n)=real(dz(n)/z)
  end do
end do

if (cls_chk.eqv..FALSE.) then
!Quantum treatment of bath modes
  p=zi
  r=ri

  do n=1,ndof
    call qfit1D(nbt, n, qpas(1:nbt,:,:), c, p(1:nbt,n), r(1:nbt,n), dp(1:nbt,n), dr(1:nbt,n))
  end do

!  call qfit2D(nbt, qpas(1:nbt,:,:), p(1:nbt,:), r(1:nbt,:), dp(1:nbt,:), dr(1:nbt,:))

else if (cls_chk.eqv..TRUE.) then
!Classical treatment of bath modes.
  p(:,1)=zi(:,1)
  r(:,1)=ri(:,1)

  call qfit1D(nbt, 1, qpas, c, p(1:nbt,1), r(1:nbt,1), dp(1:nbt,1), dr(1:nbt,1))
  call cfitbath(nbt, qpas, p(1:nbt,1:ndof), dp(1:nbt,1:ndof))

  do n=2,ndof
   r(:,n)=0d0
   dr(:,n)=0d0
  end do
end if

return
end subroutine

subroutine mom_conv(nbt, qpas, c, p, r, dp, dr, iconv)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
real*8, dimension(1:nbt,1:ndof), intent(inout) :: p, r, dp, dr
character(len=1), intent(in) :: iconv

integer :: i, j, n
real*8 :: dq, fac, cpl1, cpl2
complex*16 :: z, prod1

real*8, dimension(1:nbt) :: wdws
real*8, dimension(1:nbt,1:nbt) :: wdw
real*8, dimension(1:nbt,1:ndof) :: zi, ri, ppl, rpl, pmn, rmn, po
real*8, dimension(1:ndof) :: delta
complex*16, dimension(1:ndof) :: dz

if (cls_chk.eqv..TRUE.) then
  write(*,*) "MOMENTUM CONVOLUTION DOESN'T YET SUPPORT CLASSICAL BATH!"
  stop
end if

po=p
!For convolution of momentum.
do i=1,nbt
  z=(0d0,0d0)
  dz=(0d0,0d0)
  do j=1,nbt
    fac=1d0
    prod1=(1d0,0d0)
    do n=1,ndof
      dq=qpas(i,n,1)-qpas(j,n,1)
      fac=fac*(qpas(j,n,3)/pi)**0.25
      prod1=prod1*exp(-0.5*qpas(j,n,3)*dq**2+im*qpas(j,n,2)*dq+im*qpas(j,n,4))
    end do
    z=z+c(j)*fac*prod1
    do n=1,ndof
      dq=qpas(i,n,1)-qpas(j,n,1)
      dz(n)=dz(n)-(qpas(j,n,3)*dq-im*qpas(j,n,2))*c(j)*fac*prod1
    end do
  end do
  do n=1,ndof
    zi(i,n)=aimag(dz(n)/z)
    ri(i,n)=real(dz(n)/z)
  end do
end do

if (iconv.eq.'T') then
  wdw=1d0
  do j=1, nbt
    do i=1, nbt
      do n=1, ndof
        wdw(j,i)=wdw(j,i)*exp(-beta*(qpas(i,n,1)-qpas(j,n,1))**2)
      end do
    end do
  end do

  p=matmul(wdw,zi)
  r=matmul(wdw,ri)
  wdws=sum(wdw,DIM=2)
  do j=1,nbt
    do n=1,ndof
      p(j,n)=p(j,n)/wdws(j)
      r(j,n)=r(j,n)/wdws(j)
    end do
  end do
else if (iconv.eq.'F') then
  p=zi
  r=ri
end if

!Classical bath modes (only for frozen-real basis).
do i=1,nbt
  do n=2,ndof
!    cpl1=dvcpl(qpas(i,n-1,1),qpas(i,n,1),1)+dvcpl(qpas(i,n-1,1),qpas(i,n,1),2)
    cpl1=vcp*(qpas(i,n-1,1))
    if (n.lt.ndof) then
!      cpl2=dvcpl(qpas(i,n,1),qpas(i,n+1,1),1)+dvcpl(qpas(i,n,1),qpas(i,n+1,1),2)
      cpl2=vcp*(qpas(i,n+1,1))
    else
      cpl2=0d0
    end if
    p(i,n)=po(i,n)-(2d0*0.5d0*qpas(i,n,1)+cpl1+cpl2)*dt
  end do
end do
!End classical bath.

if (avar.eq.'ADPT'.or.svar.eq.'ADPT') then

do n=1,ndof
  delta(n)=(xdim(n,2)-xdim(n,1))/(npx)
end do

do i=1,nbt
  z=(0d0,0d0)
  dz=(0d0,0d0)
  do j=1,nbt
    fac=1d0
    prod1=(1d0,0d0)
    do n=1,ndof
      dq=qpas(i,n,1)+delta(n)-qpas(j,n,1)
      fac=fac*(qpas(j,n,3)/pi)**0.25
      prod1=prod1*exp(-0.5*qpas(j,n,3)*dq**2+im*qpas(j,n,2)*dq+im*qpas(j,n,4))
    end do
    z=z+c(j)*fac*prod1
    do n=1,ndof
      dq=qpas(i,n,1)+delta(n)-qpas(j,n,1)
      dz(n)=dz(n)-(qpas(j,n,3)*dq-im*qpas(j,n,2))*c(j)*fac*prod1
    end do
  end do
  do n=1,ndof
    zi(i,n)=aimag(dz(n)/z)
    ri(i,n)=real(dz(n)/z)
  end do
end do

if (iconv.eq.'T') then
  wdw=1d0
  do j=1, nbt
    do i=1, nbt
      do n=1, ndof
        wdw(j,i)=wdw(j,i)*exp(-beta*(qpas(i,n,1)-(qpas(j,n,1)+delta(n)))**2)
      end do
    end do
  end do

  ppl=matmul(wdw,zi)
  rpl=matmul(wdw,ri)
  wdws=sum(wdw,DIM=2)
  do j=1,nbt
    do n=1,ndof
      ppl(j,n)=ppl(j,n)/wdws(j)
      rpl(j,n)=rpl(j,n)/wdws(j)
    end do
  end do
else if (iconv.eq.'F') then
  ppl=zi
  rpl=ri
end if

do i=1,nbt
  z=(0d0,0d0)
  dz=(0d0,0d0)
  do j=1,nbt
    fac=1d0
    prod1=(1d0,0d0)
    do n=1,ndof
      dq=qpas(i,n,1)-delta(n)-qpas(j,n,1)
      fac=fac*(qpas(j,n,3)/pi)**0.25
      prod1=prod1*exp(-0.5*qpas(j,n,3)*dq**2+im*qpas(j,n,2)*dq+im*qpas(j,n,4))
    end do
    z=z+c(j)*fac*prod1
    do n=1,ndof
      dq=qpas(i,n,1)-delta(n)-qpas(j,n,1)
      dz(n)=dz(n)-(qpas(j,n,3)*dq-im*qpas(j,n,2))*c(j)*fac*prod1
    end do
  end do
  do n=1,ndof
    zi(i,n)=aimag(dz(n)/z)
    ri(i,n)=real(dz(n)/z)
  end do
end do

if (iconv.eq.'T') then
  wdw=1d0
  do j=1, nbt
    do i=1, nbt
      do n=1, ndof
        wdw(j,i)=wdw(j,i)*exp(-beta*(qpas(i,n,1)-(qpas(j,n,1)-delta(n)))**2)
      end do
    end do
  end do

  pmn=matmul(wdw,zi)
  rmn=matmul(wdw,ri)
  wdws=sum(wdw,DIM=2)
  do j=1,nbt
    do n=1,ndof
      pmn(j,n)=pmn(j,n)/wdws(j)
      rmn(j,n)=rmn(j,n)/wdws(j)
    end do
  end do
else if (iconv.eq.'F') then
  pmn=zi
  rmn=ri
end if

do n=1,ndof
  dp(:,n)=(ppl(:,n)-pmn(:,n))/(2d0*delta(n))
  dr(:,n)=(rpl(:,n)-rmn(:,n))/(2d0*delta(n))
end do

else

dp=0d0
dr=0d0

end if

return
end subroutine

subroutine psi_conv(nbt, qpas, c, p, r, dp, dr)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
real*8, dimension(1:nbt,1:ndof), intent(out) :: p, r, dp, dr

integer :: i, j, n
real*8 :: dq, dab, b1, fac
complex*16 :: z, z0, prod1

real*8, dimension(1:nbt,1:ndof) :: zi, ri, ppl, rpl, pmn, rmn
real*8, dimension(1:ndof) :: delta
complex*16, dimension(1:ndof) :: dz

if (cls_chk.eqv..TRUE.) then
  write(*,*) "WAVEFUNCTION CONVOLUTION DOESN'T YET SUPPORT CLASSICAL BATH!"
  stop
end if

b1=0d0
!For convolution of the wavefunction, psi_b.
!Wavefunction convolution for fully adaptable basis.
do i=1, nbt
  z = (0d0,0d0)
  dz = (0d0,0d0)
  do j=1, nbt
    prod1=(1d0,0d0)
    fac=1d0
    do n=1,ndof
      dab=qpas(j,n,3)+beta
      fac=fac*(qpas(j,n,3)/pi)**0.25*dsqrt(beta/dab)
      dq=qpas(j,n,1)-qpas(i,n,1)

      prod1=prod1*exp(-0.5/dab*(qpas(i,n,2)**2+qpas(i,n,3)*beta*dq**2+2d0*im*(qpas(i,n,2)*beta*dq &    
            -qpas(i,n,4)*dab)))
!      prod1=prod1*exp(-0.5/dab*(qpas(j,n,3)*beta*dq**2+2d0*im*dq*(qpas(j,n,3)*b1+qpas(j,n,2)*beta) &
!            -2d0*im*qpas(j,n,4)*(qpas(j,n,3)+beta)+(b1-qpas(j,n,2))**2))
    end do
    z0=c(j)*fac*prod1
    z=z+z0
    do n=1,ndof
      dq=qpas(j,n,1)-qpas(i,n,1)
      dab=qpas(j,n,3)+beta
      dz(n)=dz(n)-(beta/dab)*(qpas(i,n,3)*dq-im*qpas(i,n,2))*z0
!      dz(n)=dz(n)+(qpas(j,n,3)*beta*dq+im*(b1*qpas(j,n,3)+beta*qpas(j,n,2)))/dab*z0
    end do
  end do
  do n=1, ndof
    p(i,n) = aimag(dz(n)/z)
    r(i,n) = real(dz(n)/z)
  end do
end do

if (avar.eq.'ADPT'.or.svar.eq.'ADPT') then

do n=1,ndof
  delta(n)=(xdim(n,2)-xdim(n,1))/(npx)
end do

do i=1, nbt
  z = (0d0,0d0)
  dz = (0d0,0d0)
  do j=1, nbt
    prod1=(1d0,0d0)
    fac=1d0
    do n=1,ndof
      dab=qpas(j,n,3)+beta
      fac=fac*(qpas(j,n,3)/pi)**0.25*dsqrt(beta/dab)
      dq=qpas(j,n,1)-qpas(i,n,1)-delta(n)

      prod1=prod1*exp(-0.5/dab*(qpas(i,n,2)**2+qpas(i,n,3)*beta*dq**2+2d0*im*(qpas(i,n,2)*beta*dq &
            -qpas(i,n,4)*dab)))
!      prod1=prod1*exp(-0.5/dab*(qpas(j,n,3)*beta*dq**2+2d0*im*dq*(qpas(j,n,3)*b1+qpas(j,n,2)*beta) &
!            -2d0*im*qpas(j,n,4)*(qpas(j,n,3)+beta)+(b1-qpas(j,n,2))**2))
    end do
    z0=c(j)*fac*prod1
    z=z+z0
    do n=1,ndof
      dq=qpas(j,n,1)-qpas(i,n,1)-delta(n)
      dab=qpas(j,n,3)+beta
      dz(n)=dz(n)-(beta/dab)*(qpas(i,n,3)*dq-im*qpas(i,n,2))*z0
!      dz(n)=dz(n)+(qpas(j,n,3)*beta*dq+im*(b1*qpas(j,n,3)+beta*qpas(j,n,2)))/dab*z0
    end do
  end do
  do n=1, ndof
    ppl(i,n) = aimag(dz(n)/z)
    rpl(i,n) = real(dz(n)/z)
  end do
end do

do i=1, nbt
  z = (0d0,0d0)
  dz = (0d0,0d0)
  do j=1, nbt
    prod1=(1d0,0d0)
    fac=1d0
    do n=1,ndof
      dab=qpas(j,n,3)+beta
      fac=fac*(qpas(j,n,3)/pi)**0.25*dsqrt(beta/dab)
      dq=qpas(j,n,1)-qpas(i,n,1)+delta(n)

      prod1=prod1*exp(-0.5/dab*(qpas(i,n,2)**2+qpas(i,n,3)*beta*dq**2+2d0*im*(qpas(i,n,2)*beta*dq &
            -qpas(i,n,4)*dab)))
!      prod1=prod1*exp(-0.5/dab*(qpas(j,n,3)*beta*dq**2+2d0*im*dq*(qpas(j,n,3)*b1+qpas(j,n,2)*beta) &
!            -2d0*im*qpas(j,n,4)*(qpas(j,n,3)+beta)+(b1-qpas(j,n,2))**2))
    end do
    z0=c(j)*fac*prod1
    z=z+z0
    do n=1,ndof
      dq=qpas(j,n,1)-qpas(i,n,1)+delta(n)
      dab=qpas(j,n,3)+beta
      dz(n)=dz(n)-(beta/dab)*(qpas(i,n,3)*dq-im*qpas(i,n,2))*z0
!      dz(n)=dz(n)+(qpas(j,n,3)*beta*dq+im*(b1*qpas(j,n,3)+beta*qpas(j,n,2)))/dab*z0
    end do
  end do
  do n=1, ndof
    pmn(i,n) = aimag(dz(n)/z)
    rmn(i,n) = real(dz(n)/z)
  end do
end do

do n=1,ndof
  dp(:,n)=(ppl(:,n)-pmn(:,n))/(2d0*delta(n))
  dr(:,n)=(rpl(:,n)-rmn(:,n))/(2d0*delta(n))
end do

else

dp=0d0
dr=0d0

end if

return
end subroutine

subroutine basis_repulsion(nbt, qpas, p)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
real*8, dimension(1:nbt,1:ndof), intent(inout) :: p
integer :: j, jj, n
real*8 :: dq1, dq2, dq, aavg, fq
real*8, dimension(1:ndof) :: pmax

!Artificial basis repulsion for arbitrary basis ordering.
pmax=maxval(p,DIM=1)
do j=1,nbt
  do jj=1,nbt
    if (jj.ne.j) then
      do n=1,ndof
        dq1=qpas(jj,n,1)-qpas(j,n,1)
        if (n.lt.ndof) then
          dq2=qpas(jj,ndof,1)-qpas(j,ndof,1)
        else
          if (ndof.eq.1) then
            dq2=0d0
          else
            dq2=qpas(jj,ndof-1,1)-qpas(j,ndof-1,1)
          end if
        end if
        dq=abs(sqrt(dq1**2+dq2**2)*cos(atan(dq2/dq1)))
        if (dq.gt.1e-6) then
          aavg=0.5*(qpas(j,n,3)+qpas(jj,n,3))
          fq=sign(pmax(1)*exp(-aavg*dq),-dq1)
        else
          fq=0.0
        end if
        p(j,n)=p(j,n)+fq
      end do
    end if
  end do
end do
!End artificial basis repulsion.

return
end subroutine

subroutine qfit1D(nbt, sn, qpas, c, sp, sr, sdp, sdr)
!Fitting for fully quantum DoFs, where px and py are independent.

integer, intent(in) :: nbt, sn
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
real*8, dimension(1:nbt), intent(inout) :: sp, sr, sdp, sdr
complex*16, dimension(1:nbt), intent(in) :: c

integer :: i, j, k, INFO
complex*16 :: z

real*8, dimension(1:2,1:2) :: b, ov

b=0d0
ov=0d0
do i=1,2
  b(i,:) = 0d0
  do k=1, nbt
    call psi(qpas(k,:,1), nbt, c, qpas, z)
!    b(i,1) = b(i,1)+sp(k)*qpas(k,sn,1)**(i-1) !not weighted
!    b(i,2) = b(i,1)+sr(k)*qpas(k,sn,1)**(i-1) !not weighted
    b(i,1) = b(i,1)+conjg(z)*z*sp(k)*qpas(k,sn,1)**(i-1) !weighted by psi**2
    b(i,2) = b(i,2)+conjg(z)*z*sr(k)*qpas(k,sn,1)**(i-1) !weighted by psi**2
  end do

  do j=1,2
    ov(j,i)=0d0
    do k=1,nbt
      call psi(qpas(k,:,1), nbt, c, qpas, z)
!      ov(j,i)=ov(j,i)+qpas(k,sn,1)**(i-1+j-1) !not weighted
      ov(j,i)=ov(j,i)+qpas(k,sn,1)**(i-1+j-1)*conjg(z)*z !weighted by psi**2
    end do
  end do
end do

call dposv('U',2,2,ov,2,b,2,INFO)

do k=1,nbt
  sp(k)=b(1,1)+b(2,1)*qpas(k,sn,1)
  sr(k)=b(1,2)+b(2,2)*qpas(k,sn,1)
  sdp(k)=b(2,1)
  sdr(k)=b(2,2)
end do

return
end subroutine

subroutine cfitbath(nbt, qpas, p, dp)
!Update of p and dp for classical bath; no fitting.

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
real*8, dimension(1:nbt,1:ndof), intent(inout) :: p, dp

integer :: i, n
real*8 :: cpl1, cpl2

if (potvar.ne.'BM2') then
  do i=1,nbt
    do n=2,ndof
!      cpl1=dvcpl(qpas(i,n-1,1),qpas(i,n,1),1)+dvcpl(qpas(i,n-1,1),qpas(i,n,1),2)
      cpl1=vcp*(qpas(i,n-1,1))
      if (n.lt.ndof) then
!        cpl2=dvcpl(qpas(i,n,1),qpas(i,n+1,1),1)+dvcpl(qpas(i,n,1),qpas(i,n+1,1),2)
        cpl2=vcp*(qpas(i,n+1,1))
      else
        cpl2=0d0
      end if
      p(i,n)=p(i,n)-(2d0*0.5d0*qpas(i,n,1)+cpl1+cpl2)*dt
      dp(i,n)=0d0 ! Frozen b/c classical DoFs are coherent.
    end do
  end do
else
  do i=1,nbt
    do n=2,ndof
!      cpl1=dvcpl(qpas(i,n-1,1),qpas(i,n,1),1)+dvcpl(qpas(i,n-1,1),qpas(i,n,1),2)
      cpl1=2d0*vcp*qpas(i,1,1)*qpas(i,n,1)
      p(i,n)=p(i,n)-(2d0*0.5d0*qpas(i,n,1)+cpl1)*dt
      dp(i,n)=0d0 ! Frozen b/c classical DoFs are coherent.
    end do
  end do
end if

return
end subroutine

subroutine qfit2D(nbt, qpas, p, r, dp, dr)
!Fitting for fully quantum DoFs, where px and py are codependent.

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
real*8, dimension(1:nbt,1:ndof), intent(inout) :: p, r, dp, dr

integer :: i, k, n, INFO
real*8 :: qq

real*8, dimension(1:3,1:3) :: b2
real*8, dimension(1:3,1:4) :: ov2

if (ndof.ne.2) then
  write(*,*) '2D fitting invoked for N!=2'
  stop
end if

b2=0d0
do i=1,ndof+1
  do n=1,ndof
    do k=1, nbt
!      call psi(qpas(k,:,1), nbt, c, qpas, z)
      if (i.lt.ndof+1) then
        qq=qpas(k,i,1)
      else
        qq=1d0
      end if
      b2(i,2*n-1) = b2(i,2*n-1)+p(k,n)*qq
      b2(i,2*n) = b2(i,2*n)+r(k,n)*qq
    end do
  end do
end do

ov2=0d0
do k=1,nbt
!  call psi(qpas(k,:,1), nbt, c, qpas, z)
  ov2(1,1)=ov2(1,1)+qpas(k,1,1)**2
  ov2(1,2)=ov2(1,2)+qpas(k,1,1)*qpas(k,2,1)
  ov2(1,3)=ov2(1,3)+qpas(k,1,1)
  ov2(2,2)=ov2(2,2)+qpas(k,2,1)**2
  ov2(2,3)=ov2(2,3)+qpas(k,2,1)
  ov2(3,3)=ov2(3,3)+1
end do
ov2(2,1)=ov2(1,2)
ov2(3,1)=ov2(1,3)
ov2(3,2)=ov2(2,3)

call dposv('U',3,4,ov2,3,b2,3,INFO)

do n=1,ndof
  do k=1,nbt
    p(k,n)=b2(1,2*n-1)*qpas(k,1,1)+b2(2,2*n-1)*qpas(k,2,1)+b2(3,2*n-1)
    r(k,n)=b2(1,2*n)*qpas(k,1,1)+b2(2,2*n)*qpas(k,2,1)+b2(3,2*n)
    dp(k,n)=b2(n,2*n-1)
    dr(k,n)=b2(n,2*n)
  end do
end do

return
end subroutine

subroutine psi_update(nbt,qpas,c,t,fnum)

integer, intent(in) :: nbt, fnum
real*8, intent(in) :: t
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c

integer :: i, j
real*8 :: dq1, dq2, x1, x2
complex*16 :: wf
real*8, dimension(1:ndof) :: xi

if (ndof.eq.1) then
  dq1=(xdim(1,2)-xdim(1,1))/(npx-1)
  do i=1, npx
    xi(1)=xdim(1,1)+dq1*(i-1)
    call psi(xi(1:ndof),nbt,c,qpas,wf)
    write(fnum,*) xi(1), abs(wf), abs(wf)**2
  end do
else if (ndof.eq.2) then
  dq1=(xdim(1,2)-xdim(1,1))/(npx-1)
  dq2=(xdim(2,2)-xdim(2,1))/(npx-1)
  do i=1, npx
    xi(1)=xdim(1,1)+dq1*(i-1)
    do j=1,npx
      xi(2)=xdim(2,1)+dq2*(j-1)
      call psi(xi(1:ndof),nbt,c,qpas,wf)
      write(fnum,*) xi(1), xi(2), abs(wf)
!, abs(wf)**2
    end do
    write(fnum,*)
  end do
else
  write(*,*) "WF output not available"
end if

!if (fnum.ne.999) then
!  do i=1,nbt
!    call psi(qpas(i,:,1), nbt, c, qpas, wf)
!    write(205,form1) t, (qpas(i,j,1),j=1,ndof), abs(wf)**2, 0.0
!  end do
!  write(205,*)
!else
!  open(209,file='wft.dat', status='replace', action='write')
!  if (ndof.eq.1) then
!!FOR 1D CASE:
!    dq1=(xdim(1,2)-xdim(1,1))/(npx-1)
!    do i=1, npx
!      x1 = xdim(1,1)+dq1*(i-1)
!      xi(1)=x1
!      call psi(xi(1:ndof),nbt,c,gbc,wf)
!      write(209,form1) x1, abs(wf)**2
!    end do
!  else if (ndof.eq.2) then
!!FOR 2D CASE:
!    dq1=(xdim(1,2)-xdim(1,1))/(npx-1)
!    dq2=(xdim(2,2)-xdim(2,1))/(npx-1)
!    do i=1, npx
!      x1 = xdim(1,1)+dq1*(i-1)
!      xi(1)=x1
!      do j=1,npx
!        x2 = xdim(2,1)+dq2*(j-1)
!        xi(2)=x2
!        call psi(xi(1:ndof),nbt,c,gbc,wf)
!        write(209,form1) x1, x2, abs(wf)**2
!      end do
!      write(209,form1) ' '
!    end do
!  end if
!  close (209)
!end if

end subroutine

end module
