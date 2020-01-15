!This module contains reexpansion algorithms for the matrix-inversion
!propagation method. They may be of reference for the diagonalization method, but
!that algorithm has its reexpansion built into the main code file.
module reexpand

use constants
use tier1
use tier2

implicit none

contains

subroutine reex_diag(nbt, gbc, p, c, co, sobolfile)

integer, intent(inout) :: nbt
real*8, dimension(1:nbt_max,1:ndof), intent(inout) :: gbc, p
complex*16, dimension(1:nbt_max), intent(inout) :: c
complex*16, dimension(1:nbt_max), intent(out) :: co
character(*), intent(in) :: sobolfile

integer :: nbtn, i, j
real*8, dimension(1:nbt_max,1:ndof) :: gbcn

if (reex.eq.'GRID') then
  call reex_grid(nbt, gbc, c, nbtn, gbcn, co)
else if (reex.eq.'SOBL') then
  call reex_sobl(nbt, gbc, c, nbtn, gbcn, co, sobolfile)
else if (reex.eq.'MDGD') then
  call reex_mdgd(nbt, gbc, c, nbtn, gbcn, co)
else if (reex.eq.'SCAN') then
  call reex_scan(nbt, gbc, c, nbtn, gbcn, co, sobolfile)
else if (reex.eq.'NONE') then
  write(*,*) "WARNING: REEXPANSIONS MAY BE NECESSARY, BUT NO METHOD SELECTED IN PARAMETER FILE"
else
  write(*,*) "UNRECOGNIZED REEXPANSION KEYWORD!"
  stop
end if

nbt=nbtn
do j=1,ndof
  do i=1,nbt_max
    if (i.le.nbt) then
      qpas(i,j,1)=qpasn(i,j,1)
    else
      qpas(i,j,1)=0d0
    end if
  end do
end do
p=0d0
call momentum(nbt, qpas(1:nbt,:,:), c(1:nbt), p(1:nbt,:), r(1:nbt,:), dp(1:nbt,:), dr(1:nbt,:))

return
end subroutine

subroutine reex_grid(nbt, gbc, c, nbtn, gbcn, co)

integer, parameter :: lwork = 1000

integer, intent(inout) :: nbt
real*8, dimension(1:nbt_max,1:ndof), intent(in) :: gbc
complex*16, dimension(1:nbt_max), intent(inout) :: c
integer, intent(out) :: nbtn
real*8, dimension(1:nbt_max,1:ndof), intent(out) :: gbcn
complex*16, dimension(1:nbt_max), intent(out) :: co

integer, dimension(1:nbt_max) :: ipiv
complex*16, dimension(1:lwork) :: work
complex*16, dimension(1:nbt_max,1:nbt_max) :: ovn
complex*16, dimension(1:nbt_max) :: ctemp

integer :: i, k, n, INFO
real*8 :: x1, x2, dx1, q1, dq1, diff1
complex*16 :: wf
real*8, dimension(1:ndof) :: xi
real*8, dimension(ndof,2) :: bound
logical :: lb

!GRID is really only good in 1D - 6/24/19
!Updated to include multi-D in MDGD - 7/1/19
gbcn=0d0
x1=qwf(1)-dsqrt(-0.5d0/awf(1)*log(tol))
x2=qwf(1)+dsqrt(-0.5d0/awf(1)*log(tol))
dq1=(x2-x1)/((nb(1)-1))

dx1 = (xdim(1,2)-xdim(1,1))/dble(npx-1)
lb=.FALSE.

do n=1,npx
  xi(1)=xdim(1,1)+(n-1)*dx1
  call psi(xi, nbt, c(1:nbt), gbc(1:nbt,:), wf)
  if (abs(wf)**2.gt.tol) then
    if (lb.eqv..FALSE.) then
      lb=.TRUE.
      bound(1,1)=xi(1)
    end if
    bound(1,2)=xi(1)
  end if
end do

i=0
k=0
do !10 begin
  i=i+1
  q1 = bound(1,1)+(i-1)*dq1
  xi(1) = q1
  call psi(xi, nbt, c(1:nbt), gbc(1:nbt,:), wf)
  if (abs(wf)**2.gt.tol) then
    k=k+1
    gbcn(k,1) = q1
  end if
  if (q1.gt.bound(1,2)) then
    diff1 = abs((q1-dq1)-bound(1,2))
    exit
  end if
end do !10 end

gbcn = gbcn+0.5d0*diff1
nbtn = k

call b_eval_reex(nbt, nbtn, gbc, gbcn, c, co)
ctemp=co
call overlap(nbtn,gbcn(1:nbtn,:),ovn(1:nbtn,1:nbtn))
call zhesv('U',nbtn,1,ovn(1:nbtn,1:nbtn),nbtn,ipiv(1:nbtn),ctemp(1:nbtn),nbtn,work,lwork,INFO)
c=ctemp

return
end subroutine

subroutine reex_sobl(nbt, qpas, c, nbtn, qpasn, co, sobolfile)

!integer, parameter :: lwork = 2500

integer, intent(inout) :: nbt
real*8, dimension(1:nbt_max,1:ndof,4), intent(inout) :: qpas
complex*16, dimension(1:nbt_max), intent(inout) :: c
integer, intent(out) :: nbtn
real*8, dimension(1:nbt_max,1:ndof,4), intent(out) :: qpasn
complex*16, dimension(1:nbt_max), intent(out) :: co
character(*), intent(in) :: sobolfile

integer :: j, k, n, d, nbt0, nmax, lwork, INFO
real*8 :: a_norm, zn
complex*16 :: wf, kold, vold, kerx, vrx, ztot
real*8, dimension(1:ndof) :: xi
real*8, dimension(ndof,2) :: bds
real*8, dimension(1:nbt_max,1:nbt) :: ov_no

integer, dimension(1:nbt_max) :: ipiv
complex*16, dimension(1:nbt_max) :: work
complex*16, dimension(1:nbt_max,1:nbt_max) :: ovn, ovtemp
complex*16, dimension(1:nbt_max) :: ctemp
logical :: idist

lwork=nbt_max
qpasn=0d0
call energy_comp(nbt, qpas(1:nbt,:,:), c(1:nbt), kold, vold)

nbt0=0
do n=1,ndof
  nbt0=nbt0+nb(n)
end do

bds(:,1)=minval(qpas(1:nbt,:,1),DIM=1)
bds(:,2)=maxval(qpas(1:nbt,:,1),DIM=1)

do d=1,ndof
  bds(d,1)=bds(d,1)-0.05*(bds(d,2)-bds(d,1))
  bds(d,2)=bds(d,2)+0.05*(bds(d,2)-bds(d,1))
end do

co=(0d0,0d0)
ovn=(0d0,0d0)

k=0
open(541,file=sobolfile)
nmax=100000
do n=1,nmax
  idist=.FALSE.
!Read the Sobol coordinates as q for the new basis function.
  read(541,*) xi
  do d=1,ndof
    xi(d)=xi(d)*(bds(d,2)-bds(d,1))+bds(d,1)
  end do
!End Sobol read.
!Check if the basis function is being placed in an area of appropriate wf density. Continue if check passes.
  call psi(xi, nbt, c(1:nbt), gbc(1:nbt,:), wf)
  if (abs(wf)**2.gt.tol) then
!If k<nbt_max, place the new basis function.
    if (.not.idist) then
      k=k+1
      if (k.gt.nbt_max) then
        write(*,*) 'REEXPANDED BASIS EXCEEDS NBT_MAX'
        stop
      end if
      do d=1,ndof
        gbcn(k,d) = xi(d)
      end do
!End placement.
!Calculate co = ov_no x c for the k-th new basis function.
      do j=1,nbt
        ztot = (1d0,0d0)
        do d=1,ndof
          zn = (gbc(j,d)-gbcn(k,d))/2.0d0
          ztot = ztot*exp(-a(d)*zn**2)
        end do
        ov_no(k,j)=ztot
      end do
      do j=1, nbt
        co(k)=co(k)+ov_no(k,j)*c(j)
      end do
!Calculate  ovn = <g_k|g_j> for the k-th new basis function.
      do j=1,k
        ztot=(1d0,0d0)
        do d=1,ndof
          zn = (gbcn(k,d)-gbcn(j,d))/2d0
          ztot = ztot*exp(-a(d)*zn**2)
        end do
        ovn(k,j) = ztot
      end do
      if (k.gt.1) then
        do j=1,k-1
          ovn(j,k) = conjg(ovn(k,j))
        end do
      end if
!If the number of new basis functions exceeds some minimum (nbt0), check to see if we
!can stop placement.
      if (mod(k,nbt0).eq.0) then
!...but first compute the new coefficients c and store them in ctemp.
        ctemp = co
        ovtemp = ovn
        call zhesv('U',k,1,ovtemp(1:k,1:k),k,ipiv(1:k),ctemp(1:k),k,work,lwork,INFO)
        if(INFO.ne.0) then
          write(*,*) 'Reex: Matrix fails, INFO =',INFO
          stop
        end if
!Check the normalization...
        call normalize(k, gbcn(1:k,:), ctemp(1:k), a_norm)
        if (abs(a_norm-1d0).lt.0.005) then
!If norm check passes, check the kinetic energy matching...
          call energy_comp(k, gbcn(1:k,:), ctemp(1:k), kerx, vrx)
!If both norm and KE matching pass, we're done placing new basis functions.
          if (abs(real(kerx)-real(kold)).lt.0.005*real(kold)) exit
        end if
      end if
    end if
  end if
  if (n.eq.nmax) write(*,*) 'SOBOL REEXPANSION END ENCOUNTERED PREMATURELY'
end do
close(541)
rewind(541)

!Record the total number of basis functions and their coefficients.
nbtn=k
c=ctemp

return
end subroutine

subroutine reex_scan(nbt, gbc, c, nbtn, gbcn, co, sobolfile)

integer, intent(inout) :: nbt
real*8, dimension(1:nbt_max,1:ndof), intent(inout) :: gbc
complex*16, dimension(1:nbt_max), intent(inout) :: c
integer, intent(out) :: nbtn
real*8, dimension(1:nbt_max,1:ndof), intent(out) :: gbcn
complex*16, dimension(1:nbt_max), intent(out) :: co
character(*), intent(in) :: sobolfile

integer :: i, j, k, n, d, nbt0, nmax, lwork, INFO
integer :: isgn, dir, na, n1, d1
real*8 :: a_norm, zn, kerr, okerr, nerr, onerr
real*8 :: kerr2, nerr2, da
complex*16 :: wf, kold, vold, kerx, vrx, ztot
real*8, dimension(1:ndof) :: xi, atmp, abest, abest2
real*8, dimension(ndof,2) :: bds
real*8, dimension(1:nbt_max,1:nbt) :: ov_no

integer, dimension(1:nbt_max) :: ipiv
complex*16, dimension(1:nbt_max) :: work
complex*16, dimension(1:nbt_max,1:nbt_max) :: ovn, ovtemp
complex*16, dimension(1:nbt_max) :: ctemp, cbest, cbest2, cobest, cobest2
logical :: idist

lwork=nbt_max
gbcn=0d0
call energy_comp(nbt, gbc(1:nbt,:), c(1:nbt), kold, vold)

nbt0=0
do n=1,ndof
  nbt0=nbt0+nb(n)
end do

bds(:,1)=minval(gbc(1:nbt,:),DIM=1)
bds(:,2)=maxval(gbc(1:nbt,:),DIM=1)

do d=1,ndof
  bds(d,1)=bds(d,1)-0.05*(bds(d,2)-bds(d,1))
  bds(d,2)=bds(d,2)+0.05*(bds(d,2)-bds(d,1))
end do

co=(0d0,0d0)
ovn=(0d0,0d0)

!Start the error terms.
okerr=1000
onerr=1000
kerr2=1000
nerr2=1000

atmp=a
da=0.05

k=0
open(541,file=sobolfile)
nmax=100000
do n=1,nmax
!Read the Sobol coordinates as q for the new basis function.
  read(541,*) xi
  do d=1,ndof
    xi(d)=xi(d)*(bds(d,2)-bds(d,1))+bds(d,1)
  end do
!End Sobol read.
!Check if the basis function is being placed in an area of appropriate wf
!density. Continue if check passes.
  call psi(xi, nbt, c(1:nbt), gbc(1:nbt,:), wf)
  if (abs(wf)**2.gt.tol) then
!If k<nbt_max, place the new basis function.
    k=k+1
    if (k.gt.nbt_max) then
      write(*,*) 'REEXPANDED BASIS EXCEEDS NBT_MAX'
      stop
    end if
    do d=1,ndof
      gbcn(k,d) = xi(d)
    end do
!End placement.
!If the number of new basis functions exceeds some minimum (nbt0), check to see
!if we can stop placement.
    if (k.eq.nbt0) then
      do d1=1, ndof
        na=0
        dir=0
        isgn=1
        do !Begin alfa loop.
          atmp(d1) = a(d1)+isgn*na*da
!Calculate co = ov_no x c for the k-th new basis function.
          do j=1,nbt
            do i=1,k
              ztot = (1d0,0d0)
              do d=1,ndof
                zn = (gbc(j,d)-gbcn(i,d))/2.0d0
                ztot = ztot*exp(-atmp(d)*zn**2)
              end do
              ov_no(i,j)=ztot
            end do
          end do
          co=(0d0,0d0)
          do i=1,k
            do j=1, nbt
              co(i)=co(i)+ov_no(i,j)*c(j)
            end do
          end do
!Calculate  ovn = <g_k|g_j> for the k-th new basis function.
          do i=1,k
            do j=1,i
              ztot=(1d0,0d0)
              do d=1,ndof
                zn = (gbcn(i,d)-gbcn(j,d))/2d0
                ztot = ztot*exp(-atmp(d)*zn**2)
              end do
              ovn(i,j) = ztot
            end do
          end do
          if (k.gt.1) then
            do i=2,k
              do j=1,i-1
                ovn(j,i) = conjg(ovn(i,j))
              end do
            end do
          end if
!Compute the new coefficients c and store them in ctemp.
          ctemp = co
          ovtemp = ovn
          call zhesv('U',k,1,ovtemp(1:k,1:k),k,ipiv(1:k),ctemp(1:k),k,work,lwork,INFO)
          if(INFO.ne.0) then
            write(*,*) 'Reex: Matrix fails, INFO =',INFO
            stop
          end if
!Check the normalization...
          call normalize(k, gbcn(1:k,:), ctemp(1:k), a_norm)
!Check the kinetic energy matching...
          call energy_comp(k, gbcn(1:k,:), ctemp(1:k), kerx, vrx)
!Record the error in both terms...
          nerr=abs(a_norm-1d0)
          kerr=abs(real(kerx)-real(kold))
!Update basis functions.
          if (kerr.lt.okerr) then
            cbest=ctemp
            abest(d1)=atmp(d1)
            cobest=co
            na=na+1
            okerr=kerr
            onerr=nerr
!Go back and push width the other way. 
          else
            isgn=-isgn
            dir=dir+1
            na=1
          end if
          if (dir.eq.2) exit
        end do !End alfa loop.
      end do
      if (kerr.lt.1e-6) exit
      do n1=1,ndof
        nbt0=nbt0+nb(n1)
      end do
    end if
  end if
  if (n.eq.nmax) write(*,*) 'SOBOL REEXPANSION END ENCOUNTERED PREMATURELY'
end do
close(541)
rewind(541)

!Record the total number of basis functions and their coefficients.

nbtn=k
c=cbest
co=cobest
a=abest

return
end subroutine

subroutine reex_mdgd(nbt, gbc, c, nbtn, gbcn, co)

integer, parameter :: lwork = 1000

integer, intent(inout) :: nbt
real*8, dimension(1:nbt_max,1:ndof), intent(in) :: gbc
complex*16, dimension(1:nbt_max), intent(inout) :: c
integer, intent(out) :: nbtn
real*8, dimension(1:nbt_max,1:ndof), intent(out) :: gbcn
complex*16, dimension(1:nbt_max), intent(out) :: co

integer :: k, n, INFO
real*8 :: x1, x2
complex*16 :: wf
integer, dimension(1:ndof) :: grid, gmax
real*8, dimension(1:ndof) :: xi, dqs

integer, dimension(1:nbt_max) :: ipiv
complex*16, dimension(1:lwork) :: work
complex*16, dimension(1:nbt_max,1:nbt_max) :: ovn
complex*16, dimension(1:nbt_max) :: ctemp

do n=1, ndof
  x1=qwf(n)-dsqrt(-0.5/awf(n)*log(tol))
  x2=qwf(n)+dsqrt(-0.5/awf(n)*log(tol))
  dqs(n)=(x2-x1)/((nb(n)-1))
  gmax(n)=ceiling((xdim(n,2)-xdim(n,1))/dqs(n))
  grid(n)=1
end do

k=0
do !11 begin
  do n=1,ndof
    xi(n)=xdim(n,1)+(grid(n)-1)*dqs(n)
  end do
  call psi(xi, nbt, c(1:nbt), gbc(1:nbt,:), wf)
  if (abs(wf)**2.gt.tol*1e5) then
    k=k+1
    do n=1,ndof
      gbcn(k,n)=xi(n)
    end do
  end if

 grid(1)=grid(1)+1
  do n=1,ndof-1
    if (grid(n).eq.gmax(n)) then
      grid(n)=1
      grid(n+1)=grid(n+1)+1
    end if
  end do
  if (grid(ndof).eq.gmax(ndof)) exit
end do !11 end


nbtn = k
call b_eval_reex(nbt, nbtn, gbc, gbcn, c, co)
ctemp = co
call overlap(nbtn,gbcn(1:nbtn,:),ovn(1:nbtn,1:nbtn))
call zhesv('U',nbtn,1,ovn(1:nbtn,1:nbtn),nbtn,ipiv(1:nbtn),ctemp(1:nbtn),nbtn,work,lwork,INFO)
c=ctemp

return
end subroutine

subroutine b_eval_reex(nbt, nbtn, gbc, gbcn, c, co)
use constants

integer, intent(in) :: nbt, nbtn
real*8, dimension(1:nbt_max,1:ndof), intent(in) :: gbc, gbcn
complex*16, dimension(1:nbt_max), intent(in) :: c
complex*16, dimension(1:nbt_max), intent(out) :: co

integer :: i, j, k
real*8 :: dq
complex*16 :: ovtot

real*8, allocatable :: ov_no(:,:)

allocate(ov_no(nbtn,nbt))

do j=1,nbt
  do i=1,nbtn
    ovtot = (1d0,0d0)
    do k=1,ndof
      dq = (gbc(j,k)-gbcn(i,k))/2.0d0
      ovtot = ovtot*exp(-a(k)*dq**2)
    end do
    ov_no(i,j)=ovtot
  end do
end do

co=(0d0,0d0)
do i=1,nbtn
  do j=1, nbt
    co(i)=co(i)+ov_no(i,j)*c(j)
  end do
end do

deallocate(ov_no)

return
end subroutine

end module
