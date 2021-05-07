PROGRAM qtGBdiag

use constants
use initialization
!use reexpand
use tier1
use tier2
use propagators

IMPLICIT NONE

integer :: i0, i30,i60
complex*16 :: wf0, wf30, wf60

integer :: i, j, k, n, nbt, iter, ichk, idedq, nrx
real*8 :: norm, t, einit, ediag, etest, eold
complex*16 :: ke0, v0, ket, vt, et, zcorr, tcorr, KE, V
integer, allocatable :: ntrk(:,:)
real*8, allocatable :: es(:), qsl(:,:), qsu(:,:)
!, qs(:,:,:)
real*8, allocatable :: dedql(:), dedqu(:), rxtols(:)
!, dedq(:,:)
real*8, allocatable :: qpas(:,:,:), evals(:), p(:,:), xavg(:)
!pn(:,:)
complex*16, allocatable :: c(:), co(:), ov_no(:,:), edim(:), ov(:,:), H(:,:)
logical :: ireex
character(len=255) :: sobolfile
character(len=4) :: handle

call read_params()
call classical_chk()

if (reex.eq.'SOBL'.or.init.eq.'SOBL') then
  call sobol_dataset(ndof,1000,skip,sobolfile)
  write(*,*) 'Random reexpansions selected: Sobol dataset generated'
else
  sobolfile='arbitrary'
end if
write(handle,1337) ihandle
1337 format(I4.4)

allocate(qpas(nbt_max,ndof,4))
allocate(p(nbt_max,ndof))
allocate(c(nbt_max))
allocate(co(nbt_max))
allocate(ntrk(ndof,2))
allocate(rxtols(ndof))
allocate(dedql(ndof))
allocate(dedqu(ndof))
!allocate(dedq(nbt_max,ndof))
allocate(es(4))
allocate(qsl(ndof,4))
allocate(qsu(ndof,4))
!allocate(qs(nbt_max,ndof,4))
allocate(edim(ndof))
allocate(evals(nbt_max))
allocate(ov_no(nbt_max,nbt_max))
allocate(ov(nbt_max,nbt_max))
allocate(H(nbt_max,nbt_max))
allocate(xavg(ndof))

open(401,file='energy'//trim(handle)//'.txt', status='replace', action='write')
open(402,file='evals'//trim(handle)//'.txt', status='replace', action='write')
open(403,file='norm'//trim(handle)//'.txt', status='replace', action='write')
!open(404,file='dedq'//trim(handle)//'.txt', status='replace', action='write')
if (ndof.eq.1) open(202,file='wf'//trim(handle)//'.txt', status='replace', action='write')
if (ndof.gt.1) then
  open(405,file='edim'//trim(handle)//'.txt', status='replace', action='write')
end if

open(406,file='autocorr'//trim(handle)//'.txt', status='replace', action='write')
open(407,file='ovmax'//trim(handle)//'.txt', status='replace', action='write')

open(201,file='nbasis'//trim(handle)//'.txt', status='replace', action='write')
open(202,file='wfx'//trim(handle)//'.txt', status='replace', action='write')
open(205,file='gbc'//trim(handle)//'.txt', status='replace', action='write')
open(208,file='mom'//trim(handle)//'.txt', status='replace', action='write')
open(206,file='at'//trim(handle)//'.txt', status='replace', action='write')
open(207,file='avg_pos'//trim(handle)//'.txt', status='replace', action='write')

call initialize(nbt, qpas, c, co, sobolfile)
!call reex_traj(nbt, gbc(1:nbt,1:ndof), c(1:nbt), ntrk)

open(500,file='wf0_'//trim(handle)//'.txt')
call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),0d0,500)
close(500)
!c=(0d0,0d0)

open(505,file='gbc0_'//trim(handle)//'.txt')
do i0=1,nbt
call psi(qpas(i0,:,:),nbt,c(1:nbt),qpas(1:nbt,:,:),wf0)
write(505,form1) (qpas(i0,n,1), n=1,ndof), real(wf0), aimag(wf0)
end do
close(505)

open(501,file='wf30_'//trim(handle)//'.txt')
open(502,file='wf60_'//trim(handle)//'.txt')
open(503,file='gbc30_'//trim(handle)//'.txt')
open(504,file='gbc60_'//trim(handle)//'.txt')

ireex=.FALSE.
nrx=0
!1e-3 for double well
rxtols(1)=1e-3
!rxtols(2)=1e-3

do k=1,nbt
  do n=1,ndof
    p(k,n)=pwf(n)
  end do
end do

t=0.0
!idedq=0
write(201,*) 0, t, nbt

call autocorr(nbt, qpas(1:nbt,:,:), c(1:nbt), zcorr, tcorr)
write(406,form1) t, zcorr, tcorr

do iter=1, iter_max

!idedq = idedq+1

!Create overlap and hamiltonian matrices
call overlap(nbt,qpas(1:nbt,:,:),ov(1:nbt,1:nbt))

do i=1,nbt
  do j=1,nbt
    KE=(0d0,0d0)
    V=(0d0,0d0)
    call hamiltonian(qpas(i,:,:),qpas(j,:,:),KE,V)
    H(i,j) = ov(i,j)*(KE+V)
  end do
end do

if (iter.eq.1) then
  call energy_calc(nbt, c(1:nbt), H(1:nbt,1:nbt), einit)
  call normalize(nbt, qpas(1:nbt,:,:), c(1:nbt), ov(1:nbt,1:nbt), norm)
  write(*,*) "INITIAL BASIS NUMBER = ", nbt
  write(*,*) "INITIAL ENERGY = ", einit
  write(*,*) "INITIAL NORM = ", norm
  eold=einit
  c=(0d0,0d0)
end if

call basis_diag(nbt, qpas(1:nbt,:,:), ov(1:nbt,1:nbt), H(1:nbt,1:nbt), co(1:nbt), c(1:nbt), evals(1:nbt))
call energy_calc(nbt, c(1:nbt), H(1:nbt,1:nbt), ediag)

!DYNAMIC DT
!if (mod(iter,10).eq.1) then
!  if (abs(ediag-eold).gt.5e-4) then
!    dt=dt/5d0
!  else if (abs(ediag-eold).lt.1e-4) then
!    dt=dt*2d0
!  end if
!  eold=ediag
!  write(*,*) dt
!end if

!call energy_comp(nbt, qpas(1:nbt,:,:), c(1:nbt), ov(1:nbt,1:nbt), KE, V)
call normalize(nbt, qpas(1:nbt,:,:), c(1:nbt), ov(1:nbt,1:nbt), norm)
!call dedq(nbt,gbc(1:nbt,1:ndof),ntrk,ediag,idedq,es,qsl,qsu,dedql,dedqu)
!call dedq_calc(nbt,gbc,ediag,iter,es,qs,dedq)

write(*,*) 'ITER = ', iter, 'ENERGY = ', ediag, 'NORM = ', norm

write(403,*) t, norm

if (norm.lt.0.99) then
  write(*,*) 'NORM FALLEN BELOW THRESHOLD VALUE 0.99, STOPPING...'
!  open(502,file='wft_'//trim(handle)//'.txt')
!    call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),t,502)
!  close(502)
  exit
end if

!t=t+dt

!ENERGY_DIM STILL NEEDS TO BE UPDATED TO QPAS 10/8/19.
!if (ndof.gt.1) then
!  call energy_dim(nbt, qpas(1:nbt,:,:), c(1:nbt), etest, edim)
!end if

!For E reexpansions:
!if (abs(ediag-einit).gt.abs(0.01*einit)) ireex=.TRUE.
!if (abs(norm-1d0).gt.0.001) ireex=.TRUE.

!For KE reexpansions:
!call energy_comp(nbt, qpas(1:nbt,:,:), c(1:nbt), ov(1:nbt,1:nbt), ket, vt)
!if (abs(real(ket)-real(ke0)).gt.0.005*einit) ireex=.TRUE.

!For periodic reexpansions:
!if (mod(iter,150).eq.0) ireex=.TRUE.

!FOR 2D reexpansions:
!if (mod(iter,100).eq.0) ireex=.TRUE.

!dedq reexpansions:
!do n=1, ndof
!  if (abs(dedql(n)).gt.rxtols(n)) ireex=.TRUE.
!  if (abs(dedqu(n)).gt.rxtols(n)) ireex=.TRUE.
!  if (abs(dedq(1,n)).gt.(rxtols(n))) ireex=.TRUE.
!  if (abs(dedq(nbt,n)).gt.(rxtols(n))) ireex=.TRUE.
!end do
!write(404,*) t, sum(dedq,DIM=1)/nbt, (dedq(1,k), dedq(nbt,k), k=1,ndof)

!Invoke reexpansion...
!if (ireex.and.iter.ne.iter_max) then
!!  write(*,*) 'Reexpanding...'
!  nrx=nrx+1
!!  call reex_diag(nbt, gbc(1:nbt_max,1:ndof), p(1:nbt_max,1:ndof), c(1:nbt_max), co(1:nbt_max), sobolfile)
!  call reex_sobl(nbt, qpas, c, co, sobolfile)
!  write(*,*) 'Reexpanded basis =', nbt
!!  write(*,*) 'NEW ALPHA =', a
!  write(201,*) iter, t, nbt
!!  call reex_traj(nbt, gbc(1:nbt,1:ndof), c(1:nbt), ntrk)
!!  call energy_comp(nbt, gbc(1:nbt,1:ndof), c(1:nbt), ov(1:nbt,1:nbt), ke0, v0)

!!Include call(overlap) and call(hamiltonian) loops here. 1/16/20
!  call basis_diag(nbt, qpas(1:nbt,:,:), ov(1:nbt,1:nbt), H(1:nbt,1:nbt), co(1:nbt), c(1:nbt), evals(1:nbt))
!  call energy_calc(nbt, c(1:nbt), H(1:nbt,1:nbt), ediag)
!  ireex=.FALSE.
!!  idedq=0
!end if  
!End reexpansion

!ichk=mod(iter,100)
!!if (ichk.eq.0) then !periodic output if
!call autocorr(nbt, qpas(1:nbt,:,:), c(1:nbt), zcorr, tcorr)
!call overlap_max(nbt,t,ov(1:nbt,1:nbt))
!call avg_pos(nbt, qpas(1:nbt,:,:), c(1:nbt), ov(1:nbt,1:nbt), xavg)
!!call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),t,205)

!do i=1,nbt
!  write(205,form1) t, (qpas(i,n,1), n=1,ndof), real(c(i)), aimag(c(i))
!!  write(206,form1) t, (qpas(i,n,3), n=1,ndof)
!end do

!write(406,form1) t, zcorr, tcorr
!write(207,form1) t, (xavg(n), n=1,ndof)

!!end if !end periodic output if

!ichk=mod(iter,1000)
!if (ichk.eq.0.and.ndof.eq.1) then
!  call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),t,202) 
!end if

!write(401,form1) t, ediag/norm, evals(1), evals(2)
!if (ndof.gt.1) then
!  write(405,form1) t, etest/norm, (real(edim(n))/norm, n=1,ndof)
!end if
!!write(404,form1) t, (dedql(k), dedqu(k), k=1,ndof)

if (iter.lt.iter_max) then
call prop_diag(nbt, qpas(1:nbt,:,:), p(1:nbt,1:ndof), c(1:nbt), ov_no(1:nbt,1:nbt))
co=(0d0,0d0)
do i=1,nbt
  do j=1, nbt
    co(i)=co(i)+ov_no(i,j)*c(j)
  end do
end do
end if

t=t+dt
!if (t.gt.55) qvar="FRZN"

ichk=mod(iter,nout)
if (ichk.eq.0) then !periodic output if
call autocorr(nbt, qpas(1:nbt,:,:), c(1:nbt), zcorr, tcorr)
call overlap_max(nbt,t,ov(1:nbt,1:nbt))
call avg_pos(nbt, qpas(1:nbt,:,:), c(1:nbt), ov(1:nbt,1:nbt), xavg)
!call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),t,205)

do i=1,nbt
  write(205,form1) t, (qpas(i,n,1), n=1,ndof), real(c(i)), aimag(c(i))
  write(206,form1) t, (qpas(i,n,3), n=1,ndof)
end do

write(406,form1) t, zcorr, tcorr
write(207,form1) t, (xavg(n), n=1,ndof)

end if !end periodic output if

!---TOC OUTPUT
!if (iter.eq.30000) then
!call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),0d0,501)
!do i30=1,nbt
!call psi(qpas(i30,:,:),nbt,c(1:nbt),qpas(1:nbt,:,:),wf30)
!write(503,form1) (qpas(i30,n,1), n=1,ndof), real(wf30), aimag(wf30)
!end do
!else if (iter.eq.60000) then
!call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),0d0,502)
!do i60=1,nbt
!call psi(qpas(i60,:,:),nbt,c(1:nbt),qpas(1:nbt,:,:),wf60)
!write(504,form1) (qpas(i60,n,1), n=1,ndof), real(wf60), aimag(wf60)
!end do
!end if
!---END TOC OUTPUT

!ichk=mod(iter,1000)
!if (ichk.eq.0.and.ndof.eq.1) then
!  call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),t,202)
!end if

ichk=mod(iter,50000)
if (ichk.eq.0) call proj1D(nbt, qpas(1:nbt,:,:), c(1:nbt))

write(401,form1) t, ediag/norm, evals(1), evals(2)
if (ndof.gt.1) then
  write(405,form1) t, etest/norm, (real(edim(n))/norm, n=1,ndof)
end if
!write(404,form1) t, (dedql(k), dedqu(k), k=1,ndof)

end do

write(*,*) "t = ", t
write(*,*) "Number of reexpansions = ", nrx

call psi_update(nbt,qpas(1:nbt,:,:),c(1:nbt),t,202)

close (201)
close (202)
close (205)
close (206)
close (401)
close (402)
close (403)
!close (404)
!close (405)
close (406)
close(407)

close(501)
close(502)
close(503)
close(504)

deallocate(qpas)
deallocate(p)
deallocate(c)
deallocate(co)
deallocate(ov)
deallocate(H)
deallocate(dedql)
deallocate(dedqu)
!deallocate(dedq)
deallocate(es)
deallocate(qsl)
deallocate(qsu)
!deallocate(qs)
deallocate(edim)

END PROGRAM
