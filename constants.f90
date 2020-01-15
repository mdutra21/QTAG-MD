module constants

implicit none

real*8, public, parameter :: pi = 4.d0*atan(1.0)
complex*16, public, parameter :: im = (0d0,1d0)
character(len=16), public, parameter :: form1 = "(2000(e14.7,1x))"

integer, public :: npx, iter_max, nbt_max, ndof, npar, nreex, skip, ihandle
real*8, public :: beta, dt, vcp, tol
character(len=4), public :: init, reex, pcalc
character(len=4), public :: qvar, pvar, avar, svar
character(len=5), public :: potflag
character(len=16), public :: infile, potfile

integer, allocatable, public :: nb(:)
real*8, allocatable, public :: a0(:), qwf(:), pwf(:), awf(:), ms(:)
real*8, allocatable, public :: xdim(:,:)
real*8, allocatable, public :: omegas(:)
character(len=5), allocatable, public :: pottype(:)
contains

subroutine read_params()

integer :: n

infile = 'input_params.txt'

open (1, file=infile, status='old', action='read')
  read(1,*) ihandle
  read(1,*) init, reex, pcalc
  read(1,*) qvar, pvar, avar, svar
  read(1,*) ndof, nbt_max, skip
  read(1,*) npx, beta, tol
  read(1,*) dt, iter_max
  read(1,*) vcp
  read(1,*)
  allocate(nb(ndof))
  allocate(a0(ndof))
  allocate(qwf(ndof))
  allocate(pwf(ndof))
  allocate(awf(ndof))
  allocate(xdim(ndof,2))
  allocate(ms(ndof))
  allocate(pottype(ndof))
  do n=1, ndof
    read(1,*) xdim(n,1), xdim(n,2)
  end do
  read(1,*)
  do n=1, ndof
    read(1,*) nb(n), a0(n)
  end do
  read(1,*)
  do n=1, ndof
    read(1,*) qwf(n), pwf(n), awf(n)
  end do
  read(1,*)
  do n=1,ndof
    read(1,*) ms(n)
  end do
close (1)
rewind(1)

write(*,*) 'MASSES = ', (ms(n), n=1,ndof)
do n=1, ndof
  a0(n)=a0(n)*awf(n)
end do

end subroutine

end module
