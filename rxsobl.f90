subroutine reex_sobl(nbt, qpas, c, co, sobolfile)
use constants
use tier1

!integer, parameter :: lwork = 2500

integer, intent(inout) :: nbt
real*8, dimension(1:nbt_max,1:ndof,4), intent(inout) :: qpas
complex*16, dimension(1:nbt_max), intent(inout) :: c
complex*16, dimension(1:nbt_max), intent(out) :: co
character(*), intent(in) :: sobolfile

integer :: nbtn, i, j, k, n, n1, d, nbt0, nmax, lwork, INFO
real*8 :: a_norm, fac, dq, dpb, ds, dsqerr
complex*16 :: wf, kold, vold, kerx, vrx, zn, ztot, wfo, wfn
real*8, dimension(1:ndof) :: xi
real*8, dimension(ndof,2) :: bds
real*8, dimension(1:nbt_max,1:ndof,4) :: qpasn
complex*16, dimension(1:nbt_max,1:nbt) :: ov_no
!real*8, dimension(1:nbt_max,1:ndof) :: r, dp, dr 

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
!Check if the basis function is being placed in an area of appropriate wf
!density. Continue if check passes.
  call psi(xi, nbt, c(1:nbt), qpas(1:nbt,:,:), wf)
  if (abs(wf)**2.gt.tol) then
!If k<nbt_max, place the new basis function.
    if (.not.idist) then
      k=k+1
      if (k.gt.nbt_max) then
        write(*,*) 'REEXPANDED BASIS EXCEEDS NBT_MAX'
        stop
      end if
      do d=1,ndof
        qpasn(k,d,1) = xi(d)
        qpasn(k,d,2) = 0d0
        qpasn(k,d,3) = a0(d)
        qpasn(k,d,4) = 0d0
      end do
!End placement.
!Calculate co = ov_no x c for the k-th new basis function.
      do j=1,nbt
        ztot = (1d0,0d0)
        do d=1,ndof
          fac=2d0*(qpasn(k,d,3)*qpas(j,d,3))**0.25/dsqrt(2d0*qpasn(k,d,3)+2d0*(qpas(j,d,3)))
          dq=qpasn(k,d,1)-qpas(j,d,1)
          dpb=qpasn(k,d,2)-qpas(j,d,2)
          ds=qpasn(k,d,4)-qpas(j,d,4)
          zn=dpb**2+qpasn(k,d,3)*qpas(j,d,3)*dq**2+2d0*im*(qpasn(k,d,3)*(ds-qpas(j,d,2)*dq)+qpas(j,d,3)*(ds-qpasn(k,d,2)*dq))
          ztot = fac*ztot*exp(-zn/(2d0*(qpasn(k,d,3)+qpas(j,d,3))))
        end do
        ov_no(k,j)=ztot
      end do
      do j=1,nbt
        co(k)=co(k)+ov_no(k,j)*c(j)
      end do
!Calculate  ovn = <g_k|g_j> for the k-th new basis function.
      do j=1,k
        ztot=(1d0,0d0)
        do d=1,ndof
          fac=2d0*(qpasn(k,d,3)*qpasn(j,d,3))**0.25/dsqrt(2d0*qpasn(k,d,3)+2d0*(qpasn(j,d,3)))
          dq=qpasn(k,d,1)-qpasn(j,d,1)
          dpb=qpasn(k,d,2)-qpasn(j,d,2)
          ds=qpasn(k,d,4)-qpasn(j,d,4)
          zn=dpb**2+qpasn(k,d,3)*qpasn(j,d,3)*dq**2+2d0*im*(qpasn(k,d,3)*(ds-qpasn(j,d,2)*dq)+qpasn(j,d,3)*(ds-qpasn(k,d,2)*dq))
          ztot = fac*ztot*exp(-zn/(2d0*(qpasn(k,d,3)+qpasn(j,d,3))))
        end do
        ovn(k,j) = ztot
      end do
      if (k.gt.1) then
        do j=1,k-1
          ovn(j,k) = conjg(ovn(k,j))
        end do
      end if

!If the number of new basis functions exceeds some minimum (nbt0), check to see
!if we
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
!        call normalize(k, qpasn(1:k,:,:), ctemp(1:k), a_norm)
!        if (abs(a_norm-1d0).lt.0.005) then
!If norm check passes, check the kinetic energy matching...
!          call energy_comp(k, qpasn(1:k,:,:), ctemp(1:k), kerx, vrx)
!If both norm and KE matching pass, we're done placing new basis functions.
!          if (abs(real(kerx)-real(kold)).lt.abs(1e-5*real(kold))) exit
!        end if

!||wfn-wfo||^2 check.
        dsqerr=0d0
        do n1=1,npx
          xi(1)=xdim(1,1)+(n1-1)*(xdim(1,2)-xdim(1,1))/npx
          call psi(xi, nbt, c(1:nbt), qpas(1:nbt,1:ndof,1:4), wfo)
          call psi(xi, k, ctemp(1:k), qpasn(1:k,1:ndof,1:4), wfn)
          dsqerr=dsqerr+abs(wfn-wfo)**2
        end do
        if (dsqerr.lt.1e-5) exit

      end if
    end if
  end if
  if (n.eq.nmax) then
    write(*,*) 'SOBOL REEXPANSION END ENCOUNTERED PREMATURELY'
    stop
  end if
end do
close(541)
rewind(541)

!Record the total number of basis functions and their coefficients.
nbtn=k
c=ctemp

nbt=nbtn
do n=1,4
  do d=1,ndof
    do i=1,nbt_max
      if (i.le.nbt) then
        qpas(i,d,n)=qpasn(i,d,n)
      else
        qpas(i,d,n)=0d0
      end if
    end do
  end do
end do

!p=0d0
!call momentum(nbt, qpas(1:nbt,:,:), c(1:nbt), p(1:nbt,:), r(1:nbt,:), dp(1:nbt,:), dr(1:nbt,:))

return
end subroutine


