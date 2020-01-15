!------------------------------------------------------------------------------
!Module containing subroutines needed during trajectory calculations and
!updates. Included subroutines are as follows:
!overlap
!overlap_max
!hamiltonian
!psi
!energy_calc
!energy_comp
!avg_pos
!normalize
!autocorr
!------------------------------------------------------------------------------

module tier1

use constants
implicit none

contains

!------------------------------------------------------------------------------
!Compute the overlap integral <g_i|g_j> for all i,j pairs, and store the
!result in the matrix ov_mat. Note that |g_i> and |g_j> must be
!Gaussian.
!------------------------------------------------------------------------------
subroutine overlap(nbt,qpas,ov_mat)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(in) :: qpas
complex*16, dimension(1:nbt,1:nbt), intent(out) :: ov_mat

integer :: i, j, n
real*8 :: fac, dq, dp, ds
complex*16 :: zn, ztot

do i=1,nbt
  do j=1,nbt
    ztot=(1d0,0d0)
    do n=1,ndof
      fac=2d0*(qpas(i,n,3)*qpas(j,n,3))**0.25/dsqrt(2d0*qpas(i,n,3)+2d0*(qpas(j,n,3)))
      dq=qpas(i,n,1)-qpas(j,n,1)
      dp=qpas(i,n,2)-qpas(j,n,2)
      ds=qpas(i,n,4)-qpas(j,n,4)
      zn=dp**2+qpas(i,n,3)*qpas(j,n,3)*dq**2+2d0*im*(qpas(i,n,3)*(ds-qpas(j,n,2)*dq)+qpas(j,n,3)*(ds-qpas(i,n,2)*dq))
      ztot = fac*ztot*exp(-zn/(2d0*(qpas(i,n,3)+qpas(j,n,3))))
    end do
    ov_mat(i,j) = ztot
  end do
end do

!do i=2,nbt
!  do j=1,i-1
!    ov_mat(j,i) = conjg(ov_mat(i,j))
!  end do
!end do

return
end subroutine
!------------------------------------------------------------------------------
subroutine overlap_max(nbt,ov)

integer, intent(in) :: nbt
complex*16, dimension(1:nbt,1:nbt), intent(in) :: ov

integer :: i, j
real*8, dimension(1:nbt) :: ovmax
real*8, dimension(1:nbt,1:nbt) :: ov1

ov1=(0d0)
do j=1,nbt
  do i=1,j
    if (i.ne.j) then
      ov1(i,j)=dsqrt(real(conjg(ov(i,j))*ov(i,j)))
    end if
  end do
end do

ovmax=maxval(ov1,dim=1)
write(407,*) (ovmax(i), i=1,nbt)

return
end subroutine
!------------------------------------------------------------------------------
subroutine overlap_T(nbt,qpas,qpasn,ov_no)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt,1:nbt), intent(out) :: ov_no

integer :: i, j, k
real*8 :: fac, dq, dp, ds
complex*16 :: zn, ztot

real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpasn

do j=1,nbt
  do i=1,nbt
    ztot = (1d0,0d0)
    do k=1,ndof
      fac=2d0*(qpasn(i,k,3)*qpas(j,k,3))**0.25/dsqrt(2d0*qpasn(i,k,3)+2d0*(qpas(j,k,3)))
      dq=qpasn(i,k,1)-qpas(j,k,1)
      dp=qpasn(i,k,2)-qpas(j,k,2)
      ds=qpasn(i,k,4)-qpas(j,k,4)
      zn=dp**2+qpasn(i,k,3)*qpas(j,k,3)*dq**2+2d0*im*(qpasn(i,k,3)*(ds-qpas(j,k,2)*dq)+qpas(j,k,3)*(ds-qpasn(i,k,2)*dq))
      ztot = fac*ztot*exp(-zn/(2d0*(qpasn(i,k,3)+qpas(j,k,3))))
    end do
    ov_no(i,j)=ztot
  end do
end do

return
end subroutine
!------------------------------------------------------------------------------
subroutine hamiltonian(qpasi,qpasj,KE,V)

real*8, dimension(1,1:ndof,1:4), intent(in) :: qpasi, qpasj
complex*16, intent(out) :: KE, V

integer :: n
real*8 :: dq, vv0, vv1, vv2, vx, dvx, d2vx, vcpl, dvcpl, d2vcpl
real*8 :: qi, qj, pbi, pbj, ai, aj, as, ap, dp, aqs, aqp
real*8 :: xid1, xid2, xjd1, xjd2, pbid1, pbid2, vc1, vc2
real*8 :: pbjd1, pbjd2, aid1, aid2, ajd1, ajd2, qq, pp
complex*16 :: vi, vj, z, z1, z2, zz

vi=(0d0,0d0)
vj=(0d0,0d0)
KE=(0d0,0d0)
V=(0d0,0d0)

do n=1,ndof
  qi=qpasi(1,n,1)
  qj=qpasj(1,n,1)
  pbi=qpasi(1,n,2)
  pbj=qpasj(1,n,2)
  ai=qpasi(1,n,3)
  aj=qpasj(1,n,3)

  dq=qi-qj
  as=ai+aj
  ap=ai*aj
  dp=ai*pbj+aj*pbi
  aqp=(ai*qi+aj*qj)

  KE=KE+(1d0/(2d0*ms(n)*as**2))*(2d0*im*ap*dq*dp-(ap*dq)**2+dp**2+ap*as)

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


!Exact integration for quartic potential with QTAG bases.
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
!End QTAG.

!LHA for any potential
!  z=(ai*qi+aj*qj+im*(pbj-pbi))/(ai+aj)

!  vv0=vx(qj,n)-dvx(qj,n)*qj+d2vx(qj,n)/2d0*qj**2
!  vv1=-d2vx(qj,n)*qj+dvx(qj,n)
!  vv2=d2vx(qj,n)/2d0
!  vj=vj+vv0+vv1*z+vv2*(z**2+1d0/(ai+aj))

!  vv0=vx(qi,n)-dvx(qi,n)*qi+d2vx(qi,n)/2d0*qi**2
!  vv1=-d2vx(qi,n)*qi+dvx(qi,n)
!  vv2=d2vx(qi,n)/2d0
!  vi=vi+vv0+vv1*z+vv2*(z**2+1d0/(ai+aj))
!END LHA
end do

if (ndof.gt.1) then
  do n=1,ndof-1
!BATISTA MODEL 1
    xid1=qpasi(1,n,1)
    xjd1=qpasj(1,n,1)
    pbid1=qpasi(1,n,2)
    pbjd1=qpasj(1,n,2)
    aid1=qpasi(1,n,3)
    ajd1=qpasj(1,n,3)
!END MODEL 1
!BATISTA MODEL 2
!    xid1=qpasi(1,1,1)
!    xjd1=qpasj(1,1,1)
!    pbid1=qpasi(1,1,2)
!    pbjd1=qpasj(1,1,2)
!    aid1=qpasi(1,1,3)
!    ajd1=qpasj(1,1,3)
!END MODEL 2
    xid2=qpasi(1,n+1,1)
    xjd2=qpasj(1,n+1,1)
    pbid2=qpasi(1,n+1,2)
    pbjd2=qpasj(1,n+1,2)
    aid2=qpasi(1,n+1,3)
    ajd2=qpasj(1,n+1,3)

    z1=(aid1*xid1+ajd1*xjd1-im*(pbid1-pbjd1))/(aid1+ajd1)
    z2=(aid2*xid2+ajd2*xjd2-im*(pbid2-pbjd2))/(aid2+ajd2)
    zz=z1*z2

!FOR EXACT INTEGRAL
    V=V+vcp*zz

!FOR LHA OF COUPLED POTENTIAL
!    vv0=vcpl(xid1,xid2)
!    vv1=-d2vcpl(xid1,xid2)*z1*xid2-d2vcpl(xid1,xid2)*z2*xid1
!    vv2=d2vcpl(xid1,xid2)*zz

!    vi=vi+2d0*vv0+vv1+vv2

!    vv0=vcpl(xjd1,xjd2)
!    vv1=-d2vcpl(xid1,xid2)*z1*xjd2-d2vcpl(xid1,xid2)*z2*xjd1
!    vv2=d2vcpl(xjd1,xjd2)*zz

!    vj=vj+2d0*vv0+vv1+vv2
!END LHA
  end do
end if

!FOR LHA POTENTIAL
!V=(vi+vj)/2d0

return
end subroutine
!------------------------------------------------------------------------------
subroutine psi(xi, nbt, c, qpas, wf)

use constants

integer, intent(in) :: nbt
real*8, dimension(1:ndof), intent(in) :: xi
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
complex*16, intent(out) :: wf

integer :: k, n
complex*16 :: prod

wf = (0d0, 0d0)
do k=1,nbt
  prod=1.0d0
  do n =1,ndof
    prod = prod*(qpas(k,n,3)/pi)**0.25*exp(-qpas(k,n,3)/2.0*(xi(n)-qpas(k,n,1))**2 &
           +im*qpas(k,n,2)*(xi(n)-qpas(k,n,1))+im*qpas(k,n,4))
  end do
  wf = wf + c(k)*prod
end do

return
end subroutine
!------------------------------------------------------------------------------
subroutine energy_calc(nbt, qpas, c, etot)

use constants

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,1:4), intent(inout) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
real*8, intent(out) :: etot

integer :: i, j, n
real*8 :: fac, dq, dp, ds
complex*16 :: zn, M, KE, KE1, V, V1

complex*16, dimension(1:nbt,1:nbt) :: H

do i=1,nbt
  do j=1,nbt
    M=(1d0,0d0)
    KE=(0d0,0d0)
    V=(0d0,0d0)
    do n=1,ndof
      fac=2d0*(qpas(i,n,3)*qpas(j,n,3))**0.25/dsqrt(2d0*qpas(i,n,3)+2d0*(qpas(j,n,3)))
      dq=qpas(i,n,1)-qpas(j,n,1)
      dp=qpas(i,n,2)-qpas(j,n,2)
      ds=qpas(i,n,4)-qpas(j,n,4)
      zn=dp**2+qpas(i,n,3)*qpas(j,n,3)*dq**2+2d0*im*(qpas(i,n,3)*(ds-qpas(j,n,2)*dq)+qpas(j,n,3)*(ds-qpas(i,n,2)*dq))
      M = fac*M*exp(-zn/(2d0*(qpas(i,n,3)+qpas(j,n,3))))
!      KE=KE+KE1
!      V=V+V1
    end do
    call hamiltonian(qpas(i,:,:),qpas(j,:,:),KE,V)
    H(i,j) = M*(KE+V)
  end do
end do

etot = real(sum(matmul(reshape(conjg(c),(/1,nbt/)),matmul(H,c))))

return
end subroutine
!----------------------------------------------------------------------------
subroutine energy_comp(nbt, qpas, c, KE, V)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
complex*16, intent(out) :: KE, V

integer :: i, j, n
real*8 :: dq, dp, ds, fac
complex*16, dimension(1:nbt,1:nbt) :: HKE, HV
complex*16 :: zn, M, KE1, V1, H1

do i=1,nbt
  do j=1,nbt
    M=(1d0,0d0)
    KE=(0d0,0d0)
    V=(0d0,0d0)
    do n=1,ndof
      fac=2d0*(qpas(i,n,3)*qpas(j,n,3))**0.25/dsqrt(2d0*qpas(i,n,3)+2d0*(qpas(j,n,3)))
      dq=qpas(i,n,1)-qpas(j,n,1)
      dp=qpas(i,n,2)-qpas(j,n,2)
      ds=qpas(i,n,4)-qpas(j,n,4)
      zn=dp**2+qpas(i,n,3)*qpas(j,n,3)*dq**2+2d0*im*(qpas(i,n,3)*(ds-qpas(j,n,2)*dq)+qpas(j,n,3)*(ds-qpas(i,n,2)*dq))
      M = fac*M*exp(-zn/(2d0*(qpas(i,n,3)+qpas(j,n,3))))

!      call hamiltonian(qpas(i,n,:),qpas(j,n,:),KE1,V1)
!      KE=KE+KE1
!      V=V+V1
    end do
    call hamiltonian(qpas(i,:,:),qpas(j,:,:),KE,V)
    HKE(i,j) = M*KE
    HV(i,j) = M*V
  end do
end do

KE=(0d0,0d0)
V=(0d0,0d0)
KE=sum(matmul(reshape(conjg(c),(/1,nbt/)),matmul(HKE,c)))
V=sum(matmul(reshape(conjg(c),(/1,nbt/)),matmul(HV,c)))

!KE = (0d0,0d0)
!V = (0d0,0d0)
!do i = 1, nbt
!  do j = 1, nbt
!    KE = KE+conjg(c(j))*HKE(i,j)*c(i)
!    V = V+conjg(c(j))*HV(i,j)*c(i)
!  end do
!end do

return
end subroutine
!----------------------------------------------------------------------------
!subroutine energy_dim(nbt, qpas, c, etot, edim)

!integer, intent(in) :: nbt
!real*8, dimension(1:nbt,1:ndof,1:4), intent(inout) :: qpas
!complex*16, dimension(1:nbt), intent(in) :: c
!complex*16, intent(out) :: etot
!complex*16, dimension(1:ndof), intent(out) :: edim

!integer :: i, j, n, n1
!real*8 :: dx, dx2
!real*8, dimension(1:ndof) :: xc, x0
!complex*16, dimension(1:nbt,1:nbt) :: H
!complex*16 :: M, KE, V, pV

!if (potflag.eq."MORSE") then
!  do n=1,ndof
!    x0(n)=fit(n,3)
!  end do
!else
!  do n=1,ndof
!    x0(n)=0d0
!  end if
!end if

!edim =(0d0,0d0)
!etot = (0d0,0d0)

!do n=1,ndof
!do i=1,nbt
!  do j=1,nbt
!    pV = (0d0,0d0)
!    M = (1d0,0d0)
!    dx = gbc(j,n)-gbc(i,n)
!    do n1=1,ndof
!      dx2 = gbc(j,n1)-gbc(i,n1)
!      xc(n1) = (gbc(j,n1)+gbc(i,n1))/2.0d0
!      M = M*exp((-a(n1)*dx2**2)/4d0)
!    end do
!    KE = a(n)/(4d0*ms(n))*(1d0-a(n)/2d0*dx**2)
!    call pot_select(n, xc, V)
!    do n1=1,ndof
!      if(n.ne.n1) then
!        pV=pV+1d0/2d0*(xc(n1)-x0(n1))*(xc(n)-x0(n))
!      end if
!    end do
!    H(i,j) = M*(KE+V+cp*pV)
!  end do
!end do


!do i = 1, nbt
!  do j = 1, nbt
!    edim(n) = edim(n) + conjg(c(j))*H(i,j)*c(i)
!  end do
!end do
!etot = etot+edim(n)
!end do

!return
!end subroutine
!----------------------------------------------------------------------------
subroutine avg_pos(nbt, qpas, c, xavg)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
real*8, dimension(1:ndof), intent(out) :: xavg

integer :: i, j, n
complex*16 :: mat
complex*16, dimension(1:nbt,1:nbt) :: mat0

call overlap(nbt,qpas,mat0)

xavg=0d0
do n=1,ndof
  do i=1,nbt
    do j=1,nbt
      mat=mat0(i,j)*(qpas(i,n,1)+qpas(j,n,1))/2d0
      xavg(n)=xavg(n)+real(conjg(c(i))*mat*c(j))
    end do
  end do
end do

return
end subroutine
!----------------------------------------------------------------------------
subroutine normalize(nbt, qpas, c, nrm)

integer, intent(in) :: nbt
real*8,intent(in),dimension(1:nbt,1:ndof,4) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
real*8, intent(out) :: nrm

complex*16 :: z
complex*16, dimension(1:nbt,1:nbt) :: mat

call overlap(nbt, qpas, mat)
z=sum(matmul(reshape(conjg(c),(/1,nbt/)),matmul(mat,c)))
nrm = real(z)

return
end subroutine

!------------------------------------------------------------------------------
subroutine autocorr(nbt, qpas, c, zcorr, tcorr)

integer, intent(in) :: nbt
real*8, dimension(1:nbt,1:ndof,4), intent(in) :: qpas
complex*16, dimension(1:nbt), intent(in) :: c
complex*16, intent(out) :: zcorr, tcorr

integer :: i, j
real*8 :: dqz, dqt, awf1, asum, anorm, z0, z1, z2, t0, t1, t2
complex*16, dimension(1:nbt) :: bz, bt

zcorr = (0d0,0d0)
tcorr = (0d0,0d0)
!call overlap(nbt, qpas, mat)

!do i=1, nbt
!  do j=1, nbt
!    zcorr = zcorr + c(i)*mat(i,j)*c(j)
!  end do
!end do

do j=1,nbt
  z0=0d0
  z1=0d0
  z2=0d0
  t0=0d0
  t1=0d0
  t2=0d0
  anorm=1d0
  do i=1,ndof
    dqz=qpas(j,i,1)-qwf(i)
    dqt=qpas(j,i,1)+qwf(i)
    awf1=2d0*awf(i)
    asum=qpas(j,i,3)+awf1
    anorm=anorm*(dsqrt(2d0)*(qpas(j,i,3)*awf1)**0.25/dsqrt(asum))
    z0=z0-2d0*qpas(j,i,4)*asum
    z1=z1+2d0*dqz*(qpas(j,i,3)*pwf(i)+awf1*qpas(j,i,2))
    z2=z2+awf1*qpas(j,i,3)*dqz**2+(qpas(j,i,2)-pwf(i))**2
    t0=t0-2d0*qpas(j,i,4)*asum
    t1=t1+2d0*dqt*(qpas(j,i,3)*pwf(i)+awf1*qpas(j,i,2))
    t2=t2+awf1*qpas(j,i,3)*dqt**2+(qpas(j,i,2)-pwf(i))**2
  end do
  bz(j)=c(j)*anorm*exp(-(z2+im*(z0+z1))/(2d0*asum))
  bt(j)=c(j)*anorm*exp(-(t2+im*(t0+t1))/(2d0*asum))
end do
zcorr=sum(bz)
tcorr=sum(bt)

return
end subroutine
!-----------------------------------------------------------------------------
subroutine dedq(nbt,gbc,ntrk,e,idedq,es,qsl,qsu,dedql,dedqu)
integer, intent(in) :: nbt, idedq
integer, dimension(1:ndof,2), intent(in) :: ntrk
real*8, dimension(1:nbt,1:ndof), intent(in) :: gbc
real*8, dimension(1:4), intent(inout) :: es
real*8, dimension(1:ndof,4), intent(inout) :: qsl, qsu
!real*8, dimension(1:nbt_max,1:ndof,1:4), intent(inout) :: qs
complex*16, intent(in) :: e
real*8, dimension(1:ndof), intent(out) :: dedql, dedqu
!real*8, dimension(1:nbt_max,1:ndof), intent(out) :: dedq

integer :: n, n1, k
real*8 :: dedt, dqldt, dqudt
!real*8 :: dedt, dqdt

!if (idedq.gt.1) then
!  do k=4,2,-1
!    es(k)=es(k-1)
!    do n1=1,ndof
!      do n=1,nbt
!        qs(n,n1,k)=qs(n,n1,k-1)
!      end do
!    end do
!  end do
!  es(1)=real(e)
!  do n1=1,ndof
!    do n=1,nbt
!      qs(n,n1,1)=gbc(n,n1)
!    end do
!  end do
!  if (idedq.ge.4) then
!    dedt=(11.0/6.0*es(1)-3.0*es(2)+3.0/2.0*es(3)-1.0/3.0*es(4))/dt
!    do n1=1,ndof
!      do n=1,nbt
!        dqdt=(11.0/6.0*qs(n,n1,1)-3.0*qs(n,n1,2)+3.0/2.0*qs(n,n1,3)-1.0/3.0*qs(n,n1,4))/dt
!        dedq(n,n1)=dedt/dqdt
!      end do
!    end do
!  else
!    dedq=0d0
!  end if
!else
!  es=0d0
!  es(1)=e
!  qs=0d0
!  do n1=1,ndof
!    do n=1,nbt
!      qs(n,n1,1)=gbc(n,n1)
!    end do
!  end do
!  dedq=0d0
!end if

!------------------------------------------------------------
if (idedq.gt.1) then
  do k=4,2,-1
    es(k)=es(k-1)
    do n=1,ndof
      qsl(n,k)=qsl(n,k-1)
      qsu(n,k)=qsu(n,k-1)
    end do
  end do
  es(1)=real(e)
  do n=1,ndof
    qsl(n,1)=gbc(ntrk(n,1),n)
  end do
  do n=1,ndof
    qsu(n,1)=gbc(ntrk(n,2),n)
  end do
  if (idedq.ge.4) then
    dedt=(11.0/6.0*es(1)-3.0*es(2)+3.0/2.0*es(3)-1.0/3.0*es(4))/dt
    do n=1,ndof
      dqldt=(11.0/6.0*qsl(n,1)-3.0*qsl(n,2)+3.0/2.0*qsl(n,3)-1.0/3.0*qsl(n,4))/dt
      dqudt=(11.0/6.0*qsu(n,1)-3.0*qsu(n,2)+3.0/2.0*qsu(n,3)-1.0/3.0*qsu(n,4))/dt 
      dedql(n)=dedt/dqldt
      dedqu(n)=dedt/dqudt
    end do
  else
    dedql=0d0
    dedqu=0d0
  end if
else
  es=0d0
  es(1)=real(e)
  qsl=0d0
  do n=1,ndof
    qsl(n,1)=gbc(ntrk(n,1),n)
  end do
  qsu=0d0
  do n=1,ndof
    qsu(n,1)=gbc(ntrk(n,2),n)
  end do
  dedql=0d0
  dedqu=0d0
end if

return
end subroutine

end module
