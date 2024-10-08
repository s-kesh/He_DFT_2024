SUBROUTINE STEPRK(deltat)
!
!      Runge-Kutta-Gill method
!       (Ralston & Wilf Vol I, pag. 117)
!
use Para_derivnD
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use impur  !
use quantal_imp
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc rutines)

implicit none

real (kind=8) :: arun(4),brun(4),crun(4)

integer (kind=4) :: ix,iy,iz,jrun,i_imp
real    (kind=8) :: deltat,h2o2mimp
real    (kind=8) :: aux2r(3)
complex (kind=8) :: auxc(6)
complex (kind=8) :: aux1c,aux2c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)

!include 'interface_derivnD.include'  ! Per fer servir les derivades generiques

arun(1)=0.5d0
arun(2)=1.0d0-1.d0/dsqrt(2.d0)
arun(3)=1.0d0+1.d0/dsqrt(2.d0)
arun(4)=1.d0/6.d0

brun(1)=2.0d0
brun(2)=1.0d0
brun(3)=1.0d0
brun(4)=2.0d0

 crun(1)=0.5d0
 crun(2)=1.0d0-1.d0/dsqrt(2.d0)
 crun(3)=1.0d0+1.d0/dsqrt(2.d0)
 crun(4)=0.5d0


!........................................
!.. Laplacian of H^{iorder-1}�Psi (0) ...
!........................................
!

do jrun=1,4

!Write(6,'("Steprk: Abans de calcular les derivades....")')
!Call Flush(6)

  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

!Write(6,'("Steprk: Despres de calcular les derivades....")')
!Call Flush(6)

!
!   We compute H�Psi
!
  Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 - (pot4-mu4)*psi)
!  Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 - pot4*psi) - ci*uimp*psi
  Sto1c = arun(jrun)*(Sto4c - brun(jrun)*q)
  q = q + 3.*Sto1c - crun(jrun)*Sto4c
  if(jrun.eq.1)then
    hpsiold(:,:,:,2) = hpsiold(:,:,:,1)
    hpsiold(:,:,:,1) = Sto4c
    psiold(:,:,:,3) = psiold(:,:,:,2)
    psiold(:,:,:,2) = psiold(:,:,:,1)
    psiold(:,:,:,1) = psi
  endif
  psi = psi + deltat*Sto1c
  den = Abs(psi)**2
!
! Impurity evolution if it is necessary
!
  denx = 0.d0
  Do i_imp=1,N_Imp
    h2o2mimp = h2o2m4*4.002603d0/m_imp_u(i_imp)
    Call derivnD(2,nn,hx,1,psix(:,:,:,i_imp),sto1c,Icon)
    Call derivnD(2,nn,hy,2,psix(:,:,:,i_imp),sto2c,Icon)
    Call derivnD(2,nn,hz,3,psix(:,:,:,i_imp),sto3c,Icon)

!
!   We compute H�Psix
!
    Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2mimp - upotx(:,:,:,i_imp)*psix(:,:,:,i_imp))
    Sto1c = arun(jrun)*(Sto4c - brun(jrun)*qx(:,:,:,i_imp))
    qx(:,:,:,i_imp) = qx(:,:,:,i_imp) + 3.*Sto1c - crun(jrun)*Sto4c
    if(jrun.eq.1)then
      hpsixold(:,:,:,2,i_imp) = hpsixold(:,:,:,1,i_imp)
      hpsixold(:,:,:,1,i_imp) = Sto4c
      psixold(:,:,:,3,i_imp) = psixold(:,:,:,2,i_imp)
      psixold(:,:,:,2,i_imp) = psixold(:,:,:,1,i_imp)
      psixold(:,:,:,1,i_imp) = psix(:,:,:,i_imp)
    endif
    psix(:,:,:,i_imp) = psix(:,:,:,i_imp) + deltat*Sto1c
    denx = denx + Conjg(psix(:,:,:,i_imp))*psix(:,:,:,i_imp)
  EndDo
  if(jrun.le.3)then
    call potenimp()
    call poten()
  endif
EndDo

do ix=1,3
  ioldp(ix)=ix
  ioldpx(ix,:)=ix
enddo

do ix=1,2
  ioldh(ix)=ix
  ioldhx(ix,:)=ix
enddo

return
end
