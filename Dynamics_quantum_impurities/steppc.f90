SUBROUTINE STEPPC(deltat,errHe,errpsix)
!
!      Predictor-Modifier-Corrector method
!       (Ralston & Wilf Vol I, pag. 99)
!
use Para_derivnD
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use impur  !
use Quantal_imp
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc)

implicit none

real (kind=8) :: c112=112.d0/121.d0
real (kind=8) :: c9  =9.d0/121.d0
real (kind=8) :: c1o3=1.d0/3.d0
real (kind=8) :: c4o3=4.d0/3.d0
real (kind=8) :: c5o3=5.d0/3.d0

integer (kind=4) :: ix,iy,iz,iaux,i_imp
real    (kind=8) :: deltat,h2o2mimp
real    (kind=8) :: errHe, errpsix
!real    (kind=8) :: auxr(3)
complex (kind=8) :: auxc(6)
complex (kind=8) :: aux1c,aux2c,aux3c,aux4c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)


  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

!
!   We compute H·Psi
!
!
!      Predictor:
!
Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 - (pot4-mu4)*psi)
!Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 - pot4*psi) - ci*uimp*psi
Sto1c = psiold(:,:,:,ioldp(3)) + c4o3*deltat*(2.d0*Sto4c-hpsiold(:,:,:,ioldh(1))           &
        + 2.d0*hpsiold(:,:,:,ioldh(2)))
!
!      Modificador:
!
psiold(:,:,:,ioldp(3)) = psi
psi = Sto1c - c112*pc
pc  = Sto1c
hpsiold(:,:,:,ioldh(2)) = Sto4c
den = Abs(psi)**2
!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
iaux=ioldh(2)
ioldh(2)=ioldh(1)
ioldh(1)=iaux
denx = 0.d0
Do i_imp=1, N_imp
  h2o2mimp = h2o2m4*4.002603d0/m_imp_u(i_imp)
  Call derivnD(2,nn,hx,1,psix(:,:,:,i_imp),sto1c,Icon)
  Call derivnD(2,nn,hy,2,psix(:,:,:,i_imp),sto2c,Icon)
  Call derivnD(2,nn,hz,3,psix(:,:,:,i_imp),sto3c,Icon)

!
!   We compute H·Psix
!
!
!      Predictor:
!
  Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2mimp - upotx(:,:,:,i_imp)*psix(:,:,:,i_imp))
  Sto1c = psixold(:,:,:,ioldpx(3,i_imp),i_imp)                                        & 
        + c4o3*deltat*(2.d0*Sto4c-hpsixold(:,:,:,ioldhx(1,i_imp),i_imp)               &
        + 2.d0*hpsixold(:,:,:,ioldhx(2,i_imp),i_imp))
!
!      Modificador:
!
  psixold(:,:,:,ioldpx(3,i_imp),i_imp) = psix(:,:,:,i_imp)
  psix(:,:,:,i_imp) = Sto1c - c112*pcx(:,:,:,i_imp)
  pcx(:,:,:,i_imp)  = Sto1c
  hpsixold(:,:,:,ioldhx(2,i_imp),I_imp) = Sto4c
  denx = denx + Conjg(psix(:,:,:,i_imp))*Psix(:,:,:,i_imp)
!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
  iaux=ioldhx(2,i_imp)
  ioldhx(2,i_imp)=ioldhx(1,i_imp)
  ioldhx(1,i_imp)=iaux
EndDo
    call potenimp()
    call poten()
!........................................


  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4  - (pot4-mu4)*psi) 
!Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4  - pot4*psi) - ci*uimp*psi
Sto5c = 0.125d0*( 9.d0*psiold(:,:,:,ioldp(3)) - psiold(:,:,:,ioldp(2))   &
      +3.d0*deltat*(Sto4c + 2.d0*hpsiold(:,:,:,ioldh(1)) - hpsiold(:,:,:,ioldh(2))  ))
pc = pc - Sto5c
!
!     Valor final:
!
psi = Sto5c + c9*pc
errHe = c9*Sum(Abs(pc))
den = Abs(psi)**2
errPsix = 0.d0
denx = 0.d0
Do i_imp=1, N_imp
  h2o2mimp = h2o2m4*4.002603d0/m_imp_u(i_imp)
!
!   We compute H·Psix
!
  Call derivnD(2,nn,hx,1,psix(:,:,:,i_imp),sto1c,Icon)
  Call derivnD(2,nn,hy,2,psix(:,:,:,i_imp),sto2c,Icon)
  Call derivnD(2,nn,hz,3,psix(:,:,:,i_imp),sto3c,Icon)

  Sto4c = timec*((sto1c+sto2c+sto3c)*h2o2mimp  - upotx(:,:,:,i_imp)*psix(:,:,:,i_imp))
  Sto5c = 0.125d0*( 9.d0*psixold(:,:,:,ioldpx(3,i_imp),i_imp) - psixold(:,:,:,ioldpx(2,i_imp),i_imp)  &
        + 3.d0*deltat*(Sto4c + 2.d0*hpsixold(:,:,:,ioldhx(1,i_imp),i_imp)                             & 
                             -      hpsixold(:,:,:,ioldhx(2,i_imp),i_imp)))
  pcx(:,:,:,i_imp) = pcx(:,:,:,i_imp) - Sto5c
!
!     Valor final:
!
  psix(:,:,:,i_imp) = Sto5c + c9*pcx(:,:,:,i_imp)
  errPsix = errPsix + c9*Sum(Abs(pcx(:,:,:,i_imp)))
  denx = denx + Conjg(psix(:,:,:,i_imp))*psix(:,:,:,i_imp)
EndDo
errPsix=errPsix/nxyz
errHe=errHe/nxyz

!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
iaux=ioldp(3)
ioldp(3)=ioldp(2)
ioldp(2)=ioldp(1)
ioldp(1)=iaux

Do i_imp=1, N_imp
  iaux=ioldpx(3,i_imp)
  ioldpx(3,i_imp)=ioldpx(2,i_imp)
  ioldpx(2,i_imp)=ioldpx(1,i_imp)
  ioldpx(1,i_imp)=iaux
EndDo
return
end
