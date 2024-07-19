      SUBROUTINE STEPPCR(deltat,mu4,mu4err,n4)
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

integer (kind=4) :: ix,iy,iz,iaux,n4
real    (kind=8) :: deltat,hmean,mu4,mu4err,cnorm


  Call derivnD(2,nn,hx,1,psi,sto1,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3,Icon)

!
!   We compute H·Psi
!
!
!      Predictor:
!

hpsi   = -(sto1+sto2+sto3)*h2o2m4 + pot4*psi

hmean  = sum(psi*hpsi)*dxyz/n4    ! Average over all the mu4's

timec=deltat/(1.0d0+dexp((Abs(hpsi/psi)-Dmumax)/sdmu))+dtmin

Sto4   = -(hpsi-hmean*psi)

!Sto1 = psiold(:,:,:,ioldp(3)) + c4o3*deltat*(2.d0*Sto4-hpsiold(:,:,:,ioldh(1))           &
Sto1 = psiold(:,:,:,ioldp(3)) + c4o3*timec*(2.d0*Sto4-hpsiold(:,:,:,ioldh(1))           &
        + 2.d0*hpsiold(:,:,:,ioldh(2)))
!
!      Modificador:
!
psiold(:,:,:,ioldp(3)) = psi
psi = Sto1 - c112*pc
pc  = Sto1
hpsiold(:,:,:,ioldh(2)) = Sto4
cnorm = sqrt(n4/(sum(psi**2)*dxyz))
psi = Abs(psi)*cnorm
If(psimin.Gt.0.0d0)  psi = Max(psimin,psi )
den = Abs(psi)**2

!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
iaux=ioldh(2)
ioldh(2)=ioldh(1)
ioldh(1)=iaux

!........................................
    call poten()
!........................................

  Call derivnD(2,nn,hx,1,psi,sto1,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3,Icon)

Sto4   = -(sto1+sto2+sto3)*h2o2m4  + pot4*psi
mu4err  = sum(psi*Sto4)*dxyz/n4    ! Average over all the mu4's
Sto4   = -(Sto4-mu4err*psi)

Sto5 = 0.125d0*( 9.d0*psiold(:,:,:,ioldp(3)) - psiold(:,:,:,ioldp(2))   &
!      +3.d0*deltat*(Sto4 + 2.d0*hpsiold(:,:,:,ioldh(1)) - hpsiold(:,:,:,ioldh(2))  ))
      +3.d0*timec*(Sto4 + 2.d0*hpsiold(:,:,:,ioldh(1)) - hpsiold(:,:,:,ioldh(2))  ))
pc = pc - Sto5
!
!     Valor final:
!
psi = Sto5 + c9*pc
errPC = c9*Sum(Abs(pc))
cnorm = sqrt(n4/(sum(psi**2)*dxyz))
psi = Abs(psi)*cnorm
If(psimin.Gt.0.0d0)  psi = Max(psimin,psi )
den = Abs(psi)**2

errPC=errPC/nxyz

!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
      iaux=ioldp(3)
      ioldp(3)=ioldp(2)
      ioldp(2)=ioldp(1)
      ioldp(1)=iaux
mu4err = abs(1.0d0-mu4/hmean)     ! 'Error' in mu4
mu4    = hmean                    ! New chemical potential
return
end
