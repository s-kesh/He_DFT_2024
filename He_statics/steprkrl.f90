      SUBROUTINE STEPRKRL(deltat,mu4,mu4err,n4)
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
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc rutines)

implicit none

real (kind=8) :: arun(4),brun(4),crun(4)

integer (kind=4) :: ix,iy,iz,jrun,n4
real    (kind=8) :: deltat,mu4,mu4err
real    (kind=8) :: hmean,cnorm

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



do jrun=1,4

  Call derivnD(2,nn,hx,1,psil,sto1,Icon_local)
  Call derivnD(2,nn,hy,2,psil,sto2,Icon_local)
  Call derivnD(2,nn,hz,3,psil,sto3,Icon_local)

  Call derivnD(1,nn,hx,1,psil,sto4,Icon_local)
  Call derivnD(1,nn,hy,2,psil,sto5,Icon_local)
  Call derivnD(1,nn,hz,3,psil,sto6,Icon_local)

!
!   We compute H·Psi
!

hpsi   = -(sto1+sto2+sto3+sto4**2+sto5**2+sto6**2)*h2o2m4 + pot4
mu4err  = sum(den*hpsi)*dxyz/n4    ! Average over all the mu4's
If(jrun.Eq.1)Then
  hmean  = mu4err
Endif
Sto4   = -(hpsi-mu4err)

Sto1 = arun(jrun)*(Sto4 - brun(jrun)*q)
q = q + 3.*Sto1 - crun(jrun)*Sto4

if(jrun.eq.1)then
  hpsiold(:,:,:,2) = hpsiold(:,:,:,1)
  hpsiold(:,:,:,1) = Sto4
  psiold(:,:,:,3) = psiold(:,:,:,2)
  psiold(:,:,:,2) = psiold(:,:,:,1)
  psiold(:,:,:,1) = psil
endif

psil = psil + deltat*Sto1
Write(6,'("From Steprkrl: jrun, Max & Min of psil..:",I3,1p,2E15.6)')jrun,Maxval(psil),MinVal(psil)
psil(:,:,:)=Max(psil(:,:,:),underpsil)
psi = Exp(psil)
cnorm = sqrt(n4/(sum(psi**2)*dxyz))
psil = psil + Log(cnorm)
psi = Abs(psi)*cnorm
den = Abs(psi)**2

  if(jrun.le.3)then
    call poten()
  endif
enddo

do ix=1,3
  ioldp(ix)=ix
enddo

do ix=1,2
  ioldh(ix)=ix
enddo
mu4err = abs(1.0d0-mu4/hmean)     ! 'Error' in mu4
mu4    = hmean                    ! New chemical potential
return
end
