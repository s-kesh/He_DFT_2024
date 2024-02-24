!-----------------------------------------------------------------------------
!--                          Subroutine evolo                              ---
!-----------------------------------------------------------------------------
!
subroutine evolo(deltat,mu4,mu4err,n4)

Use Para_DerivnD
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use impur  !
use rho    ! (*psi*)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)

implicit none

integer (kind=4) :: n4
real    (kind=8) :: deltat,mu,muerr,n4real
real    (kind=8) :: mu4,mu4err

integer (kind=4) :: ix,iy,iz
real    (kind=8) :: hmean,a1
real    (kind=8) :: cnorm,sto,aux,dlaux


dlaux=dlog(2.0d0)


!.......................
!.. Laplacian of Psi ...
!.......................

!
!   icon   =  0 ! Take the derivative.
!   icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!   icon   =  8 ! Take the derivative. Use periodic conditions.
   
  Call derivnD(2,nn,hx,1,psi,sto1,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3,Icon)

!.................................................................. hpsi  = (T+U)*Psi_4 (With impurity)

   sto4 = -(sto1+sto2+sto3)*h2o2m4 + pot4*psi
   hpsi = sto4
  
   If(Lbulk)Then
     hmean  = sum(psi*hpsi)/sum(den)    ! Average over all the mu4's
     mu4err = abs(1.0d0-mu4/hmean)           ! Relative error
   Else
     hmean  = sum(psi*hpsi)*dxyz/n4          ! Average over all the mu4's
     mu4err = abs(1.0d0-mu4/hmean)           ! Relative error
     mu4    = hmean
   Endif

   a1=1.0d0-deltat*(hmean-mu4)
   a1=deltat/a1

   hpsi =a1*(hpsi-mu4*psi)

   if(iron.ne.0) then
     do ix=1,iron
       call ironing(hpsi,nx,ny,nz)       ! Smoothing  (H-mu)*Psi
     end do
   end if

   forall(ix=1:nx,iy=1:ny,iz=1:nz)
     psinw(ix,iy,iz) = psi(ix,iy,iz)-  hpsi(ix,iy,iz)            &
                     + vdt(1) * (psi(ix,iy,iz)-psi1(ix,iy,iz))   &
                     + vdt(2) * (psi1(ix,iy,iz)-psi2(ix,iy,iz))
     psi2(ix,iy,iz)  = psi1(ix,iy,iz)
     psi1(ix,iy,iz)  = psi(ix,iy,iz)
   end forall

   hpsi = sto4

!   psinw=max(psimin,Abs(psinw))

   If(Lbulk)Then
     cnorm=1.0d0
     n4real=sum(den)*dxyz
     n4=n4real+0.5
   Else
     cnorm = sqrt(n4/(sum(psinw**2)*dxyz))
   Endif

   forall(ix=1:nx,iy=1:ny,iz=1:nz)
     psi(ix,iy,iz) = Abs(psinw(ix,iy,iz))*cnorm
     den(ix,iy,iz) = psi(ix,iy,iz)**2
   end forall

return
end
