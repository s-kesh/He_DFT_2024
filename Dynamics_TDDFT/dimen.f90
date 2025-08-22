!..............................................................
!...                     Subroutine dimen                   ...
!..............................................................
!
! This subroutine allocates almost all the arrays....
!
!
subroutine dimen()

use deriva
use field
use grid
use gridk
use he4
use classicimp
! use impur
use lenard4
use rho
use util1
use work1
use rkpc

implicit none

!.........................................
!.. Arrays for real and momentum grids ...
!.........................................
allocate (x(nx))      ; x = 0.0d0
allocate (px(nx))     ; px = 0.0d0
allocate (y(ny))      ; y = 0.0d0
allocate (py(ny))     ; py = 0.0d0
allocate (z(nz))      ; z = 0.0d0
allocate (pz(nz))     ; pz = 0.0d0
allocate (pmod(nx/2+1,ny,nz)) ; pmod = 0.0d0

!............................................
!.. Arrays for Lennard-Jones calculations ...
!............................................
allocate (fvlj4(nx/2+1,ny,nz)) ; fvlj4 = 0.0d0
allocate (delj4(nx,ny,nz))     ; delj4 = 0.0d0

!..................................................
!.. Arrays for partial derivatives of densities ...
!..................................................
allocate(dxden(nx,ny,nz)) ; dxden = 0.0d0
allocate(dyden(nx,ny,nz)) ; dyden = 0.0d0
allocate(dzden(nx,ny,nz)) ; dzden = 0.0d0

!..........................
!.. Arrays for Helium 4 ...
!..........................
allocate (pot4(nx,ny,nz))     ; pot4 = 0.0d0
allocate (hpsi(nx,ny,nz))     ; hpsi = 0.0d0
allocate (den(nx,ny,nz))      ; den = 0.0d0
allocate (psi(nx,ny,nz))      ; psi = 0.0d0
allocate (psiold(nx,ny,nz,3)) ; psiold = 0.0d0
allocate (hpsiold(nx,ny,nz,2)); hpsiold = 0.0d0
allocate (dencg(nx,ny,nz))    ; dencg = 0.0d0
allocate (fden(nx/2+1,ny,nz)) ; fden = 0.0d0
!allocate (fpsi(nx/2+1,ny,nz))
allocate (wcgk(nx/2+1,ny,nz)) ; wcgk = 0.0d0

!..............................
!.. Arrays for the impurity ...
!..............................
allocate(uext(nx,ny,nz)) ; uext = 0.0d0
if (Lsolid) then
   allocate(penalty(nx,ny,nz))
   penalty = 0.0d0
end if
allocate(uimp(nx,ny,nz))            ; uimp = 0.0d0
allocate(uimp_k(N_imp,nx,ny,nz))    ; uimp_k = 0.0d0

allocate(qr(N_imp,3))    ; qr    = 0.0d0
allocate(qv(N_imp,3))    ; qv    = 0.0d0
allocate(Stor(N_imp,3))  ; Stor  = 0.0d0
allocate(pcr(N_imp,3))   ; pcr   = 0.0d0
allocate(pcv(N_imp,3))   ; pcv   = 0.0d0
allocate(rimp(N_imp,3))  ; rimp  = 0.0d0
allocate(vimp(N_imp,3))  ; vimp  = 0.0d0
allocate(aimp(N_imp,3))  ; aimp  = 0.0d0
allocate(F(N_imp,3))     ; F     = 0.0d0
allocate(F_ij(N_imp,N_imp,3)) ; F_ij = 0.0d0
allocate(rimpold(N_imp,3,3))  ; rimpold = 0.0d0
allocate(vimpold(N_imp,3,3))  ; vimpold = 0.0d0
allocate(aimpold(N_imp,3,2))  ; aimpold = 0.0d0

allocate(filerimp_k(N_imp))
allocate(filevimp_k(N_imp))
allocate(fileaimp_k(N_imp))
allocate(m_imp_u(N_imp))      ; m_imp_u    = 0.0d0
allocate(m_imp(N_imp))        ; m_imp      = 0.0d0
allocate(selec_gs_k(N_imp))
allocate(selec_gs_k_k(N_imp,N_imp))
allocate(drselec_gs_k_k(N_imp,N_imp))
allocate(r_cutoff_gs_k(N_imp))        ; r_cutoff_gs_k = 0.0d0
allocate(r_cutoff_gs_k_k(N_imp,N_imp)); r_cutoff_gs_k_k = 0.0d0
allocate(drr_cutoff_gs_k_k(N_imp,N_imp)) ; drr_cutoff_gs_k_k = 0.0d0
allocate(umax_gs_k(N_imp))            ; umax_gs_k = 0.0d0
allocate(umax_gs_k_k(N_imp,N_imp))    ; umax_gs_k_k = 0.0d0
allocate(drumax_gs_k_k(N_imp,N_imp))  ; drumax_gs_k_k = 0.0d0

!....................................................
!.. Arrays for temporal storage and working areas ...
!....................................................
allocate(sto1(nx,ny,nz)) ; sto1 = 0.0d0
allocate(sto2(nx,ny,nz)) ; sto2 = 0.0d0
allocate(sto3(nx,ny,nz)) ; sto3 = 0.0d0
allocate(sto4(nx,ny,nz)) ; sto4 = 0.0d0
allocate(sto5(nx,ny,nz)) ; sto5 = 0.0d0
allocate(sto6(nx,ny,nz)) ; sto6 = 0.0d0

allocate(wk1(nx/2+1,ny,nz)) ; wk1 = (0.0d0,0.0d0)
allocate(wk2(nx/2+1,ny,nz)) ; wk2 = (0.0d0,0.0d0)
allocate(wk3(nx/2+1,ny,nz)) ; wk3 = (0.0d0,0.0d0)

allocate(sto1c(nx,ny,nz)) ; sto1c = (0.0d0,0.0d0)
allocate(sto2c(nx,ny,nz)) ; sto2c = (0.0d0,0.0d0)
allocate(sto3c(nx,ny,nz)) ; sto3c = (0.0d0,0.0d0)
allocate(sto4c(nx,ny,nz)) ; sto4c = (0.0d0,0.0d0)
allocate(sto5c(nx,ny,nz)) ; sto5c = (0.0d0,0.0d0)
allocate(sto6c(nx,ny,nz)) ; sto6c = (0.0d0,0.0d0)
allocate(sto7c(nx,ny,nz)) ; sto7c = (0.0d0,0.0d0)
allocate(sto8c(nx,ny,nz)) ; sto8c = (0.0d0,0.0d0)

!...................................................................
!.. Array for Runge-Kutta-Gill & Predictor-Corrector-Modificator ...
!...................................................................
allocate(q(nx,ny,nz))  ; q  = 0.0d0
allocate(pc(nx,ny,nz)) ; pc = 0.0d0

allocate(timec(nx,ny,nz)) ; timec = 0.0d0

return
end
