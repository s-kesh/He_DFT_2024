subroutine potimpx()

!-------------------------------------------------------------------
!
!.. This subroutine computes the fourier transform of
!   the impurity potential.
!
!
! INPUT variables:
!
! (nx,ny,nz)          ----> grid points
! (pmaxx,pmaxy,pmaxz) ----> Maximum frequency (Fourier Transforms)
! (pmod(:,:,:))    -------> Array with the p-values for Fourier
! (leepot)     -----------> 'yes'/'no ' (Reads the fourier transform
!                            of patil potential from a file.
! (vxpot) ----------------> Name of the external file ....
! (nq)    ----------------> Number of points for potential
!                           using Romberg quadrature
! (elem) -----------------> Chemical symbol of alkali impurity
! (tol)  -----------------> Tolerance for Romberg quadrature
! (rmin) -----------------> Core for Patil Potential, Below
!                             this value we put the same value
!
!( OUTPUT variables:
!
! (vq(:,:,:) ------------> Fourier Transform for Patil Potential
!-------------------------------------------------------------------
Use seleccio_de_potencial
use grid  ! (nx,ny,nz)
use gridk ! (pmaxx,pmaxy,pmaxz,pmod)
use impur ! (Selec, r_cutoff, umax)
use util1 ! irespar

implicit none

real      (kind=8)              :: qmax, qmax0=100.d0, FT_V_spline
real      (kind=8)              :: FT_S_P_V_spline, V1p, V2p, Vp_Pi
real      (kind=8)              :: hq,p, twopi2, px2, py2, pz2
integer   (kind=4)              :: ix,iy,iz, Nq0=10001
!real      (kind=8)              :: Pi
   Pi     = 4.d0*Datan(1.d0)
   twopi2 = 4.d0*Pi**2
   qmax = Sqrt(pmaxx**2+pmaxy**2+pmaxz**2)
   qmax=Max(qmax0,qmax)
   Nq = Max(Nq0, Nq)
   Hq   = qmax/(Nq-1)
   V1p   = FT_S_P_V_spline(1.d0,qmax,Hq,Selec_Sigma,r_cutoff_Sigma,umax_Sigma,Selec_Pi,r_cutoff_Pi,umax_Pi,0) 
   Do iz=1, nz
     pz2 = pz(iz)**2
     Do iy=1, ny
       py2 = py(iy)**2 
       Do ix=1, nx/2+1
         px2 = px(ix)**2
!         p=pmod(ix,iy,iz)
         p=sqrt(px2+py2+pz2)
         Vp_Pi = FT_V_spline(p,qmax,Hq,Selec_Pi,r_cutoff_Pi,umax_Pi)
         V1p   = FT_S_P_V_spline(p,qmax,Hq,Selec_Sigma,r_cutoff_Sigma,umax_Sigma,Selec_Pi,r_cutoff_Pi,umax_Pi,1) 
         V2p   = FT_S_P_V_spline(p,qmax,Hq,Selec_Sigma,r_cutoff_Sigma,umax_Sigma,Selec_Pi,r_cutoff_Pi,umax_Pi,2) 
         Vp_x(ix,iy,iz)  = Vp_Pi - (V1p + V2p*px2)/twopi2
         Vp_y(ix,iy,iz)  = Vp_Pi - (V1p + V2p*py2)/twopi2
         Vp_z(ix,iy,iz)  = Vp_Pi - (V1p + V2p*pz2)/twopi2
         Vp_xy(ix,iy,iz) =       - V2p*px(ix)*py(iy)/twopi2
         Vp_xz(ix,iy,iz) =       - V2p*px(ix)*pz(iz)/twopi2
         Vp_yz(ix,iy,iz) =       - V2p*py(iy)*pz(iz)/twopi2
       EndDo
     EndDo
   EndDo
return

end 
