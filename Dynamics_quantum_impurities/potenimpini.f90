subroutine potenimpini
!use interpol !, only:DelInter,potpi,potdel,npot,vpi,delta
use grid
use gridk
use impur
use quantal_imp
implicit none
real      (kind=8)              :: qmax, qmax0=100.d0, FT_V_spline
real      (kind=8)              :: hq,p, twopi2, px2, py2, pz2
integer   (kind=4)              :: i_imp,j_imp,ix,iy,iz, Nq0=10001
real      (kind=8)              :: Pi
   Pi     = 4.d0*Datan(1.d0)
   twopi2 = 4.d0*Pi**2
   qmax = Sqrt(pmaxx**2+pmaxy**2+pmaxz**2)
   qmax=Max(qmax0,qmax)
   Nq = Max(Nq0, Nq)
   Hq   = qmax/(Nq-1)
   Write(6,'("Començem el calcul de la T.F. X-He....")')
   Do i_imp=1,N_imp
     Write(6,'("Potencial que fem servir...:"I2,2x,A)')i_imp,Selec_gs_k(i_imp)
     Do iz=1, nz
       Do iy=1, ny
         Do ix=1, nx/2 + 1
           p=pmod(ix,iy,iz)
           uimp_k(ix,iy,iz,i_imp) =                                                      &
cmplx(FT_V_spline(p,qmax,Hq,Selec_gs_k(i_imp),r_cutoff_gs_k(i_imp),umax_gs_k(i_imp)))
         EndDo
       EndDo
     EndDo
     Write(6,'("Començem el calcul de la T.F. X-X.....")')
     Do j_imp=1,N_imp
       uimp_k_k(:,:,:,i_imp,j_imp) = (0.d0, 0.d0)
       Write(6,'("Potencial que fem servir...:"I2,I2,2x,A)')i_imp,j_imp,Selec_gs_k_k(i_imp,j_imp)
       If(Trim(Adjustl(Selec_gs_k_k(i_imp,j_imp))).Eq."null".Or.i_imp.Eq.j_imp)Cycle
       Do iz=1, nz
         Do iy=1, ny
           Do ix=1, nx/2 + 1
             p=pmod(ix,iy,iz)
             uimp_k_k(ix,iy,iz,i_imp,j_imp) =                                             &
cmplx(FT_V_spline(p,qmax,Hq,Selec_gs_k_k(i_imp,j_imp),r_cutoff_gs_k_k(i_imp,j_imp),umax_gs_k_k(i_imp,j_imp)))
           EndDo
         EndDo
       EndDo
     EndDo
   EndDo
   Write(6,'("Hem acabat les T.F.: X-He i X-X.....")')

end subroutine potenimpini

