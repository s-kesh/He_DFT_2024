subroutine potenimp()
!
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
!
use quantal_imp !only: uimp!,Also2!,pairpot
use impur
use grid
use rho
use work1
!use interpol !, only: potdel,potpi,DelInter
implicit none
integer (kind=4) :: ix,iy,iz,i_imp,j_imp
integer (kind=4) :: ir
real    (kind=8) :: zt,yt,r
real    (kind=8) :: xx,yy,zz
real    (kind=8) :: rmod
!Write(6,'("Entrem a Potenimp...")')
 call fftfw_den()
potx4 = 0.d0
Do i_imp=1, N_imp
  Sto1(:,:,:) = psix(:,:,:,i_imp)*Conjg(psix(:,:,:,i_imp))
  call fftfw_1()
  wk2(:,:,:) = wk1(:,:,:)*uimp_k(:,:,:,i_imp)
  call fftbk_2()
  potx4 = potx4 + sto2 ! Get  ( int{ Psi_x**2 * V_X dr'} ) Potential due to the i_imp impurity for the 4He
  wk2(:,:,:) = fden(:,:,:)*uimp_k(:,:,:,i_imp)
  call fftbk_2()
  upotx(:,:,:,i_imp) = sto2(:,:,:) ! Get  ( int{ rho_4    * V_X dr'} ) Potential due to 4He for the impurity
  Do j_imp=1, N_imp
    If(j_imp.Eq.i_imp.Or.Trim(Adjustl(Selec_gs_k_k(i_imp,j_imp))).Eq.'null')Cycle
    Sto1(:,:,:) = psix(:,:,:,j_imp)*Conjg(psix(:,:,:,j_imp))
    call fftfw_1()
    wk2(:,:,:) = wk1(:,:,:)*uimp_k_k(:,:,:,i_imp,j_imp)
    call fftbk_2()
    upotx(:,:,:,i_imp) = upotx(:,:,:,i_imp)            &
                       + sto2(:,:,:) !Get(int{|psi_j|**2*V_X dr'})Potential due to j_imp for the i_imp
  EndDo  
EndDo  
!Write(6,'("Sortim de Potenimp...")')
!Call Flush(6)
end subroutine potenimp


