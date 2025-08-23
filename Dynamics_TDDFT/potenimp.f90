subroutine potenimp()
!
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
!
use classicimp , only: uimp!,Also2!,pairpot
use grid
use interpol !, only: potdel,potpi,DelInter
implicit none
integer (kind=4) :: ix,iy,iz
integer (kind=4) :: ir
real    (kind=8) :: zt,yt,r
real    (kind=8) :: xx,yy,zz
real    (kind=8) :: rmod

call updatepoten()

end subroutine potenimp

!---------------------------------------------------------------------------!

subroutine forceimp()
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
use classicimp!, only: uimp_k, rimp, F, F_ij, N_imp
use deriva
use grid
use rho
implicit none
integer (kind=4) :: ix,iy,iz,k,m
real    (kind=8) :: zt,yt,zt2,yt2,r,d,Select_pot
real    (kind=8) :: aux1,aux2,aux3,drV_HeXor,dzV_XTiO,dV

real (kind=8) :: sum_1, sum_2, sum_3

!.................................!
!... First, the term due to He ...!
!.................................!

F_ij = 0
do k=1,N_imp
  do m=k+1,N_imp
  aux1 = dsqrt((rimp(k,1)-rimp(m,1))**2 + (rimp(k,2)-rimp(m,2))**2 + (rimp(k,3)-rimp(m,3))**2)
  F_ij(k,m,:) = -Select_pot(drselec_gs_k_k(k,m),aux1,drr_cutoff_gs_k_k(k,m),drumax_gs_k_k(k,m)) * (rimp(k,:)-rimp(m,:))/aux1
  F_ij(m,k,:) = -F_ij(k,m,:)
  enddo
enddo

do k=1,N_imp

    sum_1 = 0.d0; sum_2 = 0.d0; sum_3 = 0.d0;

    !$omp parallel do private(ix,iy,iz) collapse(2) reduction(+:sum_1,sum_2,sum_3)
    do iz=1,nz
        do iy=1,ny
            do ix=1,nx
                sum_1 = sum_1 + uimp_k(k, ix,iy,iz) * dxden(ix,iy,iz)
                sum_2 = sum_2 + uimp_k(k, ix,iy,iz) * dyden(ix,iy,iz)
                sum_3 = sum_3 + uimp_k(k, ix,iy,iz) * dzden(ix,iy,iz)
            enddo
        enddo
    enddo
    !$omp end parallel do

    F(k,1) = -sum_1*dxyz
    F(k,2) = -sum_2*dxyz
    F(k,3) = -sum_3*dxyz

    do m=1,N_imp
        F(k,1) = F(k,1) + F_ij(k,m,1)
        F(k,2) = F(k,2) + F_ij(k,m,2)
        F(k,3) = F(k,3) + F_ij(k,m,3)
    enddo
enddo

!...........................................!
!... Second, the term due to the surface ...!
!...........................................!

!F(3) = F(3) - dzV_XTiO(rimp(3))


end subroutine forceimp


subroutine updatepoten()
use classicimp , only: uimp_k, uimp, rimp, N_imp!,Also2!,pairpot
use grid
use interpol
implicit none
real    (kind=8)              :: r,xt,yt,zt,rmod
integer (kind=4)              :: ix,iy,iz,ir,k,m
!save (lgridnoout)

!Write(*,*) rmaxinterpol

! For some reason just parallization increases uimp by a tinybit amount
!$omp parallel do private(ix,iy,iz,k,r,xt,yt,zt,ir,rmod,lstopimp,lgridnoout) default(shared) collapse(3)
do iz=1,nz; do iy=1,ny; do ix=1,nx
    xt = x(ix); yt = y(iy); zt = z(iz)
    do k = 1, N_imp
        r = dsqrt((xt-rimp(k,1))**2 &
            + (yt-rimp(k,2))**2 &
            + (zt-rimp(k,3))**2)
        ir = int(r/DelInter)+1
        if(r.gt.rmaxinterpol .and. lgridnoout)then
            lstopimp=.true.
            lgridnoout=.false.
            print *,'>>> WARNING in updatepoten, k,ix,iy,iz= ',k, ix, iy, iz,' r = ',r,' greater than rmax = ',rmaxinterpol
            r=rmaxinterpol
        endif
        rmod = mod(r,DelInter)/DelInter
        uimp_k(k,ix,iy,iz) =  potion(k,ir)*(1.d0-rmod) +  potion(k,ir+1)*rmod
    enddo
enddo; enddo; enddo
!$omp end parallel do


!$omp parallel do private(ix,iy,iz,k) collapse(3)
do iz=1,nz; do iy=1,ny; do ix=1,nx
    uimp(ix, iy, iz) = 0.d0
    do k = 1, N_imp
        uimp(ix, iy, iz) = uimp(ix, iy, iz) + uimp_k(k,ix,iy,iz)
    enddo
end do; end do; end do
!$omp end parallel do

end subroutine updatepoten
