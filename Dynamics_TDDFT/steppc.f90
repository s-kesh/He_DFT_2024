SUBROUTINE STEPPC(deltat,errHe,errimp,errvimp)
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
use classicimp
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc)
use coalescence

implicit none

real (kind=8) :: c112=112.d0/121.d0
real (kind=8) :: c9  =9.d0/121.d0
real (kind=8) :: c1o3=1.d0/3.d0
real (kind=8) :: c4o3=4.d0/3.d0
real (kind=8) :: c5o3=5.d0/3.d0

integer (kind=4) :: ix,iy,iz,iaux
real    (kind=8) :: deltat
real    (kind=8) :: errHe, errimp,errvimp
!real    (kind=8) :: auxr(3)
complex (kind=8) :: auxc(6)
complex (kind=8) :: aux1c,aux2c,aux3c,aux4c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)

complex (kind=8) :: tmp_ke, tmp_pot, tmp_tot
real (kind=8), external :: zdotc
integer:: i


if(.not. ldroplet_frozen)then
    ! Computing gradients
    Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
    Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
    Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

    ! Computing HPsi
    ! Predictor
    !$omp parallel do private(ix,iy,iz, tmp_ke, tmp_pot, tmp_tot) schedule(static) collapse(2)
    do iz=1,nz; do iy=1,ny
        !$omp simd
        do ix=1,nx
            tmp_ke = h2o2m4*(sto1c(ix,iy,iz)+sto2c(ix,iy,iz)+sto3c(ix,iy,iz))
            tmp_pot = pot4(ix,iy,iz)*psi(ix,iy,iz)
            sto4c(ix,iy,iz) = timec(ix,iy,iz)*(tmp_ke - tmp_pot) - ci*uimp(ix,iy,iz)*psi(ix,iy,iz)
        end do
        !$omp end simd
    enddo; enddo
    !$omp end parallel do

    !$omp parallel do private(ix,iy,iz) default(shared) schedule(static) collapse(2)
    do iz=1,nz; do iy=1,ny
        !$omp simd
        do ix=1,nx
            Sto1c(ix,iy,iz) = psiold(ix,iy,iz,ioldp(3)) &
            + c4o3*deltat*(2.d0*Sto4c(ix,iy,iz)-hpsiold(ix,iy,iz,ioldh(1)) &
            + 2.d0*hpsiold(ix,iy,iz,ioldh(2)))
        enddo
        !$omp end simd
    enddo; enddo
    !$omp end parallel do

    ! Modifier
    call zcopy(nx*ny*nz, psi, 1, psiold(:,:,:,ioldp(3)), 1)
    call zcopy(nx*ny*nz, sto4c, 1, hpsiold(:,:,:,ioldh(2)), 1)
    !$omp parallel do private(ix,iy,iz) schedule(static) collapse(2)
    do iz=1,nz; do iy=1,ny
        !$omp simd
        do ix=1,nx
            psi(ix,iy,iz) = Sto1c(ix,iy,iz) - c112*pc(ix,iy,iz)
            pc(ix,iy,iz) = Sto1c(ix,iy,iz)
            den(ix,iy,iz) = real(psi(ix,iy,iz))**2 + aimag(psi(ix,iy,iz))**2
        enddo
        !$omp end simd
    enddo; enddo
    !$omp end parallel do

    ! Aqui reubicamos los indices para no tener que mover las fuciones
    iaux=ioldh(2)
    ioldh(2)=ioldh(1)
    ioldh(1)=iaux
endif


if(.not. Lcoalescence ) then
    !$omp simd
    do i=1,N_imp
        ! Predictor
        stor(i,1) = rimpold(i,1,ioldr(3)) + c4o3*deltat*(2.d0*vimp(i,1) - vimpold(i,1,ioldv(1)) + 2.d0*vimpold(i,1,ioldv(2)))
        stor(i,2) = rimpold(i,2,ioldr(3)) + c4o3*deltat*(2.d0*vimp(i,2) - vimpold(i,2,ioldv(1)) + 2.d0*vimpold(i,2,ioldv(2)))
        stor(i,3) = rimpold(i,3,ioldr(3)) + c4o3*deltat*(2.d0*vimp(i,3) - vimpold(i,3,ioldv(1)) + 2.d0*vimpold(i,3,ioldv(2)))

        ! Modificador
        rimpold(i,1,ioldr(3)) = rimp(i,1)
        rimpold(i,2,ioldr(3)) = rimp(i,2)
        rimpold(i,3,ioldr(3)) = rimp(i,3)
    end do
    !$omp end simd

    !$omp simd
    do i=1,N_imp
        ! Modificador
        rimp(i,1) = stor(i,1) + c112*pcr(i,1)
        rimp(i,2) = stor(i,2) + c112*pcr(i,2)
        rimp(i,3) = stor(i,3) + c112*pcr(i,3)
        pcr(i,1) = stor(i,1)
        pcr(i,2) = stor(i,2)
        pcr(i,3) = stor(i,3)
    enddo
    !$omp end simd

    !... velocities ...!
    !$omp simd
    do i=1,N_imp
        !..................!
        ! Predictor
        Stor(i,1) = vimpold(i,1,ioldv(3)) + c4o3*deltat*(2.d0*aimp(i,1) - aimpold(i,1,iolda(1)) + 2.d0*aimpold(i,1,iolda(2)))
        Stor(i,2) = vimpold(i,2,ioldv(3)) + c4o3*deltat*(2.d0*aimp(i,2) - aimpold(i,2,iolda(1)) + 2.d0*aimpold(i,2,iolda(2)))
        Stor(i,3) = vimpold(i,3,ioldv(3)) + c4o3*deltat*(2.d0*aimp(i,3) - aimpold(i,3,iolda(1)) + 2.d0*aimpold(i,3,iolda(2)))

        ! Modificador
        vimpold(i,1,ioldv(3)) = stor(i,1) - c112*pcv(i,1)
        vimpold(i,2,ioldv(3)) = stor(i,2) - c112*pcv(i,2)
        vimpold(i,3,ioldv(3)) = stor(i,3) - c112*pcv(i,3)
        pcv(i,1) = stor(i,1)
        pcv(i,2) = stor(i,2)
        pcv(i,3) = stor(i,3)


        !... accelerations ...!
        aimpold(i,1,iolda(2)) = aimp(i,1)
        aimpold(i,2,iolda(2)) = aimp(i,2)
        aimpold(i,3,iolda(2)) = aimp(i,3)
    end do
    !$omp end simd

    ! Reubicacion indices
    iaux=iolda(2)  ; iolda(2)=iolda(1)   ; iolda(1)=iaux

    call potenimp()
    call poten()
    call forceimp()
endif


if(.not. Ldroplet_frozen)then
    Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
    Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
    Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

    !$omp parallel do private(ix,iy,iz, tmp_ke, tmp_pot, tmp_tot) schedule(static) collapse(2)
    do iz=1,nz; do iy=1,ny
        !$omp simd
        do ix=1,nx
            tmp_ke = h2o2m4*timec(ix,iy,iz)*(sto1c(ix,iy,iz)+sto2c(ix,iy,iz)+sto3c(ix,iy,iz))
            tmp_pot = (-timec(ix,iy,iz)*pot4(ix,iy,iz) - ci*uimp(ix,iy,iz))*psi(ix,iy,iz)
            tmp_tot = tmp_ke + tmp_pot
            sto4c(ix,iy,iz) = tmp_tot
        enddo
        !$omp end simd
    enddo; enddo
    !$omp end parallel do

    !$omp parallel do private(ix,iy,iz) schedule(static) collapse(2)
    do iz=1,nz; do iy=1,ny
        !$omp simd
        do ix=1,nx
            ! Corrector
            Sto5c(ix,iy,iz) = 0.125d0*( 9.d0*psiold(ix,iy,iz,ioldp(3)) - psiold(ix,iy,iz,ioldp(2))   &
                +3.d0*deltat*(Sto4c(ix,iy,iz) + 2.d0*hpsiold(ix,iy,iz,ioldh(1)) - hpsiold(ix,iy,iz,ioldh(2))  ))

            pc(ix,iy,iz) = pc(ix,iy,iz) - sto5c(ix,iy,iz)
            ! Valor final
            psi(ix,iy,iz) = Sto5c(ix,iy,iz) + c9*pc(ix,iy,iz)
            den(ix,iy,iz) = real(psi(ix,iy,iz))**2 + aimag(psi(ix,iy,iz))**2
        enddo
        !$omp end simd
    enddo; enddo
    !$omp end parallel do

    !$omp parallel do private(ix,iy,iz) collapse(2) reduction(+:errHe)
    do iz=1,nz; do iy=1,ny
        !$omp simd reduction(+:errHe)
        do ix=1,nx
            errHe = errHe + abs(pc(ix,iy,iz))
        enddo
        !$omp end simd
    enddo; enddo
    !$omp end parallel do
    errHe=errHe*c9/nxyz

    !
    ! Aqui reubicamos los indices para no tener que mover las fuciones
    !
    iaux=ioldp(3)
    ioldp(3)=ioldp(2)
    ioldp(2)=ioldp(1)
    ioldp(1)=iaux
else
    errHe=0d0  ! To avoid numerical issues
endif


if(.not. Lcoalescence ) then
    !... positions ...!
    !$omp simd
    do i=1,N_imp
        ! Corrector:
        stor(i,1) = 0.125d0*(9.d0*rimpold(i,1,ioldr(3)) - rimpold(i,1,ioldr(2))     &
                    + 3.d0*deltat*(vimpold(i,1,ioldv(3)) &
                    + 2.d0*vimp(i,i) - vimpold(i,1,ioldv(1))))
        stor(i,2) = 0.125d0*(9.d0*rimpold(i,2,ioldr(3)) - rimpold(i,2,ioldr(2))     &
                    + 3.d0*deltat*(vimpold(i,2,ioldv(3)) &
                    + 2.d0*vimp(i,i) - vimpold(i,2,ioldv(1))))
        stor(i,3) = 0.125d0*(9.d0*rimpold(i,3,ioldr(3)) - rimpold(i,3,ioldr(2))     &
                    + 3.d0*deltat*(vimpold(i,3,ioldv(3)) &
                    + 2.d0*vimp(i,i) - vimpold(i,3,ioldv(1))))
    enddo
    !$omp end simd

    !$omp simd
    do i=1,N_imp
        ! Corrector
        pcr(i,1) = pcr(i,1) - stor(i,1)
        pcr(i,2) = pcr(i,2) - stor(i,2)
        pcr(i,3) = pcr(i,3) - stor(i,3)

        ! Valor final:
        rimp(i,1) = stor(i,1) + c9*pcr(i,1)
        rimp(i,2) = stor(i,2) + c9*pcr(i,2)
        rimp(i,3) = stor(i,3) + c9*pcr(i,3)
    enddo
    !$omp end simd

    !$omp simd reduction(+:errimp)
    do i=1,N_imp
        errimp = errimp + Abs(c9*pcr(i,1)) + Abs(c9*pcr(i,2)) + Abs(c9*pcr(i,3))
    enddo
    !$omp end simd
    errimp = errimp*0.3333333333d0/N_imp

    ! Reubicacion
    iaux=ioldr(3) ; ioldr(3)=ioldr(2) ; ioldr(2)=ioldr(1) ; ioldr(1)=iaux

    !... velocities ...!
    !$omp simd
    do i=1,N_imp
        ! Corrector:
        stor(i,1) = 0.125d0*(9.0d0*vimp(i,1) - vimpold(i,1,ioldv(2))) &
                    + 3.d0*deltat*(aimp(i,1) &
                    + 2.d0*aimpold(i,1,iolda(1)) - aimpold(i,1,iolda(2)))
        stor(i,2) = 0.125d0*(9.0d0*vimp(i,2) - vimpold(i,2,ioldv(2))) &
                    + 3.d0*deltat*(aimp(i,2) &
                    + 2.d0*aimpold(i,2,iolda(1)) - aimpold(i,2,iolda(2)))
        stor(i,3) = 0.125d0*(9.0d0*vimp(i,3) - vimpold(i,3,ioldv(2))) &
                    + 3.d0*deltat*(aimp(i,3) &
                    + 2.d0*aimpold(i,3,iolda(1)) - aimpold(i,3,iolda(2)))

        pcv(i,1) = pcv(i,1) - stor(i,1)
        pcv(i,2) = pcv(i,2) - stor(i,2)
        pcv(i,3) = pcv(i,3) - stor(i,3)
        vimpold(i,1,ioldv(3)) = vimp(i,1)
        vimpold(i,2,ioldv(3)) = vimp(i,2)
        vimpold(i,3,ioldv(3)) = vimp(i,3)

        ! Valor final
        vimp(i,1) = stor(i,1) + c9*pcv(i,1)
        vimp(i,2) = stor(i,2) + c9*pcv(i,2)
        vimp(i,3) = stor(i,3) + c9*pcv(i,3)
    enddo
    !$omp end simd

    !$omp simd reduction(+:errvimp)
    do i=1,N_imp
        errvimp = errvimp + Abs(c9*pcv(i,1)) + Abs(c9*pcv(i,2)) + Abs(c9*pcv(i,3))
    enddo
    !$omp end simd
    errvimp = errvimp*0.3333333333d0/N_imp

    ! Reubicacion
    iaux=ioldv(3) ; ioldv(3)=ioldv(2) ; ioldv(2)=ioldv(1) ; ioldv(1)=iaux
endif


return
end
