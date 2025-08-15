!---------------------------------------------------------------------
!---                             Subroutine Poten                  ---
!---------------------------------------------------------------------

subroutine poten()

use alphasterm ! (ualphas)
use grid     ! (nxyz)
use he4      ! (cp4,cpp4)
! use impur    ! (vq,)
use classicimp
use lenard4  ! (wk2,pelj4,fvlj4,delj4,core4,lalphas)
use field    ! (pot4,limp)
use rho      ! (dencg,den,wcgk,)
use work1    ! (sto1,sto2,sto3,sto4)
use coalescence

implicit none

real    (kind=8) :: a0,a1, dtemp
integer (kind=4) :: ix,iy,iz


if (ldroplet_frozen) then
    ! Write(6,*) "Lfrozen_first_iteration=",Lfrozen_first_iteration
    if(Lfrozen_first_iteration)then
        !...............................................
        !... Calculate of Fourier Transforms for 4He ...
        !...............................................

        ! call fftfw(den,fden)
        call fftfw_den()

        !..............................
        !.. Coarse-graining density ...
        !.. Alfa_s density          ...
        !....................................

        !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
        do iz=1, nz; do iy=1, ny
            do ix=1, nx/2 + 1
                wk1(ix, iy, iz) = fden(ix, iy, iz)*wcgk(ix,iy,iz)
            end do
        end do; end do
        !$omp end parallel do
        ! call fftbk(wk1,dencg)  ! get Coarse graining density
        call fftbk_cg()  ! get Coarse graining density

        call derden() ! Calculate derivatives of the density

        !..............................
        ! Lennard-Jones contribution  .
        !..............................

        !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
        do iz=1, nz; do iy=1, ny
            do ix=1, nx/2 + 1
                wk1(ix, iy, iz) = fden(ix, iy, iz)*fvlj4(ix,iy,iz)
            end do
        end do; end do
        !$omp end parallel do
        ! call fftbk(wk1,delj4) ! Get delj4 -> (   int{ rho_4*V_4 dr'}  )
        call fftbk_lj() ! Get delj4 -> (   int{ rho_4*V_4 dr'}  )

        !........................
        ! Correlation terms.  ...
        !........................

        a0 = cp4 /2.d0               ! Auxiliar variable useful for saving operations
        a1 = cpp4/3.d0               ! Auxiliar variable useful for saving operations

        !   The rest of the correlation contribution is calculated as
        !   a convolution product.

        !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
        do iz=1, nz; do iy=1, ny
            do ix=1, nx
                sto1(ix,iy,iz) = den(ix,iy,iz)*dencg(ix,iy,iz)*(cp4+cpp4*dencg(ix,iy,iz))
            end do
        end do; end do
        !$omp end parallel do

        ! call fftfw(sto1,wk1)
        call fftfw_1()

        !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
        do iz=1, nz; do iy=1, ny
            do ix=1, nx/2 + 1
                wk1(ix,iy,iz) = wk1(ix,iy,iz)*wcgk(ix,iy,iz)
            end do
        end do; end do
        !$omp end parallel do

        ! call fftbk(wk1,sto1)
        call fftbk_1()


        !..........................
        !... Solid penalty term ...
        !..........................
        If(lsolid)Then
            !$omp parallel do private(ix,iy,iz, dtemp) collapse(2)
            do iz=1, nz; do iy=1, ny
                do ix=1, nx
                    dtemp = dtanh(beta*(den(ix,iy,iz)-den_m))
                    penalty(ix,iy,iz) = C*(1.d0 + dtemp + beta*den(ix,iy,iz)*(1.d0 - dtemp**2) )
                end do
            end do; end do
            !$omp end parallel do
        Endif

        !.......................
        !.. Final 'Potential' ...
        !........................


        !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
        do iz=1, nz; do iy=1, ny
            do ix=1, nx
                pot4(ix,iy,iz) = delj4(ix,iy,iz) +                                 &   ! Lennard-Jones
                                dencg(ix,iy,iz)**2*(a0+a1*dencg(ix,iy,iz)) +      &   ! Correlation
                                sto1(ix,iy,iz)                                        ! Correlation
            end do
        end do; end do
        !$omp end parallel do

        if(core4.eq.'OTC') then
            !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
            do iz=1, nz; do iy=1, ny
                do ix=1, nx/2 + 1
                    wk1(ix,iy,iz)   = fden(ix,iy,iz)*kalfs(ix,iy,iz)
                end do
            end do; end do
            !$omp end parallel do
            !     call fftbk(wk1,denalf) ! Get Alfa_s density
            call fftbk_as() ! Get Alfa_s density
            call term_alfa()    ! Calculates alpha_s contribution to the field

            call daxpy(nx*ny*nz, 1.0d0, ualphas, 1, pot4, 1)
        end if

        If(lsolid)Then
            call daxpy(nx*ny*nz, 1.0d0, penalty, 1, pot4, 1)
        Endif
    endif ! (If Lfrozen_first_iteration)

else ! (If ldroplet_frozen)

    !...............................................
    !... Calculate of Fourier Transforms for 4He ...
    !...............................................

    ! call fftfw(den,fden)
    call fftfw_den()

    !..............................
    !.. Coarse-graining density ...
    !.. Alfa_s density          ...
    !....................................

    !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
    do iz=1, nz; do iy=1, ny
        do ix=1, nx/2 + 1
            wk1(ix, iy, iz) = fden(ix, iy, iz)*wcgk(ix,iy,iz)
        end do
    end do; end do
    !$omp end parallel do
    ! call fftbk(wk1,dencg)  ! get Coarse graining density
    call fftbk_cg()  ! get Coarse graining density

    call derden() ! Calculate derivatives of the density

    !..............................
    ! Lennard-Jones contribution  .
    !..............................

    !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
    do iz=1, nz; do iy=1, ny
        do ix=1, nx/2 + 1
            wk1(ix, iy, iz) = fden(ix, iy, iz)*fvlj4(ix,iy,iz)
        end do
    end do; end do
    !$omp end parallel do
    ! call fftbk(wk1,delj4) ! Get delj4 -> (   int{ rho_4*V_4 dr'}  )
    call fftbk_lj() ! Get delj4 -> (   int{ rho_4*V_4 dr'}  )

    !........................
    ! Correlation terms.  ...
    !........................

    a0 = cp4 /2.d0               ! Auxiliar variable useful for saving operations
    a1 = cpp4/3.d0               ! Auxiliar variable useful for saving operations

    !   The rest of the correlation contribution is calculated as
    !   a convolution product.

    !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
    do iz=1, nz; do iy=1, ny
        do ix=1, nx
            sto1(ix,iy,iz) = den(ix,iy,iz)*dencg(ix,iy,iz)*(cp4+cpp4*dencg(ix,iy,iz))
        end do
    end do; end do
    !$omp end parallel do

    ! call fftfw(sto1,wk1)
    call fftfw_1()

    !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
    do iz=1, nz; do iy=1, ny
        do ix=1, nx/2 + 1
            wk1(ix,iy,iz) = wk1(ix,iy,iz)*wcgk(ix,iy,iz)
        end do
    end do; end do
    !$omp end parallel do

    ! call fftbk(wk1,sto1)
    call fftbk_1()


    !..........................
    !... Solid penalty term ...
    !..........................
    If(lsolid)Then
        !$omp parallel do default(shared) private(ix,iy,iz, dtemp) collapse(2)
        do iz=1, nz; do iy=1, ny
            do ix=1, nx
                dtemp = dtanh(beta*(den(ix,iy,iz)-den_m))
                penalty(ix,iy,iz) = C*(1.d0 + dtemp + beta*den(ix,iy,iz)*(1.d0 - dtemp**2) )
            end do
        end do; end do
        !$omp end parallel do
    Endif

    !.......................
    !.. Final 'Potential' ...
    !........................

    !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
    do iz=1, nz; do iy=1, ny
        do ix=1, nx
            pot4(ix,iy,iz) = delj4(ix,iy,iz) +                                 &   ! Lennard-Jones
                            dencg(ix,iy,iz)**2*(a0+a1*dencg(ix,iy,iz)) +      &   ! Correlation
                            sto1(ix,iy,iz)                                        ! Correlation
        end do
    end do; end do
    !$omp end parallel do

    if(core4.eq.'OTC') then
        !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
        do iz=1, nz; do iy=1, ny
            do ix=1, nx/2 + 1
                wk1(ix,iy,iz)   = fden(ix,iy,iz)*kalfs(ix,iy,iz)
            end do
        end do; end do
        !$omp end parallel do
    !     call fftbk(wk1,denalf) ! Get Alfa_s density
        call fftbk_as() ! Get Alfa_s density
        call term_alfa()    ! Calculates alpha_s contribution to the field

        call daxpy(nx*ny*nz, 1.0d0, ualphas, 1, pot4, 1)
    end if


    If(lsolid)Then
        call daxpy(nx*ny*nz, 1.0d0, penalty, 1, pot4, 1)
    Endif
endif ! If ldroplet_frozen

return

end
