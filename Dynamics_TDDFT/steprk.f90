SUBROUTINE STEPRK(deltat)
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
use impur  !
use classicimp
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc rutines)
use coalescence

implicit none

real (kind=8) :: arun(4),brun(4),crun(4)

integer (kind=4) :: i, ix,iy,iz,jrun
real    (kind=8) :: deltat
real    (kind=8) :: aux2r(3)
complex (kind=8) :: auxc(6)
complex (kind=8) :: aux1c,aux2c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)

real (kind=8) :: ar, br, cr
complex (kind=8) :: tmp_ke, tmp_pot, tmp_tot

!include 'interface_derivnD.include'  ! Per fer servir les derivades generiques

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


!........................................
!.. Laplacian of HPSi
!........................................
!

do jrun=1,4
    ar = arun(jrun)
    br = brun(jrun)
    cr = crun(jrun)

    if(ldroplet_frozen)then
        !Nothing to do
    else
        Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
        Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
        Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

        !
        !   We compute H·Psi
        !

        !$omp parallel do private(ix,iy,iz, tmp_ke, tmp_pot) default(shared) collapse(2) schedule(static)
        do iz=1,nz; do iy=1,ny
            do ix=1,nx
                tmp_ke = h2o2m4*(sto1c(ix,iy,iz) + sto2c(ix,iy,iz) + sto3c(ix,iy,iz))
                tmp_pot = pot4(ix,iy,iz)*psi(ix,iy,iz)
                sto4c(ix,iy,iz) = timec(ix,iy,iz)*(tmp_ke - tmp_pot) - ci*uimp(ix,iy,iz)*psi(ix,iy,iz)
                sto1c(ix,iy,iz) = ar*(sto4c(ix,iy,iz) - br*q(ix,iy,iz))
                q(ix,iy,iz) = q(ix,iy,iz) + 3.0*sto1c(ix,iy,iz) - cr*sto4c(ix,iy,iz)
            enddo
        enddo; enddo
        !$omp end parallel do

        if (jrun.eq.1) then
            !$omp parallel do private(ix,iy,iz) collapse(2) default(shared) schedule(static)
            do iz=1,nz; do iy=1,ny
                do ix=1,nx
                    hpsiold(ix,iy,iz,2) = hpsiold(ix,iy,iz,1)
                    hpsiold(ix,iy,iz,1) = sto4c(ix,iy,iz)
                    psiold(ix,iy,iz,3) = psiold(ix,iy,iz,2)
                    psiold(ix,iy,iz,2) = psiold(ix,iy,iz,1)
                    psiold(ix,iy,iz,1) = psi(ix,iy,iz)
                enddo
            enddo; enddo
            !$omp end parallel do
        endif


        !$omp parallel do private(ix,iy,iz) default(shared) schedule(static) collapse(2)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    psi(ix,iy,iz) = psi(ix,iy,iz) + deltat*sto1c(ix,iy,iz)
                    den(ix,iy,iz) = real(psi(ix,iy,iz))**2 + aimag(psi(ix,iy,iz))**2
                enddo
            enddo
        enddo
        !$omp end parallel do
    endif

    if(Lcoalescence ) then
    ! Write(*,*) "Coalescence between droplets, no impurit/ies"
    else
        !
        ! Impurity evolution if it is necessary
        !
        ! Positions
        stor = ar*(vimp - br*qr)
        qr = qr + 3.*stor - cr*vimp
        if(jrun.eq.1)then
            rimpold(:,:,3) = rimpold(:,:,2)
            rimpold(:,:,2) = rimpold(:,:,1)
            rimpold(:,:,1) = rimp
        endif
        rimp = rimp + deltat*stor

        ! Velocities
        Stor = ar*(aimp - br*qv)
         qv = qv + 3.*stor - cr*aimp
         if(jrun.eq.1)then
            aimpold(:,:,2) = aimpold(:,:,1)
            aimpold(:,:,1) =    aimp
            vimpold(:,:,3) = vimpold(:,:,2)
            vimpold(:,:,2) = vimpold(:,:,1)
            vimpold(:,:,1) =    vimp
         endif
         vimp = vimp + deltat*stor

        if(jrun.le.3)then
            call potenimp()
            call poten()
            call forceimp()
        endif
    endif
enddo



do ix=1,3
    ioldp(ix)=ix
    ioldr(ix)=ix
    ioldv(ix)=ix
enddo

do ix=1,2
    ioldh(ix)=ix
    iolda(ix)=ix
enddo

return
end
