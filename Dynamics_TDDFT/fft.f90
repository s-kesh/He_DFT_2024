!...................................................................
!...                Subroutine fftini                            ...
!...................................................................

subroutine fftini(nx,ny,nz,nthread)

use rho
use lenard4
use alphasterm
use work1
use fftmodule  ! fin,fout,fftwplan,pfftfw,pfftbk,nthread,renor,npx,
               ! npy, npz
use utils
implicit none

integer   (kind=4) :: nx,ny,nz    ! Size of the grid along X, Y and Z axis
integer   (kind=4) :: ierr_init_threads
! added because needed to multi-thread FFTW and var nthread is no more in module
integer (kind=4)              :: nthread   ! Number of threads
!optim fft
integer iret
!fin optim
include 'fftw3.f.include'

! allocate(fin(nx,ny,nz)  )
! allocate(fout(nx/2+1,ny,nz) )

!optim fft
call dfftw_init_threads(iret)
Print*, "fftw init threads ",iret
call dfftw_plan_with_nthreads(nthread)
Print*, "FFTW Multi-threads with : ",nthread," threads"

npx   = nx
npy   = ny
npz   = nz
renor = 1.0d0/(nx*ny*nz)

!call dfftw_init_threads(ierr_init_threads)
!If(ierr_init_threads.Eq.0)Then
!  call dfftw_plan_with_nthreads(nthread)
!Else
!  Write(6,'("Problems with dfftw_init_threads")')
!Endif


! call dfftw_plan_dft_r2c_3d(pfftfw,nx,ny,nz,fin ,fout,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk,nx,ny,nz,fout,fin ,fftwplan)

! Transformadas forward:
call dfftw_plan_dft_r2c_3d(pfftfw_den,nx,ny,nz,den ,fden,fftwplan)
call dfftw_plan_dft_r2c_3d(pfftfw_1  ,nx,ny,nz,sto1,wk1 ,fftwplan)
call dfftw_plan_dft_r2c_3d(pfftfw_2  ,nx,ny,nz,sto2,wk2 ,fftwplan)
call dfftw_plan_dft_r2c_3d(pfftfw_3  ,nx,ny,nz,sto3,wk3 ,fftwplan)
!rutinas: una para den, otra para 1, y otra para 1+2+3.

! Transformadas backward:
! call dfftw_plan_dft_c2r_3d(pfftbk_cg,nx,ny,nz,wk1 ,dencg  ,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk_lj,nx,ny,nz,wk1 ,delj4  ,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk_as,nx,ny,nz,wk1 ,denalf ,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk_ua,nx,ny,nz,wk1 ,ualphas,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk_1 ,nx,ny,nz,wk1 ,sto1   ,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk_1x,nx,ny,nz,wk1 ,intxalf,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk_2y,nx,ny,nz,wk2 ,intyalf,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk_3z,nx,ny,nz,wk3 ,intzalf,fftwplan)

call dfftw_plan_dft_c2r_3d(pfftbk_cg,nx,ny,nz,wk1 ,dencg  ,fftwplan)
call dfftw_plan_dft_c2r_3d(pfftbk_lj,nx,ny,nz,wk1 ,delj4  ,fftwplan)
if (core4.eq.'OTC') then
call dfftw_plan_dft_c2r_3d(pfftbk_as,nx,ny,nz,wk1 ,denalf ,fftwplan)
call dfftw_plan_dft_c2r_3d(pfftbk_ua,nx,ny,nz,wk1 ,ualphas,fftwplan)
call dfftw_plan_dft_c2r_3d(pfftbk_1x,nx,ny,nz,wk1 ,intxalf,fftwplan)
call dfftw_plan_dft_c2r_3d(pfftbk_2y,nx,ny,nz,wk2 ,intyalf,fftwplan)
call dfftw_plan_dft_c2r_3d(pfftbk_3z,nx,ny,nz,wk3 ,intzalf,fftwplan)
endif
call dfftw_plan_dft_c2r_3d(pfftbk_1 ,nx,ny,nz,wk1 ,sto1   ,fftwplan)

!call dfftw_plan_dft_c2r_3d(pfftbk_1 ,nx,ny,nz,wk1 ,sto1   ,plan)
!call dfftw_plan_dft_c2r_3d(pfftbk_2 ,nx,ny,nz,wk2 ,sto2   ,plan)


! ............................................................................
! call dfftw_plan_dft_r2c_3d(pfftfw1  ,nx,ny,nz,sto1,wk1 ,fftwplan)
! !
! call dfftw_plan_dft_c2r_3d(pfftbkdencg,nx,ny,nz,wk1 ,dencg,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbkdelj4,nx,ny,nz,wk1 ,delj4,fftwplan)
! call dfftw_plan_dft_c2r_3d(pfftbk1    ,nx,ny,nz,wk1 ,sto1 ,fftwplan)


return

end

! Rutinas forward:
subroutine fftfw_den()
use fftmodule
implicit none
 call dfftw_execute(pfftfw_den)
end subroutine

subroutine fftfw_1()
use fftmodule
implicit none
 call dfftw_execute(pfftfw_1)
end subroutine

subroutine fftfw_123()
use fftmodule
implicit none
 call dfftw_execute(pfftfw_1)
 call dfftw_execute(pfftfw_2)
 call dfftw_execute(pfftfw_3)
end subroutine


! Rutinas backward:
subroutine fftbk_cg()
use rho , only: dencg
use fftmodule
use utils, only: dscal
implicit none
 call dfftw_execute(pfftbk_cg)
 call dscal(npx,npy,npz, renor, dencg)
end subroutine

subroutine fftbk_lj()
use lenard4, only: delj4
use fftmodule
use utils, only: dscal
implicit none
 call dfftw_execute(pfftbk_lj)
 call dscal(npx,npy,npz, renor, delj4)
end subroutine

subroutine fftbk_1()
use work1 , only:sto1
use fftmodule
use utils, only: dscal
implicit none
 call dfftw_execute(pfftbk_1)
 call dscal(npx,npy,npz, renor, sto1)
end subroutine

subroutine fftbk_as()
use alphasterm, only:denalf
use fftmodule
use utils, only: dscal
implicit none
 call dfftw_execute(pfftbk_as)
 call dscal(npx,npy,npz, renor, denalf)
end subroutine

subroutine fftbk_ua()
use alphasterm, only:ualphas
use fftmodule
use utils, only: dscal
implicit none
 call dfftw_execute(pfftbk_ua)
 call dscal(npx,npy,npz, renor, ualphas)
end subroutine

subroutine fftbk_xyz()
use alphasterm
use fftmodule
use utils, only: dscal
implicit none
 call dfftw_execute(pfftbk_1x)
 call dfftw_execute(pfftbk_2y)
 call dfftw_execute(pfftbk_3z)
 call dscal(npx,npy,npz, renor, intxalf)
 call dscal(npx,npy,npz, renor, intyalf)
 call dscal(npx,npy,npz, renor, intzalf)
end subroutine
