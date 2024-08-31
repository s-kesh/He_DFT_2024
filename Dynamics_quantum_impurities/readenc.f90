!------------------------------------------------------------------
!---                    Subroutine readenc                      ---
!------------------------------------------------------------------

subroutine readenc(npart,densat,fileden,fileimp,mode)

! INPUT QUANTITIES:
!
! From list:
!
! npart   ----> Number of particles
! densat  ----> Density of saturation (useful for fermi distributions)
! fileden ----> Name of the file with the density
! mode    ----> Select the density:
!                       0 = continue a calculation
!                       1 = Build from scratch a new density
!                       2 = Build from scratch a new density with impurity
! rimpur  ----> Radius of the impurity
!
!
! OUTPUT QUANTITIES:
!
!  psi    ----> He Wave function)
!  den    ----> Density (module rho)
!  r_clust ---> Size of the cluster
!  psix ------> Wave function for the impurity
!
!
! NOTE:
! This subroutine check the consistency of namelist input data when an 
! external density is used as a input. The quantities to check consistency
! came from modules
!-------------------------------------------------------------------------------

use rho
use field
use grid
use impur
use quantal_imp
use util1

implicit none
integer   (kind=4), intent(in)    :: npart
real      (kind=8)                :: nph     ! Num. part in hole
real      (kind=8), intent(in)    :: densat
character (len=60), intent(in)    :: fileden
character (len=80), intent(in)    :: fileimp(*)
integer   (kind=4), intent(in)    :: mode

real      (kind=8) :: xmaxp,ymaxp,zmaxp,xcp,ycp,zcp,hxp,hyp,hzp,ximp,yimp,zimp
real      (kind=8) :: aux1, aux2, aux3,ximp0,yimp0,zimp0
real      (kind=8) :: aux1b,aux2b,aux3b
real      (kind=8),Allocatable :: dendump(:,:,:)
! real      (kind=8) :: ximp,yimp,zimp

integer   (kind=4) :: nxp,nyp,nzp
integer   (kind=4) :: ix,iy,iz,isalto,N_imp0
!logical            :: limp,Ldensity=.true.
!logical            :: limp
real      (kind=8) :: rr

!...........................................
!With 'mode' select the kind of density...
!...........................................
!

select case(mode)
!-------------------------------------------------------------------
  case(0)   ! Continue a calculation (read a previous density or wave function (vortex))
!-------------------------------------------------------------------

     open(unit=1,file=fileden,status='old')
     Go to 20
10   Continue
       Ldensity=.false.
20   Continue     
     Rewind(1)
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp0,yimp0,zimp0
     If(Ldensity)Then
       read(1,*,Err=10) den
       Write(6,'("From Readenc: We have read a density")')
     Else
       read(1,*)psi
       Write(6,'("From Readenc: We have read a complex w.f.")')
       den=Abs(psi)**2
     Endif      
     close(1)
     Allocate(dendump(nxp,nyp,nzp))
     denx = 0.d0
     Do ix=1, N_imp 
       open(unit=1,file=fileimp(ix),status='old')
       call titols(1,cchar,isalto)
       read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp0,yimp0,zimp0
       read(1,*)dendump
       psix(:,:,:,ix) = Cmplx(dendump)
       Write(6,'("W.f.(",I2,") normalization...:",1p,E15.6)')ix, (Sum(dendump**2)*hxp*hyp*hzp)
       denx(:,:,:) = denx(:,:,:) + Conjg(psix(:,:,:,ix))*psix(:,:,:,ix)
       close(1)
     EndDo
     Deallocate(dendump)
!---------------------------------------------------------------
 case(2)   ! Continue a calculation (read a previous wave functon)
!---------------------------------------------------------------
   open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp
     read(1,*)psi
     den=Conjg(psi)*psi
     close(1)
     denx = 0.d0
     Do ix=1, N_imp 
       open(unit=1,file=fileimp(ix),status='old')
       call titols(1,cchar,isalto)
       read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp
       read(1,*)psix(:,:,:,ix)
       denx(:,:,:) = denx(:,:,:) + Conjg(psix(:,:,:,ix))*psix(:,:,:,ix)
       close(1)
     EndDo

   !
!-------------------------------------------------
   case default ! Non programmed situation
!-------------------------------------------------
!
      stop 'ERROR READEN: This mode is still not programmed'
end select

return
end
