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
! use impur
use classicimp
use util1
use coalescence

implicit none
integer   (kind=4), intent(in)    :: npart
real      (kind=8)                :: nph     ! Num. part in hole
real      (kind=8), intent(in)    :: densat
character (len=60), intent(in)    :: fileden,fileimp
integer   (kind=4), intent(in)    :: mode

real      (kind=8) :: xmaxp,ymaxp,zmaxp,xcp,ycp,zcp,hxp,hyp,hzp,ximp,yimp,zimp
real      (kind=8) :: aux1, aux2, aux3
real      (kind=8) :: aux1b,aux2b,aux3b
! real      (kind=8) :: ximp,yimp,zimp

integer   (kind=4) :: nxp,nyp,nzp
integer   (kind=4) :: ix,iy,iz,isalto,N_imp0
!logical            :: limp,Ldensity=.true.
logical            :: limp
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
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
     If(Ldensity)Then
       read(1,*,Err=10) den
       Write(6,'("From Readenc: We have read a density")')
     Else
       read(1,*)psi
       Write(6,'("From Readenc: We have read a complex w.f.")')
       den=Abs(psi)**2
     Endif      
     close(1)

!---------------------------------------------------------------
 case(2)   ! Continue a calculation (read a previous wave functon)
!---------------------------------------------------------------
   open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp

  if(Lcoalescence ) then
  ! Write(*,*) "Coalescence between droplets, no impurit/ies"
  else
	 read(1,*) N_imp0
	 if (N_imp0.ne.N_imp) then
		 write(6,'("STOP: Number of impurities are inconsistent:",2I8)') N_imp, N_imp0
		 stop
	 endif
     read(1,*) rimp
     read(1,*) vimp
   endif

     read(1,*) psi
     den=Conjg(psi)*psi
   close(1)

   !

 !---------------------------------------------------------------
   case(3)   ! Coalescence  (read two previous densities)
!---------------------------------------------------------------
    open(unit=1,file=filedenin_coales_1,status='old')
    call titols(1,cchar,isalto)
    read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
    Allocate(den_coales1(nxp,nyp,nzp))
    read(1,*) den_coales1   ! Read density
    close(1)


    open(unit=2,file=filedenin_coales_2,status='old')
    call titols(2,cchar,isalto)
    read(2,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
    Allocate(den_coales2(nxp,nyp,nzp))
    read(2,*) den_coales2   ! Read density
    close(2)

    Allocate(den_coales_total(nxp,nyp,nzp))
    Allocate(xcoa(nxp))
    Allocate(ycoa(nyp))
    Allocate(zcoa(nzp))

   do ix=1,nxp  !.............Grid X
      xcoa(ix) = -xmaxp+hxp*(ix-1)
     end do
   do iy=1,nyp  !............ Grid Y
     ycoa(iy) = -ymaxp+hyp*(iy-1)
     end do
   do iz=1,nzp  !............ Grid  Z
     zcoa(iz) = -zmaxp+hzp*(iz-1)
    end do

   do iz=1,nzp
     do iy=1,nyp
       do ix=1,nxp
   psi(ix,iy,iz)=CMPLX(cos(zcoa(iz)*k_coales)*dsqrt(den_coales1(ix,iy,iz))+cos(zcoa(iz)*k_coales)*dsqrt(den_coales2(ix,iy,iz) ) , &
sin(zcoa(iz)*k_coales)*dsqrt(den_coales1(ix,iy,iz))-sin(zcoa(iz)*k_coales)*dsqrt(den_coales2(ix,iy,iz))  )


         enddo
        enddo
       enddo

     den=Conjg(psi)*psi

     Open(Unit=11,File="Den-Total_at_the_begining",status="old")
     Write(11,*)xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
     Write(11,*)den
     close(11)

!---------------------------------------------------------------
 case(7)   ! Continue a calculation (Frozen-dynamics)
!---------------------------------------------------------------
   open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp

  if(Lcoalescence ) then
  ! Write(*,*) "Coalescence between droplets, no impurit/ies"
  else
   read(1,*) N_imp0
   if (N_imp0.ne.N_imp) then
     write(6,'("STOP: Number of impurities are inconsistent:",2I8)') N_imp, N_imp0
     stop
   endif
     read(1,*) rimp
     read(1,*) vimp
   endif
 
   close(1)

 open(unit=2,file=file_density_input_frozen,status="old")
 call titols(2,cchar,isalto)
 read(2,*)xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
 read(2,*)den

 close(2)


   !


!-------------------------------------------------
   case default ! Non programmed situation
!-------------------------------------------------
!
      stop 'ERROR READEN: This mode is still not programmed'
end select

return
end
