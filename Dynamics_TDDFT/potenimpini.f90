subroutine potenimpini
use interpol !, only:DelInter,potpi,potdel,npot,vpi,delta
use grid
use classicimp!, only: rimp, lselection, N_imp
implicit none
integer (kind=4) :: i,k,j
real    (kind=8) :: r,rmax
real    (kind=8) :: r1,r2,Select_pot
character(len=8) :: fmt ! format descriptor
character(len=15),Allocatable ::poten_name(:)
Character (len=11) :: vpoten_name,number_char
character(len=80),Allocatable ::poten_name_imp(:)
Character (len=10) :: vpoten_name_imp
Character (len=3)  :: number_char_imp

!****************************************************************************************************
vpoten_name="Vpotential_"
Allocate(poten_name(N_imp))
if(N_imp.lt.10) fmt = '(I1.1)' ! an integer of width 1 with zeros at the left
if(N_imp.ge.10 .AND. N_imp.le.100 ) fmt = '(I2.1)' ! an integer of width 2 with zeros at the left
if(N_imp.ge.100 .AND. N_imp.le.999 ) fmt = '(I3.1)' ! an integer of width 3 with zeros at the left

do i=1,N_imp
  write(number_char,fmt) i  
  poten_name(i)=vpoten_name//trim(number_char)
  !!!Write(6,*)poten_name(i)
enddo
!*****************************************************************************************************


npot =100*max(nx,ny,nz)+1
rmax = 4.d0*max(xmax,ymax,zmax)
rmaxinterpol=rmax

DelInter = rmax/dfloat(npot-1)

!.....................................
! Interpolation for V_gs
!.....................................

allocate(potion(N_imp,npot))
do k=1,N_imp
  do i=1,npot
    r = dfloat(i-1)*DelInter
	potion(k,i) = Select_pot(selec_gs_k(k),r,r_cutoff_gs_k(k),umax_gs_k(k)) 
  enddo
enddo

do k=1,N_imp
  Open(unit=1229+k,File=poten_name(k))

 do i=1,npot
   r = dfloat(i-1)*DelInter
   Write(1229+k,8971) r, potion(k,i)
 enddo

  Close(1229+k)
enddo
call updatepoten()

!****************************************************************************************************
vpoten_name_imp="V_Imp-Imp_"
Allocate(poten_name_imp(N_imp))
do i=1,N_imp
  write(number_char_imp,fmt) i  
  poten_name_imp(i)=vpoten_name_imp//trim(number_char_imp)
  ! Write(6,*)poten_name_imp(i),number_char_imp
enddo
!*****************************************************************************************************

 do k=1,N_imp
   Open(unit=1230+k,File=poten_name_imp(k))
  do j=1,N_imp
      r=0.1d0
     Write(1230+k,*) "#Potential interaction",k,"->",j
   do i=1,2000
   if(k.ne.j) Write(1230+k,8971) r, Select_pot(selec_gs_k_k(k,j),r,r_cutoff_gs_k_k(k,j),umax_gs_k_k(k,j))
  
     r=r+0.01d0
 enddo
   Write(1230+k,*) " "
   Write(1230+k,*) " "
enddo
   Close(1230+k)
enddo

!*************************




!!!Format
8971 FORMAT(1X,2(E15.8,1X))  

end subroutine potenimpini

