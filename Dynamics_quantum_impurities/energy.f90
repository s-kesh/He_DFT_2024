!............................................................
!...                      Subroutine energy               ...
!............................................................

!
! In order to have actualized arrays. It is very convient
! that before to call this subroutine calls the POTEN
! subroutine (for delj4)  and derden (for derivadores of 
! the helium4 density)subroutines


subroutine energy()

use alphasterm
use Para_DerivnD
use deriva
use energies
use field
use impur
use quantal_imp
use lenard4
use grid
use gridk
use he4
use rho
use util1
use work1

implicit none

real (kind=8)    ::  aux1,aux2,aux3,aux4,aux5 ! Auxiliar variables
complex (kind=8) :: invars(6)
integer (kind=4) :: ix,iy,iz,k,m,i_imp
real (kind=8) :: h2o2mimp

!Write(6,'("Entrem a energy")')
!Call Flush(6)

  Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(1,nn,hz,3,psi,sto3c,Icon)

!Write(6,'("Energy: Despres de calcular derivades")')
!Call Flush(6)

aux1   = 0.5d0*cp4
aux2   = (1.0d0/3.0d0)*cpp4

!..................................................................
!... Calculate the density of energy (kinetic, and correlation) ...
!..................................................................

ekin4 = 0.0d0
elj4  = 0.0d0
ecor4 = 0.0d0

do iz=1,nz
   do iy=1,ny
      do ix=1,nx
          aux3  = den(ix,iy,iz)
          aux4  = dencg(ix,iy,iz)
          aux5  = aux3*aux4**2*(aux1+aux2*aux4)
          ekin4 = ekin4 + abs(sto1c(ix,iy,iz))**2+abs(sto2c(ix,iy,iz))**2+abs(sto3c(ix,iy,iz))**2
          elj4  = elj4  + delj4(ix,iy,iz)*aux3   
          ecor4 = ecor4 + aux5
      end do
   end do
end do

!......................................................
!... Calculate the density of energy (alpha_s term) ...
!......................................................

select case(core4)
   case('OTE','OTC')
     ealphas = sum(falfs*( dxden*intxalf + dyden*intyalf + dzden*intzalf) )
     ealphas = -h2o2m4*0.5d0*alphas*ealphas*dxyz
   case default
     continue
end select

esolid=0.0d0
if(lsolid)esolid = C*sum(den*(1.d0+dtanh(beta*(den-den_m))))*dxyz

ekin4 = h2o2m4*ekin4*dxyz           ! TOTAL Kinetic energy for 4He
elj4  = 0.5d0 *elj4 *dxyz           ! TOTAL Lennard-Jones energy
ecor4 =        ecor4*dxyz           ! TOTAL Correlation energy for 4He
etot4 = ekin4+elj4+ecor4 + esolid   ! TOTAL ENERGY without impurity

select case(core4)
   case('OTE','OTC')
     etot4 = etot4+ealphas          ! TOTAL ENERGY including Alpha_s term
   case default
     continue
end select

!etot4    =  etot4 + sum(uext*den)*dxyz
!
! We will compute w.f. the kinetic energy and Sum(He-X_imp) emergies
!

 ekinx = 0d0
 wk2 = (0.d0,0.d0)
 do k=1, N_imp
   Call derivnD(1,nn,hx,1,psix(:,:,:,k),sto1c,Icon)
   Call derivnD(1,nn,hy,2,psix(:,:,:,k),sto2c,Icon)
   Call derivnD(1,nn,hz,3,psix(:,:,:,k),sto3c,Icon)
   h2o2mimp = h2o2m4*4.002603d0/m_imp_u(k)
   ekinx  = ekinx + h2o2mimp*Sum(Conjg(sto1c)*sto1c + Conjg(sto2c)*sto2c + Conjg(sto3c)*sto3c)

   Sto1(:,:,:) = psix(:,:,:,k)*Conjg(psix(:,:,:,k))
   call fftfw_1()
   wk2(:,:,:) = wk2(:,:,:) + wk1(:,:,:)*uimp_k(:,:,:,k)

 enddo

 ekinx = ekinx*dxyz
 call fftbk_2()
 eHeX = sum(sto2*den)*dxyz

 eimpu_impu = 0.d0
 do k=1,N_imp
   Sto1(:,:,:) = psix(:,:,:,k)*Conjg(psix(:,:,:,k))
   call fftfw_1()
   do m=k+1,N_imp
     wk2(:,:,:) = wk1(:,:,:)*uimp_k_k(:,:,:,k,m)
     call fftbk_2()
     eimpu_impu = eimpu_impu + Sum(Sto2(:,:,:)*psix(:,:,:,m)*Conjg(psix(:,:,:,m)))
   enddo
 enddo
 eimpu_impu = eimpu_impu*dxyz
 eimpu = ekinx + eHeX + eimpu_impu
 etot   =  etot4 + eimpu

!Write(6,'("Sortim de energy")')
!Call Flush(6)
return

end
