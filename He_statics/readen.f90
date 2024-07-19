!------------------------------------------------------------------
!---                    Subroutine readen                       ---
!------------------------------------------------------------------

subroutine readen(npart,densat,fileden,fileimp,mode,rimpur,r_clust)

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

use he4
use rho
use field
use grid
use impur
use util1

implicit none
integer   (kind=4), intent(in)    :: npart
real      (kind=8)                :: nph     ! Num. part in hole
real      (kind=8), intent(in)    :: densat
character (len=60), intent(in)    :: fileden,fileimp
integer   (kind=4), intent(in)    :: mode
real      (kind=8), intent(out)   :: r_clust
real      (kind=8), intent(in)    :: rimpur

real      (kind=8) :: xmaxp,ymaxp,zmaxp,xcp,ycp,zcp,hxp,hyp,hzp
real      (kind=8) :: ximpp,yimpp,zimpp
real      (kind=8) :: aux1, aux2, aux3,aux3_cy=1d0
real      (kind=8) :: aux1b,aux2b,aux3b,sto,r_shell2,sum_to_check=0d0
real      (kind=8) , allocatable :: dendum(:,:,:)
real      (kind=8)  , allocatable :: density_initial(:,:,:)
real      (kind=8) , allocatable :: aux3_array(:), aux3b_array(:)
integer   (kind=4) :: nxp,nyp,nzp
integer   (kind=4) :: ix,iy,iz,isalto,in


allocate(dendum(nx,ny,nz))
allocate (density_initial(nx,ny,nz))

!...........................................
!With 'mode' select the kind of density...
!...........................................
!

select case(mode)
!-------------------------------------------------------------------
  case(0)   ! Continue a calculation (read a previous density)
!-------------------------------------------------------------------

     open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
!
     dendum=0.0d0         ! Be sure that all the components are initialized
!                         ! to zero.
!
!... Read the density and the wave function if it is the case
!
     read(1,*) dendum
     close(1)
!... Read the wave function if (limp=true) (calculation with impurity)
     if(limp) then 
        open(unit=1,file=fileimp,status='old')
        call titols(1,cchar,isalto)
        read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
        read(1,*) psix
        denx=psix**2
        close(1)
      end if
!-------------------------------------------------------------
   case(1)      ! Build a Fermi density for a pure substance
!-------------------------------------------------------------
      If(.Not.Lbulk)Then
        r_clust = ((3.0d0*npart)/(fourpi*densat))**0.3333333333333d0
        rfermi  = r_clust
      Endif
      do ix=1,nx
        aux1  =  (x(ix)-xc)**2
        do iy=1,ny
          aux2  = ( y(iy)-yc)**2
          do iz=1,nz
            aux3  = sqrt(aux1+aux2+(z(iz)-zc)**2)-rfermi
            dendum(ix,iy,iz)= densat/(dexp(aux3/afermi)+1)
          !write(6,*)afermi,dexp(aux3/afermi)+1,dendum(ix,iy,iz)
          end do
        end do
      end do
      If(.Not.Lbulk)Then
        aux1   = npart/(sum(dendum)*dxyz)
        dendum = aux1*dendum
      Endif
!---------------------------------------------------------------
   case(2)      ! Build a Fermi density for a dopped substance
!---------------------------------------------------------------
      If(.Not.Lbulk)Then
        nph     = (fourpi*densat/3.0d0)*rimpur**3
        r_clust = ((3.0d0*(npart+nph))/(fourpi*densat))**0.3333333333333d0
        rfermi  = r_clust
      Endif

      r_shell2=r_shell**2
      do iz=1,nz
        aux1  = defz*(z(iz)-zc)**2
        aux1b = ((z(iz)-zimp)**2)
        do iy=1,ny
          aux2  = defy*(y(iy)-yc)**2
          aux2b = ((y(iy)-yimp)**2)
          do ix=1,nx
            aux3             = sqrt(aux1+aux2+defx*(x(ix)-xc)**2)-rfermi
            aux3b            = aux1b+aux2b+(x(ix)-ximp)**2
            dendum(ix,iy,iz) = densat/(dexp(aux3/afermi)+1)
!            psix(ix,iy,iz)   = max(psimin,exp(-gwf*Abs(aux3b-r_shell2)))
            psix(ix,iy,iz)   = exp( -gwf*(Dsqrt(aux3b)-r_shell)**2 )
          end do
        end do
      end do
      aux1   = sqrt(1.0d0/(sum(psix*psix)*dxyz)) ! Normalization constant (Imp)
      psix   = aux1*psix                ! First wave function for the impurity
      denx   = psix**2
      If(.Not.Lbulk)Then
        aux1   = npart/(sum(dendum)*dxyz)          ! Normalization constant (He)
        dendum = aux1*dendum              ! First helium density
      Endif
!---------------------------------------------------------------
   case(3)      ! Build a Fermi density for a dopped substance
!---------------------------------------------------------------
      If(.Not.Lbulk)Then
        nph     = (fourpi*densat/3.0d0)*rimpur**3
        r_clust = ((3.0d0*(npart+nph))/(fourpi*densat))**0.3333333333333d0
        rfermi  = r_clust
      Endif

      If(L_anell.And.L_esfera)Then
        Write(6,'("From readen: aquestes dues opcions son autoexcluyents....")')
        Stop 'From readen   001'
      Endif
      do iz=1,nz
        aux1  = defz*(z(iz)-zc)**2
        aux1b = ((z(iz)-zimp)**2)*defz
        do iy=1,ny
          aux2  = defy*(y(iy)-yc)**2
          aux2b = ((y(iy)-yimp)**2)*defy
          do ix=1,nx
            aux3             = sqrt(aux1+aux2+defx*(x(ix)-xc)**2)-rfermi
            aux3b            = Sqrt(aux1b+aux2b+defx*(x(ix)-ximp)**2)-rimpur
            Sto=1.0d0
            If(lrandom.And.Uext(ix,iy,iz).Lt.10.d0)Call Random_Number(sto)
            aux3 = densat/(dexp(aux3/afermi)+1)*Sto*(1.0d0-1.0d0/(1.0d0+Dexp(aux3b/afermi)))
            If(L_anell)Then
              aux3b  = ( Dsqrt( defx*(x(ix)-ximp)**2 + aux2b ) - r_anell)**2
              Sto = densat*Dexp(-aux1b/a_anell)*Dexp(-aux3b/a_anell)
              aux3 = aux3 + Sto
            EndIf
            If(L_esfera)Then
              aux3b  = ( Dsqrt( defx*(x(ix)-ximp)**2 + aux1b + aux2b ) - r_esfera)**2
              Sto = densat*Dexp(-aux3b/a_esfera)
              aux3 = aux3 + Sto
            EndIf
            dendum(ix,iy,iz) = aux3
          end do
        end do
      end do
      If(.Not.Lbulk)Then
        aux2   = npart/(sum(dendum)*dxyz)          ! Normalization constant (He)
        dendum = aux2*dendum              ! First helium density
      Endif
!-------------------------------------------------------------------
  case(4)   ! Continue a calculation but ignore ximp,yimp,zimp
!-------------------------------------------------------------------

     open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximpp,yimpp,zimpp
!
     dendum=0.0d0         ! Be sure that all the components are initialized
!                         ! to zero.
!
!... Read the density and the wave function if it is the case
!
     read(1,*) dendum
     close(1)
!... Read the wave function if (limp=true) (calculation with impurity)
     if(limp) then 
        open(unit=1,file=fileimp,status='old')
        call titols(1,cchar,isalto)
        read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
        read(1,*) psix
        denx=psix**2
        close(1)
      end if

!---------------------------------------------------------------
   case(5)      ! Build a Fermi density for a dopped substance EGA and NH 2022
!---------------------------------------------------------------
      If(.Not.Lbulk)Then
        nph     = (fourpi*densat/3.0d0)*rimpur**3
        r_clust = ((3.0d0*(npart+nph))/(fourpi*densat))**0.3333333333333d0
        rfermi  = r_clust
      Endif

      If(L_anell.And.L_esfera)Then
        Write(6,'("From readen: aquestes dues opcions son autoexcluyents....")')
        Stop 'From readen   001'
      Endif
!**************************************************************************************     
      write(6,*)"reading mode 5" 
      Open(Unit=56,File=filepure) !!! EGA and NH
      call titols(56,cchar,isalto)
      Read(56,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximpp,yimpp,zimpp
     ! Write(6,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximpp,yimpp,zimpp
      Read(56,*) density_initial
      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
           sum_to_check=sum_to_check+density_initial(ix,iy,iz)         
          enddo
        enddo
      enddo
      sum_to_check=  sum_to_check*hx*hyp*hzp
      Write(6,*) "Number of particles for Pure droplet converged",sum_to_check

 
!**************************************************************************************
      do iz=1,nz
        aux1  = defz*(z(iz)-zc)**2
        aux1b = ((z(iz)-zimp)**2)*defz
        do iy=1,ny
          aux2  = defy*(y(iy)-yc)**2
          aux2b = ((y(iy)-yimp)**2)*defy
          do ix=1,nx
            aux3             = sqrt(aux1+aux2+defx*(x(ix)-xc)**2)-rfermi
            aux3b            = Sqrt(aux1b+aux2b+defx*(x(ix)-ximp)**2)-rimpur
            Sto=1.0d0
            If(lrandom.And.Uext(ix,iy,iz).Lt.10.d0)Call Random_Number(sto)
            aux3 = density_initial(ix,iy,iz)*Sto*(1.0d0-1.0d0/(1.0d0+Dexp(aux3b/afermi)))  !!! EGA and NH
            If(L_anell)Then
              aux3b  = ( Dsqrt( defx*(x(ix)-ximp)**2 + aux2b ) - r_anell)**2
              Sto = densat*Dexp(-aux1b/a_anell)*Dexp(-aux3b/a_anell)
              aux3 = aux3 + Sto
            EndIf
            If(L_esfera)Then
              aux3b  = ( Dsqrt( defx*(x(ix)-ximp)**2 + aux1b + aux2b ) - r_esfera)**2
              Sto = densat*Dexp(-aux3b/a_esfera)
              aux3 = aux3 + Sto
            EndIf
            dendum(ix,iy,iz) = aux3
          end do
        end do
      end do
      If(.Not.Lbulk)Then
        aux2   = npart/(sum(dendum)*dxyz)          ! Normalization constant (He)
        dendum = aux2*dendum              ! First helium density
      Endif

!  write(6,*)"finishing mode 5" 
!

!---------------------------------------------------------------
   case(6)      ! Build a Fermi density for multi-impurity dopped substance EGA and NH 2022
!---------------------------------------------------------------
      If(.Not.Lbulk)Then
        nph     = (fourpi*densat/3.0d0)*rimpur**3
        r_clust = ((3.0d0*(npart+nph))/(fourpi*densat))**0.3333333333333d0
        rfermi  = r_clust
      Endif

      If(L_anell.And.L_esfera)Then
        Write(6,'("From readen: aquestes dues opcions son autoexcluyents....")')
        Stop 'From readen   001'
      Endif
!**************************************************************************************     
      write(6,*)"reading mode 6" 
      Open(Unit=56,File=filepure) !!! EGA and NH
      call titols(56,cchar,isalto)
      Read(56,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximpp,yimpp,zimpp
     ! Write(6,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximpp,yimpp,zimpp
      Read(56,*) density_initial
      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
           sum_to_check=sum_to_check+density_initial(ix,iy,iz)         
          enddo
        enddo
      enddo
      sum_to_check=  sum_to_check*hx*hyp*hzp
      Write(6,*) "Number of particles for Pure droplet converged",sum_to_check

     Allocate(aux3_array(N_imp))
     Allocate(aux3b_array(N_imp))
!**************************************************************************************
      do iz=1,nz
        !aux1  = defz*(z(iz)-zc)**2
        !aux1b = ((z(iz)-zimp)**2)*defz
        do iy=1,ny
         ! aux2  = defy*(y(iy)-yc)**2
         ! aux2b = ((y(iy)-yimp)**2)*defy
          do ix=1,nx
             do in=1,N_imp
            aux1  = defz*(z(iz)-zc)**2
            aux1b = ((z(iz)-rimp(in,3))**2)*defz
             aux2  = defy*(y(iy)-yc)**2
            aux2b = ((y(iy)-rimp(in,2))**2)*defy
           
          
            aux3             = sqrt(aux1+aux2+defx*(x(ix)-xc)**2)-rfermi
            aux3b            = Sqrt(aux1b+aux2b+defx*(x(ix)-rimp(in,1))**2)-rimpur

               aux3_array(in)= aux3
               aux3b_array(in)= aux3b
  
            enddo

            Sto=1.0d0
            If(lrandom.And.Uext(ix,iy,iz).Lt.10.d0)Call Random_Number(sto)

               !!! EGA and NH
            aux3_cy=1d0
            do in=1,N_imp

              aux3_cy=aux3_cy*(1.0d0-1.0d0/(1.0d0+Dexp(aux3b_array(in)/afermi)))
             enddo
             
              aux3 = density_initial(ix,iy,iz)*Sto*aux3_cy
          
            If(L_anell)Then
              aux3b  = ( Dsqrt( defx*(x(ix)-ximp)**2 + aux2b ) - r_anell)**2
              Sto = densat*Dexp(-aux1b/a_anell)*Dexp(-aux3b/a_anell)
              aux3 = aux3 + Sto
            EndIf
            If(L_esfera)Then
              aux3b  = ( Dsqrt( defx*(x(ix)-ximp)**2 + aux1b + aux2b ) - r_esfera)**2
              Sto = densat*Dexp(-aux3b/a_esfera)
              aux3 = aux3 + Sto
            EndIf
            dendum(ix,iy,iz) = aux3
          end do
        end do
      end do
      If(.Not.Lbulk)Then
        aux2   = npart/(sum(dendum)*dxyz)          ! Normalization constant (He)
        dendum = aux2*dendum              ! First helium density
      Endif

!  write(6,*)"finishing mode 5" 
!

!-------------------------------------------------
   case default ! Non programmed situation
!-------------------------------------------------
!
      stop 'ERROR READEN: This mode is still not programmed'
end select

den = Abs(dendum)
deallocate(dendum)
return
end
