      Use seleccio_de_potencial 
      Implicit Real*8(A-H,O-Z)
      Character  (Len=80), Allocatable :: selec_gs_k(:)
      Real       (Kind=8), Allocatable :: r_cutoff_gs_k(:)
      Real       (Kind=8), Allocatable :: umax_gs_k(:)
      Real       (Kind=8), Allocatable :: rimp(:,:)
      Real       (Kind=8), Allocatable :: x(:),y(:),z(:)
      Real       (Kind=8), Allocatable :: uext(:,:,:)
      Character  (Len=80) :: selec_Ar_Ar='Ar_Ar'
      Real       (Kind=8) :: r_cutoff_Ar_Ar=2.5d0
      Real       (Kind=8) :: umax_Ar_Ar=1.224418d4
      Character  (Len=80) :: File4='imp.input'
      Character  (Len=80) :: File5='potential_multi_impurity.input'
      Character  (Len=80) :: File6='potential_multi_impurity.res'
      Character  (Len=80) :: File7='potential_multi_impurity.out'
      logical             :: limp=.false.
!      Data Fac_x_y_z/0.529d0/    ! Bohr radius
      Data Fac_x_y_z/1.d00/
      Namelist/Input/nx,ny,nz,xmax,ymax,zmax,N_imp,File4,File5,File6,File7,selec_Ar_Ar
      Namelist/imp/rimp,selec_gs_k,umax_gs_k,r_cutoff_gs_k,Fac_x_y_z,limp
      Open(Unit=5,File=File5)
      Read(5,Input)
      Open(Unit=4,File=File4)
      Open(Unit=6,File=File6)
      Open(Unit=7,File=File7)
      Write(6,Input)
      Allocate(rimp(N_imp,3))
      Allocate(r_cutoff_gs_k(N_imp))
      Allocate(umax_gs_k(N_imp))
      Allocate(selec_gs_k(N_imp))
      Allocate(uext(nx,ny,nz))
      Allocate(x(nx))
      Allocate(y(ny))
      Allocate(z(nz))
      Read(4,imp)
      Write(6,imp)
      hx=2.d0*xmax/nx
      hy=2.d0*ymax/ny
      hz=2.d0*zmax/nz
      rimp = rimp*fac_x_y_z
      Do ix=1, nx
        x(ix) = -xmax + (ix-1)*hx
      EndDo  
      Do iy=1, ny
        y(iy) = -ymax + (iy-1)*hy
      EndDo  
      Do iz=1, nz
        z(iz) = -zmax + (iz-1)*hz
      EndDo
!  
      Write(6,'("#")')
      Write(6,'("# Multi impurity potential, the following parameters had been used:")')
      Write(6,'("# N_imp ...........:",I4)')N_imp
      Do k=1,N_imp
        Write(6,'("# rimp(",I1,").............:",1p,3E15.6)')k,(rimp(k,i),i=1,3)
        Write(6,'("# selec_gs_k(",I1,").......:",A)')k,selec_gs_k(k)
        Write(6,'("# umax_gs_k(",I1,"),........:",1p,E15.6)')k,umax_gs_k(k)
        Write(6,'("# r_cutoff_gs_k(",I1,")....:",1p,E15.6)')k,r_cutoff_gs_k(k)
      EndDo  
      Write(6,'("#")')
      Write(6,*)nx,ny,nz,hx,hy,hz
      Call Flush(6)
!      
      Write(7,'("#")')
      Write(7,'("# Multi impurity potential, the following parameters had been used:")')
      Write(7,'("# N_imp............:",I4)')N_imp
      Do k=1,N_imp
        Write(7,'("# rimp(",I1,").............:",1p,3E15.6)')k,(rimp(k,i),i=1,3)
        Write(7,'("# selec_gs_k(",I1,").......:",A)')k,selec_gs_k(k)
        Write(7,'("# umax_gs_k(",I1,"),........:",1p,E15.6)')k,umax_gs_k(k)
        Write(7,'("# r_cutoff_gs_k(",I1,")....:",1p,E15.6)')k,r_cutoff_gs_k(k)
      EndDo  
      Write(7,'("#")')
      Write(7,*)xmax,ymax,zmax
      Write(7,*)hx,hy,hz
      Write(7,*)nx,ny,nz,limp
      Write(7,*)rimp(1,1),rimp(1,2),rimp(1,3)
      Do iz=1, nz
        Do iy=1, ny
          Do ix=1, nx
            Sto = 0.d0
            Do k=1, N_imp
              r=Dsqrt( (x(ix)-rimp(k,1))**2 + &
                       (y(iy)-rimp(k,2))**2 + &
                       (z(iz)-rimp(k,3))**2    )
!              Sto = Sto +  Select_pot(selec_gs_k(k),r,r_cutoff_gs_k(k),umax_gs_k(k))
              selec_gs=selec_gs_k(k)
              r_cutoff_gs=r_cutoff_gs_k(k)
              umax_gs=umax_gs_k(k)
              Sto = Sto +  V_gs(r)
            EndDo
            uext(ix,iy,iz)=Sto
          EndDo  
        EndDo  
      EndDo  
      Write(7,*)uext
      Close(7)
!           
!     We will printout the potential values along the main axis
!           
      Open(Unit=1,File="Uext-x.dat")
      Do ix=1, nx
        Write(1,'(F10.2,1p,E15.6)')x(ix),uext(ix,ny/2+1,nz/2+1)
      EndDo  
      Close(1)
      Open(Unit=1,File="Uext-y.dat")
      Do iy=1, ny
        Write(1,'(F10.2,1p,E15.6)')y(iy),uext(nx/2+1,iy,nz/2+1)
      EndDo  
      Close(1)
      Open(Unit=1,File="Uext-z.dat")
      Do iz=1, nz
        Write(1,'(F10.2,1p,E15.6)')z(iz),uext(nx/2+1,ny/2+1,iz)
      EndDo  
      Close(1)
!
!     We will compute the interditance between atoms, and the energy of
!     the cluster.
!
      Sto = 0.d0
      Do i=1, N_imp-1
        Do j=i+1, N_imp
          r=Dsqrt( (rimp(i,1)-rimp(j,1))**2 + &
                   (rimp(i,2)-rimp(j,2))**2 + &
                   (rimp(i,3)-rimp(j,3))**2    )
          Sto = Sto + Select_pot(selec_Ar_Ar,r,r_cutoff_Ar_Ar,umax_Ar_Ar)
          Write(6,'("d(",I1,",",I1,")........:",1p,E15.6)')i,j,r
        EndDo
      EndDo
      Write(6,'("Cluster energy....:",1p,E15.6)')Sto
      Stop
      End


