module seleccio_de_potencial
      Character (Len=80) :: selec_gs='Ba_plus_gs'
      Real      (Kind=8) :: r_cutoff_gs=1.05d0
      Real      (Kind=8) :: umax_gs=6.836177d3
      Character (Len=80) :: selec_pi='Ba_plus_2D1'
      Real      (Kind=8) :: r_cutoff_pi=1.05d0
      Real      (Kind=8) :: umax_pi=5.216077d4
      Character (Len=80) :: selec_sigma='Ba_plus_2D0'
      Real      (Kind=8) :: r_cutoff_sigma=1.05d0
      Real      (Kind=8) :: umax_sigma=9.674836d4
      Character (Len=80) :: selec_delta='Ba_plus_2D2'
      Real      (Kind=8) :: r_cutoff_delta=1.05d0
      Real      (Kind=8) :: umax_delta=6.260015d4
end module seleccio_de_potencial

Module Modificacio_De_Select_Pot

Integer*4 Npg_Ba_plus_gs_fixC4, Npg_Ba_plus_pi_fixC4, Npg_Ba_plus_sigma_fixC4

Parameter(Npg_Ba_plus_gs_fixC4=13)
Parameter(Npg_Ba_plus_pi_fixC4=13)
Parameter(Npg_Ba_plus_sigma_fixC4=13)

Real*16 Pg_Ba_plus_gs_fixC4
Real*16 Pg_Ba_plus_pi_fixC4
Real*16 Pg_Ba_plus_sigma_fixC4

Integer*4 k0_Ba_plus_gs_fixC4
Integer*4 k0_Ba_plus_pi_fixC4
Integer*4 k0_Ba_plus_sigma_fixC4

Real*8 rcutoff_Ba_plus_gs_fixC4, fcutoff_Ba_plus_gs_fixC4,r0_Ba_plus_gs_fixC4,aa_Ba_plus_gs_fixC4,  &
       bb_Ba_plus_gs_fixC4,cc_Ba_plus_gs_fixC4
Real*8 rcutoff_Ba_plus_pi_fixC4, fcutoff_Ba_plus_pi_fixC4,r0_Ba_plus_pi_fixC4,aa_Ba_plus_pi_fixC4,  &
       bb_Ba_plus_pi_fixC4,cc_Ba_plus_pi_fixC4
Real*8 rcutoff_Ba_plus_sigma_fixC4, fcutoff_Ba_plus_sigma_fixC4,r0_Ba_plus_sigma_fixC4,aa_Ba_plus_sigma_fixC4, &
       bb_Ba_plus_sigma_fixC4,cc_Ba_plus_sigma_fixC4
Common/Param_LJ_V_Ba_plus_gs_fixC4/rcutoff_Ba_plus_gs_fixC4,fcutoff_Ba_plus_gs_fixC4,               &
       r0_Ba_plus_gs_fixC4,aa_Ba_plus_gs_fixC4,bb_Ba_plus_gs_fixC4,cc_Ba_plus_gs_fixC4,             &
       Pg_Ba_plus_gs_fixC4(Npg_Ba_plus_gs_fixC4),k0_Ba_plus_gs_fixC4
Common/Param_LJ_V_Ba_plus_pi_fixC4/rcutoff_Ba_plus_pi_fixC4,fcutoff_Ba_plus_pi_fixC4,               &
       r0_Ba_plus_pi_fixC4,aa_Ba_plus_pi_fixC4,bb_Ba_plus_pi_fixC4,cc_Ba_plus_pi_fixC4,             &
       Pg_Ba_plus_pi_fixC4(Npg_Ba_plus_pi_fixC4),k0_Ba_plus_pi_fixC4
Common/Param_LJ_V_Ba_plus_sigma_fixC4/rcutoff_Ba_plus_sigma_fixC4,fcutoff_Ba_plus_sigma_fixC4,      &
       r0_Ba_plus_sigma_fixC4,aa_Ba_plus_sigma_fixC4,bb_Ba_plus_sigma_fixC4,cc_Ba_plus_sigma_fixC4, &
       Pg_Ba_plus_sigma_fixC4(Npg_Ba_plus_sigma_fixC4),k0_Ba_plus_sigma_fixC4

  Real (Kind=8) :: Quita_C4_Ba_plus_gs_fix_C4=0.0d0
  Real (Kind=8) :: Quita_C4_Ba_plus_pi_fix_C4=0.0d0
  Real (Kind=8) :: Quita_C4_Ba_plus_sigma_fix_C4=0.0d0

  Save

End Module Modificacio_De_Select_Pot

Double Precision Function  V_gs(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_gs=Select_pot(selec_gs,r,r_cutoff_gs,umax_gs)
  End Function V_gs
  Double Precision Function V_pi(r)
  use seleccio_de_potencial
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_pi=Select_pot(selec_pi,r,r_cutoff_pi,umax_pi)
  End Function V_pi
  Double Precision Function V_sigma(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_sigma=Select_pot(selec_sigma,r,r_cutoff_sigma,umax_sigma)
  End Function V_sigma
  Double Precision Function V_delta(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_delta=Select_pot(selec_delta,r,r_cutoff_delta,umax_delta)
  End Function V_delta
  Double Precision Function V_Sigma_Pi(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: r, V_Sigma, V_Pi
      V_Sigma_Pi=(V_Sigma(r) - V_Pi(r))/r**2
  End Function V_Sigma_Pi

Double Precision Function Select_Pot(selec,r,r_cutoff,umax)

Implicit None

Logical       :: Lcontrol=.false.        ! Per verificar el funcionament correcte dels nous potencials
Character (Len=80) :: selec
Real (Kind=8) :: r, r_cutoff, umax
Real (Kind=8) :: V_null=0.d0             ! selec='null'    Potencial de interacció nul
Real (Kind=8) :: V_Ar_Ar                 ! selec='Ar_Ar'          Potencial Xe-He Tang & Toennies J. Chem. Phys. 118, 4976 (2003)
Real (Kind=8) :: V_Xe_Xe                 ! selec='Xe_Xe'          Potencial Xe-He Tang & Toennies Z. Phys. D 1, 91-101 (1986)
Real (Kind=8) :: drV_Ar_Ar               ! selec='dr_Ar_Ar'       Potencial Xe-He Tang & Toennies J. Chem. Phys. 118, 4976 (2003)
Real (Kind=8) :: drV_Xe_Xe               ! selec='dr_Xe_Xe'       Potencial Xe-He Tang & Toennies Z. Phys. D 1, 91-101 (1986)
Real (Kind=8) :: V_Ar_He                 ! selec='Ar_He'          Potencial Xe-He Tang & Toennies J. Chem. Phys. 118, 4976 (2003)
Real (Kind=8) :: V_Sr_He                 ! selec='Sr_He'          Potencial Sr-He Lovallo et al. Jcp 120, (2004) 246, Version SR
Real (Kind=8) :: Au_Graphene             ! selec='Au_Graphene'    Potencial de interacció Au-Graphe
Real (Kind=8) :: Ag_Graphene             ! selec='Ag_Graphene'    Potencial de interacció Ag-Graphe
Real (Kind=8) :: V_Ag_Graphene_Amorf     ! selec='Ag_Graphene_Amorf'    Potencial de interacció Ag-Graphe
Real (Kind=8) :: dz_Au_Graphene          ! selec='dz_Au_Graphene' Gradient del Potencial de interacció Au-Graphe
Real (Kind=8) :: dz_Ag_Graphene          ! selec='dz_Ag_Graphene' Gradient del Potencial de interacció Ag-Graphe
Real (Kind=8) :: dzV_Ag_Graphene_Amorf   ! selec='dz_Ag_Graphene_Amorf' Gradient del Potencial de interacció Ag-Graphe
Real (Kind=8) :: V_LJ_OT                 ! selec='LJ_OT'   Potencial de Lenard-Jones capat a la Orsay-Trento
Real (Kind=8) :: V_Aziz_He               ! selec='Aziz_He' Potencial de Aziz pel He
Real (Kind=8) :: V_alka                  ! selec='LI', 'NA', 'K', 'RB' o 'CS' Potencials de Patil
Real (Kind=8) :: V_Na_Sigma              ! selec='Na_Sigma'          Pascale potential for Na(3Sigma)-He
Real (Kind=8) :: V_Na_Pi                 ! selec='Na_Pi'             Pascale potential for Na(3Pi)-He
Real (Kind=8) :: V_Cs_Cs                 ! selec='Cs_Cs'             Points obtained by Romain Vexiau, 2021
Real (Kind=8) :: V_Cs_Cs_tsu             ! selec='Cs_Cs_tsu'         Points obtained by Spis, triplet sigma u
Real (Kind=8) :: V_Rb_Rb                 ! selec='Rb_Rb'             Gregoire Guillon, Alexandra Viel, and Jean-Michel Launay !!!Marti send me this points (11/2021) 
Real (Kind=8) :: V_Cs_ssg                ! selec='Cs_ssg'      Spies's Potential   Singlet Sigma g 6S-6S
Real (Kind=8) :: V_Cs_tsu                ! selec='Cs_tsu'      Spies's Potencial   Triplet Sigma u
Real (Kind=8) :: V_Cs_ssg_ps             ! selec='Cs_ssg_ps'      Spies's Potential   Singlet Sigma g 6S-6P
Real (Kind=8) :: V_Cs_ssu_ps             ! selec='Cs_ssu_ps'      Spies's Potential   Singlet Sigma u 6S-6P
Real (Kind=8) :: V_Cs_spg_ps             ! selec='Cs_spg_ps'      Spies's Potential   Singlet Pi g 6S-6P
Real (Kind=8) :: V_Cs_spu_ps             ! selec='Cs_spu_ps'      Spies's Potential     Singlet Pi u 6S-6P
Real (Kind=8) :: V_Cs_tpg_ps             ! selec='Cs_tpg_ps'      Spies's Potential     Triplet Pi g 6S-6P
Real (Kind=8) :: V_Cs_tpu_ps             ! selec='Cs_tpu_ps'      Spies's Potential     Triplet Pi u 6S-6P
Real (Kind=8) :: V_Cs_tsg_ps             ! selec='Cs_tsg_ps'      Spies's Potential     Triplet Sigma g 6S-6P
Real (Kind=8) :: V_Cs_tsu_ps             ! selec='Cs_tsu_ps'      Spies's Potential     Triplet Sigma u 6S-6P
Real (Kind=8) :: V_Kerkines_Li_gs        ! selec='Kerkines_Li_gs'    Kerkines-Mavidris Potential for Li-He gs
Real (Kind=8) :: V_Kerkines_Li_Pi        ! selec='Kerkines_Li_Pi'    Kerkines-Mavidris Potential for Li-He 2Pi
Real (Kind=8) :: V_Kerkines_Li_Sigma     ! selec='Kerkines_Li_Sigma' Kerkines-Mavidris Potential for Li-He 2Sigma
Real (Kind=8) :: V_He2s_fci              ! selec='He2S_fci', Potencial He(2S) del Eloranta (fci) (5-11-2015), Nist: He(2S): 329179.7623 cm-1
Real (Kind=8) :: V_Eloranta_He2_plus_gs  ! selec='He_plus', Potencial He2+ de Eloranta (6-6-2016)
Real (Kind=8) :: V_Koutselos_Cs_plus_gs  ! selec='Cs_plus_Koutselos'
Real (Kind=8) :: V_Koutselos_K_plus_gs   ! selec='K_plus_Koutselos'
Real (Kind=8) :: V_Koutselos_Rb_plus_gs  ! selec='Rb_plus_Koutselos'
Real (Kind=8) :: V_Koutselos_Na_plus_gs  ! selec='Na_plus_Koutselos'
Real (Kind=8) :: V_Fausto_Rb_plus_gs     ! selec='Rb_plus_Fausto'    New: computed by Fausto
Real (Kind=8) :: V_Ca                    ! selec='Ca'          Mauschick and Meyer 1989
Real (Kind=8) :: V_Ba_gs                 ! selec='Ba_gs'          Lovallo potential
Real (Kind=8) :: V_Ba_plus_old_gs        ! selec='Ba_plus_old_gs'
Real (Kind=8) :: V_Ba_plus_gs_fixC4      ! selec='Ba_plus_gs_fix_C4'
Real (Kind=8) :: V_Ba_plus_gs            ! selec='Ba_plus_gs'
Real (Kind=8) :: V_Ba_plus_pi            ! selec='Ba_plus_pi'
Real (Kind=8) :: V_Ba_plus_pi_fixC4      ! selec='Ba_plus_pi_fix_C4'
Real (Kind=8) :: V_Ba_plus_sigma         ! selec='Ba_plus_sigma'
Real (Kind=8) :: V_Ba_plus_sigma_fixC4   ! selec='Ba_plus_sigma_fix_C4'
Real (Kind=8) :: V_Ba_plus_2D0           ! selec='Ba_plus_2D0'
Real (Kind=8) :: V_Ba_plus_2D1           ! selec='Ba_plus_2D1'
Real (Kind=8) :: V_Ba_plus_2D2           ! selec='Ba_plus_2D2'
Real (Kind=8) :: V_CH3I_0                ! selec='CH3I_0'
Real (Kind=8) :: V_CH3I_1                ! selec='CH3I_1'
Real (Kind=8) :: V_CH3I_2                ! selec='CH3I_2'
Real (Kind=8) :: V_CH3I_3                ! selec='CH3I_3'
Real (Kind=8) :: V_CH3I_4                ! selec='CH3I_4'
Real (Kind=8) :: V_Ag_gs                 ! selec='Ag_gs'
Real (Kind=8) :: V_Ag_Pi                 ! selec='Ag_Pi'
Real (Kind=8) :: V_Ag_Sig                ! selec='Ag_Sig'
Real (Kind=8) :: V_K_4p_Sigma            ! selec='K_4p_Sigma'     Potencial de Pascale
Real (Kind=8) :: V_K_4p_Pi               ! selec='K_4p_Pi'            "     "     "       
Real (Kind=8) :: V_K_5s                  ! selec='K_5s'               "     "     "       
Real (Kind=8) :: V_Cs_7s                 ! selec='Cs_7s'              "     "     "
Real (Kind=8) :: V_Rb_6s                 ! selec='Rb_6s'              "     "     "
Real (Kind=8) :: V_Rb_6p_Sigma           ! selec='Rb_6p_sigma'        "     "     "
Real (Kind=8) :: V_Rb_6p_Pi              ! selec='Rb_6p_pi'           "     "     "
Real (Kind=8) :: V_Cs_sigma              ! selec='Cs_sigma'           "     "     "
Real (Kind=8) :: V_Cs_pi                 ! selec='Cs_pi'              "     "     "
Real (Kind=8) :: V_Fausto_Cs_gs          ! selec='Cs_gs_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Cs_6p0         ! selec='Cs_6p_Sigma_Fausto' "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Cs_6p1         ! selec='Cs_6p_Pi_Fausto'    "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Cs_7s          ! selec='Cs_7s_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_Rb_5p_sigma           ! selec='Rb_sigma'           "     "  Pascale
Real (Kind=8) :: V_Rb_5p_pi              ! selec='Rb_pi'              "     "      "
Real (Kind=8) :: V_Fausto_Rb_gs          ! selec='Rb_gs_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Rb_5p0         ! selec='Rb_5p_Sigma_Fausto' "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Rb_5p1         ! selec='Rb_5p_Pi_Fausto'    "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Rb_6s          ! selec='Rb_6s_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_grafeno                   ! selec='Grafeno'
Real (Kind=8) :: V_grafeno_sin_dispersion    ! selec='Grafeno_sin_dispersion'
Real (Kind=8) :: V_TiO2                  ! selec='Tio2'
Real (Kind=8) :: V_Au                    ! selec='Au'
Real (Kind=8) :: V_Au_TiO2               ! selec='Au_TiO2'                     potencial Au-TiO2( Pilar )
Real (Kind=8) :: V_dAu_TiO2              ! selec='dAu_TiO2'       Derivada del potencial Au-TiO2( Pilar )
Real (Kind=8) :: V_pH2_He                ! selec='pH2_He'         Potencial pH2-He ajutadode Gianturco
Real (Kind=8) :: V_Xe_He                 ! selec='Xe_He'          Potencial Xe-He Tang & Toennies Z. Phys. D 1, 91-101 (1986)
Character (Len=3)  :: elem
!
! Aqui comença la seleccio del potencial
!
If(Trim(selec).Eq.'Aziz_He')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Aziz_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Aziz_He")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'null')Then
   If(r.lt.r_cutoff)Then
     Select_pot = V_null
   Else
     Select_pot = V_null
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_null")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'He2S_fci')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_He2s_fci(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_He2s_fci")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'He_plus')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Eloranta_He2_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Eloranta_He2_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'LJ_OT')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_LJ_OT(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_LJ_OT")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'LI'.Or.   &
       Trim(selec).Eq.'NA'.Or.   &
       Trim(selec).Eq.'K' .Or.   &
       Trim(selec).Eq.'RB'.Or.   &
       Trim(selec).Eq.'CS')Then
       If(r.lt.r_cutoff)Then
       Select_pot = umax
   Else
     If(Len_Trim(selec).Eq.1)elem=(Trim(selec)//'  ')
     If(Len_Trim(selec).Eq.2)elem=(Trim(selec)//' ')
     Select_pot = V_alka(r,elem)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_alka amb el parametre:",A3)')elem
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_plus_Koutselos')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_Cs_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_Cs_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Na_plus_Koutselos')Then
!
!  r_cutoff=1.0d0
!  umax=9.012825d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_Na_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_Na_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_plus_Koutselos')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_Rb_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_Rb_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_plus_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'K_plus_Koutselos')Then
!
!  r_cutoff=1.0d0, 2.0d0
!  umax=1.8628633412d5, 3.2762532479d3 
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_K_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_K_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_plus_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_old_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_old_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_old_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_gs_fix_C4')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_gs_fixC4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_gs_fixC4")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_pi_fix_C4')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_pi_fixC4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_pi_fixC4")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_sigma_fix_C4')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_sigma_fixC4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_sigma_fixC4")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_2D0')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_2D0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_2D0")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_2D1')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_2D1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_2D1")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_2D2')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_2D2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_2D2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ag_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ag_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ag_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ag_Pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ag_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ag_Pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ag_Sig')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ag_Sig(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ag_Sig")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_7s')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_7s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_7s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6s')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_6s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_6s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_gs_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_6p_Sigma_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_6p0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_6p0")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_6p_Pi_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_6p1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_6p1")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_7s_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_7s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_7s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_5p_sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_5p_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_5p_pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_5p_pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6p_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_6p_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_6p_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6p_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_6p_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_6p_Pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_gs_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_5p_Sigma_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_5p0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_5p0")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_5p_Pi_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_5p1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_5p1")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6s_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_6s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_6s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Grafeno')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_grafeno(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_grafeno")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Grafeno_sin_dispersion')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_grafeno_sin_dispersion(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_grafeno_sin_dispersion")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'TiO2')Then
!
!  r_cutoff=2.5d0
!  umax=1.887964E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_TiO2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Tio2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Au')Then
!
!  r_cutoff=1.d0
!  umax=6.645232E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Au(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_TAu")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Au_TiO2')Then
!
!  r_cutoff=2.d0
!  umax=21450.4867067095d0
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Au_TiO2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Au_TiO2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'dAu_TiO2')Then
!
!  r_cutoff=2.d0
!  umax=-6.2842194d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_dAu_TiO2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_dAu_TiO2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'pH2_He')Then
!
!  r_cutoff=
!  umax=-6.2842194d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_pH2_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_pH2_He")')
   Endif
ElseIf(Trim(selec).Eq.'Ag_Graphene')Then
!
!  r_cutoff= 1.5
!  umax= 6.1120297153d5
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = Ag_Graphene(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a Ag_Graphene")')
   Endif
ElseIf(Trim(selec).Eq.'dz_Ag_Graphene')Then
!
!  r_cutoff= 1.5
!  umax= -1.0533234197d6
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = dz_Ag_Graphene(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a dz_Ag_Graphene")')
   Endif
ElseIf(Trim(selec).Eq.'Xe_He')Then
!
!  r_cutoff= 2.0d0
!  umax= 21450.4867067095d0
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Xe_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Xe_He")')
   Endif
ElseIf(Trim(selec).Eq.'Ar_He')Then
!
!  r_cutoff= 2.0d0
!  umax= 8.0808095081d3 
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ar_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ar_He")')
   Endif
ElseIf(Trim(selec).Eq.'Au_Graphene')Then
!
!  r_cutoff= 2.3d0
!  umax= 4.9029319765d3 
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = Au_Graphene(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a Au_Grapehene")')
   Endif
ElseIf(Trim(selec).Eq.'dz_Au_Graphene')Then
!
!  r_cutoff= 2.3d0
!  umax= -2.3647589877E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = dz_Au_Graphene(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a dz_Au_Grapehene")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_0')Then
!
!  r_cutoff= 3.6d0
!  umax    = 3.1648655845d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_0(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_1')Then
!
!  r_cutoff= 3.6d0
!  umax    = -5.1115736616d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_1(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_2')Then
!
!  r_cutoff= 3.6d0
!  umax    = 5.5581502684d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_2(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_3')Then
!
!  r_cutoff= 3.6d0
!  umax    = -2.7398540818d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_3(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_3(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_4')Then
!
!  r_cutoff= 3.6d0
!  umax    = 2.4025334047d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_4(r)")')
   Endif
ElseIf(Trim(selec).Eq.'K_4p_Sigma')Then
!
!  r_cutoff= 2.2d0
!  umax    = 5.494421d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_K_4p_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_K_4p_Sigma(r)")')
   Endif
ElseIf(Trim(selec).Eq.'K_4p_Pi')Then
!
!  r_cutoff= 1.9d0
!  umax    = 5.348304d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_K_4p_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_K_4p_Pi(r)")')
   Endif
ElseIf(Trim(selec).Eq.'K_5s')Then
!
!  r_cutoff= 2.d0
!  umax    = 2.729808d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_K_5s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_K_5s(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Ag_Graphene_Amorf')Then
!
!  r_cutoff= 2.d0
!  umax    = 2.243286d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ag_Graphene_Amorf(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ag_Graphene_Amorf(r)")')
   Endif
ElseIf(Trim(selec).Eq.'dz_Ag_Graphene_Amorf')Then
!
!  r_cutoff= 2.2d0
!  umax    = -6.976782d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = dzV_Ag_Graphene_Amorf(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a dzV_Ag_Graphene_Amorf(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Ca')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ca(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ca(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Ar_Ar')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ar_Ar(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ar_Ar(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Xe_Xe')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Xe_Xe(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Xe_Xe(r)")')
   Endif
ElseIf(Trim(selec).Eq.'dr_Xe_Xe')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = drV_Xe_Xe(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a drV_Xe_Xe(r)")')
   Endif
ElseIf(Trim(selec).Eq.'dr_Ar_Ar')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = drV_Ar_Ar(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a drV_Ar_Ar(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Sr_He')Then
!
!  r_cutoff= 2.5d0
!  umax    = 6.423137d03
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Sr_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Sr_He(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Kerkines_Li_gs')Then
!
!  r_cutoff= 1.d0
!  umax    = 8.962011d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Kerkines_Li_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Kerkines_Li_gs(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Kerkines_Li_Pi')Then
!
!  r_cutoff= 1.d0
!  umax    = 3.453693d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Kerkines_Li_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Kerkines_Li_Pi(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Kerkines_Li_Sigma')Then
!
!  r_cutoff= 1.d0
!  umax    = 1.501675d5
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Kerkines_Li_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Kerkines_Li_Sigma(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Na_Sigma')Then
!
!  r_cutoff= 1.7d0
!  umax    = 1.0654129292E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Na_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Na_Sigma(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Na_Pi')Then
!
!  r_cutoff= 1.65d0
!  umax    = 1.012716E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Na_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Na_Pi(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Cs_Cs')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_Cs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_Cs(r)")')
   Endif
   ElseIf(Trim(selec).Eq.'Cs_Cs_tsu')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_Cs_tsu(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_Cs_tsu(r)")')
   Endif
    ElseIf(Trim(selec).Eq.'Rb_Rb')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_Rb(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_Rb(r)")')
   Endif
     ElseIf(Trim(selec).Eq.'Cs_ssg')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_ssg(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_ssg(r)")')
   Endif
   ElseIf(Trim(selec).Eq.'Cs_tsu')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_tsu(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_tsu(r)")')
   Endif
    ElseIf(Trim(selec).Eq.'Cs_ssg_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_ssg_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_ssg_ps(r)")')
   Endif
   ElseIf(Trim(selec).Eq.'Cs_ssu_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_ssu_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_ssu_ps(r)")')
   Endif
    ElseIf(Trim(selec).Eq.'Cs_spg_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_spg_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_spg_ps(r)")')
   Endif
    ElseIf(Trim(selec).Eq.'Cs_spu_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_spu_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_spu_ps(r)")')
   Endif
   ElseIf(Trim(selec).Eq.'Cs_tpg_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_tpg_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_tpg_ps(r)")')
   Endif
    ElseIf(Trim(selec).Eq.'Cs_tpu_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_tpu_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_tpu_ps(r)")')
   Endif
    ElseIf(Trim(selec).Eq.'Cs_tsg_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_tsg_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_tsg_ps(r)")')
   Endif
    ElseIf(Trim(selec).Eq.'Cs_tsu_ps')Then   
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_tsu_ps(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_tsu_ps(r)")')
   Endif
Else
   Write(6,'(A,": aquest potencial no està definit, de part de: Selec_Pot")')Trim(Selec)
   Stop 'From V_impur: 0001'
Endif
End Function Select_Pot
!
! Potencial Ag-Grapheno
! Mínimo del potencial: (3.4 Ang, -2.2424456319E+04 K)
!                       (2.0 Ang,  1.8372848084E+05 K)
!                       (1.5 Ang,  6.1120297153E+05 K)
!
Double Precision Function Ag_Graphene(z)
Implicit Real*8(A-H,O-Z)
Data A/11557.3D0/, beta/2.33671D0/, C4/1143.0D0/
Data Ckcal_to_K/5.03216592455d2/


Vz = A*dexp(-beta*z) - C4/z**4

Ag_Graphene = Vz * Ckcal_to_K
Return
End
!
! Gradient of Ag-Graphene potential
!
!   Some values:   (2   Ang, -5.5036249176E+05)
!                  (1.5  " , -1.0533234197E+06)
!
Double Precision Function dz_Ag_Graphene(z)
Implicit Real*8(A-H,O-Z)
Data A/11557.3D0/, beta/2.33671D0/, C4/1143.0D0/
Data Ckcal_to_K/5.03216592455d2/


dVz = -A*beta*dexp(-beta*z) + 4.d0*C4/z**5

dz_Ag_Graphene = dVz * Ckcal_to_K
Return
End
!
!  Potencial de Aziz pel He
!
double precision function V_Aziz_He(rr)
implicit none
real      (kind=8) :: Eps    = 10.948d0
real      (kind=8) :: A      = 1.8443101d5
real      (kind=8) :: alpha  = 10.43329537d0
real      (kind=8) :: beta   = -2.27965105d0
real      (kind=8) :: D      = 1.4826d0
real      (kind=8) :: C6     = 1.36745214d0
real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C10    = 0.17473318d0
real      (kind=8) :: rm     = 2.963d0
real (kind=8) :: rr,r,ff
 ff=1.d0
 r=rr/rm
 if(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 if(r==0)then
   V_Aziz_He= A*Eps
 else
   V_Aziz_He=(A*dexp( - alpha*r + beta*r**2) - ff*(C6/r**6 + C8/r**8 + C10/r**10 ))*Eps
 endif
end function

!
!  Potencial de Lenard-Jones capat a la Orsay-Trento
!
Block Data Inicio_LJ_V_LJ_OT
Implicit Real*8(A-H,O-Z)
Parameter(Npg=7)
!Real*16 Pg
Real*8 Pg
Integer*4 k0
Common/Param_LJ_V_LJ_OT/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
!Data Pg/ -1.139924227407542Q4, 0., 0., 0., 0., 0., 3.178638072971337Q6/
Data Pg/ -1.139924227407542D4, 0., 0., 0., 0., 0., 3.178638072971337D6/
Data rcutoff/ 2.190323254008688d0/
Data fcutoff /1.574659520268255d2/
Data r0/ 0.d0/
Data aa/ 0.d0/
Data bb/ 0.d0/
Data cc/ 0.d0/
Data k0/   5/
End
Double Precision Function V_LJ_OT(x)
Parameter(Npg=7)
Implicit Real*8(A-H,O-Z)
!Real*16 Pg,Sto,xq
Real*8 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_LJ_OT/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_LJ_OT=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_LJ_OT=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_LJ_OT=Sto
Return
End
!
!  Potencial de la Plata pel g.s.
!
double precision function V_Ag_gs(r)
implicit none
real      (kind=8) :: Ags    = 78741.3877
real      (kind=8) :: alphags= 0.65299
real      (kind=8) :: betags = 0.36367
real      (kind=8) :: Dgs    = 12.093
real      (kind=8) :: C6gs   = 137042.07d0
real      (kind=8) :: C8gs   = 189386.217d0
real      (kind=8) :: C10gs  = 111358462.d0
real      (kind=8) :: C12gs  = 1.17670377d10
real (kind=8) :: r,ff
 ff=1.d0
 if(r.le.Dgs)ff=dexp(-(Dgs/r-1.d0)**2)
 if(r==0)then
   V_Ag_gs= Ags
 else
   V_Ag_gs=Ags*dexp(-alphags*r-betags*r**2) - ff*(C6gs/r**6 + C8gs/r**8 + C10gs/r**10 + C12gs/r**12 )
 endif
end function

double precision function V_Ag_Pi(r)
implicit none
real      (kind=8) :: Api    = 0.7095092d6
real      (kind=8) :: alphapi= 1.9928
real      (kind=8) :: betapi = 0.45500
real      (kind=8) :: Dpi    = 4.0562
real      (kind=8) :: C6pi   = 0.21661134d6
real      (kind=8) :: C8pi   = 123454.7952
real      (kind=8) :: C10pi  = 0.307155024d6
real      (kind=8) :: C12pi  = 0.740464032d-10
real (kind=8) :: r,ff
 ff=1.d0
 if(r.le.Dpi)ff=dexp(-(Dpi/r-1.d0)**2)
 if(r==0)then
   V_Ag_Pi= Api
 else
   V_Ag_Pi=Api*dexp(-alphapi*r-betapi*r**2) - ff*(C6pi/r**6 + C8pi/r**8 + C10pi/r**10 + C12pi/r**12 )
 endif
end function

double precision function V_Ag_Sig(r)
implicit none
real      (kind=8) :: Asi    = 40005.2
real      (kind=8) :: alphasi= 0.41529
real      (kind=8) :: betasi = 0.14888
real      (kind=8) :: Dsi    = 25.730
real      (kind=8) :: C6si   = 1.97443d-9
real      (kind=8) :: C8si   = 3.37393d-6
real      (kind=8) :: C10si  = 0.00788363
real      (kind=8) :: C12si  = 7.52624d12
real (kind=8) :: r,ff
 ff=1.d0
 if(r.le.Dsi)ff=dexp(-(Dsi/r-1.d0)**2)
 if(r==0)then
   V_Ag_Sig= Asi
 else
   V_Ag_Sig=Asi*dexp(-alphasi*r-betasi*r**2) - ff*(C6si/r**6 + C8si/r**8 + C10si/r**10 + C12si/r**12 )
 endif
end function
double precision function V_alka(dist,elem)

implicit none

real (kind=8),    parameter :: evk     = 11604.448    ! 1eV   = 11604.448  K
real (kind=8),    parameter :: bohr    = 0.529177249  ! 1Bohr = 0.529 \AA
real (kind=8),    parameter :: hartree = 27.211608    ! 1H    = 27.2  eV
real (kind=8),    parameter :: factor  = evk*hartree

real      (kind=8), intent(IN) :: dist
character (len=3) , intent(IN) :: elem

real (kind=8) :: abar,u,v,roa
real (kind=8) :: a ,b ,c6 ,c8 ,agran
real (kind=8) :: av(5),bv(5),c6v(5),c8v(5),agranv(5)
real (kind=8) :: r2v(5),r4v(5),r6v(5)
real (kind=8) :: d,r2,r4,r6
real (kind=8) :: s
real (kind=8) :: vimp
real (kind=8) :: pi
real (kind=8) :: f(2)

character*3 celem(5)

integer (kind=4) ::  l,l7
integer (kind=4) ::  n,ielem

data celem /  'LI ' ,  'NA ' ,  'K  ' , 'RB '  , 'CS ' /
data av    / 1.588d0, 1.627d0, 1.771d0, 1.805d0, 1.869d0/
data bv    / 0.744d0, 0.744d0, 0.744d0, 0.744d0, 0.744d0/
data c6v   /  22.5d0,  24.7d0,  38.9d0,  44.6d0,  51.2d0/
data c8v   /  1.06d3,  1.29d3,  2.66d3,  3.18d3,  4.34d3/
data r2v   /  2.37d0,  2.37d0,  2.37d0,  2.37d0,  2.37d0/
data r4v   /  7.78d0,  7.78d0,  7.78d0,  7.78d0,  7.78d0/
data r6v   /  50.5d0,  50.5d0,  50.5d0,  50.5d0,  50.5d0/
data agranv/0.0519d0,0.0450d0,0.0259d0,0.0226d0,0.0173d0/


do ielem=1,7
  if(ielem.eq.7) STOP 'From V_impur: This impurity does not exist'
  if(elem.ne.celem(ielem)) cycle
  a     = av(ielem)
  b     = bv(ielem)
  c6    = c6v(ielem)
  c8    = c8v(ielem)
  r2    = r2v(ielem)
  r4    = r4v(ielem)
  r6    = r6v(ielem)
  agran = agranv(ielem)
  exit
end do

d=dist/bohr



pi    = 4.0d0*datan(1.0d0)
abar  = (a+b)/4.d0
u     = a-1.d0
v     = -0.5d0*a**2*u
roa   = agran*d**(2.d0*u)*(1.d0+v/d)**2*dexp(-2.d0*d/a)

do l=1,2
  s  = 0.d0
  l7 = 2*l+7
  do n=0,l7
    s = s+(d/abar)**n/fac(n)
  end do
  f(l)=1.d0-s*dexp(-d/abar)
end do
vimp = 4.d0*pi/3.d0*roa*(r2+37.d0/90.d0*r4/a**2 +            &
       28.d0/525.d0*r6/a**4)-c6*f(1)/d**6-c8*f(2)/d**8

V_alka=vimp*factor

return

contains

!....................................................................
!. Function fac 
!....................................................................

! Small version of factorial

double precision function fac(n)
implicit none
integer :: n,i
fac=1.d0
if (n.eq.0.or.n.eq.1) return
do i=2,n
  fac=fac*i
end do
return
end function fac
end function V_alka


! Everything in Angstrom Kelvin 

function V_Ba_gs(x)
!..........................................................!
! Ba-He potential. Fit with Mathematica in:
! /u/carraca/dmateo/barium/lovallo_Ba/Fit_potential.nb
!..........................................................!
implicit real*8(a-h,r-z)
    if(x<=3.7d0) then
        V_Ba_gs = 1297.1d0
    else
        V_Ba_gs = -6.146426949008034d14/x**18 + 1.6797251118695925d14/x**16 - 1.6513058293510484d13/x**14 + &
                7.015334514979032d11/x**12 - 1.1066880774828484d10/x**10 + 4.900282758445466d7/x**8 - &
                450294.6878721156d0/x**6
    endif
return
end
!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------
!
! Potencial 2D Sigma del Ba+
!
Block Data Inicio_LJ_V_Ba_plus_2D0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1087334.00470567619271105807407299Q0,  84828727.5548886783972206269121162Q0, -2706837566.34944393727530187448557Q0,  &
  47482174965.3049172156448636792347Q0, -519244658310.373320149026960630547Q0,  3776794752968.39098091679878465104Q0,  &
 -18931184609496.2743530917730707405Q0,  66467512849641.5867602022345173376Q0, -163349259737503.439380947298152502Q0,  &
  275471652871895.368050105366543138Q0, -303799513108731.163572027673779711Q0,  197342887711895.715873097734435815Q0,  &
 -57247315878299.6181404889603871296Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 2.3097950000000001D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/ 2.4148340053884458D+04/
Data bb/-1.4167500343195186D+05/
Data cc/ 2.0985459566015401D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_2D0(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_2D0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_2D0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_2D0=Sto
Return
End
!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------
Block Data Inicio_LJ_V_Ba_plus_2D1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17806.3555611251407400764850972719Q0,  1823098.97448459689302902398930546Q0, -85830730.9944556694819916079148077Q0,  &
  1725318458.71278993648593521565983Q0, -19436563702.5972378924691700268134Q0,  135815448941.441168085385272965158Q0,  &
 -614470726947.401636604960682655298Q0,  1813794100230.25415271043202962245Q0, -3408845825677.80602502264479414645Q0,  &
  3767306350864.97856677416953497814Q0, -1870543050026.55759258430089080476Q0, -263830499030.803962487458419711050Q0,  &
  490468253633.822352576865392543355Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 8.9667530000000006D+03/
Data r0/ 2.0000000000000000D+00/
Data aa/ 1.6993055762720283D+04/
Data bb/-9.4455519278364838D+04/
Data cc/ 1.2990556820533771D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_2D1(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_2D1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_2D1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_2D1=Sto
Return
End

Block Data Inicio_LJ_V_Ba_plus_2D2
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D2/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  3438.08194066326485411722455335026Q0,  0.00000000000000000000000000000000Q0, -22271935.4558369081774443426658153Q0,  &
  509815799.148361882435180461889140Q0, -5103275255.84719311793286237036308Q0,  24347543963.7285486850837178167031Q0,  &
 -22098022739.0755793944143793174099Q0, -371505231363.903388585217865272476Q0,  2188591221460.46490291662922516704Q0,  &
 -6008459393013.55312031551791210851Q0,  9242835510174.84636285804976040549Q0, -7681295358440.52888750989831902749Q0,  &
  2696390604822.62910359945307774040Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 1.1687330000000000D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/ 1.9539148260089081D+04/
Data bb/-1.0968122164089163D+05/
Data cc/ 1.5289317408495379D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_2D2(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D2/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_2D2=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_2D2=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_2D2=Sto
Return
End
!
! Nou ajust del Ba+-He amb mes punts a curtes distancies(4-06-15)
!
Block Data Inicio_LJ_V_Ba_plus_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1233025.18431181221429176595319949Q0,  84965045.0899243014427808833458580Q0, -2329830901.26373970562255233678155Q0,  &
  34487035482.8661381911727426332743Q0, -312315418217.722824241631231266832Q0,  1837177830634.86715823984693486481Q0,  &
 -7222749564186.92692954840041234856Q0,  19113448483679.6067074196915946358Q0, -33537885793985.8415962556220185604Q0,  &
  37188922413954.1348677597526953193Q0, -23134511849820.2719590286530591524Q0,  5477820267454.80225927552742559374Q0,  &
  657017642818.351074515977126633909Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 2.7557869999999999D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/-5.1918208117033464D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/ 3.2749691986506852D+04/
Data k0/   3/
End
Double Precision Function V_Ba_plus_gs(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_gs=Sto
Return
End

Block Data Inicio_LJ_V_Ba_plus_gs_fixC4
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  1608302.51173741722899670327302374Q0, -79502305.7654667384046894884587665Q0,  &
  1903942512.84206481338934071885394Q0, -30134081435.0207804471353184171407Q0,  340887729307.024694401976841468399Q0,  &
 -2749060057173.10100779144500705256Q0,  15421639364129.7422511135464734090Q0, -58747373560398.8645333932507612952Q0,  &
  147892600745628.790351767793542405Q0, -234781354096552.082032333269632683Q0,  212676835684307.885943839636875957Q0,  &
 -83791034270656.6190293274527463146Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 6.9240280000000002D+03/
Data r0/ 2.3999999999999999D+00/
Data aa/-8.5702212115060604D+02/
Data bb/ 0.0000000000000000D+00/
Data cc/ 7.7810505607043342D+03/
Data k0/   3/
End
Double Precision Function V_Ba_plus_gs_fixC4(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_gs_fixC4=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_gs_fixC4=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_gs_fixC4=Sto
Return
End
      Double Precision Function V_Ba_plus_old_gs(x)
      Implicit Real*8(A-H,O-Z)
      V_Ba_plus_old_gs=Func_Ba_plus_gs(x)
      Return
      End

      Block Data Inicio_LJ_Ba_plus_gs
      Implicit Real*8(A-H,O-Z)
      Parameter(Npg=16)
      Common/Param_LJ_Ba_plus_gs/r0,aa,bb,cc,Pg(Npg)
      Data Pg/                                                           &
      -5.75372193971343D+06, 7.15787621605274D+08,-4.96564765588082D+10, &
       1.98124418444111D+12,-4.71450937592845D+13, 6.52366641357938D+14, &
      -3.61296282269969D+15,-3.34341761554747D+16, 6.27946406049051D+17, &
      -1.25469318152625D+18,-2.47420248203689D+19,-1.56040616253124D+20, &
       4.94566594588086D+21,-2.19511834000164D+22,-6.01807193986940D+22, &
       4.58579988466807D+23                                              &
      /
      Data r0/ 3.000000000000000D+00/
      Data aa/ 1.406590885854125D+06/
      Data bb/-8.567381411545813D+06/
      Data cc/ 1.304908921806432D+07/
      End
      Double Precision Function Func_Ba_plus_gs(x)
      Parameter(Npg=16)
      Implicit Real*8(A-H,O-Z)
      Common/Param_LJ_Ba_plus_gs/r0,aa,bb,cc,Pg(Npg)
      If(x.le.r0)Then
        Func_Ba_plus_gs=aa*x**2+bb*x+cc
        Return
      Endif
      Sto=0.0d0
      Do k=1,Npg
        Sto=Sto+pg(k)/x**(2*(k-1)+6)
      EndDo
      Func_Ba_plus_gs=Sto
      Return
      End
Block Data Inicio_LJ_V_Ba_plus_pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -9826.50263754921899133897618912485Q0,  0.00000000000000000000000000000000Q0,  7905587.64166108193848904251702271Q0,  &
 -486504109.501964824574972292520630Q0,  10201523709.2410213874782756339362Q0, -115537603753.397220352368797609468Q0,  &
  811129915505.566723750323994437416Q0, -3741881843963.65413253722118668059Q0,  11556765251738.7569992644472524436Q0,  &
 -23666525048784.0598736494386530051Q0,  30832393059680.8435265222087894345Q0, -23137225457502.0957616538305554638Q0,  &
  7615700353700.85986113395012178657Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 7.953992E+03/
Data r0/ 2.1250000000000000D+00/
Data aa/ 1.1631685135869398D+04/
Data bb/-6.7226852121137432D+04/
Data cc/ 9.5880955855918975D+04/
Data k0/   3/
End
Double Precision Function V_Ba_plus_pi(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_pi=Sto
Return
End
!
! Comentari sobre la funcio
!
Block Data Inicio_LJ_V_Ba_plus_pi_fixC4
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  1480585.38629277157051414278594678Q0, -66524410.6117297064088341969433940Q0,  &
  1248284301.42647870704024354840855Q0, -13069642195.0685028205930753960471Q0,  83381807453.7017144683188590909548Q0,  &
 -330934385581.979371931451499683451Q0,  772746327472.086796494677428176487Q0, -781686189760.313352024174457834920Q0,  &
 -739084882072.026661092648940188249Q0,  3162914410401.90196849006392816034Q0, -3570205775497.43736747130271625283Q0,  &
  1460004701781.63332377962082886275Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 1.0377790000000001D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/ 1.9282248948487551D+04/
Data bb/-1.0746305336734858D+05/
Data cc/ 1.4817490133941581D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_pi_fixC4(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_pi_fixC4=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_pi_fixC4=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_pi_fixC4=Sto
Return
End
!
! Adjust of Ba+ 2P sigma (from Fausto calculations)
!
Block Data Inicio_LJ_V_Ba_plus_sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1946800.56568338068309116827100120Q0,  160853137.584288414167976754320470Q0, -5443802098.75045627541231247577964Q0,  &
  101026984958.806322366443362695122Q0, -1156732968598.98543549109440671187Q0,  8665994393409.78419887735371608207Q0,  &
 -43862571678823.8062435100392360317Q0,  152374645410796.393876923253474426Q0, -363469581548844.057191969295126660Q0,  &
  584580735231448.559942865666177553Q0, -605026417604047.637639655063714044Q0,  363284554181176.601081977921602535Q0,  &
 -95975010650140.9085581326839809273Q0/ 
Data rcutoff/ 1.5000000000000000D+00/
Data fcutoff/ 2.2626640000000000D+07/
Data r0/ 1.0000000000000000D+00/
Data aa/ 1.5500784610027036D+04/
Data bb/-9.5332029461583734D+04/
Data cc/ 1.5198462687793490D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_sigma(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_sigma=Sto
Return
End
!
! Comentari sobre la funcio
!
Block Data Inicio_LJ_V_Ba_plus_sigma_fixC4
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  17447561.0137612249712543003457326Q0, -1721357727.60473744645852747379753Q0,  &
  71971618157.3904453097848751064117Q0, -1705493726202.32479039775818481073Q0,  25631441812788.7446995165917324296Q0,  &
 -257596356972754.910546506340985159Q0,  1771762993803105.76996486810206154Q0, -8363689636005407.15424774634164701Q0,  &
  26640527998718191.1824700127470625Q0, -54721138485295325.7984174767437262Q0,  65456519239493414.9665950639016363Q0,  &
 -34651280867863688.9220625829424204Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 2.0907770000000000D+04/
Data r0/ 3.5000000000000000D+00/
Data aa/ 1.4415808830195720D+03/
Data bb/-1.3729666852098217D+04/
Data cc/ 3.3158457362188419D+04/
Data k0/   3/
End
Double Precision Function V_Ba_plus_sigma_fixC4(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_sigma_fixC4=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_sigma_fixC4=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_sigma_fixC4=Sto
Return
End
  Block Data Inicio_AP_Fit_APG_Cs_7s0                                                                  
      Implicit Real*8(A-H,O-Z)
      Parameter(Npg=12)
      Common/Param_AP_Fit_APG_Cs_7s0/r0,aa,bb,cc,Pg(Npg),Ps(Npg)
      Data Pg/								&
     -6.99807816089155D+07, 1.23439533336424D+04, 7.41947367149584D+06, &
      5.14036819939642D+07,-2.54666356840017D+13, 4.53268712277415D+09, &
      4.66829432137196D+12, 1.23144121825973D+14, 3.89994535976293D+12, &
      2.25670962111522D+16,-2.36027013243917D-04,-5.02039698946873D-04  &
     /
      Data Ps/								&
      1.29148499223708D+00, 7.42273770077802D-01, 1.83956945157573D+00, &
      3.01474078952469D+00, 6.17493897637321D+00, 4.97136805912161D+00, &
      6.59152034486467D+00, 8.41614866661145D+00, 7.36125612301946D+00, &
      1.34956539961155D+01, 1.18675333603722D+00, 1.59183920812119D+00  &
     /
      Data r0/ 2.100000000000000D+00/
      Data aa/ 2.420229505061327D+07/
      Data bb/-1.036963443170083D+08/
      Data cc/ 1.110675394653385D+08/
      End
      Double Precision Function V_Cs_7s(xx)
      Parameter(Npg=12)
      Implicit Real*8(A-H,O-Z)
      Common/Param_AP_Fit_APG_Cs_7s0/r0,aa,bb,cc,Pg(Npg),Ps(Npg)
      data x0/2.5d0/
      If(xx.lt.x0)then
        x=x0
      Else
         x=xx
      Endif   
      If(x.le.r0)Then
        V_Cs_7s=aa*x**2+bb*x+cc
        Return
      Endif
      Sto=pg(1)*Dexp(-Ps(1)*x)/x
      Do k=2,Npg
        Sto=Sto+pg(k)*x**(k-1)*Dexp(-Ps(k)*x)
      EndDo
      V_Cs_7s=Sto
      Return
      End

      Double Precision Function V_Koutselos_Cs_plus_gs(xx)
!
!      Potencial He-Cs+ de Koutselos et al, JCP 93, 7125 (1990)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Cs_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Cs+
      Real  (Kind=8)  :: v_0=0.5092d0, rho=0.9241d0, R_m=6.38d0
      Real  (Kind=8)  :: C4=0.0d0, C6=13.69d0, C8=236.4d0
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0
! Acaban los parametros para el Cs+

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats atòmiques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_Cs_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_Cs_plus_gs
      Double Precision Function V_Koutselos_Na_plus_gs(xx)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Na_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Na+
      Real  (Kind=8)  :: v_0=0.2063d0, rho=0.7640d0, R_m=4.55d0
      Real  (Kind=8)  :: C4=0.0d0, C6=1.424d0, C8=13.17d0
! Acaban los parametros para el Rb+
!
! Parametres propis del He
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats atòmiques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_Na_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_Na_plus_gs
      Double Precision Function V_Koutselos_Rb_plus_gs(xx)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Rb_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Rb+
      Real  (Kind=8)  :: v_0=0.3474d0, rho=0.8868d0, R_m=5.85d0
      Real  (Kind=8)  :: C4=0.0d0, C6=8.8d0, C8=126.2d0
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0
! Acaban los parametros para el Rb+

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats atòmiques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_Rb_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_Rb_plus_gs
      Double Precision Function V_Koutselos_K_plus_gs(xx)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Rb_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Rb+
      Real  (Kind=8)  :: v_0=0.2690d0, rho=0.8589d0, R_m=5.49d0
      Real  (Kind=8)  :: C4=0.0d0, C6=5.83d0, C8=72.05d0
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0
! Acaban los parametros para el Rb+

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats atòmiques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_K_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_K_plus_gs
Block Data Inicio_LJ_Rb_6s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Common/Param_LJ_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
Data Pg/                           &
 -926758671.131109233668080149984260    ,  17751542326.9239124078481901178123    ,  2294371208218.10069268028221768684    ,  &
 -127580931840930.107534925325246956    ,  3010274903906904.53291694447191940    , -40628265880993902.4497466283119501    ,  &
  341766354107911552.266136623510917    , -1829660300635868428.84732815469319    ,  6079663336799721325.17593352865994    ,  &
 -11458268191609309405.5213761979227    ,  9376679057615973867.05134498742729    /
Data rcutoff/ 4.5000000000000000D+00/
Data fcutoff/ 8.2666990000000005D+03/
Data r0/ 4.5000000000000000D+00/
Data aa/ 8.8941971711569586D+03/
Data bb/-8.8209867027930217D+04/
Data cc/ 2.1983426802106350D+05/
End
Double Precision Function V_Rb_6s(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto
Common/Param_LJ_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
If(x.le.rcutoff)Then
  V_Rb_6s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_6s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
Do k=1,Npg
  Sto=Sto+pg(k)/x**(k+5)
EndDo
V_Rb_6s=Sto
Return
End
Block Data Inicio_LJ_V_Cs_sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Common/Param_LJ_V_Cs_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
Data Pg/                           &
  218972446.498381935851370116102201    , -20328460725.5889253542688904081392    ,  792064635008.609454121016549101098    ,  &
 -16988413450311.4536905549403585190    ,  220781047872801.503342285464821515    , -1818628789150376.23619182940382982    ,  &
  9689518147992622.13883260672871893    , -33234418435959134.0581023528704222    ,  70598735214769031.0197085111104892    ,  &
 -83973545896479856.9435608557993711    ,  42292199911291753.0362255745407528    /
Data rcutoff/ 3.9688300000000001D+00/
Data fcutoff/ 1.4295530000000001D+03/
Data r0/ 3.7999999999999998D+00/
Data aa/-5.2840987940870477D+01/
Data bb/ 0.0000000000000000D+00/
Data cc/ 2.2618840928664013D+03/
End
Double Precision Function V_Cs_sigma(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto
Common/Param_LJ_V_Cs_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
If(x.le.rcutoff)Then
  V_Cs_sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Cs_sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
Do k=1,Npg
  Sto=Sto+pg(k)/x**(k+5)
EndDo
V_Cs_sigma=Sto
Return
End
Block Data Inicio_LJ_Cs_pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Common/Param_LJ_Cs_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
Data Pg/                           &
 -271926.174454643323364441958407885    , -18233848.3121256520794825424693032    ,  953937957.545172931044438821385762    ,  &
 -27223560403.3593039055504509753318    ,  364166058669.852409015532291696960    , -2714770540249.03637313691657921121    ,  &
  12344772907874.9689365059244245234    , -35232833044851.8149649455556894484    ,  61848915415723.4235994491545114415    ,  &
 -61192851763498.3879914475300334623    ,  26164252463999.8193908959197343657    /
Data rcutoff/ 2.3813000000000000D+00/
Data fcutoff/ 2.5438229999999999D+03/
Data r0/ 2.3799999999999999D+00/
Data aa/-4.6845761002283029D+02/
Data bb/ 0.0000000000000000D+00/
Data cc/ 5.2002536629712176D+03/
End
Double Precision Function V_Cs_pi(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto
Common/Param_LJ_Cs_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
If(x.le.rcutoff)Then
  V_Cs_pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Cs_pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
Do k=1,Npg
  Sto=Sto+pg(k)/x**(k+5)
EndDo
V_Cs_pi=Sto
Return
End
!
! Function for He-Cs gs(6s) potential, from Fausto Cargnoni
!
Block Data Inicio_LJ_V_Fausto_Cs_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -692107812.206136401684998300995763Q0,  53612999836.2373067596697398126653Q0, -1870657670551.24152231789523420133Q0,  &
  38854872489697.0224741123416213730Q0, -535664068036490.962168104702289226Q0,  5175823532151103.79719017070590004Q0,  &
 -36081508593491482.7030173645408564Q0,  183930078114942817.654215375874979Q0, -686481003395905469.369448400855558Q0,  &
  1855128368497788766.28473667982076Q0, -3532497504702993458.11254607147916Q0,  4493344993463040582.19605074943107Q0,  &
 -3425510460926691233.35152150484476Q0,  1183127973060179869.81759271465610Q0/ 
Data rcutoff/ 2.7500000000000000D+00/
Data fcutoff/ 2.1223260000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/ 2.1884620798613876D+04/
Data bb/-1.2336001555678876D+05/
Data cc/ 1.7538255947080589D+05/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_gs(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_gs=Sto
Return
End
!
! Function for He-Cs(6p_Sigma), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Cs_6p0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  983012682.924529743263329958061894Q0, -84137127560.7181859398054466255159Q0,  3147363940329.22685865151335515968Q0,  &
 -68037224947910.0508055503708771544Q0,  948500179815737.587517335900477079Q0, -9034983578856463.08414720204540198Q0,  &
  60838778286068358.5669137120545770Q0, -294962476058276873.288394571504344Q0,  1035102455473919644.13005948263137Q0,  &
 -2608241792508316463.19876570809066Q0,  4603440277322316076.27999884747413Q0, -5404901585348550051.78209124683525Q0,  &
  3792737249528805624.44606122782020Q0, -1203728882114200188.77827500640065Q0/ 
Data rcutoff/ 2.6000000000000001D+00/
Data fcutoff/ 2.3170690000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/ 8.0635334129850098D+03/
Data bb/-4.6993574624164961D+04/
Data cc/ 7.0021714221092683D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_6p0(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_6p0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_6p0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_6p0=Sto
Return
End
!
! Function for He-Cs(7s), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Cs_7s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_7s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  3825265036.06471340558006213331397Q0, -106391699401.629756164352751452134Q0, -2415651814425.17303291100897069416Q0,  &
  160947890567832.110284939811381405Q0, -3494992921526184.29345177011789542Q0,  43853249263146779.6148549547055530Q0,  &
 -362634101182696589.346055675040595Q0,  2082630471929281529.33806179721184Q0, -8475927125141034661.39181490017403Q0,  &
  24430137804313288176.7891778391113Q0, -48844600104781726185.3454398788720Q0,  64493166634665917457.9924699139579Q0,  &
 -50603788390020951995.5383689991595Q0,  17874696930230973941.5742344069639Q0/ 
Data rcutoff/ 2.8999999999999999D+00/
Data fcutoff/ 1.3262239999999999D+03/
Data r0/ 2.7000000000000002D+00/
Data aa/ 2.7945811089763697D+03/
Data bb/-1.9344908344361000D+04/
Data cc/ 3.3924031337637469D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_7s(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_7s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_7s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_7s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_7s=Sto
Return
End
!
! Function for He-Cs(6p_Pi), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Cs_6p1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -28916389.8646917185914187350216402Q0,  1698212759.47458891331191572250302Q0, -46334669587.0885737998848348571058Q0,  &
  757625158504.461142433946867564510Q0, -8235739328760.26112437184337510191Q0,  62619687246020.4804186731767136627Q0,  &
 -342217765797573.666436871785982512Q0,  1361896256113801.28031727108845724Q0, -3953071532136512.29016059435905051Q0,  &
  8282683882411200.31626276523433664Q0, -12202657183519187.5983436477674604Q0,  11995345991890722.2461560807048455Q0,  &
 -7064895065092973.46664535404252966Q0,  1885957227405596.86900447464232992Q0/ 
Data rcutoff/ 2.3999999999999999D+00/
Data fcutoff/ 2.6430850000000000D+03/
Data r0/ 2.1000000000000001D+00/
Data aa/ 3.6295506204203572D+04/
Data bb/-1.8721605892180325D+05/
Data cc/ 2.4278672572368380D+05/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_6p1(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_6p1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_6p1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_6p1=Sto
Return
End
!
! Function for He-Rb gs(5s), from Fausto Cargnoni
!
Block Data Inicio_LJ_V_Fausto_Rb_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  19100771.2956833483968680974522071Q0, -1172294215.25983631779261524084372Q0,  33435877600.2318788619551266974640Q0,  &
 -626736755694.868731078365328255351Q0,  8889903709250.47980290206377580454Q0, -98809797408403.8989098633620172560Q0,  &
  838535230285386.179156473882388379Q0, -5250871933035878.90801737880175866Q0,  23702922123542559.7263712283775712Q0,  &
 -75653235518269118.4852893687161765Q0,  166080556623113322.510632192031701Q0, -238324245940614178.157128370536680Q0,  &
  201208460378793600.224928931873042Q0, -75780945879534690.8963660066902431Q0/ 
Data rcutoff/ 2.7500000000000000D+00/
Data fcutoff/ 1.7084310000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/ 2.5322093281985394D+03/
Data bb/-1.6991936870778129D+04/
Data cc/ 2.9338073941076589D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_gs(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_gs=Sto
Return
End
!
! Function for He-Rb(5p_Sigma), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Rb_5p0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  755939876.954161676057804724548838Q0, -61892314758.9221075767972745682921Q0,  2217888483872.67461384962575591948Q0,  &
 -45968257647250.0207867967385818922Q0,  614701763499008.868032004191219577Q0, -5618036039629895.48401381180406123Q0,  &
  36301923962954331.4203297556619919Q0, -168897353040772297.040299997362214Q0,  568755905574416789.719347690523984Q0,  &
 -1375110600761493082.00782733743110Q0,  2328499818364179948.01741052654689Q0, -2622699174021865946.26117767851036Q0,  &
  1765450061925162442.37646295715979Q0, -537479750788251198.126286046944468Q0/ 
Data rcutoff/ 2.3999999999999999D+00/
Data fcutoff/ 3.3863560000000002D+03/
Data r0/ 2.3999999999999999D+00/
Data aa/ 5.1234896077967460D+03/
Data bb/-2.8817394513216812D+04/
Data cc/ 4.3036802950585261D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_5p0(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_5p0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_5p0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_5p0=Sto
Return
End
!
! Function for He-Rb(6p_Pi), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Rb_5p1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  21065388.4588612130967750242379954Q0, -1356769725.00046943557636990593514Q0,  38031190615.3371756235401161001838Q0,  &
 -628459250285.123546375190971037063Q0,  6849776778387.97809323521638065669Q0, -52146078017554.2500446984422353406Q0,  &
  285792287469017.676532214877262703Q0, -1143063591357831.95828165866077991Q0,  3339589997615048.00774407054849081Q0,  &
 -7045930377343015.98657569250566483Q0,  10445469637201842.5898884645426793Q0, -10315641071334019.0899871943233809Q0,  &
  6090035673546989.48603755367194306Q0, -1625247766723261.80485107267140556Q0/ 
Data rcutoff/ 2.2000000000000002D+00/
Data fcutoff/ 2.7312600000000002D+03/
Data r0/ 1.0000000000000000D+00/
Data aa/ 3.2376196740691164D+04/
Data bb/-1.5866800201408935D+05/
Data cc/ 1.9510007184216869D+05/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_5p1(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_5p1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_5p1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_5p1=Sto
Return
End
!
! Potencial del Rb+ calculat per en Fasuto et al. amb long range del tipus 1/r**4
!
Block Data Inicio_LJ_V_Fausto_Rb_plus_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -406622.375170394176501817185213738Q0,  21332931.6604765915055670610776133Q0, -515840595.362288581156823966609460Q0,  &
  7286963557.49977496017331947024218Q0, -67050629086.2414609732865947739718Q0,  423545364741.341210425847019599572Q0,  &
 -1884804263958.89912314686654731090Q0,  5961212677076.96608868548002263051Q0, -13327238149555.0874553559671429585Q0,  &
  20600169905046.3280051663330937625Q0, -20966862545953.0513329124868195299Q0,  12657369730747.1509804505111593368Q0,  &
 -3437256517661.33873687152088791033Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 6.9130180000000000D+03/
Data r0/ 2.1000000000000001D+00/
Data aa/ 2.2004608525779306D+04/
Data bb/-1.1197802949237564D+05/
Data cc/ 1.4285064271139764D+05/
Data k0/   3/
End
Double Precision Function V_Fausto_Rb_plus_gs(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_plus_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_plus_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Fausto_Rb_plus_gs=Sto
Return
End
!
! Function for He-Rb(6s), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Rb_6s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -4780502820.90050189325416969698484Q0,  401314863378.456544916850650103099Q0, -14753972047325.5303952564485971199Q0,  &
  314957742848247.648994876344500017Q0, -4375394053042189.01319139625006427Q0,  42025746888366128.5201928564392822Q0,  &
 -288764135760081719.682880070083058Q0,  1443952790541320699.20726573664655Q0, -5274837452555578768.65696176481809Q0,  &
  13946175032217141852.1624775870371Q0, -26003985198509507719.3486720078825Q0,  32447777154760679012.5074313031215Q0,  &
 -24326312738770487365.6660811096881Q0,  8287435576981232747.15789107114440Q0/ 
Data rcutoff/ 3.0000000000000000D+00/
Data fcutoff/ 2.5755680000000002D+03/
Data r0/ 2.8999999999999999D+00/
Data aa/-8.1144161152598926D+01/
Data bb/ 2.1803190565823604D+02/
Data cc/ 2.6582835372153691D+03/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_6s(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_6s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_6s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_6s=Sto
Return
End

Double Precision Function V_grafeno(z)
!
! Potencial He - Superficie de Grafeno, de Pilar de Lara
!
Implicit Real*8(A-H,O-Z)
!
! Parametres del potencial
!
!  z_cutoff=1.4, umax=3.282645d4
!
!
! A     = 554811 meV
! alpha = 3.21193 Ang.-1
! C4    = 1531.83 meV. Angs4
! C6    = 22259 meV Angs6
!
Data C_meV_to_K/11.604448d0/
Data A/554811.d0/, alpha/3.21193d0/, C4/1531.83d0/, C6/22259.d0/

V_average=A*exp(-alpha*z)-C4/(z**4)-C6/(z**6)

V_grafeno=V_average*C_meV_to_K
Return
End
Double precision Function V_grafeno_sin_dispersion(z)
!
! Potencial He - Superficie de Grafeno, de Pilar de Lara
!
Implicit Real*8(A-H,O-Z)
!
! Parametres del potencial
!
!  z_cutoff=1.4, umax=5.107576d4
!
! alpha = 1.54605 Angs.-1
! D     = 0.493191 meV
! Ze    = 4.34871 Angs.
!
Data C_meV_to_K/11.604448d0/
Data D/0.493191d0/, alpha/1.54605d0/, Ze/4.34871d0/

V_Morse=D*(1.0d0-dexp(-alpha*(z-Ze)))**2-D

V_grafeno_sin_dispersion=V_Morse*C_meV_to_K
Return
End
Double Precision Function V_TiO2(Z)
!
!  Ultim potencial He-TiO2, calculat per la Pilar
!
!     Z es la distancia
!     al Ti(5f) y estÃ¡ en AA, VZ en cm-1

IMPLICIT REAL*8(A-H,O-Z)
Data cmm1toK/1.4387770d0/
Data de/95.9583d0/,alpha/1.47089d0/,ze/4.22534d0/
Data z0/4.77765d0/,c3/8.85163d0/, aa/1.06161d0/
vz= de*(1-dexp(-alpha*(z-ze)))**2-de 
vzd= - 0.5d0 * (1.d0 + dtanh ((z-z0)/aa)*c3*(c3/z)**3)
V_TiO2=(vz+vzd)*cmm1toK
return
end
Double Precision Function V_Au(rr)
Implicit None
!Implicit Real*8(A-H,R-Z)
!
! Interaccion Oro-He: entrada en Ang, salida en K
!

Real  (Kind=8)  :: rr,a1,a2,a3,a4,a5,a6,a7,a8,a9,r,Re,De


Data a1/1.6417d0/,a2/-0.6164d0/,a3/0.3770d0/,a4/-0.0549d0/
Data a5/0.0018d0/,a6/0.0046d0/,a7/-0.0010d0/,a8/0.0001d0/
Data a9/2.7d-6/,Re/4.124d0/,De/20.301143d0/
! rms = 0.0034
! De = 14.110 ! in cm^-1

r = rr - Re
!V_Au = -De*(1.d0 + a1*r + a2*r**2 + a3*r**3 + a4*r**4 + a5*r**5 + a6*r**6 + a7*r**7 + a8*r**8 + a9*r**9)*dexp(-a1*r)

V_Au = -De*(1.d0 + r*(a1 + r*(a2 + r*(a3 + r*(a4 + r*(a5 + r*(a6 + r*(a7 + r*(a8 + r*a9)))))))))*dexp(-a1*r)

End Function V_Au
      Double Precision Function V_Au_TiO2(z)
      Implicit None
      Real (Kind=8)  :: C6=-4.18518d6, C8=6.54459D7, C10=-3.73112D8, C12=9.44942D8, C14=-8.97265D8
      Real (kind=8)  :: z, Aux, Aux2, Aux6, CmeV_to_K=11.604448d0
      If(z.Gt.2.)Then
        Aux=1.d0/z
        Aux2=Aux**2
        Aux6=Aux2**3
        V_Au_TiO2 = (((((C14*Aux2+C12)*Aux2+C10)*Aux2+C8)*Aux2+C6)*Aux6)*CmeV_to_K
      Else
        V_Au_TiO2 = 1821.44287109375d0*CmeV_to_K
      Endif
      Return
      End

      Double Precision Function V_dAu_TiO2(z)
      Implicit None
      Real (Kind=8)  :: C6=-4.18518d6, C8=6.54459D7, C10=-3.73112D8, C12=9.44942D8, C14=-8.97265D8
      Real (kind=8)  :: z, Aux, Aux2, Aux7, CmeV_to_K=11.604448d0
      If(z.Gt.2.)Then
        Aux=1.d0/z
        Aux2=Aux**2
        Aux7=Aux2**3*Aux
        V_dAu_TiO2 = (((((-14.0d0*C14*Aux2-12.0d0*C12)*Aux2-10.0d0*C10)*Aux2-8.0d0*C8)*Aux2-6.0d0*C6)*Aux7)*CmeV_to_K
      Else
        V_dAu_TiO2 = -5415.35400390625d0*CmeV_to_K
      Endif
      Return
      End
!
! Potencial pH2-He de Gianturco ajustado (nouevo promedio 29/5/2015)
!
Block Data Inicio_LJ_V_ph2_He
Implicit Real*8(A-H,O-Z)
Parameter(Npg=20)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_ph2_He/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1163439.97300733181599835248769270Q0,  47292172.4090617880619017432408668Q0, -841086014.436697107304093367838305Q0,  &
  8465900672.40128294507703823716660Q0, -53828298395.5314532464705938241342Q0,  227557134285.265792622701141972689Q0,  &
 -650430073601.165531760983033096775Q0,  1218814073699.63463676672914897444Q0, -1216214794159.68959709558137823136Q0,  &
 -567844642080.136519610799404824244Q0,  4746308391191.44071750176385675985Q0, -10143616872740.8356624374938660844Q0,  &
  14089430994274.8504311602782253396Q0, -14374323621667.8460913194491051522Q0,  11078725941332.5159330085651217405Q0,  &
 -6402732017618.46278859102735729547Q0,  2688141511623.61334745699379906995Q0, -771949591132.651962769099634129229Q0,  &
  135129364439.500429837521587722877Q0, -10838377535.0308515558156365736038Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 6.5305870000000003D+04/
Data r0/ 5.0000000000000000D-01/
Data aa/-3.6541132591800960D+06/
Data bb/ 0.0000000000000000D+00/
Data cc/ 8.5553041453631725D+06/
Data k0/   5/
End
Double Precision Function V_ph2_He(x)
Parameter(Npg=20)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_ph2_He/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_ph2_He=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_ph2_He=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_ph2_He=Sto
Return
End
!
! Function for He-Rb(6p_Sigma), From Pascale calculation (06-05-2015)
!
Block Data Inicio_LJ_V_Pascale_Rb_6p0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -27238209654.0970551997286039358492Q0,  3155212616112.62367579232721486859Q0, -159215352636228.495876363247390109Q0,  &
  4657247444544819.83119004845414435Q0, -89015024933389427.4025630902664066Q0,  1190266214999619563.79458870093428Q0,  &
 -11626907744717667292.0430483506713Q0,  85342874325414564884.2727262312123Q0, -479206558620624734149.448129236565Q0,  &
  2079235750526732316320.70782350523Q0, -6995096564286351837337.16158701126Q0,  18191735476673693147864.7261298775Q0,  &
 -36198334606432096939869.4997435395Q0,  54039097780171191560513.9501513440Q0, -58516086287882353120407.4077308311Q0,  &
  43341928857121805849873.8181966873Q0, -19613221630184108600706.5386003175Q0,  4082648274688574570356.65959774196Q0/ 
Data rcutoff/ 2.3500000000000001D+00/
Data fcutoff/ 2.2097469999999998D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 1.2723431741929555D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-6.9728651248261986D+03/
Data k0/   5/
End
Double Precision Function V_Pascale_Rb_6p0(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Pascale_Rb_6p0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Pascale_Rb_6p0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Pascale_Rb_6p0=Sto
Return
End
!
! Function for He-Rb(6p_Pi), From Pascale calculation (06-05-2015)
!
Block Data Inicio_LJ_V_Pascale_Rb_6p1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  473687428.390644638435766046481170Q0, -64347148419.7845939200849296048072Q0,  3637808075975.30386305837579355438Q0,  &
 -117799828212905.378494464152006355Q0,  2485273579964282.73699052550747524Q0, -36637637078170660.1214956082180008Q0,  &
  393964578545592969.739990894637550Q0, -3176400344113005754.50295649136964Q0,  19540934127522939346.8829848452431Q0,  &
 -92640862627641172344.0073533069276Q0,  339664614593611552754.420694144245Q0, -960564218001611305059.888936433560Q0,  &
  2074959812887202322685.47117302676Q0, -3359492547653334985410.92994317689Q0,  3944818157081899283205.45238007755Q0,  &
 -3171181660093470107341.15793633147Q0,  1560672268303756719160.88518302906Q0, -354536702973118696499.223764295971Q0/ 
Data rcutoff/ 2.1499999999999999D+00/
Data fcutoff/ 3.0895619999999999D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 1.5092579075989399D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-6.9138253346225456D+03/
Data k0/   5/
End
Double Precision Function V_Pascale_Rb_6p1(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Pascale_Rb_6p1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Pascale_Rb_6p1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Pascale_Rb_6p1=Sto
Return
End
!
!  Funcio que calcula el potencial Ar-He, entrada \AA, sortida K
!
!  r_cutoff = 2.0d0; umax = 8.0808095081d3
!
!
Double precision function V_Ar_He(r)
Implicit none
real      (kind=8) :: AA     = 124.3d0
real      (kind=8) :: BB     = 2.153d0
real      (kind=8) :: C6     = 9.538d0
real      (kind=8) :: C8     = 167.5d0
real      (kind=8) :: C10    = 3701.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or2, c1or6
real (kind=8) :: r_cutoff=1.9d0, f_cutoff=1.5101972249d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
if(r.Le.r_cutoff)Then
  V_Ar_He = f_cutoff
  Return
Endif        
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
c1or2 = 1.d0/(rr*rr)
c1or6 = c1or2*c1or2*c1or2
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = AA*Stoe - (C6*F6 + (C8*F8 + C10*F10*c1or2)*c1or2)*c1or6
V_Ar_He = Aux*Hartree_to_K
Return
End
!
!  Funcio que calcula el potencial Xe-He, entrada \AA, sortida K
!
!  r_cutoff = 2.0d0; umax = 21450.4867067095d0
!
double precision function V_Xe_He(r)
implicit none
real      (kind=8) :: AA     = 74.9002d0
real      (kind=8) :: BB     = 1.8175d0
real      (kind=8) :: C6     = 21.5674d0
real      (kind=8) :: C8     = 430.148d0
real      (kind=8) :: C10    = 10610.21d0
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, Sto1, Aux
real (kind=8) :: B(18), C(7), S(18)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i, kmax
Save     :: Lfirst, C, B
If(Lfirst)Then
!  Write(6,'("Control per veure si nomes hi pasem una vegada")')
  B(1) = BB
  Do i=2, 18
    B(i) = B(i-1)*BB/Dfloat(i)
  EndDo
  C(1) = C6; C(2) = C8; C(3) = C10
  Do i=4, 7
    C(i) = (C(i-1)/C(i-2))**3*C(i-3)
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
Stoe = Exp(-BB*rr)
Stor = rr
S(1) = B(1)*Stor
Do i=2, 18
  Stor = Stor*rr
  S(i) = S(i-1) + b(i)*Stor
Enddo
Aux = Aa*Stoe
Sto1 = 1.0d0/rr**2
Stor = Sto1**3
Do i = 1, 7
  kmax = 2*i + 4
  Aux = Aux + ( (S(kmax)+1.0d0)*Stoe - 1.0d0 )*C(i)*Stor
  Stor = Stor*Sto1
EndDo
V_Xe_He = Aux*Hartree_to_K
Return
End
!
! Ajust potencial 5p_Sigma de Pascale (29-10-2015) (unitats correctes)!!
!
Block Data Inicio_LJ_V_Rb_5p_sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_5p_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -16369756.1679657643374850977670248Q0,  358314586.951384914502183247424336Q0,  25272036141.5453163057909074234291Q0,  &
 -1211017031542.79911725659637667590Q0,  21309419454172.0588581060167965039Q0, -197508142095030.863285799109515175Q0,  &
  1082206330864640.68613317062793854Q0, -3643715772086083.49000422254676492Q0,  7432352868423233.80679404088747201Q0,  &
 -8450047598780025.77211999338779863Q0,  4117216576621054.92026930080355380Q0/ 
Data rcutoff/ 3.1750600000000002D+00/
Data fcutoff/ 2.7347030000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/-1.4702350246680987D+02/
Data bb/ 0.0000000000000000D+00/
Data cc/ 4.3142207107069107D+03/
Data k0/   5/
End
Double Precision Function V_Rb_5p_sigma(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_5p_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_5p_sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_5p_sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_5p_sigma=Sto
Return
End
!
! Ajust potencial 5_Pi de Pascale (29-10-2015) (unitats correctes)!!
!
Block Data Inicio_LJ_V_Rb_5p_pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_5p_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -8741753.92488395686001702665966285Q0,  486176414.868784278332186303239805Q0, -11591463068.8021123878989034988742Q0,  &
  148126647810.351946447172854810821Q0, -1164177194524.65854784396663864098Q0,  5992071282192.61641098759405074505Q0,  &
 -20640521486966.5879282525969265615Q0,  47219174940119.0277918081596800239Q0, -68878071746805.1901021018758721574Q0,  &
  57983399307116.7535300990213004201Q0, -21431560362273.3380432160483734331Q0/ 
Data rcutoff/ 2.1167099999999999D+00/
Data fcutoff/ 4.8593370000000004D+03/
Data r0/ 2.0000000000000000D+00/
Data aa/-3.0976037301688293D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/ 1.8738030191588245D+04/
Data k0/   5/
End
Double Precision Function V_Rb_5p_pi(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_5p_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_5p_pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_5p_pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_5p_pi=Sto
Return
End
!
! Function for He-Rb(6p_Sigma), From Pascale calculation (30-09-2015)unitats corre
!
Block Data Inicio_LJ_V_Rb_6p_Sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -39855258349.8562779111031444973651Q0,  4616917836286.01917946198638628870Q0, -233048345693666.247944477233277838Q0,  &
  6821255093289797.55038593246614875Q0, -130498399702449880.755284684685938Q0,  1747068821575881055.36281494557793Q0,  &
 -17090582433391375124.6575563934023Q0,  125655099724694531869.495505442861Q0, -706879748387175611849.223558587160Q0,  &
  3073458276584230945887.87133405050Q0, -10363697041788686021395.9401520884Q0,  27020880156532784314594.6844200979Q0,  &
 -53918859334910322896215.0051685128Q0,  80747814635804792925159.0352864619Q0, -87748853076398042247783.7300727085Q0,  &
  65257436010976109146648.5575394867Q0, -29668229502584876388311.8179117254Q0,  6209323876594826013383.91461241976Q0/ 
Data rcutoff/ 2.3500000000000001D+00/
Data fcutoff/ 3.0013080000000000D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 1.8414125504060087D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-1.0102382215546504D+04/
Data k0/   5/
End
Double Precision Function V_Rb_6p_Sigma(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_6p_Sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_6p_Sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_6p_Sigma=Sto
Return
End
!
! Function for He-Rb(6p_Pi), From Pascale calculation (30-09-2015) unitast correct
!
Block Data Inicio_LJ_V_Rb_6p_Pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  840075983.737711643987850807092476Q0, -109565252185.238302267920541824536Q0,  6040041697493.86281523568261964119Q0,  &
 -192096523030023.441530962612113581Q0,  3997548804666069.64682575835922540Q0, -58293811547018697.3031323722230219Q0,  &
  621277997758149591.165075250840335Q0, -4971929235445194262.19280905259923Q0,  30392879922429262840.1034288721564Q0,  &
 -143297441127511031342.841171428469Q0,  522870652282050925787.809971272961Q0, -1472390377555816831804.88187901110Q0,  &
  3168561475038540774092.35836032836Q0, -5112734757976406776681.55293208016Q0,  5985212425925539459439.17268289123Q0,  &
 -4798131130638083481859.59074371543Q0,  2355427891789631538140.18011927531Q0, -533852506199621938645.045228250685Q0/ 
Data rcutoff/ 2.1499999999999999D+00/
Data fcutoff/ 4.4447330000000002D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 2.1849245258600372D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-1.0021627430299464D+04/
Data k0/   5/
End
Double Precision Function V_Rb_6p_Pi(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_6p_Pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_6p_Pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_6p_Pi=Sto
Return
End
!
! Potencial He(2S) del Eloranta (fci) (5-11-2015), Nist: He(2S): 329179.7623 cm-1
!
Block Data Inicio_LJ_V_He2s_fci
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_He2s_fci/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -84641112.1038843174868366388591186Q0,  5501813411.89974475880756271682097Q0, -158511366392.752661259845175521651Q0,  &
  2656659334094.79191838532512036583Q0, -28623286583317.0595767811072611642Q0,  206395570574320.545052095736748454Q0,  &
 -1005349880685171.95006561484544998Q0,  3261504889875669.31289981948240792Q0, -6745427649591000.99921588607933893Q0,  &
  8040867895208215.67859457559937689Q0, -4202197034021318.61309730280594230Q0/ 
Data rcutoff/ 3.1750200000000000D+00/
Data fcutoff/ 5.8000649999999996D+02/
Data r0/ 0.0000000000000000D+00/
Data aa/ 4.3188872751600064D+16/
Data bb/ 0.0000000000000000D+00/
Data cc/-1.3088461815357608D+17/
Data k0/   5/
End
Double Precision Function V_He2s_fci(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_He2s_fci/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_He2s_fci=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_He2s_fci=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_He2s_fci=Sto
Return
End
!
!     Potential Au/Graphene (hollow) with Vz in K  and z in Angs. 
!
Double Precision Function Au_Graphene(z)     
Implicit Double Precision (a-h,o-z)
Data xmeV_To_K/11.604448d0/, A/448058.D0/, beta/2.12336D0/, C4/83077.3D0/      
Au_Graphene=(A*dexp(-beta*z) - C4/z**4)*xmeV_To_K
Return      
End
!
!     Gradiente del potential Au/Graphene (hollow) with dz_Vz in K  and Z in Angs. 
!
Double Precision Function dz_Au_Graphene(z)     
Implicit Double Precision (a-h,o-z)
Data xmeV_To_K/11.604448d0/, A/448058.D0/, beta/2.12336D0/, C4/83077.3D0/      
dz_Au_Graphene=(-A*beta*dexp(-beta*z) + 4.d0*C4/z**5)*xmeV_To_K
Return      
End
!
! Function for He-He+ potential, from Jussi Eloranta (6/6/2016)
!
Block Data Inicio_LJ_V_Eloranta_He2_plus_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=17)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Eloranta_He2_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  588108.342377894438906335100185956Q0, -12708070.5256971354232632580211002Q0,  &
  46593490.4538337953748571219487504Q0,  903614394.043264050274603032772751Q0, -11755047059.0947744782716618472676Q0,  &
  66544986310.3971244712772558577499Q0, -230023812497.692645654143106555337Q0,  540723165541.193441108405864043649Q0,  &
 -908192814328.760315407843038697665Q0,  1115091165095.81869466120443671210Q0, -1006768340907.62729321655030853465Q0,  &
  662532608654.360730944885641242397Q0, -309620499229.935417801050724076383Q0,  97461010879.4325864938071757129741Q0,  &
 -18542204834.6544171733101740365183Q0,  1611683836.74445477113707981227279Q0/ 
Data rcutoff/ 5.2917700000000001D-01/
Data fcutoff/ 1.1804110000000001D+04/
Data r0/ 7.3999999999999999D-01/
Data aa/-3.9258895086799792D+04/
Data bb/ 0.0000000000000000D+00/
Data cc/ 2.2797707492869155D+04/
Data k0/   3/
End
Double Precision Function V_Eloranta_He2_plus_gs(x)
Parameter(Npg=17)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Eloranta_He2_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Eloranta_He2_plus_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Eloranta_He2_plus_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Eloranta_He2_plus_gs=Sto
Return
End
!
!  Potencial de la molecula CH3I con lambda=0
! 
double precision function V_CH3I_0(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6d0,  umax=3.1648655845d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_0 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_0 = Radpot(0,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=1
! 
double precision function V_CH3I_1(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6, umax=-5.1115736616d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_1 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_1 = Radpot(1,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=2
! 
double precision function V_CH3I_2(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6d0, umax=5.5581502684d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_2 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_2 = Radpot(2,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=3
! 
double precision function V_CH3I_3(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6d0, umax=-2.7398540818d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_3 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_3 = Radpot(3,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=4
! 
double precision function V_CH3I_4(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6, umax=2.4025334047d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_4 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_4 = Radpot(4,r_au)*Hartree*eV_K
EndIf
Return
end function
      function radpot(l,r)
      implicit none
!
!     the He-CH3I interaction potential iin grd state, only m=0
!


      integer l
      real*8  radpot,r

      real*8  d,b,c


      select case (l) 

!  fort.50 2:3
      case (0)
      d               = 0.000185445   !  +/- 2.003e-06    (1.08%)
      b               = 0.782913      !  +/- 0.005996     (0.7659%)
      c               = 9.5243        !  +/- 0.01072      (0.1125%)

!  fort.50 2:4
      case (1)
      d               = -4.51657e-05  !  +/- 5.099e-07    (1.129%)
      b               = 0.837569      !  +/- 0.004355     (0.52%)
      c               = 10.3772       !  +/- 0.007319     (0.07053%)

!  fort.50 2:5
      case (2)
      d               = 2.53688e-05   !  +/- 2.933e-07    (1.156%)
      b               = 0.849398      !  +/- 0.002787     (0.3281%)
      c               = 10.6988       !  +/- 0.004548     (0.04251%)

!  fort.50 2:6
      case (3)
      d               = -6.51825e-06  !  +/- 4.938e-08    (0.7576%)
      b               = 0.779254      !  +/- 0.0061       (0.7828%)
      c               = 11.454        !  +/- 0.007429     (0.06486%)

!  fort.50 2:7
      case (4)
      d               = 1.5032e-06    !  +/- 2.002e-08    (1.332%)
      b               = 0.787553      !  +/- 0.0054       (0.6856%)
      c               = 12.236        !  +/- 0.01018      (0.08318%)

      case default
      stop 'in radpot: l out of range'
      end select

!     Morse function
      radpot=d*(1.d0-dexp(-b*(r-c)))**2-d

      return 
      end
!
! Function for He-K(4p_Sigma), From Pascale data potential
!
Block Data Inicio_LJ_V_K_4p_Sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_K_4p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  123270229.922828519094376988757267Q0, -12085480784.0001990705905344017368Q0,  493376520440.846742309241896876082Q0,  &
 -11130615786861.2195078000310223820Q0,  155081464596787.947785296266388633Q0, -1423219373370643.41282606175556080Q0,  &
  8982444221608931.65773550693292715Q0, -40017473420474543.0863116482077444Q0,  127156098740259262.973889637297298Q0,  &
 -286289630234854071.709781745781639Q0,  443029672222162513.944773691873259Q0, -435778229290691161.879891464441776Q0,  &
  216297197302294980.110743898796384Q0, -1493759170257458.23144984879814468Q0,  3521195510768782.35743497581930372Q0,  &
 -105695005528937470.204215869914495Q0,  101576258377617829.001630054279239Q0, -30021724906855482.6250169475416685Q0/ 
Data rcutoff/ 2.2000000000000002D+00/
Data fcutoff/ 5.4944210000000003D+03/
Data r0/-2.0000000000000000D+00/
Data aa/-1.1846279611402621D+13/
Data bb/ 2.7476162593191215D+13/
Data cc/-1.5810332386853908D+13/
Data k0/   5/
End
Double Precision Function V_K_4p_Sigma(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_K_4p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_K_4p_Sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_K_4p_Sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_K_4p_Sigma=Sto
Return
End
!
! Function for He-K(4p_Pi), From Pascale data potential
!
Block Data Inicio_LJ_V_K_4p_Pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_K_4p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  7463823.30565045302747409754900998Q0, -718023677.708100048067350320776803Q0,  28527955672.4973076800039235370599Q0,  &
 -649780615925.820487084894875824681Q0,  9494248151960.38025074812833414839Q0, -95423287484468.9225481965073472623Q0,  &
  688611363436857.181033088432509010Q0, -3659192698303041.57349366577052330Q0,  14504208321803535.5304968565259890Q0,  &
 -43002303720081040.5938560703278609Q0,  94613503021063793.2518086792502445Q0, -151026197668471021.305996587780281Q0,  &
  166277977852662387.539989377379199Q0, -110888302469029095.889154928172435Q0,  22770818008615233.2321789074210083Q0,  &
  27284389361289459.6643521776459617Q0, -23355843002431104.0076310950131271Q0,  5879948233732488.71454983137084160Q0/ 
Data rcutoff/ 1.8999999999999999D+00/
Data fcutoff/ 5.3483040000000001D+03/
Data r0/-2.0000000000000000D+00/
Data aa/ 9.6643483830893945D+11/
Data bb/-2.2413727796818149D+12/
Data cc/ 1.2896495964123530D+12/
Data k0/   5/
End
Double Precision Function V_K_4p_Pi(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_K_4p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_K_4p_Pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_K_4p_Pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_K_4p_Pi=Sto
Return
End
!!
!! Function for He-K(5s), From Pascale data potential (old Fit by Marti)
!!
!Block Data Inicio_LJ_V_K_5s
!Implicit Real*8(A-H,O-Z)
!Parameter(Npg=20)
!Real*16 Pg
!Integer*4 k0
!Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
!Data Pg/                           &
!  22013374032.6756992466989907471231Q0, -2607473486972.96115845928542359422Q0,  137278559465513.504318508804350585Q0,  &
! -4284085369698061.73976060584229719Q0,  89194755667933302.9197442829881481Q0, -1322401314191165810.56518310532255Q0,  &
!  14549200062430442236.7345712007885Q0, -122132199482850453201.933729576988Q0,  797237288979113873142.546850234427Q0,  &
! -4098091530865115675119.54131140732Q0,  16712660424005787663159.1277264440Q0, -54228388406852069443691.0811371931Q0,  &
!  139769338235053275092644.042876545Q0, -284264086314955716327432.688486366Q0,  450393287366663046416997.448177291Q0,  &
! -544150565669437497681315.135493803Q0,  484098752094243106107299.506791578Q0, -298854633397265157385967.162411014Q0,  &
!  114324971752633317584036.132040739Q0, -20408914798270344038969.9560381946Q0/ 
!Data rcutoff/ 2.0000000000000000D+00/
!Data fcutoff/ 2.7298080000000000D+03/
!Data r0/-2.0000000000000000D+00/
!Data aa/-1.1118034495930548D+05/
!Data bb/ 4.7153152909318462D+05/
!Data cc/-4.9833086495060596D+05/
!Data k0/   5/
!End
!Double Precision Function V_K_5s(x)
!Parameter(Npg=20)
!Implicit Real*8(A-H,O-Z)
!Real*16 Pg,Sto,xq,Aux
!Integer*4 k0
!Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
!If(x.le.rcutoff)Then
!  V_K_5s=fcutoff
!  Return
!Endif
!If(x.le.r0)Then
!  V_K_5s=Abs(aa*x**2+bb*x+cc)
!  Return
!Endif
!Sto=0.0q0
!xq=1.0q0/x
!Aux=xq**(k0+1)
!Do k=1,Npg
!  Sto=Sto+pg(k)*Aux
!  Aux=Aux*xq
!EndDo
!V_K_5s=Sto
!Return
!End
!
! Function for He-K(5s), From Pascale data potential (new fit from Maxime Martinez)
!
Block Data Inicio_LJ_V_K_5s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=20)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1168397771.45542931556701660156250Q0,  118762977694.007186889648437500000Q0, -5067632180013.47558593750000000000Q0,  &
  119820387630240.453125000000000000Q0, -1770902027637011.75000000000000000Q0,  17565225586361198.0000000000000000Q0,  &
 -122078985717169440.000000000000000Q0,  608590024952670464.000000000000000Q0, -2190107282135850496.00000000000000Q0,  &
  5613121256762224640.00000000000000Q0, -9769037125165985792.00000000000000Q0,  9996223050685026304.00000000000000Q0,  &
 -2343730974240809984.00000000000000Q0, -7106519104385017856.00000000000000Q0,  5621562826636969984.00000000000000Q0,  &
  4921614491484340224.00000000000000Q0, -6352243461664580608.00000000000000Q0, -4076025095303220224.00000000000000Q0,  &
  8675010493626412032.00000000000000Q0, -3523847995886618112.00000000000000Q0/
Data rcutoff/ 1.980000000000000D+00/
Data fcutoff/ 2.946118600000D+03/
Data r0/-2.0000000000000000D+00/
Data aa/-1.1118034495930548D+05/
Data bb/ 4.7153152909318462D+05/
Data cc/-4.9833086495060596D+05/
Data k0/   5/
End
Double Precision Function V_K_5s(x)
Parameter(Npg=20)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_K_5s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_K_5s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_K_5s=Sto
Return
End
!
! Function for Ag-Amorf_Graphene, From Pilar & Ricardo
!
Block Data Inicio_LJ_V_Ag_Graphene_Amorf
Implicit Real*8(A-H,O-Z)
Parameter(Npg=20)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ag_Graphene_Amorf/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  6520053292.61616126880185962278597Q0, -918464943394.461620370979352927232Q0,  55277361903554.6322625907166690706Q0,  &
 -1908687721418280.26519637838045051Q0,  42876752253768590.0379531867092950Q0, -673631689880263937.577513289837356Q0,  &
  7753571958855461494.03223442954254Q0, -67456829997827292426.0143033339887Q0,  453140328806614255951.042844599700Q0,  &
 -2383712935315990337409.86157129657Q0,  9903388137504520111389.34436541863Q0, -32614132404123601870878.9652482507Q0,  &
  85047981814704420869972.3572341441Q0, -174535601306065036127719.999618677Q0,  278400480074652198718256.087011167Q0,  &
 -337958051090618684547853.261359119Q0,  301590125879860950437663.888010620Q0, -186494245763497980212184.175837030Q0,  &
  71375363687513215919486.0552441231Q0, -12734873402816982253348.4841107942Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 2.2432860000000001D+03/
Data r0/-2.0000000000000000D+00/
Data aa/-8.9557209992220247D+04/
Data bb/ 3.8053908634304180D+05/
Data cc/-4.0283715595284803D+05/
Data k0/   5/
End
Double Precision Function V_Ag_Graphene_Amorf(x)
Parameter(Npg=20)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ag_Graphene_Amorf/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ag_Graphene_Amorf=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ag_Graphene_Amorf=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ag_Graphene_Amorf=Sto
Return
End
!
! Function for dzV_Ag-Amorf_Graphene, From Pilar & Ricardo
!
Block Data Inicio_LJ_dzV_Ag_Graphene_Amorf
Implicit Real*8(A-H,O-Z)
Parameter(Npg=20)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_dzV_Ag_Graphene_Amorf/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  6520053292.61616126880185962278597Q0, -918464943394.461620370979352927232Q0,  55277361903554.6322625907166690706Q0,  &
 -1908687721418280.26519637838045051Q0,  42876752253768590.0379531867092950Q0, -673631689880263937.577513289837356Q0,  &
  7753571958855461494.03223442954254Q0, -67456829997827292426.0143033339887Q0,  453140328806614255951.042844599700Q0,  &
 -2383712935315990337409.86157129657Q0,  9903388137504520111389.34436541863Q0, -32614132404123601870878.9652482507Q0,  &
  85047981814704420869972.3572341441Q0, -174535601306065036127719.999618677Q0,  278400480074652198718256.087011167Q0,  &
 -337958051090618684547853.261359119Q0,  301590125879860950437663.888010620Q0, -186494245763497980212184.175837030Q0,  &
  71375363687513215919486.0552441231Q0, -12734873402816982253348.4841107942Q0/ 
Data rcutoff/ 2.2000000000000000D+00/
Data fcutoff/ -6.976782d3/
Data r0/-2.0000000000000000D+00/
Data aa/-8.9557209992220247D+04/
Data bb/ 3.8053908634304180D+05/
Data cc/-4.0283715595284803D+05/
Data k0/   5/
End
Double Precision Function dzV_Ag_Graphene_Amorf(x)
Parameter(Npg=20)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_dzV_Ag_Graphene_Amorf/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  dzV_Ag_Graphene_Amorf=fcutoff
  Return
Endif
If(x.le.r0)Then
  dzV_Ag_Graphene_Amorf=Abs(2.*aa*x+bb)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+2)
Do k=1,Npg
  Sto=Sto-pg(k)*Aux*(k0+k)
  Aux=Aux*xq
EndDo
dzV_Ag_Graphene_Amorf=Sto
Return
End
!....................
!.. Funcion for Ca ..
!....................
double precision function V_ca(rr)

!...ground state potential for Ca-He. Mauschick and Meyer 1989

implicit none

real (kind=8),    parameter :: bohr    = 0.529177249  ! 1Bohr = 0.529 \AA
integer (kind=4), parameter :: mlr=5
integer (kind=4), parameter :: mp =6
real    (kind=8)            :: clr(mlr),p(mp)
real    (kind=8)            :: r,rr
integer (kind=4)            :: k,m

data p/ 9.397738d-6,1.233746d0,0.06802583d0,    &
        0.4441812d0,1.641748d0,4.73051d0/

data clr/36.404,2.0834d+3,1.2405d+5,7.68460D+6,4.95275D+8/ 

interface  
function gam(l,x)
   integer (kind=4) :: l
   real    (kind=8) :: x
   real    (kind=8) :: gam
end function gam
end interface
r=rr/bohr
V_ca = p(1)*(1+p(3)*r**p(4))*exp(-p(2)*(r-11)) ! Repulsion

do k=1,mlr
  m    = 4+2*k
  V_ca = V_ca-clr(k)/r**m*gam(m,p(5)*(r-p(6))) ! Attraction
end do
V_ca = V_ca/3.1667d-6
return

end

!...................
!.. Function gam ...
!...................

double precision function gam(l,x)
implicit none

integer (kind=4) :: i,l
real    (kind=8) :: x,xx


gam = 0.0d0
if(x.lt.0.d0) return
gam = 1.0d0
xx  = 1.0d0

do i=1,l
  xx  = xx*x/i
  gam = gam+xx
end do
gam = 1.0d0-gam*exp(-x)

return
end 
!
!  Funcio que calcula el potencial Ar-Ar, entrada \AA, sortida K
!
!  r_cutoff = 2.0d0; umax = 8.0808095081d3
!
!
Double precision function V_Ar_Ar(r)
Implicit none
real      (kind=8) :: AA     = 748.3d0
real      (kind=8) :: BB     = 2.031d0
real      (kind=8) :: C6     = 64.3d0
real      (kind=8) :: C8     = 1623.d0
real      (kind=8) :: C10    = 49060.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or2, c1or6
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=1.224418d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
if(r.Le.r_cutoff)Then
  V_Ar_Ar = f_cutoff
  Return
Endif        
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
c1or2 = 1.d0/(rr*rr)
c1or6 = c1or2*c1or2*c1or2
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = AA*Stoe - (C6*F6 + (C8*F8 + C10*F10*c1or2)*c1or2)*c1or6
V_Ar_Ar = Aux*Hartree_to_K
Return
End
!
!  Funcio que calcula el potencial Xe-Xe, entrada \AA, sortida K
!
!
Double precision function V_Xe_Xe(r)
Implicit none
real      (kind=8) :: AA     = 951.8d0
real      (kind=8) :: BB     = 1.681d0
real      (kind=8) :: C6     = 285.9d0
real      (kind=8) :: C8     = 12810.d0
real      (kind=8) :: C10    = 619800.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or2, c1or6
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=8.858633d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
If(r.Le.r_cutoff)Then
  V_Xe_Xe = f_cutoff
  Return
EndIf  
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
c1or2 = 1.d0/(rr*rr)
c1or6 = c1or2*c1or2*c1or2
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = AA*Stoe - (C6*F6 + (C8*F8 + C10*F10*c1or2)*c1or2)*c1or6
V_Xe_Xe = Aux*Hartree_to_K
Return
End
!
!  Funcio que calcula la derivada del potencial Ar-Ar, entrada \AA, sortida K
!
!  r_cutoff = 2.0d0; umax = 8.0808095081d3
!
!
Double precision function drV_Ar_Ar(r)
Implicit none
real      (kind=8) :: AA     = 748.3d0
real      (kind=8) :: BB     = 2.031d0
real      (kind=8) :: C6     = 64.3d0
real      (kind=8) :: C8     = 1623.d0
real      (kind=8) :: C10    = 49060.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or, c1or2, c1or7, b2, b6, b8, b10
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=-5.266478d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
If(r.Le.r_cutoff)Then
  drV_Ar_Ar = f_cutoff
  Return
EndIf  
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
b2=bb*bb
b6=b2*b2*b2
b8=b6*b2
b10=b8*b2
rr = r/rbohr
c1or = 1.d0/rr
c1or2 = c1or*c1or
c1or7 = c1or2*c1or2*c1or2*c1or
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = -bb*AA*Stoe + (6.d0*C6*F6 + (8.d0*C8*F8 + 10.d0*C10*F10*c1or2)*c1or2)*c1or7      &
    -  bb*Stoe*(b6*C6/Fact(6) + b8*C8/Fact(8) + b10*C10/fact(10))
drV_Ar_Ar = Aux*Hartree_to_K/rbohr
Return
End
!
!  Funcio que calcula la derivada del potencial Ar-Ar, entrada \AA, sortida K
!
!
Double precision function drV_Xe_Xe(r)
Implicit none
real      (kind=8) :: AA     = 951.8d0
real      (kind=8) :: BB     = 1.681d0
real      (kind=8) :: C6     = 285.9d0
real      (kind=8) :: C8     = 12810.d0
real      (kind=8) :: C10    = 619800.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or, c1or2, c1or7, b2, b6, b8, b10
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=-3.017394d5
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
If(r.Le.r_cutoff)Then
  drV_Xe_Xe = f_cutoff
  Return
EndIf  
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
b2=bb*bb
b6=b2*b2*b2
b8=b6*b2
b10=b8*b2
rr = r/rbohr
c1or = 1.d0/rr
c1or2 = c1or*c1or
c1or7 = c1or2*c1or2*c1or2*c1or
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = -bb*AA*Stoe + (6.d0*C6*F6 + (8.d0*C8*F8 + 10.d0*C10*F10*c1or2)*c1or2)*c1or7      &
    -  bb*Stoe*(b6*C6/Fact(6) + b8*C8/Fact(8) + b10*C10/fact(10))
drV_Xe_Xe = Aux*Hartree_to_K/rbohr
Return
End
!
! Function for He-Sr(5s2), From Lovallo data potential
!
Block Data Inicio_LJ_V_Sr_He
Implicit Real*8(A-H,O-Z)
Parameter(N=16)
Integer (Kind=4) :: k0
Common/Param_Spline_LJ_V_Sr_He/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
Data X/                                                                                                   &
  2.00000000000000D0,  3.00000000000000D0,  4.00000000000000D0,  5.00000000000000D0,  5.30000000000000D0, &
  6.30000000000000D0,  6.40000000000000D0,  6.50000000000000D0,  6.60000000000000D0,  6.70000000000000D0, &
  6.80000000000000D0,  7.50000000000000D0,  8.00000000000000D0,  9.00000000000000D0,  10.0000000000000D0, &
  12.0000000000000D0/ 
Data A/                                                                                                   &
  19682.2373956429D0,  135538.884184861D0, -28755.0411506975D0,  82682.1648908597D0, -1028.21398783985D0, &
  12723.4933807024D0, -3297.73585567146D0,  5511.51880624365D0, -1397.11955900432D0,  555.849036582276D0, &
  824.670054265175D0, -628.064864640832D0, -6.25389576350048D0, -82.9390390800452D0, -21.4578500114020D0, &
  0.00000000000000D0/ 
Data B/                                                                                                   &
  12722.4786274321D0, -103134.168161786D0,  20086.2758398830D0, -46776.0477850513D0,  607.185542514463D0, &
 -5941.24653774376D0,  1568.70466680649D0, -2497.10517715432D0,  643.184988867479D0, -231.278561395176D0, &
 -349.876069196455D0,  231.217898365948D0, -1.96121496305145D0,  23.6004994757968D0,  5.15614275520381D0, &
  0.00000000000000D0/ 
Data C/                                                                                                   &
 -12360.7671325954D0,  26258.1151304771D0, -4546.99586994003D0,  8825.46885504684D0, -114.763848267464D0, &
  924.669815265587D0, -248.760060445390D0,  376.749146317813D0, -99.0523939885207D0,  31.4645538118754D0, &
  48.9053637826517D0, -28.5738318923353D0, 0.573557273789606D0, -2.26663321941575D0,-0.422197547356456D0, &
  0.00000000000000D0/ 
Data D/                                                                                                   &
  2060.12785543257D0, -2230.85906268660D0,  336.233520681496D0, -555.264127650962D0,  7.01465872100645D0, &
 -47.9818314130174D0,  13.1343079469293D0, -18.9430872716965D0,  5.08729355185571D0, -1.40608693572618D0, &
 -2.26102860096032D0,  1.18249120681688D0,-3.198334177165370D0, 7.320889871743369D0, 1.172770964879045D0, &
  0.00000000000000D0/ 
Data ak/-4.6914394784477353D+05/
Data a12/ 1.2247267783468750D+10/
Data rf/ 8.0000000000000000D+00/
Data r0/ 2.0000000000000000D+00/
Data aa/-2.9997639094396982D+03/
Data cc/ 2.4164204601344718D+04/
Data k0/   6/
End
Double Precision Function V_Sr_He(xx)
Parameter(N=16)
Implicit Real*8(A-H,O-Z)
Integer*4 k0
Common/Param_Spline_LJ_V_Sr_He/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
If(xx.le.r0)Then
  V_Sr_He=aa*xx**2+cc
  Return
Endif
If(xx.gt.rf)Then
  V_Sr_He=ak/xx**k0 + a12/xx**12
  Return
Endif
Call Findi2(X,xx,N,ix)
V_Sr_He= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return
End
double precision function V_Kerkines_Li_gs(r)
implicit none
real      (kind=8) :: A      = 1.6725d4
real      (kind=8) :: alpha  = 0.42685d0
real      (kind=8) :: beta   = 0.19706d0
real      (kind=8) :: D      = 13.988d0
real      (kind=8) :: C6     = 1.5992d5
!real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C12    = 3.2761d10
real (kind=8) :: r,ff
 If(r==0)then
   V_Kerkines_Li_gs = A
   Return
 EndIf  
 ff=1.d0
 If(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 V_Kerkines_Li_gs = A*dexp( -alpha*r -beta*r**2) - ff*(C6/r**6 + C12/r**12)
 Return
end function
double precision function V_Kerkines_Li_Pi(r)
implicit none
real      (kind=8) :: A      = 6.2669d5
real      (kind=8) :: alpha  = 0.d0
real      (kind=8) :: beta   = 2.8981d0
real      (kind=8) :: D      = 4.5734d0
real      (kind=8) :: C6     = 1.9521d5
!real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C10    = 3.7221d6
!real      (kind=8) :: C12    = 1.7527d12
real (kind=8) :: r,ff
 If(r==0)then
   V_Kerkines_Li_Pi = A
   Return
 EndIf  
 ff=1.d0
 If(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 V_Kerkines_Li_Pi = A*dexp( -Beta*r**2) - ff*(C6/r**6 + C10/r**10)
 Return
end function
double precision function V_Kerkines_Li_Sigma(r)
implicit none
real      (kind=8) :: A      = 8.6936d5
real      (kind=8) :: alpha  = 1.6978d0
real      (kind=8) :: beta   = 0.d0
real      (kind=8) :: D      = 21.456d0
!real      (kind=8) :: C6     = 1.5992d5
!real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C10    = 7.0088d3
real      (kind=8) :: C12    = 1.7527d12
real (kind=8) :: r,ff
 If(r==0)then
   V_Kerkines_Li_Sigma = A
   Return
 EndIf  
 ff=1.d0
 If(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 V_Kerkines_Li_Sigma = A*dexp( -alpha*r) - ff*(C10/r**10 + C12/r**12)
 Return
end function
!
! Function for He-Na(3p0), From Pascale data potential
!
Block Data Inicio_LJ_V_Na_Sigma
Implicit Real*8(A-H,O-Z)
Parameter(N=34)
Integer (Kind=4) :: k0
Common/Param_Spline2_LJ_V_Na_Sigma/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
Data X/                                                                                                                      &
  1.058354000000000D+00,   1.322942500000000D+00,   1.587531000000000D+00,   1.852119500000000D+00,   2.116708000000000D+00, &
  2.381296500000000D+00,   2.645885000000000D+00,   2.910473500000000D+00,   3.175062000000000D+00,   3.439650500000000D+00, &
  3.704239000000000D+00,   3.968827500000000D+00,   4.233416000000000D+00,   4.762593000000000D+00,   5.291770000000000D+00, &
  5.820947000000000D+00,   6.350124000000000D+00,   6.879301000000000D+00,   7.408478000000000D+00,   7.937655000000000D+00, &
  8.466832000000000D+00,   8.996009000000001D+00,   9.525186000000000D+00,   1.005436300000000D+01,   1.058354000000000D+01, &
  1.190648250000000D+01,   1.322942500000000D+01,   1.455236750000000D+01,   1.587531000000000D+01,   1.719825250000000D+01, &
  1.852119500000000D+01,   2.116708000000000D+01,   2.381296500000000D+01,   2.645885000000000D+01  /
Data A/                                                                                                                      &
 -3.890915094145514D+05,   2.842048955360034D+06,  -9.387352154073080D+04,   4.860945168944493D+05,  -3.682190378428085D+04, &
  3.674455587895416D+04,  -1.721438688442292D+04,   2.567156632954004D+03,   1.303065248252341D+04,   1.997174092036358D+04, &
  2.697079033034076D+04,   2.768983007896994D+04,   2.813607007945593D+04,   2.551325162844106D+04,   2.110508160272032D+04, &
  1.557964097612617D+04,   1.047575045279671D+04,   6.640191937636850D+03,   4.126079405769719D+03,   2.436647851185728D+03, &
  1.417993363266505D+03,   8.060063945920606D+02,   3.925260172221136D+02,   2.248211095194046D+02,   3.369649265259418D+01, &
 -2.516352812110249D+01,  -1.056065295577870D+01,  -1.446393496628078D+01,   2.398632234516816D+00,  -7.540154642234392D+00, &
  3.404083314497583D-02,  -1.285427052711144D+00,   3.352404111107370D-01,   0.000000000000000D+00  /
Data B/                                                                                                                      &
  2.038171376789891D+06,  -5.288996201788737D+06,   2.590954138724962D+05,  -6.803170351112438D+05,   6.080987612837456D+04, &
 -3.187046795717086D+04,   2.931013070457417D+04,   8.920104630764952D+03,  -9.664700404409080D+02,  -7.020359909054157D+03, &
 -1.268877067572756D+04,  -1.323228617139665D+04,  -1.354851306746780D+04,  -1.189637622666603D+04,  -9.397305012844907D+03, &
 -6.549603104584538D+03,  -4.138363328796164D+03,  -2.465711478632961D+03,  -1.447643044657535D+03,  -8.091302012860393D+02, &
 -4.481966828510839D+02,  -2.441105257537164D+02,  -1.138830286623377D+02,  -6.384358577489398D+01,  -9.667586950263340D+00, &
  5.162994820724863D+00,   1.851541337601345D+00,   2.656210820522679D+00,  -5.303538262396706D-01,   1.203331653157357D+00, &
 -2.351069832904046D-02,   1.634968462906406D-01,  -4.067774479726644D-02,   0.000000000000000D+00  /
Data C/                                                                                                                      &
 -2.228002766666602D+06,   3.310535437812148D+06,  -1.842571776777128D+05,   3.229522377966883D+05,  -2.717962320616065D+04, &
  1.174049618491284D+04,  -1.138242815297223D+04,  -4.376686278390921D+03,  -1.262865297507440D+03,   4.971652255393781D+02, &
  2.027414965546275D+03,   2.164361077633301D+03,   2.239058885755905D+03,   1.892160286442092D+03,   1.419904082219095D+03, &
  9.306878244073738D+02,   5.509724399852153D+02,   3.078297354919530D+02,   1.704103581278810D+02,   8.996936599966477D+01, &
  4.734025531990225D+01,   2.465395552880960D+01,   1.098204433553960D+01,   6.005155905368167D+00,   8.862638499093484D-01, &
 -3.593249938980653D-01,  -1.090152878357445D-01,  -1.643100350939714D-01,   3.641452693110717D-02,  -6.439138223885225D-02, &
  1.848531086740251D-03,  -6.986298498619706D-03,   1.587795113322655D-03,   0.000000000000000D+00  /
Data D/                                                                                                                      &
  7.017194519875840D+05,  -6.937904824902311D+05,   4.000991059182749D+04,  -5.127459084457212D+04,   3.863208197281073D+03, &
 -1.584821113792885D+03,   1.328246327143395D+03,   5.258864047709351D+02,   1.989824492068674D+02,   2.841952689379911D+01, &
 -1.092830475537368D+02,  -1.207848586122947D+02,  -1.266664609339661D+02,  -1.023870680826480D+02,  -7.263922667211827D+01, &
 -4.462452632572438D+01,  -2.469234944108620D+01,  -1.291098267658196D+01,  -6.728004321060446D+00,  -3.349970878848767D+00, &
 -1.671692207449495D+00,  -8.310861197539162D-01,  -3.526383780415761D-01,  -1.876390827050994D-01,  -2.641728434466360D-02, &
  8.454163311614218D-03,   2.147270254097337D-03,   3.413839588565351D-03,  -8.007670768599447D-04,   1.153034317980596D-03, &
 -3.911182099860281D-05,   1.000166367348203D-04,  -2.000332734696399D-05,   0.000000000000000D+00  /
Data ak/-6.7791488520017383D+05/
Data a12/-1.2075208543818171D+12/
Data rf/ 1.5000000000000000D+01/
Data r0/ 1.3000000000000000D+00/
Data aa/-7.5737767294701451D+04/
Data cc/ 1.6488106749051623D+05/
Data k0/   6/
End
Double Precision Function V_Na_Sigma(xx)
Parameter(N=34)
Implicit Real*8(A-H,O-Z)
Integer*4 k0
Common/Param_Spline2_LJ_V_Na_Sigma/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
If(xx.le.r0)Then
  V_Na_Sigma=aa*xx**2+cc
  Return
Endif
If(xx.gt.rf)Then
  V_Na_Sigma=ak/xx**k0 + a12/xx**12
  Return
Endif
Call Findi2(X,xx,N,ix)
V_Na_Sigma= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return
End
!
! Function for He-Na(3p1), From Pascale data potential
!
Block Data Inicio_LJ_V_Na_Pi
Implicit Real*8(A-H,O-Z)
Parameter(N=34)
Integer (Kind=4) :: k0
Common/Param_Spline2_LJ_V_Na_Pi/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
Data X/                                                                                                                      &
  1.058354000000000D+00,   1.322942500000000D+00,   1.587531000000000D+00,   1.852119500000000D+00,   2.116708000000000D+00, &
  2.381296500000000D+00,   2.645885000000000D+00,   2.910473500000000D+00,   3.175062000000000D+00,   3.439650500000000D+00, &
  3.704239000000000D+00,   3.968827500000000D+00,   4.233416000000000D+00,   4.762593000000000D+00,   5.291770000000000D+00, &
  5.820947000000000D+00,   6.350124000000000D+00,   6.879301000000000D+00,   7.408478000000000D+00,   7.937655000000000D+00, &
  8.466832000000000D+00,   8.996009000000001D+00,   9.525186000000000D+00,   1.005436300000000D+01,   1.058354000000000D+01, &
  1.190648250000000D+01,   1.322942500000000D+01,   1.455236750000000D+01,   1.587531000000000D+01,   1.719825250000000D+01, &
  1.852119500000000D+01,   2.116708000000000D+01,   2.381296500000000D+01,   2.645885000000000D+01  /
Data A/                                                                                                                      &
 -3.829612265209504D+05,   2.858488533719470D+06,  -8.850505962450517D+04,   4.985161221791919D+05,  -9.374435653208049D+03, &
  7.506304846840330D+04,   1.048141319844631D+04,   2.389023908616203D+04,   1.761808333336876D+04,   1.864857334491393D+04, &
  1.600855163427982D+04,   1.423083304231262D+04,   1.009103672139864D+04,   5.509370743020555D+03,   2.782969659106633D+03, &
  1.139314871927326D+03,   3.791845333915003D+02,   5.516886143493618D+01,  -3.906518534900924D+01,  -6.883528596251043D+01, &
 -6.820116280662592D+01,  -4.342108704388850D+01,  -3.307391379617849D+01,  -2.971075725994607D+01,  -1.441692751124337D+01, &
 -5.637482470438435D+00,  -4.645791803909031D-01,  -6.208047560804336D+00,   3.766387287683818D+00,  -2.791993659373083D+00, &
  1.963023548370577D-01,  -6.934189615083726D-01,   1.987803926698458D-01,   0.000000000000000D+00  /
Data B/                                                                                                                      &
  2.042259675306722D+06,  -5.308286006551154D+06,   2.607269954197847D+05,  -6.901098946464834D+05,   2.972093393128923D+04, &
 -7.665492155954697D+04,  -3.429930371392179D+03,  -1.725121328743561D+04,  -1.132488263361537D+04,  -1.222365709765917D+04, &
 -1.008554853287642D+04,  -8.741787491132578D+03,  -5.808130850079128D+03,  -2.922098401970659D+03,  -1.376452379611076D+03, &
 -5.293454807614694D+02,  -1.702357979257522D+02,  -2.893553851426660D+01,   9.223735273390762D+00,   2.047520713006675D+01, &
  2.025052221041162D+01,   1.198683246885931D+01,   8.727943902995086D+00,   7.724452223940203D+00,   3.389278052905926D+00, &
  1.177177617499013D+00,   4.132691499388059D-03,   1.188160317330051D+00,  -6.967354450496733D-01,   4.472844395967870D-01, &
 -3.674955727168822D-02,   8.935021412031723D-02,  -2.305065920791614D-02,   0.000000000000000D+00  /
Data C/                                                                                                                      &
 -2.238951111328112D+06,   3.317258385197909D+06,  -1.907127989685347D+05,   3.226649231850672D+05,  -1.740599286843306D+04, &
  2.726541007985003D+04,  -4.096367147592592D+02,   4.339172341926340D+03,   2.472647954742643D+03,   2.733946148859692D+03, &
  2.156740259935474D+03,   1.818161412216886D+03,   1.125185186621894D+03,   5.192059546975057D+02,   2.271210714997131D+02, &
  8.159373739933278D+01,   2.504213575569288D+01,   4.502220521356484D+00,  -6.485370009925553D-01,  -2.066017586464712D+00, &
 -2.039480515733283D+00,  -1.120885420780364D+00,  -7.787515699700785D-01,  -6.789449806262401D-01,  -2.693302230844084D-01, &
 -8.354063867886183D-02,   5.128742341071029D-03,  -7.623448778836754D-02,   4.249678501069743D-02,  -2.402275728852516D-02, &
  2.111301386872158D-03,  -3.846051795144703D-03,   8.741026807147035D-04,   0.000000000000000D+00  /
Data D/                                                                                                                      &
  7.051676821832495D+05,  -6.947947743674628D+05,   4.177271820778849D+04,  -5.062191798309119D+04,   2.931511218212643D+03, &
 -3.321581435065812D+03,   1.649705496966011D+02,  -3.789053358692815D+02,  -1.829485968424848D+02,  -2.082707619982226D+02, &
 -1.563297390310924D+02,  -1.278932257519952D+02,  -7.332923412010020D+01,  -3.091681775193069D+01,  -1.251812927291563D+01, &
 -4.184586403208442D+00,  -1.216056358564700D+00,  -2.208043102177145D-01,   1.094628337470038D-02,   7.047186473225737D-02, &
  6.942711856937026D-02,   3.538998443752289D-02,   2.341703011025180D-02,   2.010813196179497D-02,   7.207131680797889D-03, &
  2.005775069513611D-03,  -2.283704311394874D-04,   1.635317845686456D-03,  -8.576785703219029D-04,   4.315908346020876D-04, &
 -3.875424258105890D-05,   5.506051098937314D-05,  -1.101210219787461D-05,   0.000000000000000D+00  /
Data ak/-2.9611722908468673D+05/
Data a12/ 8.8120472293789087D+11/
Data rf/ 1.5000000000000000D+01/
Data r0/ 1.3000000000000000D+00/
Data aa/-7.8389640568190327D+04/
Data cc/ 1.6988086355011899D+05/
Data k0/   6/
End
Double Precision Function V_Na_Pi(xx)
Parameter(N=34)
Implicit Real*8(A-H,O-Z)
Integer*4 k0
Common/Param_Spline2_LJ_V_Na_Pi/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
If(xx.le.r0)Then
  V_Na_Pi=aa*xx**2+cc
  Return
Endif
If(xx.gt.rf)Then
  V_Na_Pi=ak/xx**k0 + a12/xx**12
  Return
Endif
Call Findi2(X,xx,N,ix)
V_Na_Pi= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return
End


Block Data Inicio_Cs_Cs
Parameter(N=4350)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_Cs/RCs(N),ECs(N),roCs,ECsmax,rfCs
Data RCs/                                                                            &
 0.33000012D+01,  0.33096679D+01,  0.33193346D+01,  0.33290013D+01,  0.33386680D+01, &
 0.33483348D+01,  0.33580015D+01,  0.33676682D+01,  0.33773349D+01,  0.33870016D+01, &
 0.33966683D+01,  0.34063350D+01,  0.34160018D+01,  0.34256685D+01,  0.34353352D+01, &
 0.34450019D+01,  0.34546686D+01,  0.34643353D+01,  0.34740021D+01,  0.34836688D+01, &
 0.34933355D+01,  0.34944007D+01,  0.34953546D+01,  0.34962360D+01,  0.34970763D+01, &
 0.34979010D+01,  0.34987303D+01,  0.34995805D+01,  0.35004645D+01,  0.35013924D+01, &
 0.35023722D+01,  0.35034099D+01,  0.35045101D+01,  0.35056762D+01,  0.35069104D+01, &
 0.35082144D+01,  0.35095892D+01,  0.35110351D+01,  0.35125523D+01,  0.35141405D+01, &
 0.35157993D+01,  0.35175281D+01,  0.35193264D+01,  0.35211933D+01,  0.35231281D+01, &
 0.35251302D+01,  0.35271989D+01,  0.35293334D+01,  0.35315332D+01,  0.35337978D+01, &
 0.35361267D+01,  0.35385195D+01,  0.35409759D+01,  0.35434956D+01,  0.35460784D+01, &
 0.35487242D+01,  0.35514330D+01,  0.35542046D+01,  0.35570392D+01,  0.35599368D+01, &
 0.35628975D+01,  0.35659215D+01,  0.35690091D+01,  0.35721604D+01,  0.35753756D+01, &
 0.35786552D+01,  0.35819995D+01,  0.35854087D+01,  0.35888833D+01,  0.35924237D+01, &
 0.35960303D+01,  0.35997037D+01,  0.36034442D+01,  0.36072525D+01,  0.36111291D+01, &
 0.36150745D+01,  0.36190895D+01,  0.36231747D+01,  0.36273308D+01,  0.36315585D+01, &
 0.36358586D+01,  0.36402319D+01,  0.36446793D+01,  0.36492017D+01,  0.36538000D+01, &
 0.36584753D+01,  0.36632286D+01,  0.36680610D+01,  0.36729737D+01,  0.36779680D+01, &
 0.36830449D+01,  0.36882060D+01,  0.36934527D+01,  0.36987863D+01,  0.37042084D+01, &
 0.37097206D+01,  0.37153247D+01,  0.37210224D+01,  0.37268154D+01,  0.37327059D+01, &
 0.37386957D+01,  0.37447870D+01,  0.37509821D+01,  0.37572832D+01,  0.37636929D+01, &
 0.37702138D+01,  0.37768486D+01,  0.37836001D+01,  0.37904714D+01,  0.37974658D+01, &
 0.38045866D+01,  0.38118374D+01,  0.38192221D+01,  0.38267446D+01,  0.38344094D+01, &
 0.38422210D+01,  0.38501844D+01,  0.38583047D+01,  0.38665876D+01,  0.38750391D+01, &
 0.38836657D+01,  0.38924743D+01,  0.39014723D+01,  0.39106679D+01,  0.39200697D+01, &
 0.39296873D+01,  0.39395309D+01,  0.39496118D+01,  0.39599420D+01,  0.39705352D+01, &
 0.39814058D+01,  0.39925702D+01,  0.40040463D+01,  0.40158539D+01,  0.40280152D+01, &
 0.40405550D+01,  0.40535014D+01,  0.40668858D+01,  0.40807444D+01,  0.40951184D+01, &
 0.41100557D+01,  0.41256119D+01,  0.41418525D+01,  0.41588558D+01,  0.41767158D+01, &
 0.41955480D+01,  0.42154963D+01,  0.42367439D+01,  0.42595303D+01,  0.42841787D+01, &
 0.43111436D+01,  0.43411000D+01,  0.43751299D+01,  0.44151780D+01,  0.44654999D+01, &
 0.44707116D+01,  0.44757616D+01,  0.44808116D+01,  0.44858616D+01,  0.44909116D+01, &
 0.44959616D+01,  0.45010116D+01,  0.45060616D+01,  0.45111116D+01,  0.45161616D+01, &
 0.45212116D+01,  0.45262616D+01,  0.45313116D+01,  0.45363616D+01,  0.45414116D+01, &
 0.45464616D+01,  0.45515216D+01,  0.45565716D+01,  0.45616216D+01,  0.45666716D+01, &
 0.45717216D+01,  0.45767716D+01,  0.45818216D+01,  0.45868716D+01,  0.45919216D+01, &
 0.45969716D+01,  0.46020216D+01,  0.46070716D+01,  0.46121216D+01,  0.46171716D+01, &
 0.46222216D+01,  0.46272717D+01,  0.46323217D+01,  0.46373717D+01,  0.46424217D+01, &
 0.46474717D+01,  0.46525317D+01,  0.46575817D+01,  0.46626317D+01,  0.46676817D+01, &
 0.46727317D+01,  0.46777817D+01,  0.46828317D+01,  0.46878817D+01,  0.46929317D+01, &
 0.46979817D+01,  0.47030317D+01,  0.47080817D+01,  0.47131317D+01,  0.47181817D+01, &
 0.47232317D+01,  0.47282817D+01,  0.47333317D+01,  0.47383817D+01,  0.47434317D+01, &
 0.47484817D+01,  0.47535417D+01,  0.47585917D+01,  0.47636417D+01,  0.47686917D+01, &
 0.47737417D+01,  0.47787917D+01,  0.47838417D+01,  0.47888917D+01,  0.47939417D+01, &
 0.47989917D+01,  0.48040417D+01,  0.48090917D+01,  0.48141417D+01,  0.48191917D+01, &
 0.48242417D+01,  0.48292917D+01,  0.48343417D+01,  0.48393917D+01,  0.48444417D+01, &
 0.48470530D+01,  0.49085855D+01,  0.49599128D+01,  0.50052932D+01,  0.50466740D+01, &
 0.50851395D+01,  0.51213670D+01,  0.51558132D+01,  0.51888034D+01,  0.52205791D+01, &
 0.52513255D+01,  0.52811885D+01,  0.53102852D+01,  0.53387119D+01,  0.53665484D+01, &
 0.53938621D+01,  0.54207106D+01,  0.54471436D+01,  0.54732041D+01,  0.54989302D+01, &
 0.55243555D+01,  0.55495098D+01,  0.55744201D+01,  0.55991107D+01,  0.56236034D+01, &
 0.56479186D+01,  0.56720748D+01,  0.56960889D+01,  0.57199768D+01,  0.57437534D+01, &
 0.57674323D+01,  0.57910267D+01,  0.58145487D+01,  0.58380100D+01,  0.58614217D+01, &
 0.58847943D+01,  0.59081380D+01,  0.59314625D+01,  0.59547772D+01,  0.59780913D+01, &
 0.60014137D+01,  0.60247529D+01,  0.60481174D+01,  0.60715156D+01,  0.60949556D+01, &
 0.61184455D+01,  0.61419933D+01,  0.61656068D+01,  0.61892940D+01,  0.62130627D+01, &
 0.62369207D+01,  0.62608760D+01,  0.62849365D+01,  0.63091100D+01,  0.63334046D+01, &
 0.63578285D+01,  0.63823898D+01,  0.64070968D+01,  0.64319581D+01,  0.64569823D+01, &
 0.64821781D+01,  0.65075547D+01,  0.65331211D+01,  0.65588869D+01,  0.65848617D+01, &
 0.66110555D+01,  0.66374787D+01,  0.66641417D+01,  0.66910556D+01,  0.67182317D+01, &
 0.67456817D+01,  0.67734177D+01,  0.68014524D+01,  0.68297989D+01,  0.68584707D+01, &
 0.68874823D+01,  0.69168483D+01,  0.69465842D+01,  0.69767061D+01,  0.70072311D+01, &
 0.70381768D+01,  0.70695617D+01,  0.71014053D+01,  0.71337281D+01,  0.71665516D+01, &
 0.71998984D+01,  0.72337924D+01,  0.72682587D+01,  0.73033239D+01,  0.73390161D+01, &
 0.73753650D+01,  0.74124020D+01,  0.74501606D+01,  0.74886761D+01,  0.75279863D+01, &
 0.75681313D+01,  0.76091538D+01,  0.76510994D+01,  0.76940169D+01,  0.77379583D+01, &
 0.77829796D+01,  0.78291407D+01,  0.78765060D+01,  0.79251449D+01,  0.79751320D+01, &
 0.80265481D+01,  0.80794804D+01,  0.81340235D+01,  0.81902799D+01,  0.82483612D+01, &
 0.83083887D+01,  0.83704952D+01,  0.84348253D+01,  0.85015379D+01,  0.85708071D+01, &
 0.86428245D+01,  0.87178014D+01,  0.87959710D+01,  0.88775915D+01,  0.89629493D+01, &
 0.90523630D+01,  0.91461877D+01,  0.92448199D+01,  0.93487039D+01,  0.94583381D+01, &
 0.95742835D+01,  0.96971728D+01,  0.98277212D+01,  0.99667389D+01,  0.10115146D+02, &
 0.10273987D+02,  0.10444455D+02,  0.10627908D+02,  0.10825893D+02,  0.11040176D+02, &
 0.11070073D+02,  0.11099970D+02,  0.11129868D+02,  0.11159765D+02,  0.11189662D+02, &
 0.11219559D+02,  0.11249457D+02,  0.11279354D+02,  0.11309251D+02,  0.11339148D+02, &
 0.11369046D+02,  0.11398943D+02,  0.11428840D+02,  0.11458738D+02,  0.11488635D+02, &
 0.11518532D+02,  0.11548429D+02,  0.11578327D+02,  0.11608224D+02,  0.11638121D+02, &
 0.11668018D+02,  0.11697916D+02,  0.11727813D+02,  0.11757710D+02,  0.11787607D+02, &
 0.11817505D+02,  0.12000004D+02,  0.12500004D+02,  0.13000005D+02,  0.13500005D+02, &
 0.14000005D+02,  0.14500005D+02,  0.15000005D+02,  0.15500006D+02,  0.16000006D+02, &
 0.16500006D+02,  0.17000006D+02,  0.17500006D+02,  0.18000006D+02,  0.18500007D+02, &
 0.19000007D+02,  0.19500007D+02,  0.20000007D+02,  0.20500007D+02,  0.21000007D+02, &
 0.21500008D+02,  0.22000008D+02,  0.22500008D+02,  0.23000008D+02,  0.23500008D+02, &
 0.24000009D+02,  0.24500009D+02,  0.25000009D+02,  0.25665094D+02,  0.25929682D+02, &
 0.26194271D+02,  0.26458859D+02,  0.26723448D+02,  0.26988037D+02,  0.27252625D+02, &
 0.27517214D+02,  0.27781802D+02,  0.28046391D+02,  0.28310980D+02,  0.28575568D+02, &
 0.28840157D+02,  0.29104745D+02,  0.29369334D+02,  0.29633923D+02,  0.29898511D+02, &
 0.30163100D+02,  0.30427688D+02,  0.30692277D+02,  0.30956866D+02,  0.31221454D+02, &
 0.31486043D+02,  0.31750631D+02,  0.32015220D+02,  0.32279809D+02,  0.32544397D+02, &
 0.32808986D+02,  0.33073574D+02,  0.33338163D+02,  0.33602751D+02,  0.33867340D+02, &
 0.34131929D+02,  0.34396517D+02,  0.34661106D+02,  0.34925694D+02,  0.35190283D+02, &
 0.35454872D+02,  0.35719460D+02,  0.35984049D+02,  0.36248637D+02,  0.36513226D+02, &
 0.36777815D+02,  0.37042403D+02,  0.37306992D+02,  0.37571580D+02,  0.37836169D+02, &
 0.38100758D+02,  0.38365346D+02,  0.38629935D+02,  0.38894523D+02,  0.39159112D+02, &
 0.39423701D+02,  0.39688289D+02,  0.39952878D+02,  0.40217466D+02,  0.40482055D+02, &
 0.40746644D+02,  0.41011232D+02,  0.41275821D+02,  0.41540409D+02,  0.41804998D+02, &
 0.42069587D+02,  0.42334175D+02,  0.42598764D+02,  0.42863352D+02,  0.43127941D+02, &
 0.43392529D+02,  0.43657118D+02,  0.43921707D+02,  0.44186295D+02,  0.44450884D+02, &
 0.44715472D+02,  0.44980061D+02,  0.45244650D+02,  0.45509238D+02,  0.45773827D+02, &
 0.46038415D+02,  0.46303004D+02,  0.46567593D+02,  0.46832181D+02,  0.47096770D+02, &
 0.47361358D+02,  0.47625947D+02,  0.47890536D+02,  0.48155124D+02,  0.48419713D+02, &
 0.48684301D+02,  0.48948890D+02,  0.49213479D+02,  0.49478067D+02,  0.49742656D+02, &
 0.50007244D+02,  0.50271833D+02,  0.50536422D+02,  0.50801010D+02,  0.51065599D+02, &
 0.51330187D+02,  0.51594776D+02,  0.51859365D+02,  0.52123953D+02,  0.52388542D+02, &
 0.52653130D+02,  0.52917719D+02,  0.53182307D+02,  0.53446896D+02,  0.53711485D+02, &
 0.53976073D+02,  0.54240662D+02,  0.54505250D+02,  0.54769839D+02,  0.55034428D+02, &
 0.55299016D+02,  0.55563605D+02,  0.55828193D+02,  0.56092782D+02,  0.56357371D+02, &
 0.56621959D+02,  0.56886548D+02,  0.57151136D+02,  0.57415725D+02,  0.57680314D+02, &
 0.57944902D+02,  0.58209491D+02,  0.58474079D+02,  0.58738668D+02,  0.59003257D+02, &
 0.59267845D+02,  0.59532434D+02,  0.59797022D+02,  0.60061611D+02,  0.60326200D+02, &
 0.60590788D+02,  0.60855377D+02,  0.61119965D+02,  0.61384554D+02,  0.61649143D+02, &
 0.61913731D+02,  0.62178320D+02,  0.62442908D+02,  0.62707497D+02,  0.62972085D+02, &
 0.63236674D+02,  0.63501263D+02,  0.63765851D+02,  0.64030440D+02,  0.64295028D+02, &
 0.64559617D+02,  0.64824206D+02,  0.65088794D+02,  0.65353383D+02,  0.65617971D+02, &
 0.65882560D+02,  0.66147149D+02,  0.66411737D+02,  0.66676326D+02,  0.66940914D+02, &
 0.67205503D+02,  0.67470092D+02,  0.67734680D+02,  0.67999269D+02,  0.68263857D+02, &
 0.68528446D+02,  0.68793035D+02,  0.69057623D+02,  0.69322212D+02,  0.69586800D+02, &
 0.69851389D+02,  0.70115978D+02,  0.70380566D+02,  0.70645155D+02,  0.70909743D+02, &
 0.71174332D+02,  0.71438920D+02,  0.71703509D+02,  0.71968098D+02,  0.72232686D+02, &
 0.72497275D+02,  0.72761863D+02,  0.73026452D+02,  0.73291041D+02,  0.73555629D+02, &
 0.73820218D+02,  0.74084806D+02,  0.74349395D+02,  0.74613984D+02,  0.74878572D+02, &
 0.75143161D+02,  0.75407749D+02,  0.75672338D+02,  0.75936927D+02,  0.76201515D+02, &
 0.76466104D+02,  0.76730692D+02,  0.76995281D+02,  0.77259870D+02,  0.77524458D+02, &
 0.77789047D+02,  0.78053635D+02,  0.78318224D+02,  0.78582813D+02,  0.78847401D+02, &
 0.79111990D+02,  0.79376578D+02,  0.79641167D+02,  0.79905756D+02,  0.80170344D+02, &
 0.80434933D+02,  0.80699521D+02,  0.80964110D+02,  0.81228698D+02,  0.81493287D+02, &
 0.81757876D+02,  0.82022464D+02,  0.82287053D+02,  0.82551641D+02,  0.82816230D+02, &
 0.83080819D+02,  0.83345407D+02,  0.83609996D+02,  0.83874584D+02,  0.84139173D+02, &
 0.84403762D+02,  0.84668350D+02,  0.84932939D+02,  0.85197527D+02,  0.85462116D+02, &
 0.85726705D+02,  0.85991293D+02,  0.86255882D+02,  0.86520470D+02,  0.86785059D+02, &
 0.87049648D+02,  0.87314236D+02,  0.87578825D+02,  0.87843413D+02,  0.88108002D+02, &
 0.88372591D+02,  0.88637179D+02,  0.88901768D+02,  0.89166356D+02,  0.89430945D+02, &
 0.89695534D+02,  0.89960122D+02,  0.90224711D+02,  0.90489299D+02,  0.90753888D+02, &
 0.91018476D+02,  0.91283065D+02,  0.91547654D+02,  0.91812242D+02,  0.92076831D+02, &
 0.92341419D+02,  0.92606008D+02,  0.92870597D+02,  0.93135185D+02,  0.93399774D+02, &
 0.93664362D+02,  0.93928951D+02,  0.94193540D+02,  0.94458128D+02,  0.94722717D+02, &
 0.94987305D+02,  0.95251894D+02,  0.95516483D+02,  0.95781071D+02,  0.96045660D+02, &
 0.96310248D+02,  0.96574837D+02,  0.96839426D+02,  0.97104014D+02,  0.97368603D+02, &
 0.97633191D+02,  0.97897780D+02,  0.98162369D+02,  0.98426957D+02,  0.98691546D+02, &
 0.98956134D+02,  0.99220723D+02,  0.99485312D+02,  0.99749900D+02,  0.10001449D+03, &
 0.10027908D+03,  0.10054367D+03,  0.10080825D+03,  0.10107284D+03,  0.10133743D+03, &
 0.10160202D+03,  0.10186661D+03,  0.10213120D+03,  0.10239579D+03,  0.10266037D+03, &
 0.10292496D+03,  0.10318955D+03,  0.10345414D+03,  0.10371873D+03,  0.10398332D+03, &
 0.10424791D+03,  0.10451249D+03,  0.10477708D+03,  0.10504167D+03,  0.10530626D+03, &
 0.10557085D+03,  0.10583544D+03,  0.10610003D+03,  0.10636461D+03,  0.10662920D+03, &
 0.10689379D+03,  0.10715838D+03,  0.10742297D+03,  0.10768756D+03,  0.10795215D+03, &
 0.10821674D+03,  0.10848132D+03,  0.10874591D+03,  0.10901050D+03,  0.10927509D+03, &
 0.10953968D+03,  0.10980427D+03,  0.11006886D+03,  0.11033344D+03,  0.11059803D+03, &
 0.11086262D+03,  0.11112721D+03,  0.11139180D+03,  0.11165639D+03,  0.11192098D+03, &
 0.11218556D+03,  0.11245015D+03,  0.11271474D+03,  0.11297933D+03,  0.11324392D+03, &
 0.11350851D+03,  0.11377310D+03,  0.11403768D+03,  0.11430227D+03,  0.11456686D+03, &
 0.11483145D+03,  0.11509604D+03,  0.11536063D+03,  0.11562522D+03,  0.11588980D+03, &
 0.11615439D+03,  0.11641898D+03,  0.11668357D+03,  0.11694816D+03,  0.11721275D+03, &
 0.11747734D+03,  0.11774192D+03,  0.11800651D+03,  0.11827110D+03,  0.11853569D+03, &
 0.11880028D+03,  0.11906487D+03,  0.11932946D+03,  0.11959404D+03,  0.11985863D+03, &
 0.12012322D+03,  0.12038781D+03,  0.12065240D+03,  0.12091699D+03,  0.12118158D+03, &
 0.12144616D+03,  0.12171075D+03,  0.12197534D+03,  0.12223993D+03,  0.12250452D+03, &
 0.12276911D+03,  0.12303370D+03,  0.12329829D+03,  0.12356287D+03,  0.12382746D+03, &
 0.12409205D+03,  0.12435664D+03,  0.12462123D+03,  0.12488582D+03,  0.12515041D+03, &
 0.12541499D+03,  0.12567958D+03,  0.12594417D+03,  0.12620876D+03,  0.12647335D+03, &
 0.12673794D+03,  0.12700253D+03,  0.12726711D+03,  0.12753170D+03,  0.12779629D+03, &
 0.12806088D+03,  0.12832547D+03,  0.12859006D+03,  0.12885465D+03,  0.12911923D+03, &
 0.12938382D+03,  0.12964841D+03,  0.12991300D+03,  0.13017759D+03,  0.13044218D+03, &
 0.13070677D+03,  0.13097135D+03,  0.13123594D+03,  0.13150053D+03,  0.13176512D+03, &
 0.13202971D+03,  0.13229430D+03,  0.13255889D+03,  0.13282347D+03,  0.13308806D+03, &
 0.13335265D+03,  0.13361724D+03,  0.13388183D+03,  0.13414642D+03,  0.13441101D+03, &
 0.13467559D+03,  0.13494018D+03,  0.13520477D+03,  0.13546936D+03,  0.13573395D+03, &
 0.13599854D+03,  0.13626313D+03,  0.13652771D+03,  0.13679230D+03,  0.13705689D+03, &
 0.13732148D+03,  0.13758607D+03,  0.13785066D+03,  0.13811525D+03,  0.13837983D+03, &
 0.13864442D+03,  0.13890901D+03,  0.13917360D+03,  0.13943819D+03,  0.13970278D+03, &
 0.13996737D+03,  0.14023196D+03,  0.14049654D+03,  0.14076113D+03,  0.14102572D+03, &
 0.14129031D+03,  0.14155490D+03,  0.14181949D+03,  0.14208408D+03,  0.14234866D+03, &
 0.14261325D+03,  0.14287784D+03,  0.14314243D+03,  0.14340702D+03,  0.14367161D+03, &
 0.14393620D+03,  0.14420078D+03,  0.14446537D+03,  0.14472996D+03,  0.14499455D+03, &
 0.14525914D+03,  0.14552373D+03,  0.14578832D+03,  0.14605290D+03,  0.14631749D+03, &
 0.14658208D+03,  0.14684667D+03,  0.14711126D+03,  0.14737585D+03,  0.14764044D+03, &
 0.14790502D+03,  0.14816961D+03,  0.14843420D+03,  0.14869879D+03,  0.14896338D+03, &
 0.14922797D+03,  0.14949256D+03,  0.14975714D+03,  0.15002173D+03,  0.15028632D+03, &
 0.15055091D+03,  0.15081550D+03,  0.15108009D+03,  0.15134468D+03,  0.15160926D+03, &
 0.15187385D+03,  0.15213844D+03,  0.15240303D+03,  0.15266762D+03,  0.15293221D+03, &
 0.15319680D+03,  0.15346138D+03,  0.15372597D+03,  0.15399056D+03,  0.15425515D+03, &
 0.15451974D+03,  0.15478433D+03,  0.15504892D+03,  0.15531350D+03,  0.15557809D+03, &
 0.15584268D+03,  0.15610727D+03,  0.15637186D+03,  0.15663645D+03,  0.15690104D+03, &
 0.15716563D+03,  0.15743021D+03,  0.15769480D+03,  0.15795939D+03,  0.15822398D+03, &
 0.15848857D+03,  0.15875316D+03,  0.15901775D+03,  0.15928233D+03,  0.15954692D+03, &
 0.15981151D+03,  0.16007610D+03,  0.16034069D+03,  0.16060528D+03,  0.16086987D+03, &
 0.16113445D+03,  0.16139904D+03,  0.16166363D+03,  0.16192822D+03,  0.16219281D+03, &
 0.16245740D+03,  0.16272199D+03,  0.16298657D+03,  0.16325116D+03,  0.16351575D+03, &
 0.16378034D+03,  0.16404493D+03,  0.16430952D+03,  0.16457411D+03,  0.16483869D+03, &
 0.16510328D+03,  0.16536787D+03,  0.16563246D+03,  0.16589705D+03,  0.16616164D+03, &
 0.16642623D+03,  0.16669081D+03,  0.16695540D+03,  0.16721999D+03,  0.16748458D+03, &
 0.16774917D+03,  0.16801376D+03,  0.16827835D+03,  0.16854293D+03,  0.16880752D+03, &
 0.16907211D+03,  0.16933670D+03,  0.16960129D+03,  0.16986588D+03,  0.17013047D+03, &
 0.17039505D+03,  0.17065964D+03,  0.17092423D+03,  0.17118882D+03,  0.17145341D+03, &
 0.17171800D+03,  0.17198259D+03,  0.17224717D+03,  0.17251176D+03,  0.17277635D+03, &
 0.17304094D+03,  0.17330553D+03,  0.17357012D+03,  0.17383471D+03,  0.17409930D+03, &
 0.17436388D+03,  0.17462847D+03,  0.17489306D+03,  0.17515765D+03,  0.17542224D+03, &
 0.17568683D+03,  0.17595142D+03,  0.17621600D+03,  0.17648059D+03,  0.17674518D+03, &
 0.17700977D+03,  0.17727436D+03,  0.17753895D+03,  0.17780354D+03,  0.17806812D+03, &
 0.17833271D+03,  0.17859730D+03,  0.17886189D+03,  0.17912648D+03,  0.17939107D+03, &
 0.17965566D+03,  0.17992024D+03,  0.18018483D+03,  0.18044942D+03,  0.18071401D+03, &
 0.18097860D+03,  0.18124319D+03,  0.18150778D+03,  0.18177236D+03,  0.18203695D+03, &
 0.18230154D+03,  0.18256613D+03,  0.18283072D+03,  0.18309531D+03,  0.18335990D+03, &
 0.18362448D+03,  0.18388907D+03,  0.18415366D+03,  0.18441825D+03,  0.18468284D+03, &
 0.18494743D+03,  0.18521202D+03,  0.18547660D+03,  0.18574119D+03,  0.18600578D+03, &
 0.18627037D+03,  0.18653496D+03,  0.18679955D+03,  0.18706414D+03,  0.18732872D+03, &
 0.18759331D+03,  0.18785790D+03,  0.18812249D+03,  0.18838708D+03,  0.18865167D+03, &
 0.18891626D+03,  0.18918085D+03,  0.18944543D+03,  0.18971002D+03,  0.18997461D+03, &
 0.19023920D+03,  0.19050379D+03,  0.19076838D+03,  0.19103297D+03,  0.19129755D+03, &
 0.19156214D+03,  0.19182673D+03,  0.19209132D+03,  0.19235591D+03,  0.19262050D+03, &
 0.19288509D+03,  0.19314967D+03,  0.19341426D+03,  0.19367885D+03,  0.19394344D+03, &
 0.19420803D+03,  0.19447262D+03,  0.19473721D+03,  0.19500179D+03,  0.19526638D+03, &
 0.19553097D+03,  0.19579556D+03,  0.19606015D+03,  0.19632474D+03,  0.19658933D+03, &
 0.19685391D+03,  0.19711850D+03,  0.19738309D+03,  0.19764768D+03,  0.19791227D+03, &
 0.19817686D+03,  0.19844145D+03,  0.19870603D+03,  0.19897062D+03,  0.19923521D+03, &
 0.19949980D+03,  0.19976439D+03,  0.20002898D+03,  0.20029357D+03,  0.20055815D+03, &
 0.20082274D+03,  0.20108733D+03,  0.20135192D+03,  0.20161651D+03,  0.20188110D+03, &
 0.20214569D+03,  0.20241027D+03,  0.20267486D+03,  0.20293945D+03,  0.20320404D+03, &
 0.20346863D+03,  0.20373322D+03,  0.20399781D+03,  0.20426239D+03,  0.20452698D+03, &
 0.20479157D+03,  0.20505616D+03,  0.20532075D+03,  0.20558534D+03,  0.20584993D+03, &
 0.20611452D+03,  0.20637910D+03,  0.20664369D+03,  0.20690828D+03,  0.20717287D+03, &
 0.20743746D+03,  0.20770205D+03,  0.20796664D+03,  0.20823122D+03,  0.20849581D+03, &
 0.20876040D+03,  0.20902499D+03,  0.20928958D+03,  0.20955417D+03,  0.20981876D+03, &
 0.21008334D+03,  0.21034793D+03,  0.21061252D+03,  0.21087711D+03,  0.21114170D+03, &
 0.21140629D+03,  0.21167088D+03,  0.21193546D+03,  0.21220005D+03,  0.21246464D+03, &
 0.21272923D+03,  0.21299382D+03,  0.21325841D+03,  0.21352300D+03,  0.21378758D+03, &
 0.21405217D+03,  0.21431676D+03,  0.21458135D+03,  0.21484594D+03,  0.21511053D+03, &
 0.21537512D+03,  0.21563970D+03,  0.21590429D+03,  0.21616888D+03,  0.21643347D+03, &
 0.21669806D+03,  0.21696265D+03,  0.21722724D+03,  0.21749182D+03,  0.21775641D+03, &
 0.21802100D+03,  0.21828559D+03,  0.21855018D+03,  0.21881477D+03,  0.21907936D+03, &
 0.21934394D+03,  0.21960853D+03,  0.21987312D+03,  0.22013771D+03,  0.22040230D+03, &
 0.22066689D+03,  0.22093148D+03,  0.22119606D+03,  0.22146065D+03,  0.22172524D+03, &
 0.22198983D+03,  0.22225442D+03,  0.22251901D+03,  0.22278360D+03,  0.22304819D+03, &
 0.22331277D+03,  0.22357736D+03,  0.22384195D+03,  0.22410654D+03,  0.22437113D+03, &
 0.22463572D+03,  0.22490031D+03,  0.22516489D+03,  0.22542948D+03,  0.22569407D+03, &
 0.22595866D+03,  0.22622325D+03,  0.22648784D+03,  0.22675243D+03,  0.22701701D+03, &
 0.22728160D+03,  0.22754619D+03,  0.22781078D+03,  0.22807537D+03,  0.22833996D+03, &
 0.22860455D+03,  0.22886913D+03,  0.22913372D+03,  0.22939831D+03,  0.22966290D+03, &
 0.22992749D+03,  0.23019208D+03,  0.23045667D+03,  0.23072125D+03,  0.23098584D+03, &
 0.23125043D+03,  0.23151502D+03,  0.23177961D+03,  0.23204420D+03,  0.23230879D+03, &
 0.23257337D+03,  0.23283796D+03,  0.23310255D+03,  0.23336714D+03,  0.23363173D+03, &
 0.23389632D+03,  0.23416091D+03,  0.23442549D+03,  0.23469008D+03,  0.23495467D+03, &
 0.23521926D+03,  0.23548385D+03,  0.23574844D+03,  0.23601303D+03,  0.23627761D+03, &
 0.23654220D+03,  0.23680679D+03,  0.23707138D+03,  0.23733597D+03,  0.23760056D+03, &
 0.23786515D+03,  0.23812973D+03,  0.23839432D+03,  0.23865891D+03,  0.23892350D+03, &
 0.23918809D+03,  0.23945268D+03,  0.23971727D+03,  0.23998186D+03,  0.24024644D+03, &
 0.24051103D+03,  0.24077562D+03,  0.24104021D+03,  0.24130480D+03,  0.24156939D+03, &
 0.24183398D+03,  0.24209856D+03,  0.24236315D+03,  0.24262774D+03,  0.24289233D+03, &
 0.24315692D+03,  0.24342151D+03,  0.24368610D+03,  0.24395068D+03,  0.24421527D+03, &
 0.24447986D+03,  0.24474445D+03,  0.24500904D+03,  0.24527363D+03,  0.24553822D+03, &
 0.24580280D+03,  0.24606739D+03,  0.24633198D+03,  0.24659657D+03,  0.24686116D+03, &
 0.24712575D+03,  0.24739034D+03,  0.24765492D+03,  0.24791951D+03,  0.24818410D+03, &
 0.24844869D+03,  0.24871328D+03,  0.24897787D+03,  0.24924246D+03,  0.24950704D+03, &
 0.24977163D+03,  0.25003622D+03,  0.25030081D+03,  0.25056540D+03,  0.25082999D+03, &
 0.25109458D+03,  0.25135916D+03,  0.25162375D+03,  0.25188834D+03,  0.25215293D+03, &
 0.25241752D+03,  0.25268211D+03,  0.25294670D+03,  0.25321128D+03,  0.25347587D+03, &
 0.25374046D+03,  0.25400505D+03,  0.25426964D+03,  0.25453423D+03,  0.25479882D+03, &
 0.25506341D+03,  0.25532799D+03,  0.25559258D+03,  0.25585717D+03,  0.25612176D+03, &
 0.25638635D+03,  0.25665094D+03,  0.25691553D+03,  0.25718011D+03,  0.25744470D+03, &
 0.25770929D+03,  0.25797388D+03,  0.25823847D+03,  0.25850306D+03,  0.25876765D+03, &
 0.25903223D+03,  0.25929682D+03,  0.25956141D+03,  0.25982600D+03,  0.26009059D+03, &
 0.26035518D+03,  0.26061977D+03,  0.26088435D+03,  0.26114894D+03,  0.26141353D+03, &
 0.26167812D+03,  0.26194271D+03,  0.26220730D+03,  0.26247189D+03,  0.26273647D+03, &
 0.26300106D+03,  0.26326565D+03,  0.26353024D+03,  0.26379483D+03,  0.26405942D+03, &
 0.26432401D+03,  0.26458859D+03,  0.26485318D+03,  0.26511777D+03,  0.26538236D+03, &
 0.26564695D+03,  0.26591154D+03,  0.26617613D+03,  0.26644071D+03,  0.26670530D+03, &
 0.26696989D+03,  0.26723448D+03,  0.26749907D+03,  0.26776366D+03,  0.26802825D+03, &
 0.26829283D+03,  0.26855742D+03,  0.26882201D+03,  0.26908660D+03,  0.26935119D+03, &
 0.26961578D+03,  0.26988037D+03,  0.27014495D+03,  0.27040954D+03,  0.27067413D+03, &
 0.27093872D+03,  0.27120331D+03,  0.27146790D+03,  0.27173249D+03,  0.27199708D+03, &
 0.27226166D+03,  0.27252625D+03,  0.27279084D+03,  0.27305543D+03,  0.27332002D+03, &
 0.27358461D+03,  0.27384920D+03,  0.27411378D+03,  0.27437837D+03,  0.27464296D+03, &
 0.27490755D+03,  0.27517214D+03,  0.27543673D+03,  0.27570132D+03,  0.27596590D+03, &
 0.27623049D+03,  0.27649508D+03,  0.27675967D+03,  0.27702426D+03,  0.27728885D+03, &
 0.27755344D+03,  0.27781802D+03,  0.27808261D+03,  0.27834720D+03,  0.27861179D+03, &
 0.27887638D+03,  0.27914097D+03,  0.27940556D+03,  0.27967014D+03,  0.27993473D+03, &
 0.28019932D+03,  0.28046391D+03,  0.28072850D+03,  0.28099309D+03,  0.28125768D+03, &
 0.28152226D+03,  0.28178685D+03,  0.28205144D+03,  0.28231603D+03,  0.28258062D+03, &
 0.28284521D+03,  0.28310980D+03,  0.28337438D+03,  0.28363897D+03,  0.28390356D+03, &
 0.28416815D+03,  0.28443274D+03,  0.28469733D+03,  0.28496192D+03,  0.28522650D+03, &
 0.28549109D+03,  0.28575568D+03,  0.28602027D+03,  0.28628486D+03,  0.28654945D+03, &
 0.28681404D+03,  0.28707862D+03,  0.28734321D+03,  0.28760780D+03,  0.28787239D+03, &
 0.28813698D+03,  0.28840157D+03,  0.28866616D+03,  0.28893075D+03,  0.28919533D+03, &
 0.28945992D+03,  0.28972451D+03,  0.28998910D+03,  0.29025369D+03,  0.29051828D+03, &
 0.29078287D+03,  0.29104745D+03,  0.29131204D+03,  0.29157663D+03,  0.29184122D+03, &
 0.29210581D+03,  0.29237040D+03,  0.29263499D+03,  0.29289957D+03,  0.29316416D+03, &
 0.29342875D+03,  0.29369334D+03,  0.29395793D+03,  0.29422252D+03,  0.29448711D+03, &
 0.29475169D+03,  0.29501628D+03,  0.29528087D+03,  0.29554546D+03,  0.29581005D+03, &
 0.29607464D+03,  0.29633923D+03,  0.29660381D+03,  0.29686840D+03,  0.29713299D+03, &
 0.29739758D+03,  0.29766217D+03,  0.29792676D+03,  0.29819135D+03,  0.29845593D+03, &
 0.29872052D+03,  0.29898511D+03,  0.29924970D+03,  0.29951429D+03,  0.29977888D+03, &
 0.30004347D+03,  0.30030805D+03,  0.30057264D+03,  0.30083723D+03,  0.30110182D+03, &
 0.30136641D+03,  0.30163100D+03,  0.30189559D+03,  0.30216017D+03,  0.30242476D+03, &
 0.30268935D+03,  0.30295394D+03,  0.30321853D+03,  0.30348312D+03,  0.30374771D+03, &
 0.30401230D+03,  0.30427688D+03,  0.30454147D+03,  0.30480606D+03,  0.30507065D+03, &
 0.30533524D+03,  0.30559983D+03,  0.30586442D+03,  0.30612900D+03,  0.30639359D+03, &
 0.30665818D+03,  0.30692277D+03,  0.30718736D+03,  0.30745195D+03,  0.30771654D+03, &
 0.30798112D+03,  0.30824571D+03,  0.30851030D+03,  0.30877489D+03,  0.30903948D+03, &
 0.30930407D+03,  0.30956866D+03,  0.30983324D+03,  0.31009783D+03,  0.31036242D+03, &
 0.31062701D+03,  0.31089160D+03,  0.31115619D+03,  0.31142078D+03,  0.31168536D+03, &
 0.31194995D+03,  0.31221454D+03,  0.31247913D+03,  0.31274372D+03,  0.31300831D+03, &
 0.31327290D+03,  0.31353748D+03,  0.31380207D+03,  0.31406666D+03,  0.31433125D+03, &
 0.31459584D+03,  0.31486043D+03,  0.31512502D+03,  0.31538960D+03,  0.31565419D+03, &
 0.31591878D+03,  0.31618337D+03,  0.31644796D+03,  0.31671255D+03,  0.31697714D+03, &
 0.31724172D+03,  0.31750631D+03,  0.31777090D+03,  0.31803549D+03,  0.31830008D+03, &
 0.31856467D+03,  0.31882926D+03,  0.31909384D+03,  0.31935843D+03,  0.31962302D+03, &
 0.31988761D+03,  0.32015220D+03,  0.32041679D+03,  0.32068138D+03,  0.32094597D+03, &
 0.32121055D+03,  0.32147514D+03,  0.32173973D+03,  0.32200432D+03,  0.32226891D+03, &
 0.32253350D+03,  0.32279809D+03,  0.32306267D+03,  0.32332726D+03,  0.32359185D+03, &
 0.32385644D+03,  0.32412103D+03,  0.32438562D+03,  0.32465021D+03,  0.32491479D+03, &
 0.32517938D+03,  0.32544397D+03,  0.32570856D+03,  0.32597315D+03,  0.32623774D+03, &
 0.32650233D+03,  0.32676691D+03,  0.32703150D+03,  0.32729609D+03,  0.32756068D+03, &
 0.32782527D+03,  0.32808986D+03,  0.32835445D+03,  0.32861903D+03,  0.32888362D+03, &
 0.32914821D+03,  0.32941280D+03,  0.32967739D+03,  0.32994198D+03,  0.33020657D+03, &
 0.33047115D+03,  0.33073574D+03,  0.33100033D+03,  0.33126492D+03,  0.33152951D+03, &
 0.33179410D+03,  0.33205869D+03,  0.33232327D+03,  0.33258786D+03,  0.33285245D+03, &
 0.33311704D+03,  0.33338163D+03,  0.33364622D+03,  0.33391081D+03,  0.33417539D+03, &
 0.33443998D+03,  0.33470457D+03,  0.33496916D+03,  0.33523375D+03,  0.33549834D+03, &
 0.33576293D+03,  0.33602751D+03,  0.33629210D+03,  0.33655669D+03,  0.33682128D+03, &
 0.33708587D+03,  0.33735046D+03,  0.33761505D+03,  0.33787964D+03,  0.33814422D+03, &
 0.33840881D+03,  0.33867340D+03,  0.33893799D+03,  0.33920258D+03,  0.33946717D+03, &
 0.33973176D+03,  0.33999634D+03,  0.34026093D+03,  0.34052552D+03,  0.34079011D+03, &
 0.34105470D+03,  0.34131929D+03,  0.34158388D+03,  0.34184846D+03,  0.34211305D+03, &
 0.34237764D+03,  0.34264223D+03,  0.34290682D+03,  0.34317141D+03,  0.34343600D+03, &
 0.34370058D+03,  0.34396517D+03,  0.34422976D+03,  0.34449435D+03,  0.34475894D+03, &
 0.34502353D+03,  0.34528812D+03,  0.34555270D+03,  0.34581729D+03,  0.34608188D+03, &
 0.34634647D+03,  0.34661106D+03,  0.34687565D+03,  0.34714024D+03,  0.34740482D+03, &
 0.34766941D+03,  0.34793400D+03,  0.34819859D+03,  0.34846318D+03,  0.34872777D+03, &
 0.34899236D+03,  0.34925694D+03,  0.34952153D+03,  0.34978612D+03,  0.35005071D+03, &
 0.35031530D+03,  0.35057989D+03,  0.35084448D+03,  0.35110906D+03,  0.35137365D+03, &
 0.35163824D+03,  0.35190283D+03,  0.35216742D+03,  0.35243201D+03,  0.35269660D+03, &
 0.35296118D+03,  0.35322577D+03,  0.35349036D+03,  0.35375495D+03,  0.35401954D+03, &
 0.35428413D+03,  0.35454872D+03,  0.35481331D+03,  0.35507789D+03,  0.35534248D+03, &
 0.35560707D+03,  0.35587166D+03,  0.35613625D+03,  0.35640084D+03,  0.35666543D+03, &
 0.35693001D+03,  0.35719460D+03,  0.35745919D+03,  0.35772378D+03,  0.35798837D+03, &
 0.35825296D+03,  0.35851755D+03,  0.35878213D+03,  0.35904672D+03,  0.35931131D+03, &
 0.35957590D+03,  0.35984049D+03,  0.36010508D+03,  0.36036967D+03,  0.36063425D+03, &
 0.36089884D+03,  0.36116343D+03,  0.36142802D+03,  0.36169261D+03,  0.36195720D+03, &
 0.36222179D+03,  0.36248637D+03,  0.36275096D+03,  0.36301555D+03,  0.36328014D+03, &
 0.36354473D+03,  0.36380932D+03,  0.36407391D+03,  0.36433849D+03,  0.36460308D+03, &
 0.36486767D+03,  0.36513226D+03,  0.36539685D+03,  0.36566144D+03,  0.36592603D+03, &
 0.36619061D+03,  0.36645520D+03,  0.36671979D+03,  0.36698438D+03,  0.36724897D+03, &
 0.36751356D+03,  0.36777815D+03,  0.36804273D+03,  0.36830732D+03,  0.36857191D+03, &
 0.36883650D+03,  0.36910109D+03,  0.36936568D+03,  0.36963027D+03,  0.36989486D+03, &
 0.37015944D+03,  0.37042403D+03,  0.37068862D+03,  0.37095321D+03,  0.37121780D+03, &
 0.37148239D+03,  0.37174698D+03,  0.37201156D+03,  0.37227615D+03,  0.37254074D+03, &
 0.37280533D+03,  0.37306992D+03,  0.37333451D+03,  0.37359910D+03,  0.37386368D+03, &
 0.37412827D+03,  0.37439286D+03,  0.37465745D+03,  0.37492204D+03,  0.37518663D+03, &
 0.37545122D+03,  0.37571580D+03,  0.37598039D+03,  0.37624498D+03,  0.37650957D+03, &
 0.37677416D+03,  0.37703875D+03,  0.37730334D+03,  0.37756792D+03,  0.37783251D+03, &
 0.37809710D+03,  0.37836169D+03,  0.37862628D+03,  0.37889087D+03,  0.37915546D+03, &
 0.37942004D+03,  0.37968463D+03,  0.37994922D+03,  0.38021381D+03,  0.38047840D+03, &
 0.38074299D+03,  0.38100758D+03,  0.38127216D+03,  0.38153675D+03,  0.38180134D+03, &
 0.38206593D+03,  0.38233052D+03,  0.38259511D+03,  0.38285970D+03,  0.38312428D+03, &
 0.38338887D+03,  0.38365346D+03,  0.38391805D+03,  0.38418264D+03,  0.38444723D+03, &
 0.38471182D+03,  0.38497640D+03,  0.38524099D+03,  0.38550558D+03,  0.38577017D+03, &
 0.38603476D+03,  0.38629935D+03,  0.38656394D+03,  0.38682853D+03,  0.38709311D+03, &
 0.38735770D+03,  0.38762229D+03,  0.38788688D+03,  0.38815147D+03,  0.38841606D+03, &
 0.38868065D+03,  0.38894523D+03,  0.38920982D+03,  0.38947441D+03,  0.38973900D+03, &
 0.39000359D+03,  0.39026818D+03,  0.39053277D+03,  0.39079735D+03,  0.39106194D+03, &
 0.39132653D+03,  0.39159112D+03,  0.39185571D+03,  0.39212030D+03,  0.39238489D+03, &
 0.39264947D+03,  0.39291406D+03,  0.39317865D+03,  0.39344324D+03,  0.39370783D+03, &
 0.39397242D+03,  0.39423701D+03,  0.39450159D+03,  0.39476618D+03,  0.39503077D+03, &
 0.39529536D+03,  0.39555995D+03,  0.39582454D+03,  0.39608913D+03,  0.39635371D+03, &
 0.39661830D+03,  0.39688289D+03,  0.39714748D+03,  0.39741207D+03,  0.39767666D+03, &
 0.39794125D+03,  0.39820583D+03,  0.39847042D+03,  0.39873501D+03,  0.39899960D+03, &
 0.39926419D+03,  0.39952878D+03,  0.39979337D+03,  0.40005795D+03,  0.40032254D+03, &
 0.40058713D+03,  0.40085172D+03,  0.40111631D+03,  0.40138090D+03,  0.40164549D+03, &
 0.40191007D+03,  0.40217466D+03,  0.40243925D+03,  0.40270384D+03,  0.40296843D+03, &
 0.40323302D+03,  0.40349761D+03,  0.40376220D+03,  0.40402678D+03,  0.40429137D+03, &
 0.40455596D+03,  0.40482055D+03,  0.40508514D+03,  0.40534973D+03,  0.40561432D+03, &
 0.40587890D+03,  0.40614349D+03,  0.40640808D+03,  0.40667267D+03,  0.40693726D+03, &
 0.40720185D+03,  0.40746644D+03,  0.40773102D+03,  0.40799561D+03,  0.40826020D+03, &
 0.40852479D+03,  0.40878938D+03,  0.40905397D+03,  0.40931856D+03,  0.40958314D+03, &
 0.40984773D+03,  0.41011232D+03,  0.41037691D+03,  0.41064150D+03,  0.41090609D+03, &
 0.41117068D+03,  0.41143526D+03,  0.41169985D+03,  0.41196444D+03,  0.41222903D+03, &
 0.41249362D+03,  0.41275821D+03,  0.41302280D+03,  0.41328738D+03,  0.41355197D+03, &
 0.41381656D+03,  0.41408115D+03,  0.41434574D+03,  0.41461033D+03,  0.41487492D+03, &
 0.41513950D+03,  0.41540409D+03,  0.41566868D+03,  0.41593327D+03,  0.41619786D+03, &
 0.41646245D+03,  0.41672704D+03,  0.41699162D+03,  0.41725621D+03,  0.41752080D+03, &
 0.41778539D+03,  0.41804998D+03,  0.41831457D+03,  0.41857916D+03,  0.41884374D+03, &
 0.41910833D+03,  0.41937292D+03,  0.41963751D+03,  0.41990210D+03,  0.42016669D+03, &
 0.42043128D+03,  0.42069587D+03,  0.42096045D+03,  0.42122504D+03,  0.42148963D+03, &
 0.42175422D+03,  0.42201881D+03,  0.42228340D+03,  0.42254799D+03,  0.42281257D+03, &
 0.42307716D+03,  0.42334175D+03,  0.42360634D+03,  0.42387093D+03,  0.42413552D+03, &
 0.42440011D+03,  0.42466469D+03,  0.42492928D+03,  0.42519387D+03,  0.42545846D+03, &
 0.42572305D+03,  0.42598764D+03,  0.42625223D+03,  0.42651681D+03,  0.42678140D+03, &
 0.42704599D+03,  0.42731058D+03,  0.42757517D+03,  0.42783976D+03,  0.42810435D+03, &
 0.42836893D+03,  0.42863352D+03,  0.42889811D+03,  0.42916270D+03,  0.42942729D+03, &
 0.42969188D+03,  0.42995647D+03,  0.43022105D+03,  0.43048564D+03,  0.43075023D+03, &
 0.43101482D+03,  0.43127941D+03,  0.43154400D+03,  0.43180859D+03,  0.43207317D+03, &
 0.43233776D+03,  0.43260235D+03,  0.43286694D+03,  0.43313153D+03,  0.43339612D+03, &
 0.43366071D+03,  0.43392529D+03,  0.43418988D+03,  0.43445447D+03,  0.43471906D+03, &
 0.43498365D+03,  0.43524824D+03,  0.43551283D+03,  0.43577742D+03,  0.43604200D+03, &
 0.43630659D+03,  0.43657118D+03,  0.43683577D+03,  0.43710036D+03,  0.43736495D+03, &
 0.43762954D+03,  0.43789412D+03,  0.43815871D+03,  0.43842330D+03,  0.43868789D+03, &
 0.43895248D+03,  0.43921707D+03,  0.43948166D+03,  0.43974624D+03,  0.44001083D+03, &
 0.44027542D+03,  0.44054001D+03,  0.44080460D+03,  0.44106919D+03,  0.44133378D+03, &
 0.44159836D+03,  0.44186295D+03,  0.44212754D+03,  0.44239213D+03,  0.44265672D+03, &
 0.44292131D+03,  0.44318590D+03,  0.44345048D+03,  0.44371507D+03,  0.44397966D+03, &
 0.44424425D+03,  0.44450884D+03,  0.44477343D+03,  0.44503802D+03,  0.44530260D+03, &
 0.44556719D+03,  0.44583178D+03,  0.44609637D+03,  0.44636096D+03,  0.44662555D+03, &
 0.44689014D+03,  0.44715472D+03,  0.44741931D+03,  0.44768390D+03,  0.44794849D+03, &
 0.44821308D+03,  0.44847767D+03,  0.44874226D+03,  0.44900684D+03,  0.44927143D+03, &
 0.44953602D+03,  0.44980061D+03,  0.45006520D+03,  0.45032979D+03,  0.45059438D+03, &
 0.45085896D+03,  0.45112355D+03,  0.45138814D+03,  0.45165273D+03,  0.45191732D+03, &
 0.45218191D+03,  0.45244650D+03,  0.45271109D+03,  0.45297567D+03,  0.45324026D+03, &
 0.45350485D+03,  0.45376944D+03,  0.45403403D+03,  0.45429862D+03,  0.45456321D+03, &
 0.45482779D+03,  0.45509238D+03,  0.45535697D+03,  0.45562156D+03,  0.45588615D+03, &
 0.45615074D+03,  0.45641533D+03,  0.45667991D+03,  0.45694450D+03,  0.45720909D+03, &
 0.45747368D+03,  0.45773827D+03,  0.45800286D+03,  0.45826745D+03,  0.45853203D+03, &
 0.45879662D+03,  0.45906121D+03,  0.45932580D+03,  0.45959039D+03,  0.45985498D+03, &
 0.46011957D+03,  0.46038415D+03,  0.46064874D+03,  0.46091333D+03,  0.46117792D+03, &
 0.46144251D+03,  0.46170710D+03,  0.46197169D+03,  0.46223627D+03,  0.46250086D+03, &
 0.46276545D+03,  0.46303004D+03,  0.46329463D+03,  0.46355922D+03,  0.46382381D+03, &
 0.46408839D+03,  0.46435298D+03,  0.46461757D+03,  0.46488216D+03,  0.46514675D+03, &
 0.46541134D+03,  0.46567593D+03,  0.46594051D+03,  0.46620510D+03,  0.46646969D+03, &
 0.46673428D+03,  0.46699887D+03,  0.46726346D+03,  0.46752805D+03,  0.46779263D+03, &
 0.46805722D+03,  0.46832181D+03,  0.46858640D+03,  0.46885099D+03,  0.46911558D+03, &
 0.46938017D+03,  0.46964476D+03,  0.46990934D+03,  0.47017393D+03,  0.47043852D+03, &
 0.47070311D+03,  0.47096770D+03,  0.47123229D+03,  0.47149688D+03,  0.47176146D+03, &
 0.47202605D+03,  0.47229064D+03,  0.47255523D+03,  0.47281982D+03,  0.47308441D+03, &
 0.47334900D+03,  0.47361358D+03,  0.47387817D+03,  0.47414276D+03,  0.47440735D+03, &
 0.47467194D+03,  0.47493653D+03,  0.47520112D+03,  0.47546570D+03,  0.47573029D+03, &
 0.47599488D+03,  0.47625947D+03,  0.47652406D+03,  0.47678865D+03,  0.47705324D+03, &
 0.47731782D+03,  0.47758241D+03,  0.47784700D+03,  0.47811159D+03,  0.47837618D+03, &
 0.47864077D+03,  0.47890536D+03,  0.47916994D+03,  0.47943453D+03,  0.47969912D+03, &
 0.47996371D+03,  0.48022830D+03,  0.48049289D+03,  0.48075748D+03,  0.48102206D+03, &
 0.48128665D+03,  0.48155124D+03,  0.48181583D+03,  0.48208042D+03,  0.48234501D+03, &
 0.48260960D+03,  0.48287418D+03,  0.48313877D+03,  0.48340336D+03,  0.48366795D+03, &
 0.48393254D+03,  0.48419713D+03,  0.48446172D+03,  0.48472631D+03,  0.48499089D+03, &
 0.48525548D+03,  0.48552007D+03,  0.48578466D+03,  0.48604925D+03,  0.48631384D+03, &
 0.48657843D+03,  0.48684301D+03,  0.48710760D+03,  0.48737219D+03,  0.48763678D+03, &
 0.48790137D+03,  0.48816596D+03,  0.48843055D+03,  0.48869513D+03,  0.48895972D+03, &
 0.48922431D+03,  0.48948890D+03,  0.48975349D+03,  0.49001808D+03,  0.49028267D+03, &
 0.49054725D+03,  0.49081184D+03,  0.49107643D+03,  0.49134102D+03,  0.49160561D+03, &
 0.49187020D+03,  0.49213479D+03,  0.49239937D+03,  0.49266396D+03,  0.49292855D+03, &
 0.49319314D+03,  0.49345773D+03,  0.49372232D+03,  0.49398691D+03,  0.49425149D+03, &
 0.49451608D+03,  0.49478067D+03,  0.49504526D+03,  0.49530985D+03,  0.49557444D+03, &
 0.49583903D+03,  0.49610361D+03,  0.49636820D+03,  0.49663279D+03,  0.49689738D+03, &
 0.49716197D+03,  0.49742656D+03,  0.49769115D+03,  0.49795573D+03,  0.49822032D+03, &
 0.49848491D+03,  0.49874950D+03,  0.49901409D+03,  0.49927868D+03,  0.49954327D+03, &
 0.49980785D+03,  0.50007244D+03,  0.50033703D+03,  0.50060162D+03,  0.50086621D+03, &
 0.50113080D+03,  0.50139539D+03,  0.50165998D+03,  0.50192456D+03,  0.50218915D+03, &
 0.50245374D+03,  0.50271833D+03,  0.50298292D+03,  0.50324751D+03,  0.50351210D+03, &
 0.50377668D+03,  0.50404127D+03,  0.50430586D+03,  0.50457045D+03,  0.50483504D+03, &
 0.50509963D+03,  0.50536422D+03,  0.50562880D+03,  0.50589339D+03,  0.50615798D+03, &
 0.50642257D+03,  0.50668716D+03,  0.50695175D+03,  0.50721634D+03,  0.50748092D+03, &
 0.50774551D+03,  0.50801010D+03,  0.50827469D+03,  0.50853928D+03,  0.50880387D+03, &
 0.50906846D+03,  0.50933304D+03,  0.50959763D+03,  0.50986222D+03,  0.51012681D+03, &
 0.51039140D+03,  0.51065599D+03,  0.51092058D+03,  0.51118516D+03,  0.51144975D+03, &
 0.51171434D+03,  0.51197893D+03,  0.51224352D+03,  0.51250811D+03,  0.51277270D+03, &
 0.51303728D+03,  0.51330187D+03,  0.51356646D+03,  0.51383105D+03,  0.51409564D+03, &
 0.51436023D+03,  0.51462482D+03,  0.51488940D+03,  0.51515399D+03,  0.51541858D+03, &
 0.51568317D+03,  0.51594776D+03,  0.51621235D+03,  0.51647694D+03,  0.51674152D+03, &
 0.51700611D+03,  0.51727070D+03,  0.51753529D+03,  0.51779988D+03,  0.51806447D+03, &
 0.51832906D+03,  0.51859365D+03,  0.51885823D+03,  0.51912282D+03,  0.51938741D+03, &
 0.51965200D+03,  0.51991659D+03,  0.52018118D+03,  0.52044577D+03,  0.52071035D+03, &
 0.52097494D+03,  0.52123953D+03,  0.52150412D+03,  0.52176871D+03,  0.52203330D+03, &
 0.52229789D+03,  0.52256247D+03,  0.52282706D+03,  0.52309165D+03,  0.52335624D+03, &
 0.52362083D+03,  0.52388542D+03,  0.52415001D+03,  0.52441459D+03,  0.52467918D+03, &
 0.52494377D+03,  0.52520836D+03,  0.52547295D+03,  0.52573754D+03,  0.52600213D+03, &
 0.52626671D+03,  0.52653130D+03,  0.52679589D+03,  0.52706048D+03,  0.52732507D+03, &
 0.52758966D+03,  0.52785425D+03,  0.52811883D+03,  0.52838342D+03,  0.52864801D+03, &
 0.52891260D+03,  0.52917719D+03,  0.52944178D+03,  0.52970637D+03,  0.52997095D+03, &
 0.53023554D+03,  0.53050013D+03,  0.53076472D+03,  0.53102931D+03,  0.53129390D+03, &
 0.53155849D+03,  0.53182307D+03,  0.53208766D+03,  0.53235225D+03,  0.53261684D+03, &
 0.53288143D+03,  0.53314602D+03,  0.53341061D+03,  0.53367519D+03,  0.53393978D+03, &
 0.53420437D+03,  0.53446896D+03,  0.53473355D+03,  0.53499814D+03,  0.53526273D+03, &
 0.53552732D+03,  0.53579190D+03,  0.53605649D+03,  0.53632108D+03,  0.53658567D+03, &
 0.53685026D+03,  0.53711485D+03,  0.53737944D+03,  0.53764402D+03,  0.53790861D+03, &
 0.53817320D+03,  0.53843779D+03,  0.53870238D+03,  0.53896697D+03,  0.53923156D+03, &
 0.53949614D+03,  0.53976073D+03,  0.54002532D+03,  0.54028991D+03,  0.54055450D+03, &
 0.54081909D+03,  0.54108368D+03,  0.54134826D+03,  0.54161285D+03,  0.54187744D+03, &
 0.54214203D+03,  0.54240662D+03,  0.54267121D+03,  0.54293580D+03,  0.54320038D+03, &
 0.54346497D+03,  0.54372956D+03,  0.54399415D+03,  0.54425874D+03,  0.54452333D+03, &
 0.54478792D+03,  0.54505250D+03,  0.54531709D+03,  0.54558168D+03,  0.54584627D+03, &
 0.54611086D+03,  0.54637545D+03,  0.54664004D+03,  0.54690462D+03,  0.54716921D+03, &
 0.54743380D+03,  0.54769839D+03,  0.54796298D+03,  0.54822757D+03,  0.54849216D+03, &
 0.54875674D+03,  0.54902133D+03,  0.54928592D+03,  0.54955051D+03,  0.54981510D+03, &
 0.55007969D+03,  0.55034428D+03,  0.55060887D+03,  0.55087345D+03,  0.55113804D+03, &
 0.55140263D+03,  0.55166722D+03,  0.55193181D+03,  0.55219640D+03,  0.55246099D+03, &
 0.55272557D+03,  0.55299016D+03,  0.55325475D+03,  0.55351934D+03,  0.55378393D+03, &
 0.55404852D+03,  0.55431311D+03,  0.55457769D+03,  0.55484228D+03,  0.55510687D+03, &
 0.55537146D+03,  0.55563605D+03,  0.55590064D+03,  0.55616523D+03,  0.55642981D+03, &
 0.55669440D+03,  0.55695899D+03,  0.55722358D+03,  0.55748817D+03,  0.55775276D+03, &
 0.55801735D+03,  0.55828193D+03,  0.55854652D+03,  0.55881111D+03,  0.55907570D+03, &
 0.55934029D+03,  0.55960488D+03,  0.55986947D+03,  0.56013405D+03,  0.56039864D+03, &
 0.56066323D+03,  0.56092782D+03,  0.56119241D+03,  0.56145700D+03,  0.56172159D+03, &
 0.56198617D+03,  0.56225076D+03,  0.56251535D+03,  0.56277994D+03,  0.56304453D+03, &
 0.56330912D+03,  0.56357371D+03,  0.56383829D+03,  0.56410288D+03,  0.56436747D+03, &
 0.56463206D+03,  0.56489665D+03,  0.56516124D+03,  0.56542583D+03,  0.56569041D+03, &
 0.56595500D+03,  0.56621959D+03,  0.56648418D+03,  0.56674877D+03,  0.56701336D+03, &
 0.56727795D+03,  0.56754254D+03,  0.56780712D+03,  0.56807171D+03,  0.56833630D+03, &
 0.56860089D+03,  0.56886548D+03,  0.56913007D+03,  0.56939466D+03,  0.56965924D+03, &
 0.56992383D+03,  0.57018842D+03,  0.57045301D+03,  0.57071760D+03,  0.57098219D+03, &
 0.57124678D+03,  0.57151136D+03,  0.57177595D+03,  0.57204054D+03,  0.57230513D+03, &
 0.57256972D+03,  0.57283431D+03,  0.57309890D+03,  0.57336348D+03,  0.57362807D+03, &
 0.57389266D+03,  0.57415725D+03,  0.57442184D+03,  0.57468643D+03,  0.57495102D+03, &
 0.57521560D+03,  0.57548019D+03,  0.57574478D+03,  0.57600937D+03,  0.57627396D+03, &
 0.57653855D+03,  0.57680314D+03,  0.57706772D+03,  0.57733231D+03,  0.57759690D+03, &
 0.57786149D+03,  0.57812608D+03,  0.57839067D+03,  0.57865526D+03,  0.57891984D+03, &
 0.57918443D+03,  0.57944902D+03,  0.57971361D+03,  0.57997820D+03,  0.58024279D+03, &
 0.58050738D+03,  0.58077196D+03,  0.58103655D+03,  0.58130114D+03,  0.58156573D+03, &
 0.58183032D+03,  0.58209491D+03,  0.58235950D+03,  0.58262408D+03,  0.58288867D+03, &
 0.58315326D+03,  0.58341785D+03,  0.58368244D+03,  0.58394703D+03,  0.58421162D+03, &
 0.58447621D+03,  0.58474079D+03,  0.58500538D+03,  0.58526997D+03,  0.58553456D+03, &
 0.58579915D+03,  0.58606374D+03,  0.58632833D+03,  0.58659291D+03,  0.58685750D+03, &
 0.58712209D+03,  0.58738668D+03,  0.58765127D+03,  0.58791586D+03,  0.58818045D+03, &
 0.58844503D+03,  0.58870962D+03,  0.58897421D+03,  0.58923880D+03,  0.58950339D+03, &
 0.58976798D+03,  0.59003257D+03,  0.59029715D+03,  0.59056174D+03,  0.59082633D+03, &
 0.59109092D+03,  0.59135551D+03,  0.59162010D+03,  0.59188469D+03,  0.59214927D+03, &
 0.59241386D+03,  0.59267845D+03,  0.59294304D+03,  0.59320763D+03,  0.59347222D+03, &
 0.59373681D+03,  0.59400139D+03,  0.59426598D+03,  0.59453057D+03,  0.59479516D+03, &
 0.59505975D+03,  0.59532434D+03,  0.59558893D+03,  0.59585351D+03,  0.59611810D+03, &
 0.59638269D+03,  0.59664728D+03,  0.59691187D+03,  0.59717646D+03,  0.59744105D+03, &
 0.59770563D+03,  0.59797022D+03,  0.59823481D+03,  0.59849940D+03,  0.59876399D+03, &
 0.59902858D+03,  0.59929317D+03,  0.59955775D+03,  0.59982234D+03,  0.60008693D+03, &
 0.60035152D+03,  0.60061611D+03,  0.60088070D+03,  0.60114529D+03,  0.60140988D+03, &
 0.60167446D+03,  0.60193905D+03,  0.60220364D+03,  0.60246823D+03,  0.60273282D+03, &
 0.60299741D+03,  0.60326200D+03,  0.60352658D+03,  0.60379117D+03,  0.60405576D+03, &
 0.60432035D+03,  0.60458494D+03,  0.60484953D+03,  0.60511412D+03,  0.60537870D+03, &
 0.60564329D+03,  0.60590788D+03,  0.60617247D+03,  0.60643706D+03,  0.60670165D+03, &
 0.60696624D+03,  0.60723082D+03,  0.60749541D+03,  0.60776000D+03,  0.60802459D+03, &
 0.60828918D+03,  0.60855377D+03,  0.60881836D+03,  0.60908294D+03,  0.60934753D+03, &
 0.60961212D+03,  0.60987671D+03,  0.61014130D+03,  0.61040589D+03,  0.61067048D+03, &
 0.61093506D+03,  0.61119965D+03,  0.61146424D+03,  0.61172883D+03,  0.61199342D+03, &
 0.61225801D+03,  0.61252260D+03,  0.61278718D+03,  0.61305177D+03,  0.61331636D+03, &
 0.61358095D+03,  0.61384554D+03,  0.61411013D+03,  0.61437472D+03,  0.61463930D+03, &
 0.61490389D+03,  0.61516848D+03,  0.61543307D+03,  0.61569766D+03,  0.61596225D+03, &
 0.61622684D+03,  0.61649143D+03,  0.61675601D+03,  0.61702060D+03,  0.61728519D+03, &
 0.61754978D+03,  0.61781437D+03,  0.61807896D+03,  0.61834355D+03,  0.61860813D+03, &
 0.61887272D+03,  0.61913731D+03,  0.61940190D+03,  0.61966649D+03,  0.61993108D+03, &
 0.62019567D+03,  0.62046025D+03,  0.62072484D+03,  0.62098943D+03,  0.62125402D+03, &
 0.62151861D+03,  0.62178320D+03,  0.62204779D+03,  0.62231237D+03,  0.62257696D+03, &
 0.62284155D+03,  0.62310614D+03,  0.62337073D+03,  0.62363532D+03,  0.62389991D+03, &
 0.62416449D+03,  0.62442908D+03,  0.62469367D+03,  0.62495826D+03,  0.62522285D+03, &
 0.62548744D+03,  0.62575203D+03,  0.62601661D+03,  0.62628120D+03,  0.62654579D+03, &
 0.62681038D+03,  0.62707497D+03,  0.62733956D+03,  0.62760415D+03,  0.62786873D+03, &
 0.62813332D+03,  0.62839791D+03,  0.62866250D+03,  0.62892709D+03,  0.62919168D+03, &
 0.62945627D+03,  0.62972085D+03,  0.62998544D+03,  0.63025003D+03,  0.63051462D+03, &
 0.63077921D+03,  0.63104380D+03,  0.63130839D+03,  0.63157297D+03,  0.63183756D+03, &
 0.63210215D+03,  0.63236674D+03,  0.63263133D+03,  0.63289592D+03,  0.63316051D+03, &
 0.63342510D+03,  0.63368968D+03,  0.63395427D+03,  0.63421886D+03,  0.63448345D+03, &
 0.63474804D+03,  0.63501263D+03,  0.63527722D+03,  0.63554180D+03,  0.63580639D+03, &
 0.63607098D+03,  0.63633557D+03,  0.63660016D+03,  0.63686475D+03,  0.63712934D+03, &
 0.63739392D+03,  0.63765851D+03,  0.63792310D+03,  0.63818769D+03,  0.63845228D+03, &
 0.63871687D+03,  0.63898146D+03,  0.63924604D+03,  0.63951063D+03,  0.63977522D+03, &
 0.64003981D+03,  0.64030440D+03,  0.64056899D+03,  0.64083358D+03,  0.64109816D+03, &
 0.64136275D+03,  0.64162734D+03,  0.64189193D+03,  0.64215652D+03,  0.64242111D+03, &
 0.64268570D+03,  0.64295028D+03,  0.64321487D+03,  0.64347946D+03,  0.64374405D+03, &
 0.64400864D+03,  0.64427323D+03,  0.64453782D+03,  0.64480240D+03,  0.64506699D+03, &
 0.64533158D+03,  0.64559617D+03,  0.64586076D+03,  0.64612535D+03,  0.64638994D+03, &
 0.64665452D+03,  0.64691911D+03,  0.64718370D+03,  0.64744829D+03,  0.64771288D+03, &
 0.64797747D+03,  0.64824206D+03,  0.64850664D+03,  0.64877123D+03,  0.64903582D+03, &
 0.64930041D+03,  0.64956500D+03,  0.64982959D+03,  0.65009418D+03,  0.65035877D+03, &
 0.65062335D+03,  0.65088794D+03,  0.65115253D+03,  0.65141712D+03,  0.65168171D+03, &
 0.65194630D+03,  0.65221089D+03,  0.65247547D+03,  0.65274006D+03,  0.65300465D+03, &
 0.65326924D+03,  0.65353383D+03,  0.65379842D+03,  0.65406301D+03,  0.65432759D+03, &
 0.65459218D+03,  0.65485677D+03,  0.65512136D+03,  0.65538595D+03,  0.65565054D+03, &
 0.65591513D+03,  0.65617971D+03,  0.65644430D+03,  0.65670889D+03,  0.65697348D+03, &
 0.65723807D+03,  0.65750266D+03,  0.65776725D+03,  0.65803183D+03,  0.65829642D+03, &
 0.65856101D+03,  0.65882560D+03,  0.65909019D+03,  0.65935478D+03,  0.65961937D+03, &
 0.65988395D+03,  0.66014854D+03,  0.66041313D+03,  0.66067772D+03,  0.66094231D+03, &
 0.66120690D+03,  0.66147149D+03,  0.66173607D+03,  0.66200066D+03,  0.66226525D+03, &
 0.66252984D+03,  0.66279443D+03,  0.66305902D+03,  0.66332361D+03,  0.66358819D+03, &
 0.66385278D+03,  0.66411737D+03,  0.66438196D+03,  0.66464655D+03,  0.66491114D+03, &
 0.66517573D+03,  0.66544032D+03,  0.66570490D+03,  0.66596949D+03,  0.66623408D+03, &
 0.66649867D+03,  0.66676326D+03,  0.66702785D+03,  0.66729244D+03,  0.66755702D+03, &
 0.66782161D+03,  0.66808620D+03,  0.66835079D+03,  0.66861538D+03,  0.66887997D+03, &
 0.66914456D+03,  0.66940914D+03,  0.66967373D+03,  0.66993832D+03,  0.67020291D+03, &
 0.67046750D+03,  0.67073209D+03,  0.67099668D+03,  0.67126126D+03,  0.67152585D+03, &
 0.67179044D+03,  0.67205503D+03,  0.67231962D+03,  0.67258421D+03,  0.67284880D+03, &
 0.67311338D+03,  0.67337797D+03,  0.67364256D+03,  0.67390715D+03,  0.67417174D+03, &
 0.67443633D+03,  0.67470092D+03,  0.67496550D+03,  0.67523009D+03,  0.67549468D+03, &
 0.67575927D+03,  0.67602386D+03,  0.67628845D+03,  0.67655304D+03,  0.67681762D+03, &
 0.67708221D+03,  0.67734680D+03,  0.67761139D+03,  0.67787598D+03,  0.67814057D+03, &
 0.67840516D+03,  0.67866974D+03,  0.67893433D+03,  0.67919892D+03,  0.67946351D+03, &
 0.67972810D+03,  0.67999269D+03,  0.68025728D+03,  0.68052186D+03,  0.68078645D+03, &
 0.68105104D+03,  0.68131563D+03,  0.68158022D+03,  0.68184481D+03,  0.68210940D+03, &
 0.68237399D+03,  0.68263857D+03,  0.68290316D+03,  0.68316775D+03,  0.68343234D+03, &
 0.68369693D+03,  0.68396152D+03,  0.68422611D+03,  0.68449069D+03,  0.68475528D+03, &
 0.68501987D+03,  0.68528446D+03,  0.68554905D+03,  0.68581364D+03,  0.68607823D+03, &
 0.68634281D+03,  0.68660740D+03,  0.68687199D+03,  0.68713658D+03,  0.68740117D+03, &
 0.68766576D+03,  0.68793035D+03,  0.68819493D+03,  0.68845952D+03,  0.68872411D+03, &
 0.68898870D+03,  0.68925329D+03,  0.68951788D+03,  0.68978247D+03,  0.69004705D+03, &
 0.69031164D+03,  0.69057623D+03,  0.69084082D+03,  0.69110541D+03,  0.69137000D+03, &
 0.69163459D+03,  0.69189917D+03,  0.69216376D+03,  0.69242835D+03,  0.69269294D+03, &
 0.69295753D+03,  0.69322212D+03,  0.69348671D+03,  0.69375129D+03,  0.69401588D+03, &
 0.69428047D+03,  0.69454506D+03,  0.69480965D+03,  0.69507424D+03,  0.69533883D+03, &
 0.69560341D+03,  0.69586800D+03,  0.69613259D+03,  0.69639718D+03,  0.69666177D+03, &
 0.69692636D+03,  0.69719095D+03,  0.69745553D+03,  0.69772012D+03,  0.69798471D+03, &
 0.69824930D+03,  0.69851389D+03,  0.69877848D+03,  0.69904307D+03,  0.69930766D+03, &
 0.69957224D+03,  0.69983683D+03,  0.70010142D+03,  0.70036601D+03,  0.70063060D+03, &
 0.70089519D+03,  0.70115978D+03,  0.70142436D+03,  0.70168895D+03,  0.70195354D+03, &
 0.70221813D+03,  0.70248272D+03,  0.70274731D+03,  0.70301190D+03,  0.70327648D+03, &
 0.70354107D+03,  0.70380566D+03,  0.70407025D+03,  0.70433484D+03,  0.70459943D+03, &
 0.70486402D+03,  0.70512860D+03,  0.70539319D+03,  0.70565778D+03,  0.70592237D+03, &
 0.70618696D+03,  0.70645155D+03,  0.70671614D+03,  0.70698072D+03,  0.70724531D+03, &
 0.70750990D+03,  0.70777449D+03,  0.70803908D+03,  0.70830367D+03,  0.70856826D+03, &
 0.70883284D+03,  0.70909743D+03,  0.70936202D+03,  0.70962661D+03,  0.70989120D+03, &
 0.71015579D+03,  0.71042038D+03,  0.71068496D+03,  0.71094955D+03,  0.71121414D+03, &
 0.71147873D+03,  0.71174332D+03,  0.71200791D+03,  0.71227250D+03,  0.71253708D+03, &
 0.71280167D+03,  0.71306626D+03,  0.71333085D+03,  0.71359544D+03,  0.71386003D+03, &
 0.71412462D+03,  0.71438920D+03,  0.71465379D+03,  0.71491838D+03,  0.71518297D+03, &
 0.71544756D+03,  0.71571215D+03,  0.71597674D+03,  0.71624133D+03,  0.71650591D+03, &
 0.71677050D+03,  0.71703509D+03,  0.71729968D+03,  0.71756427D+03,  0.71782886D+03, &
 0.71809345D+03,  0.71835803D+03,  0.71862262D+03,  0.71888721D+03,  0.71915180D+03, &
 0.71941639D+03,  0.71968098D+03,  0.71994557D+03,  0.72021015D+03,  0.72047474D+03, &
 0.72073933D+03,  0.72100392D+03,  0.72126851D+03,  0.72153310D+03,  0.72179769D+03, &
 0.72206227D+03,  0.72232686D+03,  0.72259145D+03,  0.72285604D+03,  0.72312063D+03, &
 0.72338522D+03,  0.72364981D+03,  0.72391439D+03,  0.72417898D+03,  0.72444357D+03, &
 0.72470816D+03,  0.72497275D+03,  0.72523734D+03,  0.72550193D+03,  0.72576651D+03, &
 0.72603110D+03,  0.72629569D+03,  0.72656028D+03,  0.72682487D+03,  0.72708946D+03, &
 0.72735405D+03,  0.72761863D+03,  0.72788322D+03,  0.72814781D+03,  0.72841240D+03, &
 0.72867699D+03,  0.72894158D+03,  0.72920617D+03,  0.72947075D+03,  0.72973534D+03, &
 0.72999993D+03,  0.73026452D+03,  0.73052911D+03,  0.73079370D+03,  0.73105829D+03, &
 0.73132288D+03,  0.73158746D+03,  0.73185205D+03,  0.73211664D+03,  0.73238123D+03, &
 0.73264582D+03,  0.73291041D+03,  0.73317500D+03,  0.73343958D+03,  0.73370417D+03, &
 0.73396876D+03,  0.73423335D+03,  0.73449794D+03,  0.73476253D+03,  0.73502712D+03, &
 0.73529170D+03,  0.73555629D+03,  0.73582088D+03,  0.73608547D+03,  0.73635006D+03, &
 0.73661465D+03,  0.73687924D+03,  0.73714382D+03,  0.73740841D+03,  0.73767300D+03, &
 0.73793759D+03,  0.73820218D+03,  0.73846677D+03,  0.73873136D+03,  0.73899594D+03, &
 0.73926053D+03,  0.73952512D+03,  0.73978971D+03,  0.74005430D+03,  0.74031889D+03, &
 0.74058348D+03,  0.74084806D+03,  0.74111265D+03,  0.74137724D+03,  0.74164183D+03, &
 0.74190642D+03,  0.74217101D+03,  0.74243560D+03,  0.74270018D+03,  0.74296477D+03, &
 0.74322936D+03,  0.74349395D+03,  0.74375854D+03,  0.74402313D+03,  0.74428772D+03, &
 0.74455230D+03,  0.74481689D+03,  0.74508148D+03,  0.74534607D+03,  0.74561066D+03, &
 0.74587525D+03,  0.74613984D+03,  0.74640442D+03,  0.74666901D+03,  0.74693360D+03, &
 0.74719819D+03,  0.74746278D+03,  0.74772737D+03,  0.74799196D+03,  0.74825655D+03, &
 0.74852113D+03,  0.74878572D+03,  0.74905031D+03,  0.74931490D+03,  0.74957949D+03, &
 0.74984408D+03,  0.75010867D+03,  0.75037325D+03,  0.75063784D+03,  0.75090243D+03, &
 0.75116702D+03,  0.75143161D+03,  0.75169620D+03,  0.75196079D+03,  0.75222537D+03, &
 0.75248996D+03,  0.75275455D+03,  0.75301914D+03,  0.75328373D+03,  0.75354832D+03, &
 0.75381291D+03,  0.75407749D+03,  0.75434208D+03,  0.75460667D+03,  0.75487126D+03, &
 0.75513585D+03,  0.75540044D+03,  0.75566503D+03,  0.75592961D+03,  0.75619420D+03, &
 0.75645879D+03,  0.75672338D+03,  0.75698797D+03,  0.75725256D+03,  0.75751715D+03, &
 0.75778173D+03,  0.75804632D+03,  0.75831091D+03,  0.75857550D+03,  0.75884009D+03, &
 0.75910468D+03,  0.75936927D+03,  0.75963385D+03,  0.75989844D+03,  0.76016303D+03, &
 0.76042762D+03,  0.76069221D+03,  0.76095680D+03,  0.76122139D+03,  0.76148597D+03, &
 0.76175056D+03,  0.76201515D+03,  0.76227974D+03,  0.76254433D+03,  0.76280892D+03, &
 0.76307351D+03,  0.76333809D+03,  0.76360268D+03,  0.76386727D+03,  0.76413186D+03, &
 0.76439645D+03,  0.76466104D+03,  0.76492563D+03,  0.76519022D+03,  0.76545480D+03, &
 0.76571939D+03,  0.76598398D+03,  0.76624857D+03,  0.76651316D+03,  0.76677775D+03, &
 0.76704234D+03,  0.76730692D+03,  0.76757151D+03,  0.76783610D+03,  0.76810069D+03, &
 0.76836528D+03,  0.76862987D+03,  0.76889446D+03,  0.76915904D+03,  0.76942363D+03, &
 0.76968822D+03,  0.76995281D+03,  0.77021740D+03,  0.77048199D+03,  0.77074658D+03, &
 0.77101116D+03,  0.77127575D+03,  0.77154034D+03,  0.77180493D+03,  0.77206952D+03, &
 0.77233411D+03,  0.77259870D+03,  0.77286328D+03,  0.77312787D+03,  0.77339246D+03, &
 0.77365705D+03,  0.77392164D+03,  0.77418623D+03,  0.77445082D+03,  0.77471540D+03, &
 0.77497999D+03,  0.77524458D+03,  0.77550917D+03,  0.77577376D+03,  0.77603835D+03, &
 0.77630294D+03,  0.77656752D+03,  0.77683211D+03,  0.77709670D+03,  0.77736129D+03, &
 0.77762588D+03,  0.77789047D+03,  0.77815506D+03,  0.77841964D+03,  0.77868423D+03, &
 0.77894882D+03,  0.77921341D+03,  0.77947800D+03,  0.77974259D+03,  0.78000718D+03, &
 0.78027176D+03,  0.78053635D+03,  0.78080094D+03,  0.78106553D+03,  0.78133012D+03, &
 0.78159471D+03,  0.78185930D+03,  0.78212389D+03,  0.78238847D+03,  0.78265306D+03, &
 0.78291765D+03,  0.78318224D+03,  0.78344683D+03,  0.78371142D+03,  0.78397601D+03, &
 0.78424059D+03,  0.78450518D+03,  0.78476977D+03,  0.78503436D+03,  0.78529895D+03, &
 0.78556354D+03,  0.78582813D+03,  0.78609271D+03,  0.78635730D+03,  0.78662189D+03, &
 0.78688648D+03,  0.78715107D+03,  0.78741566D+03,  0.78768025D+03,  0.78794483D+03, &
 0.78820942D+03,  0.78847401D+03,  0.78873860D+03,  0.78900319D+03,  0.78926778D+03, &
 0.78953237D+03,  0.78979695D+03,  0.79006154D+03,  0.79032613D+03,  0.79059072D+03, &
 0.79085531D+03,  0.79111990D+03,  0.79138449D+03,  0.79164907D+03,  0.79191366D+03, &
 0.79217825D+03,  0.79244284D+03,  0.79270743D+03,  0.79297202D+03,  0.79323661D+03, &
 0.79350119D+03,  0.79376578D+03,  0.79403037D+03,  0.79429496D+03,  0.79455955D+03, &
 0.79482414D+03,  0.79508873D+03,  0.79535331D+03,  0.79561790D+03,  0.79588249D+03, &
 0.79614708D+03,  0.79641167D+03,  0.79667626D+03,  0.79694085D+03,  0.79720544D+03, &
 0.79747002D+03,  0.79773461D+03,  0.79799920D+03,  0.79826379D+03,  0.79852838D+03, &
 0.79879297D+03,  0.79905756D+03,  0.79932214D+03,  0.79958673D+03,  0.79985132D+03, &
 0.80011591D+03,  0.80038050D+03,  0.80064509D+03,  0.80090968D+03,  0.80117426D+03, &
 0.80143885D+03,  0.80170344D+03,  0.80196803D+03,  0.80223262D+03,  0.80249721D+03, &
 0.80276180D+03,  0.80302638D+03,  0.80329097D+03,  0.80355556D+03,  0.80382015D+03, &
 0.80408474D+03,  0.80434933D+03,  0.80461392D+03,  0.80487850D+03,  0.80514309D+03, &
 0.80540768D+03,  0.80567227D+03,  0.80593686D+03,  0.80620145D+03,  0.80646604D+03, &
 0.80673062D+03,  0.80699521D+03,  0.80725980D+03,  0.80752439D+03,  0.80778898D+03, &
 0.80805357D+03,  0.80831816D+03,  0.80858274D+03,  0.80884733D+03,  0.80911192D+03, &
 0.80937651D+03,  0.80964110D+03,  0.80990569D+03,  0.81017028D+03,  0.81043486D+03, &
 0.81069945D+03,  0.81096404D+03,  0.81122863D+03,  0.81149322D+03,  0.81175781D+03, &
 0.81202240D+03,  0.81228698D+03,  0.81255157D+03,  0.81281616D+03,  0.81308075D+03, &
 0.81334534D+03,  0.81360993D+03,  0.81387452D+03,  0.81413911D+03,  0.81440369D+03, &
 0.81466828D+03,  0.81493287D+03,  0.81519746D+03,  0.81546205D+03,  0.81572664D+03, &
 0.81599123D+03,  0.81625581D+03,  0.81652040D+03,  0.81678499D+03,  0.81704958D+03, &
 0.81731417D+03,  0.81757876D+03,  0.81784335D+03,  0.81810793D+03,  0.81837252D+03, &
 0.81863711D+03,  0.81890170D+03,  0.81916629D+03,  0.81943088D+03,  0.81969547D+03, &
 0.81996005D+03,  0.82022464D+03,  0.82048923D+03,  0.82075382D+03,  0.82101841D+03, &
 0.82128300D+03,  0.82154759D+03,  0.82181217D+03,  0.82207676D+03,  0.82234135D+03, &
 0.82260594D+03,  0.82287053D+03,  0.82313512D+03,  0.82339971D+03,  0.82366429D+03, &
 0.82392888D+03,  0.82419347D+03,  0.82445806D+03,  0.82472265D+03,  0.82498724D+03, &
 0.82525183D+03,  0.82551641D+03,  0.82578100D+03,  0.82604559D+03,  0.82631018D+03, &
 0.82657477D+03,  0.82683936D+03,  0.82710395D+03,  0.82736853D+03,  0.82763312D+03, &
 0.82789771D+03,  0.82816230D+03,  0.82842689D+03,  0.82869148D+03,  0.82895607D+03, &
 0.82922065D+03,  0.82948524D+03,  0.82974983D+03,  0.83001442D+03,  0.83027901D+03, &
 0.83054360D+03,  0.83080819D+03,  0.83107278D+03,  0.83133736D+03,  0.83160195D+03, &
 0.83186654D+03,  0.83213113D+03,  0.83239572D+03,  0.83266031D+03,  0.83292490D+03, &
 0.83318948D+03,  0.83345407D+03,  0.83371866D+03,  0.83398325D+03,  0.83424784D+03, &
 0.83451243D+03,  0.83477702D+03,  0.83504160D+03,  0.83530619D+03,  0.83557078D+03, &
 0.83583537D+03,  0.83609996D+03,  0.83636455D+03,  0.83662914D+03,  0.83689372D+03, &
 0.83715831D+03,  0.83742290D+03,  0.83768749D+03,  0.83795208D+03,  0.83821667D+03, &
 0.83848126D+03,  0.83874584D+03,  0.83901043D+03,  0.83927502D+03,  0.83953961D+03, &
 0.83980420D+03,  0.84006879D+03,  0.84033338D+03,  0.84059796D+03,  0.84086255D+03, &
 0.84112714D+03,  0.84139173D+03,  0.84165632D+03,  0.84192091D+03,  0.84218550D+03, &
 0.84245008D+03,  0.84271467D+03,  0.84297926D+03,  0.84324385D+03,  0.84350844D+03, &
 0.84377303D+03,  0.84403762D+03,  0.84430220D+03,  0.84456679D+03,  0.84483138D+03, &
 0.84509597D+03,  0.84536056D+03,  0.84562515D+03,  0.84588974D+03,  0.84615433D+03, &
 0.84641891D+03,  0.84668350D+03,  0.84694809D+03,  0.84721268D+03,  0.84747727D+03, &
 0.84774186D+03,  0.84800645D+03,  0.84827103D+03,  0.84853562D+03,  0.84880021D+03, &
 0.84906480D+03,  0.84932939D+03,  0.84959398D+03,  0.84985857D+03,  0.85012315D+03, &
 0.85038774D+03,  0.85065233D+03,  0.85091692D+03,  0.85118151D+03,  0.85144610D+03, &
 0.85171069D+03,  0.85197527D+03,  0.85223986D+03,  0.85250445D+03,  0.85276904D+03, &
 0.85303363D+03,  0.85329822D+03,  0.85356281D+03,  0.85382739D+03,  0.85409198D+03, &
 0.85435657D+03,  0.85462116D+03,  0.85488575D+03,  0.85515034D+03,  0.85541493D+03, &
 0.85567951D+03,  0.85594410D+03,  0.85620869D+03,  0.85647328D+03,  0.85673787D+03, &
 0.85700246D+03,  0.85726705D+03,  0.85753163D+03,  0.85779622D+03,  0.85806081D+03, &
 0.85832540D+03,  0.85858999D+03,  0.85885458D+03,  0.85911917D+03,  0.85938375D+03, &
 0.85964834D+03,  0.85991293D+03,  0.86017752D+03,  0.86044211D+03,  0.86070670D+03, &
 0.86097129D+03,  0.86123587D+03,  0.86150046D+03,  0.86176505D+03,  0.86202964D+03, &
 0.86229423D+03,  0.86255882D+03,  0.86282341D+03,  0.86308800D+03,  0.86335258D+03, &
 0.86361717D+03,  0.86388176D+03,  0.86414635D+03,  0.86441094D+03,  0.86467553D+03, &
 0.86494012D+03,  0.86520470D+03,  0.86546929D+03,  0.86573388D+03,  0.86599847D+03, &
 0.86626306D+03,  0.86652765D+03,  0.86679224D+03,  0.86705682D+03,  0.86732141D+03, &
 0.86758600D+03,  0.86785059D+03,  0.86811518D+03,  0.86837977D+03,  0.86864436D+03, &
 0.86890894D+03,  0.86917353D+03,  0.86943812D+03,  0.86970271D+03,  0.86996730D+03, &
 0.87023189D+03,  0.87049648D+03,  0.87076106D+03,  0.87102565D+03,  0.87129024D+03, &
 0.87155483D+03,  0.87181942D+03,  0.87208401D+03,  0.87234860D+03,  0.87261318D+03, &
 0.87287777D+03,  0.87314236D+03,  0.87340695D+03,  0.87367154D+03,  0.87393613D+03, &
 0.87420072D+03,  0.87446530D+03,  0.87472989D+03,  0.87499448D+03,  0.87525907D+03, &
 0.87552366D+03,  0.87578825D+03,  0.87605284D+03,  0.87631742D+03,  0.87658201D+03, &
 0.87684660D+03,  0.87711119D+03,  0.87737578D+03,  0.87764037D+03,  0.87790496D+03, &
 0.87816954D+03,  0.87843413D+03,  0.87869872D+03,  0.87896331D+03,  0.87922790D+03, &
 0.87949249D+03,  0.87975708D+03,  0.88002167D+03,  0.88028625D+03,  0.88055084D+03, &
 0.88081543D+03,  0.88108002D+03,  0.88134461D+03,  0.88160920D+03,  0.88187379D+03, &
 0.88213837D+03,  0.88240296D+03,  0.88266755D+03,  0.88293214D+03,  0.88319673D+03, &
 0.88346132D+03,  0.88372591D+03,  0.88399049D+03,  0.88425508D+03,  0.88451967D+03, &
 0.88478426D+03,  0.88504885D+03,  0.88531344D+03,  0.88557803D+03,  0.88584261D+03, &
 0.88610720D+03,  0.88637179D+03,  0.88663638D+03,  0.88690097D+03,  0.88716556D+03, &
 0.88743015D+03,  0.88769473D+03,  0.88795932D+03,  0.88822391D+03,  0.88848850D+03, &
 0.88875309D+03,  0.88901768D+03,  0.88928227D+03,  0.88954685D+03,  0.88981144D+03, &
 0.89007603D+03,  0.89034062D+03,  0.89060521D+03,  0.89086980D+03,  0.89113439D+03, &
 0.89139897D+03,  0.89166356D+03,  0.89192815D+03,  0.89219274D+03,  0.89245733D+03, &
 0.89272192D+03,  0.89298651D+03,  0.89325109D+03,  0.89351568D+03,  0.89378027D+03, &
 0.89404486D+03,  0.89430945D+03,  0.89457404D+03,  0.89483863D+03,  0.89510321D+03, &
 0.89536780D+03,  0.89563239D+03,  0.89589698D+03,  0.89616157D+03,  0.89642616D+03, &
 0.89669075D+03,  0.89695534D+03,  0.89721992D+03,  0.89748451D+03,  0.89774910D+03, &
 0.89801369D+03,  0.89827828D+03,  0.89854287D+03,  0.89880746D+03,  0.89907204D+03, &
 0.89933663D+03,  0.89960122D+03,  0.89986581D+03,  0.90013040D+03,  0.90039499D+03, &
 0.90065958D+03,  0.90092416D+03,  0.90118875D+03,  0.90145334D+03,  0.90171793D+03, &
 0.90198252D+03,  0.90224711D+03,  0.90251170D+03,  0.90277628D+03,  0.90304087D+03, &
 0.90330546D+03,  0.90357005D+03,  0.90383464D+03,  0.90409923D+03,  0.90436382D+03, &
 0.90462840D+03,  0.90489299D+03,  0.90515758D+03,  0.90542217D+03,  0.90568676D+03, &
 0.90595135D+03,  0.90621594D+03,  0.90648052D+03,  0.90674511D+03,  0.90700970D+03, &
 0.90727429D+03,  0.90753888D+03,  0.90780347D+03,  0.90806806D+03,  0.90833264D+03, &
 0.90859723D+03,  0.90886182D+03,  0.90912641D+03,  0.90939100D+03,  0.90965559D+03, &
 0.90992018D+03,  0.91018476D+03,  0.91044935D+03,  0.91071394D+03,  0.91097853D+03, &
 0.91124312D+03,  0.91150771D+03,  0.91177230D+03,  0.91203689D+03,  0.91230147D+03, &
 0.91256606D+03,  0.91283065D+03,  0.91309524D+03,  0.91335983D+03,  0.91362442D+03, &
 0.91388901D+03,  0.91415359D+03,  0.91441818D+03,  0.91468277D+03,  0.91494736D+03, &
 0.91521195D+03,  0.91547654D+03,  0.91574113D+03,  0.91600571D+03,  0.91627030D+03, &
 0.91653489D+03,  0.91679948D+03,  0.91706407D+03,  0.91732866D+03,  0.91759325D+03, &
 0.91785783D+03,  0.91812242D+03,  0.91838701D+03,  0.91865160D+03,  0.91891619D+03, &
 0.91918078D+03,  0.91944537D+03,  0.91970995D+03,  0.91997454D+03,  0.92023913D+03, &
 0.92050372D+03,  0.92076831D+03,  0.92103290D+03,  0.92129749D+03,  0.92156207D+03, &
 0.92182666D+03,  0.92209125D+03,  0.92235584D+03,  0.92262043D+03,  0.92288502D+03, &
 0.92314961D+03,  0.92341419D+03,  0.92367878D+03,  0.92394337D+03,  0.92420796D+03, &
 0.92447255D+03,  0.92473714D+03,  0.92500173D+03,  0.92526631D+03,  0.92553090D+03, &
 0.92579549D+03,  0.92606008D+03,  0.92632467D+03,  0.92658926D+03,  0.92685385D+03, &
 0.92711843D+03,  0.92738302D+03,  0.92764761D+03,  0.92791220D+03,  0.92817679D+03, &
 0.92844138D+03,  0.92870597D+03,  0.92897056D+03,  0.92923514D+03,  0.92949973D+03, &
 0.92976432D+03,  0.93002891D+03,  0.93029350D+03,  0.93055809D+03,  0.93082268D+03, &
 0.93108726D+03,  0.93135185D+03,  0.93161644D+03,  0.93188103D+03,  0.93214562D+03, &
 0.93241021D+03,  0.93267480D+03,  0.93293938D+03,  0.93320397D+03,  0.93346856D+03, &
 0.93373315D+03,  0.93399774D+03,  0.93426233D+03,  0.93452692D+03,  0.93479150D+03, &
 0.93505609D+03,  0.93532068D+03,  0.93558527D+03,  0.93584986D+03,  0.93611445D+03, &
 0.93637904D+03,  0.93664362D+03,  0.93690821D+03,  0.93717280D+03,  0.93743739D+03, &
 0.93770198D+03,  0.93796657D+03,  0.93823116D+03,  0.93849574D+03,  0.93876033D+03, &
 0.93902492D+03,  0.93928951D+03,  0.93955410D+03,  0.93981869D+03,  0.94008328D+03, &
 0.94034786D+03,  0.94061245D+03,  0.94087704D+03,  0.94114163D+03,  0.94140622D+03, &
 0.94167081D+03,  0.94193540D+03,  0.94219998D+03,  0.94246457D+03,  0.94272916D+03, &
 0.94299375D+03,  0.94325834D+03,  0.94352293D+03,  0.94378752D+03,  0.94405210D+03, &
 0.94431669D+03,  0.94458128D+03,  0.94484587D+03,  0.94511046D+03,  0.94537505D+03, &
 0.94563964D+03,  0.94590423D+03,  0.94616881D+03,  0.94643340D+03,  0.94669799D+03, &
 0.94696258D+03,  0.94722717D+03,  0.94749176D+03,  0.94775635D+03,  0.94802093D+03, &
 0.94828552D+03,  0.94855011D+03,  0.94881470D+03,  0.94907929D+03,  0.94934388D+03, &
 0.94960847D+03,  0.94987305D+03,  0.95013764D+03,  0.95040223D+03,  0.95066682D+03, &
 0.95093141D+03,  0.95119600D+03,  0.95146059D+03,  0.95172517D+03,  0.95198976D+03, &
 0.95225435D+03,  0.95251894D+03,  0.95278353D+03,  0.95304812D+03,  0.95331271D+03, &
 0.95357729D+03,  0.95384188D+03,  0.95410647D+03,  0.95437106D+03,  0.95463565D+03, &
 0.95490024D+03,  0.95516483D+03,  0.95542941D+03,  0.95569400D+03,  0.95595859D+03, &
 0.95622318D+03,  0.95648777D+03,  0.95675236D+03,  0.95701695D+03,  0.95728153D+03, &
 0.95754612D+03,  0.95781071D+03,  0.95807530D+03,  0.95833989D+03,  0.95860448D+03, &
 0.95886907D+03,  0.95913365D+03,  0.95939824D+03,  0.95966283D+03,  0.95992742D+03, &
 0.96019201D+03,  0.96045660D+03,  0.96072119D+03,  0.96098577D+03,  0.96125036D+03, &
 0.96151495D+03,  0.96177954D+03,  0.96204413D+03,  0.96230872D+03,  0.96257331D+03, &
 0.96283790D+03,  0.96310248D+03,  0.96336707D+03,  0.96363166D+03,  0.96389625D+03, &
 0.96416084D+03,  0.96442543D+03,  0.96469002D+03,  0.96495460D+03,  0.96521919D+03, &
 0.96548378D+03,  0.96574837D+03,  0.96601296D+03,  0.96627755D+03,  0.96654214D+03, &
 0.96680672D+03,  0.96707131D+03,  0.96733590D+03,  0.96760049D+03,  0.96786508D+03, &
 0.96812967D+03,  0.96839426D+03,  0.96865884D+03,  0.96892343D+03,  0.96918802D+03, &
 0.96945261D+03,  0.96971720D+03,  0.96998179D+03,  0.97024638D+03,  0.97051096D+03, &
 0.97077555D+03,  0.97104014D+03,  0.97130473D+03,  0.97156932D+03,  0.97183391D+03, &
 0.97209850D+03,  0.97236308D+03,  0.97262767D+03,  0.97289226D+03,  0.97315685D+03, &
 0.97342144D+03,  0.97368603D+03,  0.97395062D+03,  0.97421520D+03,  0.97447979D+03, &
 0.97474438D+03,  0.97500897D+03,  0.97527356D+03,  0.97553815D+03,  0.97580274D+03, &
 0.97606732D+03,  0.97633191D+03,  0.97659650D+03,  0.97686109D+03,  0.97712568D+03, &
 0.97739027D+03,  0.97765486D+03,  0.97791945D+03,  0.97818403D+03,  0.97844862D+03, &
 0.97871321D+03,  0.97897780D+03,  0.97924239D+03,  0.97950698D+03,  0.97977157D+03, &
 0.98003615D+03,  0.98030074D+03,  0.98056533D+03,  0.98082992D+03,  0.98109451D+03, &
 0.98135910D+03,  0.98162369D+03,  0.98188827D+03,  0.98215286D+03,  0.98241745D+03, &
 0.98268204D+03,  0.98294663D+03,  0.98321122D+03,  0.98347581D+03,  0.98374039D+03, &
 0.98400498D+03,  0.98426957D+03,  0.98453416D+03,  0.98479875D+03,  0.98506334D+03, &
 0.98532793D+03,  0.98559251D+03,  0.98585710D+03,  0.98612169D+03,  0.98638628D+03, &
 0.98665087D+03,  0.98691546D+03,  0.98718005D+03,  0.98744463D+03,  0.98770922D+03, &
 0.98797381D+03,  0.98823840D+03,  0.98850299D+03,  0.98876758D+03,  0.98903217D+03, &
 0.98929675D+03,  0.98956134D+03,  0.98982593D+03,  0.99009052D+03,  0.99035511D+03, &
 0.99061970D+03,  0.99088429D+03,  0.99114887D+03,  0.99141346D+03,  0.99167805D+03, &
 0.99194264D+03,  0.99220723D+03,  0.99247182D+03,  0.99273641D+03,  0.99300099D+03, &
 0.99326558D+03,  0.99353017D+03,  0.99379476D+03,  0.99405935D+03,  0.99432394D+03, &
 0.99458853D+03,  0.99485312D+03,  0.99511770D+03,  0.99538229D+03,  0.99564688D+03, &
 0.99591147D+03,  0.99617606D+03,  0.99644065D+03,  0.99670524D+03,  0.99696982D+03, &
 0.99723441D+03,  0.99749900D+03,  0.99776359D+03,  0.99802818D+03,  0.99829277D+03, &
 0.99855736D+03,  0.99882194D+03,  0.99908653D+03,  0.99935112D+03,  0.99961571D+03, &
 0.99988030D+03,  0.10001449D+04,  0.10004095D+04,  0.10006741D+04,  0.10009387D+04, &
 0.10012032D+04,  0.10014678D+04,  0.10017324D+04,  0.10019970D+04,  0.10022616D+04, &
 0.10025262D+04,  0.10027908D+04,  0.10030554D+04,  0.10033200D+04,  0.10035845D+04, &
 0.10038491D+04,  0.10041137D+04,  0.10043783D+04,  0.10046429D+04,  0.10049075D+04, &
 0.10051721D+04,  0.10054367D+04,  0.10057012D+04,  0.10059658D+04,  0.10062304D+04, &
 0.10064950D+04,  0.10067596D+04,  0.10070242D+04,  0.10072888D+04,  0.10075534D+04, &
 0.10078180D+04,  0.10080825D+04,  0.10083471D+04,  0.10086117D+04,  0.10088763D+04, &
 0.10091409D+04,  0.10094055D+04,  0.10096701D+04,  0.10099347D+04,  0.10101993D+04, &
 0.10104638D+04,  0.10107284D+04,  0.10109930D+04,  0.10112576D+04,  0.10115222D+04, &
 0.10117868D+04,  0.10120514D+04,  0.10123160D+04,  0.10125806D+04,  0.10128451D+04, &
 0.10131097D+04,  0.10133743D+04,  0.10136389D+04,  0.10139035D+04,  0.10141681D+04, &
 0.10144327D+04,  0.10146973D+04,  0.10149618D+04,  0.10152264D+04,  0.10154910D+04, &
 0.10157556D+04,  0.10160202D+04,  0.10162848D+04,  0.10165494D+04,  0.10168140D+04, &
 0.10170786D+04,  0.10173431D+04,  0.10176077D+04,  0.10178723D+04,  0.10181369D+04, &
 0.10184015D+04,  0.10186661D+04,  0.10189307D+04,  0.10191953D+04,  0.10194599D+04, &
 0.10197244D+04,  0.10199890D+04,  0.10202536D+04,  0.10205182D+04,  0.10207828D+04, &
 0.10210474D+04,  0.10213120D+04,  0.10215766D+04,  0.10218412D+04,  0.10221057D+04, &
 0.10223703D+04,  0.10226349D+04,  0.10228995D+04,  0.10231641D+04,  0.10234287D+04, &
 0.10236933D+04,  0.10239579D+04,  0.10242224D+04,  0.10244870D+04,  0.10247516D+04, &
 0.10250162D+04,  0.10252808D+04,  0.10255454D+04,  0.10258100D+04,  0.10260746D+04, &
 0.10263392D+04,  0.10266037D+04,  0.10268683D+04,  0.10271329D+04,  0.10273975D+04, &
 0.10276621D+04,  0.10279267D+04,  0.10281913D+04,  0.10284559D+04,  0.10287205D+04, &
 0.10289850D+04,  0.10292496D+04,  0.10295142D+04,  0.10297788D+04,  0.10300434D+04, &
 0.10303080D+04,  0.10305726D+04,  0.10308372D+04,  0.10311018D+04,  0.10313663D+04, &
 0.10316309D+04,  0.10318955D+04,  0.10321601D+04,  0.10324247D+04,  0.10326893D+04, &
 0.10329539D+04,  0.10332185D+04,  0.10334830D+04,  0.10337476D+04,  0.10340122D+04, &
 0.10342768D+04,  0.10345414D+04,  0.10348060D+04,  0.10350706D+04,  0.10353352D+04, &
 0.10355998D+04,  0.10358643D+04,  0.10361289D+04,  0.10363935D+04,  0.10366581D+04, &
 0.10369227D+04,  0.10371873D+04,  0.10374519D+04,  0.10377165D+04,  0.10379811D+04, &
 0.10382456D+04,  0.10385102D+04,  0.10387748D+04,  0.10390394D+04,  0.10393040D+04, &
 0.10395686D+04,  0.10398332D+04,  0.10400978D+04,  0.10403624D+04,  0.10406269D+04, &
 0.10408915D+04,  0.10411561D+04,  0.10414207D+04,  0.10416853D+04,  0.10419499D+04, &
 0.10422145D+04,  0.10424791D+04,  0.10427437D+04,  0.10430082D+04,  0.10432728D+04, &
 0.10435374D+04,  0.10438020D+04,  0.10440666D+04,  0.10443312D+04,  0.10445958D+04, &
 0.10448604D+04,  0.10451249D+04,  0.10453895D+04,  0.10456541D+04,  0.10459187D+04, &
 0.10461833D+04,  0.10464479D+04,  0.10467125D+04,  0.10469771D+04,  0.10472417D+04, &
 0.10475062D+04,  0.10477708D+04,  0.10480354D+04,  0.10483000D+04,  0.10485646D+04, &
 0.10488292D+04,  0.10490938D+04,  0.10493584D+04,  0.10496230D+04,  0.10498875D+04, &
 0.10501521D+04,  0.10504167D+04,  0.10506813D+04,  0.10509459D+04,  0.10512105D+04, &
 0.10514751D+04,  0.10517397D+04,  0.10520043D+04,  0.10522688D+04,  0.10525334D+04, &
 0.10527980D+04,  0.10530626D+04,  0.10533272D+04,  0.10535918D+04,  0.10538564D+04, &
 0.10541210D+04,  0.10543855D+04,  0.10546501D+04,  0.10549147D+04,  0.10551793D+04, &
 0.10554439D+04,  0.10557085D+04,  0.10559731D+04,  0.10562377D+04,  0.10565023D+04, &
 0.10567668D+04,  0.10570314D+04,  0.10572960D+04,  0.10575606D+04,  0.10578252D+04, &
 0.10580898D+04,  0.10583544D+04,  0.10586190D+04,  0.10588836D+04,  0.10591481D+04, &
 0.10594127D+04,  0.10596773D+04,  0.10599419D+04,  0.10602065D+04,  0.10604711D+04, &
 0.10607357D+04,  0.10610003D+04,  0.10612649D+04,  0.10615294D+04,  0.10617940D+04, &
 0.10620586D+04,  0.10623232D+04,  0.10625878D+04,  0.10628524D+04,  0.10631170D+04, &
 0.10633816D+04,  0.10636461D+04,  0.10639107D+04,  0.10641753D+04,  0.10644399D+04, &
 0.10647045D+04,  0.10649691D+04,  0.10652337D+04,  0.10654983D+04,  0.10657629D+04 /


Data ECs/                                                                            &
  0.25601842D+04,  0.24043121D+04,  0.22515502D+04,  0.21018365D+04,  0.19551101D+04, &
 0.18113115D+04,  0.16703822D+04,  0.15322650D+04,  0.13969038D+04,  0.12642435D+04, &
 0.11342303D+04,  0.10068113D+04,  0.88172419D+03,  0.75889690D+03,  0.63821433D+03, &
 0.51967428D+03,  0.40342033D+03,  0.28968449D+03,  0.17866408D+03,  0.70410549D+02, &
-0.35185240D+02, -0.46421983D+02, -0.56192196D+02, -0.64527179D+02, -0.72307670D+02, &
-0.80197851D+02, -0.89207485D+02, -0.98923790D+02, -0.10931569D+03, -0.12028342D+03, &
-0.13195240D+03, -0.14424997D+03, -0.15718642D+03, -0.17077070D+03, -0.18498973D+03, &
-0.19984922D+03, -0.21534346D+03, -0.23147097D+03, -0.24822766D+03, -0.26561061D+03, &
-0.28361569D+03, -0.30223866D+03, -0.32147561D+03, -0.34132241D+03, -0.36177424D+03, &
-0.38282654D+03, -0.40447456D+03, -0.42671330D+03, -0.44953761D+03, -0.47294205D+03, &
-0.49692101D+03, -0.52146855D+03, -0.54657854D+03, -0.57224460D+03, -0.59846007D+03, &
-0.62521808D+03, -0.65251156D+03, -0.68033326D+03, -0.70867581D+03, -0.73753171D+03, &
-0.76689345D+03, -0.79675349D+03, -0.82710437D+03, -0.85793868D+03, -0.88924921D+03, &
-0.92102889D+03, -0.95327089D+03, -0.98596863D+03, -0.10191158D+04, -0.10527064D+04, &
-0.10867346D+04, -0.11211949D+04, -0.11560823D+04, -0.11913918D+04, -0.12271188D+04, &
-0.12632587D+04, -0.12998073D+04, -0.13367606D+04, -0.13741147D+04, -0.14118656D+04, &
-0.14500098D+04, -0.14885437D+04, -0.15274636D+04, -0.15667661D+04, -0.16064477D+04, &
-0.16465050D+04, -0.16869346D+04, -0.17277332D+04, -0.17688975D+04, -0.18104241D+04, &
-0.18523097D+04, -0.18945511D+04, -0.19371451D+04, -0.19800884D+04, -0.20233780D+04, &
-0.20670107D+04, -0.21109835D+04, -0.21552932D+04, -0.21999369D+04, -0.22449116D+04, &
-0.22902144D+04, -0.23358423D+04, -0.23817925D+04, -0.24280622D+04, -0.24746485D+04, &
-0.25215488D+04, -0.25687603D+04, -0.26162803D+04, -0.26641062D+04, -0.27122355D+04, &
-0.27606655D+04, -0.28093939D+04, -0.28584182D+04, -0.29077360D+04, -0.29573451D+04, &
-0.30072431D+04, -0.30574279D+04, -0.31078974D+04, -0.31586496D+04, -0.32096824D+04, &
-0.32609938D+04, -0.33125822D+04, -0.33644456D+04, -0.34165823D+04, -0.34689907D+04, &
-0.35216691D+04, -0.35746160D+04, -0.36278298D+04, -0.36813091D+04, -0.37350524D+04, &
-0.37890585D+04, -0.38433258D+04, -0.38978533D+04, -0.39526395D+04, -0.40076833D+04, &
-0.40629834D+04, -0.41185387D+04, -0.41743480D+04, -0.42304102D+04, -0.42867241D+04, &
-0.43432887D+04, -0.44001029D+04, -0.44571656D+04, -0.45144758D+04, -0.45720323D+04, &
-0.46298342D+04, -0.46878806D+04, -0.47461703D+04, -0.48047026D+04, -0.48634767D+04, &
-0.49224917D+04, -0.49817472D+04, -0.50412428D+04, -0.51009787D+04, -0.51609554D+04, &
-0.51662373D+04, -0.51711966D+04, -0.51759971D+04, -0.51806397D+04, -0.51851255D+04, &
-0.51894553D+04, -0.51936301D+04, -0.51976507D+04, -0.52015182D+04, -0.52052333D+04, &
-0.52087970D+04, -0.52122104D+04, -0.52154741D+04, -0.52185893D+04, -0.52215568D+04, &
-0.52243775D+04, -0.52270524D+04, -0.52295823D+04, -0.52319682D+04, -0.52342110D+04, &
-0.52363116D+04, -0.52382709D+04, -0.52400898D+04, -0.52417692D+04, -0.52433101D+04, &
-0.52447133D+04, -0.52459797D+04, -0.52471102D+04, -0.52481057D+04, -0.52489671D+04, &
-0.52496954D+04, -0.52502914D+04, -0.52507560D+04, -0.52510901D+04, -0.52512946D+04, &
-0.52513705D+04, -0.52513186D+04, -0.52511397D+04, -0.52508349D+04, -0.52504049D+04, &
-0.52498506D+04, -0.52491728D+04, -0.52483726D+04, -0.52474506D+04, -0.52464078D+04, &
-0.52452451D+04, -0.52439633D+04, -0.52425634D+04, -0.52410460D+04, -0.52394123D+04, &
-0.52376629D+04, -0.52357988D+04, -0.52338209D+04, -0.52317300D+04, -0.52295269D+04, &
-0.52272125D+04, -0.52247877D+04, -0.52222532D+04, -0.52196100D+04, -0.52168589D+04, &
-0.52140007D+04, -0.52110362D+04, -0.52079663D+04, -0.52047916D+04, -0.52015132D+04, &
-0.51981317D+04, -0.51946481D+04, -0.51910631D+04, -0.51873775D+04, -0.51835922D+04, &
-0.51797080D+04, -0.51757256D+04, -0.51716460D+04, -0.51674699D+04, -0.51631982D+04, &
-0.51609562D+04, -0.51009643D+04, -0.50412107D+04, -0.49816962D+04, -0.49224217D+04, &
-0.48633879D+04, -0.48045955D+04, -0.47460454D+04, -0.46877382D+04, -0.46296748D+04, &
-0.45718562D+04, -0.45142832D+04, -0.44569569D+04, -0.43998782D+04, -0.43430483D+04, &
-0.42864682D+04, -0.42301391D+04, -0.41740622D+04, -0.41182386D+04, -0.40626696D+04, &
-0.40073564D+04, -0.39523003D+04, -0.38975027D+04, -0.38429649D+04, -0.37886882D+04, &
-0.37346740D+04, -0.36809237D+04, -0.36274389D+04, -0.35742209D+04, -0.35212713D+04, &
-0.34685917D+04, -0.34161835D+04, -0.33640484D+04, -0.33121880D+04, -0.32606040D+04, &
-0.32092981D+04, -0.31582721D+04, -0.31075277D+04, -0.30570669D+04, -0.30068913D+04, &
-0.29570031D+04, -0.29074041D+04, -0.28580964D+04, -0.28090820D+04, -0.27603631D+04, &
-0.27119418D+04, -0.26638203D+04, -0.26160010D+04, -0.25684861D+04, -0.25212781D+04, &
-0.24743794D+04, -0.24277926D+04, -0.23815202D+04, -0.23355648D+04, -0.22899292D+04, &
-0.22446161D+04, -0.21996284D+04, -0.21549690D+04, -0.21106408D+04, -0.20666469D+04, &
-0.20229903D+04, -0.19796743D+04, -0.19367021D+04, -0.18940771D+04, -0.18518025D+04, &
-0.18098819D+04, -0.17683189D+04, -0.17271169D+04, -0.16862797D+04, -0.16458110D+04, &
-0.16057147D+04, -0.15659947D+04, -0.15266549D+04, -0.14876994D+04, -0.14491322D+04, &
-0.14109577D+04, -0.13731801D+04, -0.13358037D+04, -0.12988330D+04, -0.12622725D+04, &
-0.12261267D+04, -0.11904004D+04, -0.11550983D+04, -0.11202251D+04, -0.10857859D+04, &
-0.10517856D+04, -0.10182292D+04, -0.98512193D+03, -0.95246900D+03, -0.92027572D+03, &
-0.88854748D+03, -0.85728975D+03, -0.82650808D+03, -0.79620811D+03, -0.76639556D+03, &
-0.73707621D+03, -0.70825595D+03, -0.67994072D+03, -0.65213655D+03, -0.62484953D+03, &
-0.59808585D+03, -0.57185173D+03, -0.54615346D+03, -0.52099739D+03, -0.49638991D+03, &
-0.47233746D+03, -0.44884649D+03, -0.42592345D+03, -0.40357481D+03, -0.38180699D+03, &
-0.36062637D+03, -0.34003925D+03, -0.32005182D+03, -0.30067012D+03, -0.28190001D+03, &
-0.26374715D+03, -0.24621690D+03, -0.22931433D+03, -0.21304415D+03, -0.19741067D+03, &
-0.18241776D+03, -0.16806880D+03, -0.15436659D+03, -0.14131334D+03, -0.12891053D+03, &
-0.11715878D+03, -0.10605764D+03, -0.95605313D+02, -0.85798257D+02, -0.76630744D+02, &
-0.68094406D+02, -0.60177870D+02, -0.52867036D+02, -0.46142259D+02, -0.39997955D+02, &
-0.39259914D+02, -0.38509271D+02, -0.37776739D+02, -0.37061843D+02, -0.36364117D+02, &
-0.35683112D+02, -0.35018395D+02, -0.34369546D+02, -0.33736155D+02, -0.33117827D+02, &
-0.32514178D+02, -0.31924838D+02, -0.31349444D+02, -0.30787647D+02, -0.30239109D+02, &
-0.29703499D+02, -0.29180496D+02, -0.28669792D+02, -0.28171082D+02, -0.27684075D+02, &
-0.27208484D+02, -0.26744033D+02, -0.26290450D+02, -0.25847473D+02, -0.25414847D+02, &
-0.24992323D+02, -0.22204644D+02, -0.16793404D+02, -0.12898635D+02, -0.10039785D+02, &
-0.79075179D+01, -0.62932102D+01, -0.50558621D+01, -0.40961979D+01, -0.33437175D+01, &
-0.27495027D+01, -0.22761451D+01, -0.18963080D+01, -0.15884097D+01, -0.13380625D+01, &
-0.11323174D+01, -0.96254175D+00, -0.82154161D+00, -0.70500068D+00, -0.60716385D+00, &
-0.52515357D+00, -0.45609228D+00, -0.39710243D+00, -0.34674523D+00, -0.30358193D+00, &
-0.26617373D+00, -0.23452064D+00, -0.20718387D+00, -0.17663192D+00, -0.16587662D+00, &
-0.15588234D+00, -0.14658277D+00, -0.13792738D+00, -0.12986248D+00, -0.12234388D+00, &
-0.11533051D+00, -0.10878450D+00, -0.10266794D+00, -0.96949250D-01, -0.91600021D-01, &
-0.86591829D-01, -0.81899412D-01, -0.77503824D-01, -0.73382960D-01, -0.69514715D-01, &
-0.65883303D-01, -0.62469774D-01, -0.59261500D-01, -0.56245849D-01, -0.53407031D-01, &
-0.50735574D-01, -0.48218847D-01, -0.45844219D-01, -0.43608532D-01, -0.41495997D-01, &
-0.39503456D-01, -0.37624595D-01, -0.35846782D-01, -0.34166858D-01, -0.32575352D-01, &
-0.31073526D-01, -0.29650644D-01, -0.28303548D-01, -0.27027501D-01, -0.25818082D-01, &
-0.24671187D-01, -0.23583658D-01, -0.22551705D-01, -0.21571855D-01, -0.20641582D-01, &
-0.19758044D-01, -0.18918082D-01, -0.18119803D-01, -0.17360364D-01, -0.16637871D-01, &
-0.15950428D-01, -0.15295827D-01, -0.14672171D-01, -0.14078514D-01, -0.13512329D-01, &
-0.12972670D-01, -0.12457957D-01, -0.11967242D-01, -0.11498632D-01, -0.11051495D-01, &
-0.10624251D-01, -0.10216585D-01, -0.98269190D-02, -0.94543044D-02, -0.90984260D-02, &
-0.87577047D-02, -0.84321407D-02, -0.81204707D-02, -0.78220633D-02, -0.75366027D-02, &
-0.72631415D-02, -0.70010482D-02, -0.67500071D-02, -0.65097023D-02, -0.62788707D-02, &
-0.60578282D-02, -0.58459432D-02, -0.56422683D-02, -0.54471193D-02, -0.52598647D-02, &
-0.50798729D-02, -0.49071440D-02, -0.47410463D-02, -0.45818957D-02, -0.44284290D-02, &
-0.42812779D-02, -0.41398107D-02, -0.40037116D-02, -0.38726650D-02, -0.37466707D-02, &
-0.36254131D-02, -0.35088921D-02, -0.33964762D-02, -0.32884812D-02, -0.31842754D-02, &
-0.30841116D-02, -0.29874844D-02, -0.28943624D-02, -0.28046507D-02, -0.27181283D-02, &
-0.26347321D-02, -0.25543042D-02, -0.24767183D-02, -0.24018796D-02, -0.23296619D-02, &
-0.22599703D-02, -0.21926787D-02, -0.21277553D-02, -0.20650424D-02, -0.20044767D-02, &
-0.19459636D-02, -0.18894399D-02, -0.18348424D-02, -0.17820764D-02, -0.17310471D-02, &
-0.16817231D-02, -0.16340411D-02, -0.15879379D-02, -0.15433505D-02, -0.15001840D-02, &
-0.14584386D-02, -0.14180509D-02, -0.13789580D-02, -0.13411281D-02, -0.13044982D-02, &
-0.12690051D-02, -0.12346488D-02, -0.12013977D-02, -0.11691571D-02, -0.11379269D-02, &
-0.11076757D-02, -0.10783717D-02, -0.10499520D-02, -0.10223848D-02, -0.99567025D-03, &
-0.96977670D-03, -0.94467258D-03, -0.92032633D-03, -0.89667478D-03, -0.87374951D-03, &
-0.85151895D-03, -0.82995151D-03, -0.80898405D-03, -0.78864814D-03, -0.76891220D-03, &
-0.74974466D-03, -0.73114551D-03, -0.71308318D-03, -0.69552608D-03, -0.67847423D-03, &
-0.66189604D-03, -0.64579152D-03, -0.63016065D-03, -0.61497187D-03, -0.60019360D-03, &
-0.58582584D-03, -0.57186858D-03, -0.55829026D-03, -0.54509086D-03, -0.53223882D-03, &
-0.51976570D-03, -0.50760836D-03, -0.49579838D-03, -0.48430417D-03, -0.47312573D-03, &
-0.46223149D-03, -0.45165303D-03, -0.44132718D-03, -0.43128554D-03, -0.42152809D-03, &
-0.41199168D-03, -0.40273947D-03, -0.39370831D-03, -0.38492976D-03, -0.37637226D-03, &
-0.36803580D-03, -0.35992038D-03, -0.35199443D-03, -0.34428952D-03, -0.33677407D-03, &
-0.32944809D-03, -0.32228000D-03, -0.31533295D-03, -0.30854694D-03, -0.30193145D-03, &
-0.29548333D-03, -0.28918993D-03, -0.28305758D-03, -0.27707364D-03, -0.27123812D-03, &
-0.26554470D-03, -0.25999021D-03, -0.25456836D-03, -0.24928228D-03, -0.24412252D-03, &
-0.23908591D-03, -0.23416929D-03, -0.22937267D-03, -0.22468972D-03, -0.22011730D-03, &
-0.21565540D-03, -0.21129455D-03, -0.20703790D-03, -0.20288230D-03, -0.19882459D-03, &
-0.19485846D-03, -0.19098705D-03, -0.18720407D-03, -0.18350634D-03, -0.17989388D-03, &
-0.17636667D-03, -0.17291841D-03, -0.16954909D-03, -0.16625555D-03, -0.16303781D-03, &
-0.15988953D-03, -0.15681388D-03, -0.15380770D-03, -0.15086784D-03, -0.14799428D-03, &
-0.14518073D-03, -0.14243349D-03, -0.13974624D-03, -0.13711583D-03, -0.13454543D-03, &
-0.13202870D-03, -0.12956881D-03, -0.12715945D-03, -0.12480692D-03, -0.12250177D-03, &
-0.12024713D-03, -0.11803986D-03, -0.11587996D-03, -0.11376427D-03, -0.11169594D-03, &
-0.10967183D-03, -0.10768876D-03, -0.10574990D-03, -0.10384893D-03, -0.10198902D-03, &
-0.10016700D-03, -0.98382869D-04, -0.96636633D-04, -0.94925132D-04, -0.93251524D-04, &
-0.91609494D-04, -0.90002199D-04, -0.88426482D-04, -0.86885500D-04, -0.85372937D-04, &
-0.83891953D-04, -0.82442545D-04, -0.81021558D-04, -0.79628990D-04, -0.78261684D-04, &
-0.76922798D-04, -0.75612331D-04, -0.74327127D-04, -0.73067184D-04, -0.71832504D-04, &
-0.70623086D-04, -0.69435772D-04, -0.68270562D-04, -0.67130614D-04, -0.66012770D-04, &
-0.64913873D-04, -0.63837080D-04, -0.62782392D-04, -0.61746650D-04, -0.60733012D-04, &
-0.59738320D-04, -0.58759418D-04, -0.57802619D-04, -0.56861610D-04, -0.55939547D-04, &
-0.55036430D-04, -0.54149102D-04, -0.53277563D-04, -0.52421813D-04, -0.51585009D-04, &
-0.50760836D-04, -0.49952452D-04, -0.49159857D-04, -0.48383050D-04, -0.47618875D-04, &
-0.46867330D-04, -0.46131574D-04, -0.45408449D-04, -0.44701113D-04, -0.44003251D-04, &
-0.43318019D-04, -0.42648576D-04, -0.41988606D-04, -0.41341267D-04, -0.40703402D-04, &
-0.40081325D-04, -0.39465563D-04, -0.38862433D-04, -0.38271934D-04, -0.37690908D-04, &
-0.37119355D-04, -0.36557275D-04, -0.36007827D-04, -0.35467852D-04, -0.34934192D-04, &
-0.34413163D-04, -0.33898450D-04, -0.33393210D-04, -0.32897443D-04, -0.32411149D-04, &
-0.31931171D-04, -0.31460666D-04, -0.30998056D-04, -0.30543340D-04, -0.30096518D-04, &
-0.29656959D-04, -0.29225295D-04, -0.28800578D-04, -0.28383439D-04, -0.27972931D-04, &
-0.27569686D-04, -0.27173073D-04, -0.26783091D-04, -0.26399424D-04, -0.26022389D-04, &
-0.25651669D-04, -0.25286949D-04, -0.24928228D-04, -0.24575508D-04, -0.24228787D-04, &
-0.23887434D-04, -0.23552081D-04, -0.23221780D-04, -0.22897163D-04, -0.22577599D-04, &
-0.22263403D-04, -0.21954259D-04, -0.21650168D-04, -0.21350813D-04, -0.21056511D-04, &
-0.20766629D-04, -0.20481484D-04, -0.20201076D-04, -0.19925089D-04, -0.19653206D-04, &
-0.19386061D-04, -0.19123020D-04, -0.18864085D-04, -0.18608938D-04, -0.18358213D-04, &
-0.18111277D-04, -0.17868130D-04, -0.17629088D-04, -0.17393520D-04, -0.17161426D-04, &
-0.16933436D-04, -0.16708604D-04, -0.16487246D-04, -0.16269361D-04, -0.16054950D-04, &
-0.15843696D-04, -0.15635916D-04, -0.15430978D-04, -0.15229198D-04, -0.15030576D-04, &
-0.14835111D-04, -0.14642488D-04, -0.14452707D-04, -0.14265769D-04, -0.14081988D-04, &
-0.13900733D-04, -0.13722004D-04, -0.13546117D-04, -0.13372757D-04, -0.13202238D-04, &
-0.13034246D-04, -0.12868464D-04, -0.12705208D-04, -0.12544479D-04, -0.12385960D-04, &
-0.12229967D-04, -0.12076185D-04, -0.11924613D-04, -0.11775251D-04, -0.11628100D-04, &
-0.11483159D-04, -0.11340113D-04, -0.11199277D-04, -0.11060652D-04, -0.10923606D-04, &
-0.10788770D-04, -0.10656144D-04, -0.10525098D-04, -0.10395946D-04, -0.10268688D-04, &
-0.10143326D-04, -0.10019542D-04, -0.98976526D-05, -0.97773423D-05, -0.96589267D-05, &
-0.95420899D-05, -0.94268320D-05, -0.93134688D-05, -0.92013686D-05, -0.90911631D-05, &
-0.89822208D-05, -0.88748573D-05, -0.87690726D-05, -0.86645511D-05, -0.85619242D-05, &
-0.84602446D-05, -0.83601440D-05, -0.82616221D-05, -0.81640477D-05, -0.80680520D-05, &
-0.79733195D-05, -0.78798501D-05, -0.77879596D-05, -0.76970164D-05, -0.76073363D-05, &
-0.75189193D-05, -0.74317654D-05, -0.73455588D-05, -0.72606153D-05, -0.71769349D-05, &
-0.70942019D-05, -0.70127319D-05, -0.69322093D-05, -0.68529497D-05, -0.67746375D-05, &
-0.66975884D-05, -0.66211709D-05, -0.65460164D-05, -0.64718093D-05, -0.63985495D-05, &
-0.63262370D-05, -0.62548718D-05, -0.61844540D-05, -0.61149835D-05, -0.60464603D-05, &
-0.59788844D-05, -0.59119401D-05, -0.58462589D-05, -0.57812093D-05, -0.57167912D-05, &
-0.56533204D-05, -0.55907969D-05, -0.55292208D-05, -0.54682762D-05, -0.54079632D-05, &
-0.53485975D-05, -0.52898633D-05, -0.52320765D-05, -0.51749212D-05, -0.51183975D-05, &
-0.50628211D-05, -0.50075604D-05, -0.49532471D-05, -0.48998812D-05, -0.48468310D-05, &
-0.47944123D-05, -0.47429410D-05, -0.46917854D-05, -0.46415772D-05, -0.45916847D-05, &
-0.45427396D-05, -0.44941102D-05, -0.44464282D-05, -0.43990620D-05, -0.43523273D-05, &
-0.43062241D-05, -0.42607525D-05, -0.42155967D-05, -0.41713882D-05, -0.41274954D-05, &
-0.40839185D-05, -0.40412888D-05, -0.39989750D-05, -0.39569769D-05, -0.39159262D-05, &
-0.38751912D-05, -0.38347720D-05, -0.37949843D-05, -0.37555124D-05, -0.37166721D-05, &
-0.36784633D-05, -0.36405703D-05, -0.36029931D-05, -0.35660474D-05, -0.35294175D-05, &
-0.34934192D-05, -0.34577366D-05, -0.34223698D-05, -0.33876345D-05, -0.33532151D-05, &
-0.33191114D-05, -0.32853234D-05, -0.32521670D-05, -0.32193264D-05, -0.31868016D-05, &
-0.31548452D-05, -0.31231729D-05, -0.30918481D-05, -0.30608705D-05, -0.30302719D-05, &
-0.30000207D-05, -0.29701484D-05, -0.29405918D-05, -0.29113826D-05, -0.28825208D-05, &
-0.28540063D-05, -0.28258076D-05, -0.27979247D-05, -0.27703575D-05, -0.27431377D-05, &
-0.27162021D-05, -0.26895823D-05, -0.26632466D-05, -0.26372268D-05, -0.26115227D-05, &
-0.25861028D-05, -0.25609355D-05, -0.25360840D-05, -0.25115167D-05, -0.24872336D-05, &
-0.24632031D-05, -0.24394568D-05, -0.24159632D-05, -0.23927537D-05, -0.23697653D-05, &
-0.23470611D-05, -0.23246095D-05, -0.23024105D-05, -0.22804641D-05, -0.22587388D-05, &
-0.22372661D-05, -0.22160460D-05, -0.21950470D-05, -0.21742690D-05, -0.21537436D-05, &
-0.21334393D-05, -0.21133560D-05, -0.20934621D-05, -0.20738209D-05, -0.20544008D-05, &
-0.20351701D-05, -0.20161604D-05, -0.19973402D-05, -0.19787411D-05, -0.19603314D-05, &
-0.19421428D-05, -0.19241436D-05, -0.19063339D-05, -0.18887136D-05, -0.18712828D-05, &
-0.18540415D-05, -0.18369897D-05, -0.18200957D-05, -0.18034228D-05, -0.17869077D-05, &
-0.17705506D-05, -0.17543829D-05, -0.17384047D-05, -0.17225844D-05, -0.17069219D-05, &
-0.16914490D-05, -0.16761023D-05, -0.16609451D-05, -0.16459458D-05, -0.16311043D-05, &
-0.16164208D-05, -0.16018636D-05, -0.15874958D-05, -0.15732544D-05, -0.15591708D-05, &
-0.15452451D-05, -0.15314457D-05, -0.15178043D-05, -0.15042891D-05, -0.14909318D-05, &
-0.14777008D-05, -0.14645962D-05, -0.14516494D-05, -0.14388289D-05, -0.14261348D-05, &
-0.14135669D-05, -0.14011254D-05, -0.13888417D-05, -0.13766528D-05, -0.13645902D-05, &
-0.13526539D-05, -0.13408439D-05, -0.13291603D-05, -0.13176029D-05, -0.13061403D-05, &
-0.12948039D-05, -0.12835939D-05, -0.12724786D-05, -0.12614897D-05, -0.12505954D-05, &
-0.12398275D-05, -0.12291543D-05, -0.12185758D-05, -0.12081237D-05, -0.11977663D-05, &
-0.11875352D-05, -0.11773672D-05, -0.11673256D-05, -0.11573786D-05, -0.11475265D-05, &
-0.11378006D-05, -0.11281379D-05, -0.11185699D-05, -0.11091282D-05, -0.10997497D-05, &
-0.10904659D-05, -0.10812769D-05, -0.10721825D-05, -0.10631830D-05, -0.10542781D-05, &
-0.10454364D-05, -0.10366894D-05, -0.10280372D-05, -0.10194481D-05, -0.10109538D-05, &
-0.10025541D-05, -0.99421769D-06, -0.98597596D-06, -0.97779738D-06, -0.96971354D-06, &
-0.96169286D-06, -0.95376690D-06, -0.94590411D-06, -0.93813604D-06, -0.93043113D-06, &
-0.92278937D-06, -0.91521077D-06, -0.90772690D-06, -0.90030619D-06, -0.89298021D-06, &
-0.88568581D-06, -0.87848614D-06, -0.87134962D-06, -0.86427626D-06, -0.85726605D-06, &
-0.85031900D-06, -0.84346669D-06, -0.83664595D-06, -0.82988836D-06, -0.82322551D-06, &
-0.81659423D-06, -0.81005769D-06, -0.80355272D-06, -0.79711091D-06, -0.79073226D-06, &
-0.78441676D-06, -0.77816441D-06, -0.77197522D-06, -0.76581761D-06, -0.75975473D-06, &
-0.75372342D-06, -0.74775527D-06, -0.74181870D-06, -0.73597687D-06, -0.73016660D-06, &
-0.72438792D-06, -0.71870397D-06, -0.71305160D-06, -0.70746238D-06, -0.70190474D-06, &
-0.69641025D-06, -0.69094735D-06, -0.68557917D-06, -0.68021099D-06, -0.67490597D-06, &
-0.66966411D-06, -0.66445382D-06, -0.65930669D-06, -0.65419113D-06, -0.64913873D-06, &
-0.64411791D-06, -0.63916024D-06, -0.63423415D-06, -0.62933964D-06, -0.62450828D-06, &
-0.61970850D-06, -0.61494030D-06, -0.61023525D-06, -0.60556178D-06, -0.60095146D-06, &
-0.59637272D-06, -0.59182556D-06, -0.58730998D-06, -0.58285755D-06, -0.57840512D-06, &
-0.57404743D-06, -0.56968973D-06, -0.56536362D-06, -0.56110065D-06, -0.55686927D-06, &
-0.55266946D-06, -0.54850123D-06, -0.54439615D-06, -0.54029108D-06, -0.53624916D-06, &
-0.53223882D-06, -0.52826005D-06, -0.52431286D-06, -0.52039725D-06, -0.51651322D-06, &
-0.51269234D-06, -0.50887146D-06, -0.50508216D-06, -0.50135602D-06, -0.49762987D-06, &
-0.49396688D-06, -0.49030389D-06, -0.48670406D-06, -0.48310422D-06, -0.47953596D-06, &
-0.47603086D-06, -0.47252576D-06, -0.46905223D-06, -0.46564186D-06, -0.46223149D-06, &
-0.45885270D-06, -0.45550548D-06, -0.45218984D-06, -0.44887421D-06, -0.44562172D-06, &
-0.44240082D-06, -0.43917991D-06, -0.43599059D-06, -0.43283284D-06, -0.42970666D-06, &
-0.42661207D-06, -0.42354905D-06, -0.42048603D-06, -0.41745459D-06, -0.41445473D-06, &
-0.41148644D-06, -0.40854974D-06, -0.40561303D-06, -0.40270790D-06, -0.39983434D-06, &
-0.39699237D-06, -0.39415039D-06, -0.39137157D-06, -0.38859275D-06, -0.38581393D-06, &
-0.38309827D-06, -0.38038260D-06, -0.37766694D-06, -0.37501443D-06, -0.37236192D-06, &
-0.36974098D-06, -0.36715163D-06, -0.36456227D-06, -0.36200450D-06, -0.35944672D-06, &
-0.35695210D-06, -0.35442590D-06, -0.35196285D-06, -0.34949980D-06, -0.34706834D-06, &
-0.34463687D-06, -0.34223698D-06, -0.33986867D-06, -0.33750035D-06, -0.33516362D-06, &
-0.33285846D-06, -0.33055330D-06, -0.32827972D-06, -0.32600614D-06, -0.32376414D-06, &
-0.32152214D-06, -0.31931171D-06, -0.31710129D-06, -0.31494139D-06, -0.31278148D-06, &
-0.31064053D-06, -0.30851536D-06, -0.30640599D-06, -0.30431556D-06, -0.30224091D-06, &
-0.30018522D-06, -0.29814215D-06, -0.29611804D-06, -0.29410971D-06, -0.29211717D-06, &
-0.29013726D-06, -0.28817629D-06, -0.28623112D-06, -0.28429858D-06, -0.28238182D-06, &
-0.28048086D-06, -0.27859568D-06, -0.27672313D-06, -0.27486638D-06, -0.27302541D-06, &
-0.27119707D-06, -0.26938452D-06, -0.26758460D-06, -0.26580047D-06, -0.26402898D-06, &
-0.26227011D-06, -0.26052703D-06, -0.25879658D-06, -0.25707877D-06, -0.25537674D-06, &
-0.25368419D-06, -0.25200742D-06, -0.25034329D-06, -0.24869178D-06, -0.24704975D-06, &
-0.24542351D-06, -0.24380990D-06, -0.24220892D-06, -0.24062057D-06, -0.23904170D-06, &
-0.23747545D-06, -0.23592184D-06, -0.23438086D-06, -0.23285251D-06, -0.23133363D-06, &
-0.22982738D-06, -0.22833377D-06, -0.22684962D-06, -0.22537496D-06, -0.22391292D-06, &
-0.22246351D-06, -0.22102358D-06, -0.21959627D-06, -0.21817844D-06, -0.21677009D-06, &
-0.21537436D-06, -0.21398811D-06, -0.21261133D-06, -0.21124718D-06, -0.20989251D-06, &
-0.20854730D-06, -0.20721158D-06, -0.20588848D-06, -0.20457170D-06, -0.20326755D-06, &
-0.20196971D-06, -0.20068451D-06, -0.19940877D-06, -0.19814252D-06, -0.19688257D-06, &
-0.19563526D-06, -0.19439742D-06, -0.19316590D-06, -0.19194701D-06, -0.19073443D-06, &
-0.18953133D-06, -0.18833770D-06, -0.18715039D-06, -0.18597570D-06, -0.18480734D-06, &
-0.18364844D-06, -0.18249586D-06, -0.18135276D-06, -0.18021913D-06, -0.17909181D-06, &
-0.17797397D-06, -0.17686560D-06, -0.17576354D-06, -0.17467096D-06, -0.17358469D-06, &
-0.17250790D-06, -0.17143742D-06, -0.17037326D-06, -0.16931857D-06, -0.16827336D-06, &
-0.16723446D-06, -0.16620187D-06, -0.16517560D-06, -0.16415881D-06, -0.16314833D-06, &
-0.16214732D-06, -0.16115263D-06, -0.16016425D-06, -0.15918219D-06, -0.15820645D-06, &
-0.15724018D-06, -0.15628022D-06, -0.15532658D-06, -0.15437926D-06, -0.15343825D-06, &
-0.15250671D-06, -0.15157833D-06, -0.15065943D-06, -0.14974684D-06, -0.14883740D-06, &
-0.14793744D-06, -0.14704380D-06, -0.14615332D-06, -0.14527230D-06, -0.14439761D-06, &
-0.14352607D-06, -0.14266400D-06, -0.14180509D-06, -0.14095566D-06, -0.14010938D-06, &
-0.13926942D-06, -0.13843577D-06, -0.13760844D-06, -0.13678743D-06, -0.13596957D-06, &
-0.13515803D-06, -0.13435280D-06, -0.13355389D-06, -0.13276130D-06, -0.13197186D-06, &
-0.13118874D-06, -0.13041193D-06, -0.12963828D-06, -0.12887411D-06, -0.12810993D-06, &
-0.12735523D-06, -0.12660368D-06, -0.12585845D-06, -0.12511638D-06, -0.12438063D-06, &
-0.12365119D-06, -0.12292490D-06, -0.12220494D-06, -0.12148813D-06, -0.12077763D-06, &
-0.12007030D-06, -0.11936928D-06, -0.11867457D-06, -0.11798302D-06, -0.11729464D-06, &
-0.11661256D-06, -0.11593364D-06, -0.11526104D-06, -0.11459160D-06, -0.11392847D-06, &
-0.11326850D-06, -0.11261485D-06, -0.11196435D-06, -0.11131701D-06, -0.11067599D-06, &
-0.11003813D-06, -0.10940342D-06, -0.10877503D-06, -0.10814979D-06, -0.10752771D-06, &
-0.10691195D-06, -0.10629935D-06, -0.10568990D-06, -0.10508677D-06, -0.10448680D-06, &
-0.10388999D-06, -0.10329633D-06, -0.10270899D-06, -0.10212480D-06, -0.10154378D-06, &
-0.10096907D-06, -0.10039436D-06, -0.99825961D-07, -0.99260723D-07, -0.98698644D-07, &
-0.98139722D-07, -0.97587116D-07, -0.97034509D-07, -0.96488219D-07, -0.95945085D-07, &
-0.95405110D-07, -0.94868293D-07, -0.94334633D-07, -0.93807289D-07, -0.93279944D-07, &
-0.92758915D-07, -0.92237887D-07, -0.91723173D-07, -0.91211618D-07, -0.90703220D-07, &
-0.90197980D-07, -0.89695898D-07, -0.89196973D-07, -0.88701206D-07, -0.88208597D-07, &
-0.87719146D-07, -0.87232852D-07, -0.86749717D-07, -0.86269739D-07, -0.85789760D-07, &
-0.85316098D-07, -0.84845593D-07, -0.84378246D-07, -0.83914057D-07, -0.83453025D-07, &
-0.82991994D-07, -0.82537278D-07, -0.82085719D-07, -0.81634161D-07, -0.81188918D-07, &
-0.80743675D-07, -0.80301590D-07, -0.79862663D-07, -0.79426894D-07, -0.78994282D-07, &
-0.78564828D-07, -0.78135374D-07, -0.77712235D-07, -0.77289097D-07, -0.76869116D-07, &
-0.76452293D-07, -0.76038628D-07, -0.75628120D-07, -0.75217612D-07, -0.74810263D-07, &
-0.74406071D-07, -0.74005036D-07, -0.73607160D-07, -0.73209283D-07, -0.72817722D-07, &
-0.72426161D-07, -0.72037758D-07, -0.71649355D-07, -0.71267267D-07, -0.70885179D-07, &
-0.70506249D-07, -0.70127319D-07, -0.69754704D-07, -0.69382090D-07, -0.69012633D-07, &
-0.68643176D-07, -0.68276877D-07, -0.67916894D-07, -0.67553752D-07, -0.67196927D-07, &
-0.66840101D-07, -0.66486433D-07, -0.66132765D-07, -0.65782255D-07, -0.65434902D-07, &
-0.65090707D-07, -0.64746512D-07, -0.64405475D-07, -0.64067596D-07, -0.63729717D-07, &
-0.63394995D-07, -0.63060274D-07, -0.62731868D-07, -0.62400304D-07, -0.62075056D-07, &
-0.61749807D-07, -0.61427717D-07, -0.61105626D-07, -0.60786694D-07, -0.60470918D-07, &
-0.60155143D-07, -0.59842526D-07, -0.59529909D-07, -0.59220449D-07, -0.58914148D-07, &
-0.58607846D-07, -0.58304702D-07, -0.58001558D-07, -0.57701571D-07, -0.57404743D-07, &
-0.57107914D-07, -0.56811086D-07, -0.56517415D-07, -0.56226902D-07, -0.55939547D-07, &
-0.55649034D-07, -0.55364836D-07, -0.55080639D-07, -0.54796441D-07, -0.54515401D-07, &
-0.54237519D-07, -0.53959637D-07, -0.53684913D-07, -0.53410189D-07, -0.53135465D-07, &
-0.52867056D-07, -0.52595489D-07, -0.52327080D-07, -0.52061829D-07, -0.51796578D-07, &
-0.51534485D-07, -0.51272392D-07, -0.51013456D-07, -0.50754521D-07, -0.50498743D-07, &
-0.50242965D-07, -0.49987187D-07, -0.49734567D-07, -0.49485105D-07, -0.49235643D-07, &
-0.48986181D-07, -0.48739876D-07, -0.48496729D-07, -0.48253583D-07, -0.48010436D-07, &
-0.47770447D-07, -0.47530458D-07, -0.47293626D-07, -0.47056795D-07, -0.46819964D-07, &
-0.46586290D-07, -0.46352617D-07, -0.46122101D-07, -0.45891585D-07, -0.45664227D-07, &
-0.45436869D-07, -0.45212669D-07, -0.44988469D-07, -0.44764268D-07, -0.44543226D-07, &
-0.44322183D-07, -0.44101141D-07, -0.43883256D-07, -0.43668529D-07, -0.43450644D-07, &
-0.43235917D-07, -0.43024348D-07, -0.42812779D-07, -0.42601209D-07, -0.42392798D-07, &
-0.42184386D-07, -0.41975975D-07, -0.41770721D-07, -0.41565467D-07, -0.41363371D-07, &
-0.41161275D-07, -0.40959179D-07, -0.40760241D-07, -0.40561303D-07, -0.40362364D-07, &
-0.40166584D-07, -0.39970803D-07, -0.39775023D-07, -0.39582400D-07, -0.39389777D-07, &
-0.39200312D-07, -0.39007690D-07, -0.38821382D-07, -0.38631917D-07, -0.38445610D-07, &
-0.38259303D-07, -0.38072996D-07, -0.37889846D-07, -0.37706696D-07, -0.37526705D-07, &
-0.37346713D-07, -0.37166721D-07, -0.36986729D-07, -0.36809895D-07, -0.36633061D-07, &
-0.36456227D-07, -0.36282551D-07, -0.36108875D-07, -0.35935199D-07, -0.35764680D-07, &
-0.35594162D-07, -0.35423643D-07, -0.35256282D-07, -0.35085764D-07, -0.34921561D-07, &
-0.34754200D-07, -0.34589997D-07, -0.34425794D-07, -0.34261591D-07, -0.34097388D-07, &
-0.33936343D-07, -0.33775297D-07, -0.33617410D-07, -0.33456365D-07, -0.33298477D-07, &
-0.33143747D-07, -0.32985860D-07, -0.32831130D-07, -0.32676400D-07, -0.32521670D-07, &
-0.32370098D-07, -0.32218526D-07, -0.32066954D-07, -0.31915382D-07, -0.31766968D-07, &
-0.31618554D-07, -0.31469824D-07, -0.31322673D-07, -0.31176469D-07, -0.31030897D-07, &
-0.30886272D-07, -0.30742278D-07, -0.30598916D-07, -0.30456502D-07, -0.30315035D-07, &
-0.30174199D-07, -0.30033995D-07, -0.29894738D-07, -0.29756113D-07, -0.29618435D-07, &
-0.29481388D-07, -0.29344974D-07, -0.29209190D-07, -0.29074354D-07, -0.28940150D-07, &
-0.28806893D-07, -0.28673952D-07, -0.28541958D-07, -0.28410911D-07, -0.28280180D-07, &
-0.28150397D-07, -0.28021245D-07, -0.27892724D-07, -0.27764835D-07, -0.27637578D-07, &
-0.27511268D-07, -0.27385590D-07, -0.27260227D-07, -0.27135812D-07, -0.27012028D-07, &
-0.26888875D-07, -0.26766671D-07, -0.26644781D-07, -0.26523524D-07, -0.26402898D-07, &
-0.26283219D-07, -0.26163856D-07, -0.26045125D-07, -0.25927340D-07, -0.25809872D-07, &
-0.25693035D-07, -0.25576830D-07, -0.25461572D-07, -0.25346630D-07, -0.25232320D-07, &
-0.25118641D-07, -0.25005277D-07, -0.24892861D-07, -0.24781077D-07, -0.24669609D-07, &
-0.24558772D-07, -0.24448566D-07, -0.24338992D-07, -0.24230050D-07, -0.24121423D-07, &
-0.24013744D-07, -0.23906380D-07, -0.23799648D-07, -0.23693232D-07, -0.23587763D-07, &
-0.23482610D-07, -0.23378089D-07, -0.23273883D-07, -0.23170309D-07, -0.23067366D-07, &
-0.22965055D-07, -0.22863060D-07, -0.22761696D-07, -0.22660964D-07, -0.22560547D-07, &
-0.22460762D-07, -0.22361609D-07, -0.22262771D-07, -0.22164565D-07, -0.22066675D-07, &
-0.21969416D-07, -0.21872789D-07, -0.21776478D-07, -0.21680798D-07, -0.21585434D-07, &
-0.21490701D-07, -0.21396285D-07, -0.21302499D-07, -0.21209346D-07, -0.21116508D-07, &
-0.21023986D-07, -0.20932095D-07, -0.20840521D-07, -0.20749577D-07, -0.20659266D-07, &
-0.20568954D-07, -0.20479590D-07, -0.20390225D-07, -0.20301493D-07, -0.20213391D-07, &
-0.20125606D-07, -0.20038136D-07, -0.19951298D-07, -0.19864776D-07, -0.19778885D-07, &
-0.19693310D-07, -0.19608051D-07, -0.19523423D-07, -0.19439111D-07, -0.19355115D-07, &
-0.19271750D-07, -0.19188701D-07, -0.19105968D-07, -0.19023867D-07, -0.18942081D-07, &
-0.18860611D-07, -0.18779457D-07, -0.18698934D-07, -0.18618727D-07, -0.18539152D-07, &
-0.18459577D-07, -0.18380633D-07, -0.18302005D-07, -0.18224009D-07, -0.18146012D-07, &
-0.18068647D-07, -0.17991598D-07, -0.17914865D-07, -0.17838763D-07, -0.17762977D-07, &
-0.17687191D-07, -0.17612352D-07, -0.17537514D-07, -0.17462991D-07, -0.17389099D-07, &
-0.17315208D-07, -0.17241948D-07, -0.17169004D-07, -0.17096376D-07, -0.17024379D-07, &
-0.16952383D-07, -0.16881017D-07, -0.16809652D-07, -0.16738919D-07, -0.16668501D-07, &
-0.16598399D-07, -0.16528612D-07, -0.16459142D-07, -0.16389987D-07, -0.16321464D-07, &
-0.16252941D-07, -0.16184733D-07, -0.16117158D-07, -0.16049582D-07, -0.15982637D-07, &
-0.15915693D-07, -0.15849380D-07, -0.15783383D-07, -0.15717386D-07, -0.15652021D-07, &
-0.15586971D-07, -0.15521922D-07, -0.15457504D-07, -0.15393401D-07, -0.15329299D-07, &
-0.15265828D-07, -0.15202357D-07, -0.15139518D-07, -0.15076995D-07, -0.15014471D-07, &
-0.14952264D-07, -0.14890687D-07, -0.14829111D-07, -0.14767851D-07, -0.14707222D-07, &
-0.14646593D-07, -0.14586280D-07, -0.14526283D-07, -0.14466602D-07, -0.14406920D-07, &
-0.14347870D-07, -0.14289136D-07, -0.14230402D-07, -0.14172299D-07, -0.14114197D-07, &
-0.14056410D-07, -0.13998939D-07, -0.13941783D-07, -0.13884628D-07, -0.13828104D-07, &
-0.13771581D-07, -0.13715688D-07, -0.13659796D-07, -0.13604220D-07, -0.13548959D-07, &
-0.13493699D-07, -0.13439070D-07, -0.13384440D-07, -0.13330127D-07, -0.13276130D-07, &
-0.13222448D-07, -0.13168766D-07, -0.13115400D-07, -0.13062666D-07, -0.13009616D-07, &
-0.12957197D-07, -0.12905094D-07, -0.12852991D-07, -0.12801204D-07, -0.12749733D-07, &
-0.12698261D-07, -0.12647422D-07, -0.12596582D-07, -0.12546058D-07, -0.12495534D-07, &
-0.12445326D-07, -0.12395749D-07, -0.12345856D-07, -0.12296595D-07, -0.12247335D-07, &
-0.12198389D-07, -0.12149760D-07, -0.12101131D-07, -0.12052817D-07, -0.12004819D-07, &
-0.11957137D-07, -0.11909455D-07, -0.11862089D-07, -0.11815039D-07, -0.11767988D-07, &
-0.11721253D-07, -0.11674834D-07, -0.11628416D-07, -0.11582628D-07, -0.11536525D-07, &
-0.11491053D-07, -0.11445582D-07, -0.11400426D-07, -0.11355270D-07, -0.11310430D-07, &
-0.11265906D-07, -0.11221382D-07, -0.11177173D-07, -0.11133280D-07, -0.11089388D-07, &
-0.11045811D-07, -0.11002549D-07, -0.10959288D-07, -0.10916343D-07, -0.10873397D-07, &
-0.10830768D-07, -0.10788454D-07, -0.10746140D-07, -0.10704142D-07, -0.10662460D-07, &
-0.10620777D-07, -0.10579411D-07, -0.10538044D-07, -0.10496994D-07, -0.10456259D-07, &
-0.10415524D-07, -0.10374789D-07, -0.10334685D-07, -0.10294266D-07, -0.10254478D-07, &
-0.10214691D-07, -0.10174903D-07, -0.10135747D-07, -0.10096275D-07, -0.10057435D-07, &
-0.10018279D-07, -0.99797541D-08, -0.99412295D-08, -0.99027050D-08, -0.98644962D-08, &
-0.98266032D-08, -0.97887102D-08, -0.97511330D-08, -0.97135557D-08, -0.96759785D-08, &
-0.96390328D-08, -0.96017714D-08, -0.95651415D-08, -0.95285116D-08, -0.94918817D-08, &
-0.94555675D-08, -0.94192534D-08, -0.93832551D-08, -0.93472567D-08, -0.93115741D-08, &
-0.92762073D-08, -0.92408405D-08, -0.92054737D-08, -0.91704227D-08, -0.91353717D-08, &
-0.91006364D-08, -0.90662169D-08, -0.90314817D-08, -0.89973780D-08, -0.89632743D-08, &
-0.89291706D-08, -0.88953826D-08, -0.88615947D-08, -0.88281225D-08, -0.87946504D-08, &
-0.87614940D-08, -0.87283376D-08, -0.86951813D-08, -0.86623407D-08, -0.86298158D-08, &
-0.85972910D-08, -0.85647662D-08, -0.85325571D-08, -0.85003481D-08, -0.84684548D-08, &
-0.84365615D-08, -0.84049840D-08, -0.83734065D-08, -0.83421448D-08, -0.83108830D-08, &
-0.82796213D-08, -0.82486754D-08, -0.82177294D-08, -0.81870992D-08, -0.81564691D-08, &
-0.81258389D-08, -0.80955245D-08, -0.80655258D-08, -0.80352114D-08, -0.80055286D-08, &
-0.79755300D-08, -0.79458471D-08, -0.79164800D-08, -0.78867972D-08, -0.78577459D-08, &
-0.78283788D-08, -0.77993275D-08, -0.77705920D-08, -0.77418564D-08, -0.77131209D-08, &
-0.76843854D-08, -0.76559656D-08, -0.76278617D-08, -0.75997577D-08, -0.75716537D-08, &
-0.75435497D-08, -0.75157615D-08, -0.74882891D-08, -0.74605009D-08, -0.74330285D-08, &
-0.74058718D-08, -0.73787152D-08, -0.73515585D-08, -0.73244018D-08, -0.72975610D-08, &
-0.72707201D-08, -0.72441950D-08, -0.72176699D-08, -0.71911448D-08, -0.71649355D-08, &
-0.71387261D-08, -0.71128326D-08, -0.70866233D-08, -0.70607297D-08, -0.70351519D-08, &
-0.70095741D-08, -0.69839964D-08, -0.69584186D-08, -0.69331566D-08, -0.69078946D-08, &
-0.68829484D-08, -0.68580021D-08, -0.68330559D-08, -0.68081097D-08, -0.67834792D-08, &
-0.67588488D-08, -0.67345341D-08, -0.67102194D-08, -0.66859047D-08, -0.66615901D-08, &
-0.66375912D-08, -0.66135923D-08, -0.65899091D-08, -0.65662260D-08, -0.65425429D-08, &
-0.65188597D-08, -0.64954924D-08, -0.64721250D-08, -0.64487577D-08, -0.64257061D-08, &
-0.64026545D-08, -0.63796030D-08, -0.63568672D-08, -0.63338156D-08, -0.63113955D-08, &
-0.62886597D-08, -0.62662397D-08, -0.62438197D-08, -0.62213997D-08, -0.61992954D-08, &
-0.61771912D-08, -0.61550869D-08, -0.61332984D-08, -0.61111942D-08, -0.60897215D-08, &
-0.60679330D-08, -0.60464603D-08, -0.60249876D-08, -0.60035149D-08, -0.59820422D-08, &
-0.59608853D-08, -0.59397283D-08, -0.59188872D-08, -0.58977303D-08, -0.58768891D-08, &
-0.58560480D-08, -0.58355226D-08, -0.58146814D-08, -0.57941560D-08, -0.57739464D-08, &
-0.57534211D-08, -0.57332115D-08, -0.57130019D-08, -0.56927923D-08, -0.56728984D-08, &
-0.56530046D-08, -0.56331108D-08, -0.56132170D-08, -0.55936389D-08, -0.55740609D-08, &
-0.55544828D-08, -0.55349047D-08, -0.55156425D-08, -0.54963802D-08, -0.54771179D-08, &
-0.54578556D-08, -0.54389091D-08, -0.54199626D-08, -0.54010161D-08, -0.53820696D-08, &
-0.53634389D-08, -0.53444924D-08, -0.53258617D-08, -0.53075467D-08, -0.52889160D-08, &
-0.52706010D-08, -0.52522861D-08, -0.52339711D-08, -0.52159720D-08, -0.51979728D-08, &
-0.51799736D-08, -0.51619744D-08, -0.51439753D-08, -0.51262919D-08, -0.51086085D-08, &
-0.50909251D-08, -0.50732417D-08, -0.50558740D-08, -0.50381906D-08, -0.50208230D-08, &
-0.50037711D-08, -0.49864035D-08, -0.49693517D-08, -0.49522998D-08, -0.49352480D-08, &
-0.49181961D-08, -0.49014600D-08, -0.48844082D-08, -0.48676721D-08, -0.48509360D-08, &
-0.48345157D-08, -0.48177797D-08, -0.48013594D-08, -0.47849390D-08, -0.47685187D-08, &
-0.47524142D-08, -0.47359939D-08, -0.47198894D-08, -0.47037849D-08, -0.46876803D-08, &
-0.46718916D-08, -0.46561028D-08, -0.46399983D-08, -0.46242096D-08, -0.46087366D-08, &
-0.45929478D-08, -0.45774749D-08, -0.45616861D-08, -0.45462131D-08, -0.45310559D-08, &
-0.45155829D-08, -0.45004257D-08, -0.44849528D-08, -0.44697956D-08, -0.44549541D-08, &
-0.44397969D-08, -0.44246397D-08, -0.44097983D-08, -0.43949569D-08, -0.43801155D-08, &
-0.43652740D-08, -0.43507484D-08, -0.43362227D-08, -0.43213813D-08, -0.43068557D-08, &
-0.42926458D-08, -0.42781201D-08, -0.42639102D-08, -0.42493846D-08, -0.42351747D-08, &
-0.42209648D-08, -0.42070707D-08, -0.41928609D-08, -0.41789668D-08, -0.41647569D-08, &
-0.41508628D-08, -0.41369687D-08, -0.41233904D-08, -0.41094963D-08, -0.40959179D-08, &
-0.40823396D-08, -0.40687613D-08, -0.40551830D-08, -0.40416046D-08, -0.40280263D-08, &
-0.40147637D-08, -0.40015012D-08, -0.39882386D-08, -0.39749761D-08, -0.39617135D-08, &
-0.39487668D-08, -0.39355042D-08, -0.39225574D-08, -0.39096107D-08, -0.38966639D-08, &
-0.38837171D-08, -0.38710861D-08, -0.38581393D-08, -0.38455083D-08, -0.38328773D-08, &
-0.38202463D-08, -0.38076153D-08, -0.37949843D-08, -0.37826691D-08, -0.37700381D-08, &
-0.37577229D-08, -0.37454076D-08, -0.37330924D-08, -0.37207772D-08, -0.37087777D-08, &
-0.36964625D-08, -0.36844631D-08, -0.36724636D-08, -0.36604642D-08, -0.36484647D-08, &
-0.36364653D-08, -0.36244658D-08, -0.36127821D-08, -0.36010985D-08, -0.35894148D-08, &
-0.35774153D-08, -0.35660474D-08, -0.35543638D-08, -0.35426801D-08, -0.35313122D-08, &
-0.35196285D-08, -0.35082606D-08, -0.34968927D-08, -0.34855248D-08, -0.34741569D-08, &
-0.34631048D-08, -0.34517369D-08, -0.34406847D-08, -0.34296326D-08, -0.34182647D-08, &
-0.34072126D-08, -0.33964762D-08, -0.33854241D-08, -0.33743720D-08, -0.33636356D-08, &
-0.33525835D-08, -0.33418472D-08, -0.33311108D-08, -0.33203745D-08, -0.33096381D-08, &
-0.32992175D-08, -0.32884812D-08, -0.32780606D-08, -0.32673242D-08, -0.32569037D-08, &
-0.32464831D-08, -0.32360625D-08, -0.32256419D-08, -0.32155371D-08, -0.32051166D-08, &
-0.31950118D-08, -0.31845912D-08, -0.31744864D-08, -0.31643816D-08, -0.31542768D-08, &
-0.31442351D-08, -0.31341935D-08, -0.31242150D-08, -0.31142681D-08, -0.31043528D-08, &
-0.30945006D-08, -0.30846484D-08, -0.30748594D-08, -0.30650703D-08, -0.30553445D-08, &
-0.30456502D-08, -0.30360190D-08, -0.30263879D-08, -0.30167883D-08, -0.30072519D-08, &
-0.29977471D-08, -0.29882739D-08, -0.29788006D-08, -0.29693905D-08, -0.29600436D-08, &
-0.29506966D-08, -0.29413813D-08, -0.29320975D-08, -0.29228768D-08, -0.29136562D-08, &
-0.29044987D-08, -0.28953728D-08, -0.28862469D-08, -0.28771842D-08, -0.28681530D-08, &
-0.28591534D-08, -0.28501854D-08, -0.28412490D-08, -0.28323441D-08, -0.28234709D-08, &
-0.28146292D-08, -0.28058190D-08, -0.27970405D-08, -0.27882935D-08, -0.27795781D-08, &
-0.27708943D-08, -0.27622737D-08, -0.27536530D-08, -0.27450639D-08, -0.27365064D-08, &
-0.27279805D-08, -0.27194861D-08, -0.27110234D-08, -0.27025922D-08, -0.26941926D-08, &
-0.26858245D-08, -0.26774881D-08, -0.26691832D-08, -0.26609099D-08, -0.26526366D-08, &
-0.26444264D-08, -0.26362478D-08, -0.26280693D-08, -0.26199539D-08, -0.26118384D-08, &
-0.26037862D-08, -0.25957339D-08, -0.25877132D-08, -0.25797241D-08, -0.25717666D-08, &
-0.25638406D-08, -0.25559463D-08, -0.25480835D-08, -0.25402522D-08, -0.25324210D-08, &
-0.25246214D-08, -0.25168849D-08, -0.25091484D-08, -0.25014435D-08, -0.24937702D-08, &
-0.24861284D-08, -0.24784866D-08, -0.24709080D-08, -0.24633294D-08, -0.24558140D-08, &
-0.24482985D-08, -0.24408147D-08, -0.24333308D-08, -0.24259101D-08, -0.24185210D-08, &
-0.24111318D-08, -0.24037743D-08, -0.23964483D-08, -0.23891539D-08, -0.23818595D-08, &
-0.23746282D-08, -0.23673970D-08, -0.23601973D-08, -0.23530292D-08, -0.23458927D-08, &
-0.23387562D-08, -0.23316513D-08, -0.23245779D-08, -0.23175361D-08, -0.23105259D-08, &
-0.23035157D-08, -0.22965686D-08, -0.22896216D-08, -0.22826745D-08, -0.22757907D-08, &
-0.22689068D-08, -0.22620544D-08, -0.22552337D-08, -0.22484445D-08, -0.22416554D-08, &
-0.22348978D-08, -0.22281718D-08, -0.22214773D-08, -0.22147829D-08, -0.22081201D-08, &
-0.22014888D-08, -0.21948575D-08, -0.21882894D-08, -0.21817213D-08, -0.21751847D-08, &
-0.21686482D-08, -0.21621432D-08, -0.21556698D-08, -0.21492280D-08, -0.21427862D-08, &
-0.21363760D-08, -0.21299973D-08, -0.21236502D-08, -0.21173032D-08, -0.21109877D-08, &
-0.21046722D-08, -0.20984198D-08, -0.20921675D-08, -0.20859151D-08, -0.20797259D-08, &
-0.20735367D-08, -0.20673476D-08, -0.20612215D-08, -0.20550955D-08, -0.20489694D-08, &
-0.20429066D-08, -0.20368437D-08, -0.20308124D-08, -0.20247811D-08, -0.20187814D-08, &
-0.20128132D-08, -0.20068451D-08, -0.20009085D-08, -0.19950035D-08, -0.19890985D-08, &
-0.19832251D-08, -0.19773832D-08, -0.19715414D-08, -0.19657311D-08, -0.19599209D-08, &
-0.19541422D-08, -0.19483951D-08, -0.19426796D-08, -0.19369640D-08, -0.19312485D-08, &
-0.19255646D-08, -0.19199122D-08, -0.19142914D-08, -0.19086706D-08, -0.19030814D-08, &
-0.18974922D-08, -0.18919345D-08, -0.18863769D-08, -0.18808508D-08, -0.18753563D-08, &
-0.18698618D-08, -0.18643989D-08, -0.18589676D-08, -0.18535363D-08, -0.18481049D-08, &
-0.18427368D-08, -0.18373370D-08, -0.18320004D-08, -0.18266638D-08, -0.18213272D-08, &
-0.18160538D-08, -0.18107488D-08, -0.18054753D-08, -0.18002335D-08, -0.17950232D-08, &
-0.17898129D-08, -0.17846026D-08, -0.17794239D-08, -0.17742767D-08, -0.17691296D-08, &
-0.17640141D-08, -0.17588985D-08, -0.17538145D-08, -0.17487305D-08, -0.17436781D-08, &
-0.17386257D-08, -0.17336049D-08, -0.17286157D-08, -0.17236264D-08, -0.17186688D-08, &
-0.17137111D-08, -0.17087534D-08, -0.17038273D-08, -0.16989328D-08, -0.16940383D-08, &
-0.16891754D-08, -0.16843124D-08, -0.16794811D-08, -0.16746497D-08, -0.16698499D-08, &
-0.16650502D-08, -0.16602820D-08, -0.16555453D-08, -0.16507771D-08, -0.16460721D-08, &
-0.16413355D-08, -0.16366620D-08, -0.16319885D-08, -0.16273150D-08, -0.16226732D-08, &
-0.16180313D-08, -0.16134209D-08, -0.16088106D-08, -0.16042319D-08, -0.15996532D-08, &
-0.15951060D-08, -0.15905588D-08, -0.15860117D-08, -0.15815277D-08, -0.15770121D-08, &
-0.15725281D-08, -0.15680757D-08, -0.15636232D-08, -0.15591708D-08, -0.15547499D-08, &
-0.15503607D-08, -0.15459714D-08, -0.15415821D-08, -0.15372244D-08, -0.15328667D-08, &
-0.15285406D-08, -0.15242145D-08, -0.15198884D-08, -0.15156254D-08, -0.15113309D-08, &
-0.15070679D-08, -0.15028050D-08, -0.14985736D-08, -0.14943422D-08, -0.14901424D-08, &
-0.14859426D-08, -0.14817743D-08, -0.14776061D-08, -0.14734379D-08, -0.14693012D-08, &
-0.14651961D-08, -0.14610595D-08, -0.14569860D-08, -0.14528809D-08, -0.14488074D-08, &
-0.14447655D-08, -0.14407236D-08, -0.14366817D-08, -0.14326713D-08, -0.14286610D-08, &
-0.14246506D-08, -0.14206719D-08, -0.14166931D-08, -0.14127459D-08, -0.14087987D-08, &
-0.14048831D-08, -0.14009675D-08, -0.13970519D-08, -0.13931679D-08, -0.13892838D-08, &
-0.13854314D-08, -0.13815789D-08, -0.13777265D-08, -0.13739056D-08, -0.13700847D-08, &
-0.13662954D-08, -0.13625061D-08, -0.13587168D-08, -0.13549591D-08, -0.13512014D-08, &
-0.13474436D-08, -0.13437175D-08, -0.13399913D-08, -0.13362968D-08, -0.13326022D-08, &
-0.13289392D-08, -0.13252447D-08, -0.13215817D-08, -0.13179502D-08, -0.13143188D-08, &
-0.13106874D-08, -0.13070876D-08, -0.13034878D-08, -0.12998879D-08, -0.12963197D-08, &
-0.12927514D-08, -0.12891831D-08, -0.12856465D-08, -0.12821098D-08, -0.12786047D-08, &
-0.12750996D-08, -0.12715945D-08, -0.12681209D-08, -0.12646158D-08, -0.12611739D-08, &
-0.12577004D-08, -0.12542584D-08, -0.12508481D-08, -0.12474061D-08, -0.12440273D-08, &
-0.12406169D-08, -0.12372381D-08, -0.12338594D-08, -0.12304806D-08, -0.12271333D-08, &
-0.12237861D-08, -0.12204389D-08, -0.12171233D-08, -0.12138076D-08, -0.12105236D-08, &
-0.12072079D-08, -0.12039239D-08, -0.12006714D-08, -0.11974189D-08, -0.11941664D-08, &
-0.11909140D-08, -0.11876930D-08, -0.11844721D-08, -0.11812512D-08, -0.11780619D-08, &
-0.11748726D-08, -0.11716833D-08, -0.11685255D-08, -0.11653678D-08, -0.11622100D-08, &
-0.11590838D-08, -0.11559577D-08, -0.11528315D-08, -0.11497053D-08, -0.11466107D-08, &
-0.11435161D-08, -0.11404531D-08, -0.11373901D-08, -0.11343271D-08, -0.11312640D-08, &
-0.11282326D-08, -0.11252012D-08, -0.11221697D-08, -0.11191383D-08, -0.11161384D-08, &
-0.11131701D-08, -0.11101703D-08, -0.11072020D-08, -0.11042337D-08, -0.11012654D-08, &
-0.10983287D-08, -0.10953920D-08, -0.10924553D-08, -0.10895186D-08, -0.10866135D-08, &
-0.10837083D-08, -0.10808348D-08, -0.10779296D-08, -0.10750561D-08, -0.10722141D-08, &
-0.10693406D-08, -0.10664986D-08, -0.10636566D-08, -0.10608146D-08, -0.10580042D-08, &
-0.10551938D-08, -0.10523834D-08, -0.10495730D-08, -0.10467942D-08, -0.10440154D-08, &
-0.10412366D-08, -0.10384893D-08, -0.10357421D-08, -0.10329949D-08, -0.10302476D-08, &
-0.10275320D-08, -0.10248163D-08, -0.10221006D-08, -0.10193850D-08, -0.10167009D-08, &
-0.10140168D-08, -0.10113327D-08, -0.10086802D-08, -0.10059961D-08, -0.10033436D-08, &
-0.10007227D-08, -0.99807014D-09, -0.99544921D-09, -0.99282828D-09, -0.99020734D-09, &
-0.98761799D-09, -0.98502863D-09, -0.98243928D-09, -0.97984992D-09, -0.97726057D-09, &
-0.97470279D-09, -0.97214501D-09, -0.96961881D-09, -0.96706103D-09, -0.96453483D-09, &
-0.96200863D-09, -0.95948243D-09, -0.95698781D-09, -0.95446161D-09, -0.95196699D-09, &
-0.94950394D-09, -0.94700932D-09, -0.94454627D-09, -0.94208323D-09, -0.93962018D-09, &
-0.93715714D-09, -0.93472567D-09, -0.93229420D-09, -0.92986273D-09, -0.92743127D-09, &
-0.92503138D-09, -0.92263149D-09, -0.92023160D-09, -0.91783171D-09, -0.91546339D-09, &
-0.91306350D-09, -0.91069519D-09, -0.90835845D-09, -0.90599014D-09, -0.90365341D-09, &
-0.90131667D-09, -0.89897994D-09, -0.89664320D-09, -0.89433804D-09, -0.89200131D-09, &
-0.88969615D-09, -0.88742257D-09, -0.88511741D-09, -0.88284383D-09, -0.88053867D-09, &
-0.87829667D-09, -0.87602309D-09, -0.87374951D-09, -0.87150751D-09, -0.86926551D-09, &
-0.86702350D-09, -0.86481308D-09, -0.86257108D-09, -0.86036065D-09, -0.85815022D-09, &
-0.85593980D-09, -0.85376095D-09, -0.85155053D-09, -0.84937168D-09, -0.84719283D-09, &
-0.84501398D-09, -0.84286671D-09, -0.84068787D-09, -0.83854060D-09, -0.83639333D-09, &
-0.83427763D-09, -0.83213036D-09, -0.83001467D-09, -0.82789898D-09, -0.82578328D-09, &
-0.82366759D-09, -0.82155190D-09, -0.81946778D-09, -0.81738367D-09, -0.81529955D-09, &
-0.81321544D-09, -0.81116290D-09, -0.80907879D-09, -0.80702625D-09, -0.80497371D-09, &
-0.80295275D-09, -0.80090021D-09, -0.79887925D-09, -0.79682671D-09, -0.79480575D-09, &
-0.79281637D-09, -0.79079541D-09, -0.78877445D-09, -0.78678507D-09, -0.78479569D-09, &
-0.78280630D-09, -0.78084850D-09, -0.77885911D-09, -0.77690131D-09, -0.77494350D-09, &
-0.77298570D-09, -0.77102789D-09, -0.76907009D-09, -0.76714386D-09, -0.76521763D-09, &
-0.76329141D-09, -0.76136518D-09, -0.75943895D-09, -0.75751272D-09, -0.75561807D-09, &
-0.75372342D-09, -0.75182877D-09, -0.74993412D-09, -0.74807105D-09, -0.74617640D-09, &
-0.74431333D-09, -0.74245025D-09, -0.74058718D-09, -0.73872411D-09, -0.73686104D-09, &
-0.73502954D-09, -0.73319805D-09, -0.73136655D-09, -0.72953505D-09, -0.72770356D-09, &
-0.72587206D-09, -0.72407215D-09, -0.72227223D-09, -0.72047231D-09, -0.71867239D-09, &
-0.71687248D-09, -0.71507256D-09, -0.71330422D-09, -0.71153588D-09, -0.70976754D-09, &
-0.70799920D-09, -0.70623086D-09, -0.70446252D-09, -0.70272575D-09, -0.70098899D-09, &
-0.69922065D-09, -0.69751547D-09, -0.69577870D-09, -0.69404194D-09, -0.69233676D-09, &
-0.69059999D-09, -0.68889481D-09, -0.68718962D-09, -0.68548444D-09, -0.68381083D-09, &
-0.68210565D-09, -0.68043204D-09, -0.67872685D-09, -0.67705324D-09, -0.67537964D-09, &
-0.67373761D-09, -0.67206400D-09, -0.67039039D-09, -0.66874836D-09, -0.66710633D-09, &
-0.66546430D-09, -0.66382227D-09, -0.66218024D-09, -0.66056979D-09, -0.65892776D-09, &
-0.65731731D-09, -0.65570685D-09, -0.65409640D-09, -0.65248595D-09, -0.65087549D-09, &
-0.64929662D-09, -0.64768617D-09, -0.64610729D-09, -0.64452842D-09, -0.64294954D-09, &
-0.64137067D-09, -0.63982337D-09, -0.63824449D-09, -0.63669720D-09, -0.63511832D-09, &
-0.63357102D-09, -0.63202373D-09, -0.63047643D-09, -0.62896071D-09, -0.62741341D-09, &
-0.62589769D-09, -0.62435039D-09, -0.62283467D-09, -0.62131895D-09, -0.61980323D-09, &
-0.61828751D-09, -0.61680337D-09, -0.61528765D-09, -0.61380351D-09, -0.61231936D-09, &
-0.61083522D-09, -0.60935108D-09, -0.60786694D-09, -0.60638279D-09, -0.60493023D-09, &
-0.60344608D-09, -0.60199352D-09, -0.60054095D-09, -0.59908839D-09, -0.59763582D-09, &
-0.59618326D-09, -0.59476227D-09, -0.59330971D-09, -0.59188872D-09, -0.59046773D-09, &
-0.58901517D-09, -0.58759418D-09, -0.58620477D-09, -0.58478378D-09, -0.58336279D-09, &
-0.58197338D-09, -0.58055239D-09, -0.57916298D-09, -0.57777357D-09, -0.57638416D-09, &
-0.57499475D-09, -0.57360534D-09, -0.57224751D-09, -0.57085810D-09, -0.56950027D-09, &
-0.56814244D-09, -0.56675303D-09, -0.56539519D-09, -0.56406894D-09, -0.56271111D-09, &
-0.56135327D-09, -0.56002702D-09, -0.55866919D-09, -0.55734293D-09, -0.55601668D-09, &
-0.55469042D-09, -0.55336416D-09, -0.55203791D-09, -0.55071165D-09, -0.54938540D-09, &
-0.54809072D-09, -0.54679604D-09, -0.54546979D-09, -0.54417511D-09, -0.54288043D-09, &
-0.54158576D-09, -0.54032266D-09, -0.53902798D-09, -0.53773330D-09, -0.53647020D-09, &
-0.53517552D-09, -0.53391242D-09, -0.53264932D-09, -0.53138622D-09, -0.53012312D-09, &
-0.52886002D-09, -0.52762850D-09, -0.52636540D-09, -0.52513388D-09, -0.52387078D-09, &
-0.52263925D-09, -0.52140773D-09, -0.52017621D-09, -0.51894469D-09, -0.51771316D-09, &
-0.51648164D-09, -0.51528170D-09, -0.51405017D-09, -0.51285023D-09, -0.51165028D-09, &
-0.51045034D-09, -0.50921882D-09, -0.50801887D-09, -0.50685050D-09, -0.50565056D-09, &
-0.50445061D-09, -0.50328224D-09, -0.50208230D-09, -0.50091393D-09, -0.49974556D-09, &
-0.49854562D-09, -0.49737725D-09, -0.49620888D-09, -0.49507209D-09, -0.49390373D-09, &
-0.49273536D-09, -0.49159857D-09, -0.49043020D-09, -0.48929341D-09, -0.48815662D-09, &
-0.48701983D-09, -0.48588304D-09, -0.48474625D-09, -0.48360946D-09, -0.48247267D-09, &
-0.48133588D-09, -0.48023067D-09, -0.47909388D-09, -0.47798866D-09, -0.47688345D-09, &
-0.47577824D-09, -0.47464145D-09, -0.47353624D-09, -0.47246260D-09, -0.47135739D-09, &
-0.47025218D-09, -0.46917854D-09, -0.46807333D-09, -0.46699969D-09, -0.46589448D-09, &
-0.46482085D-09, -0.46374721D-09, -0.46267358D-09, -0.46159994D-09, -0.46052631D-09, &
-0.45945267D-09, -0.45841061D-09, -0.45733698D-09, -0.45629492D-09, -0.45522129D-09, &
-0.45417923D-09, -0.45313717D-09, -0.45206353D-09, -0.45102148D-09, -0.44997942D-09, &
-0.44896894D-09, -0.44792688D-09, -0.44688482D-09, -0.44587434D-09, -0.44483229D-09, &
-0.44382181D-09, -0.44277975D-09, -0.44176927D-09, -0.44075879D-09, -0.43974831D-09, &
-0.43873783D-09, -0.43772735D-09, -0.43671687D-09, -0.43570639D-09, -0.43472749D-09, &
-0.43371701D-09, -0.43273810D-09, -0.43172762D-09, -0.43074872D-09, -0.42976982D-09, &
-0.42875934D-09, -0.42778043D-09, -0.42680153D-09, -0.42585421D-09, -0.42487530D-09, &
-0.42389640D-09, -0.42291750D-09, -0.42197017D-09, -0.42099127D-09, -0.42004395D-09, &
-0.41906504D-09, -0.41811772D-09, -0.41717039D-09, -0.41622307D-09, -0.41527574D-09, &
-0.41432842D-09, -0.41338109D-09, -0.41243377D-09, -0.41151802D-09, -0.41057070D-09, &
-0.40965495D-09, -0.40870762D-09, -0.40779188D-09, -0.40684455D-09, -0.40592880D-09, &
-0.40501306D-09, -0.40409731D-09, -0.40318156D-09, -0.40226581D-09, -0.40135006D-09, &
-0.40046589D-09, -0.39955015D-09, -0.39863440D-09, -0.39775023D-09, -0.39683448D-09, &
-0.39595031D-09, -0.39506614D-09, -0.39415039D-09, -0.39326622D-09, -0.39238205D-09, &
-0.39149788D-09, -0.39061371D-09, -0.38972954D-09, -0.38884537D-09, -0.38799278D-09, &
-0.38710861D-09, -0.38622444D-09, -0.38537185D-09, -0.38451926D-09, -0.38363509D-09, &
-0.38278249D-09, -0.38192990D-09, -0.38107731D-09, -0.38019314D-09, -0.37934055D-09, &
-0.37851953D-09, -0.37766694D-09, -0.37681434D-09, -0.37596175D-09, -0.37514074D-09, &
-0.37428814D-09, -0.37343555D-09, -0.37261454D-09, -0.37179352D-09, -0.37094093D-09, &
-0.37011991D-09, -0.36929890D-09, -0.36847788D-09, -0.36765687D-09, -0.36683585D-09, &
-0.36601484D-09, -0.36519382D-09, -0.36437281D-09, -0.36358337D-09, -0.36276236D-09, &
-0.36194134D-09, -0.36115190D-09, -0.36033089D-09, -0.35954145D-09, -0.35875201D-09, &
-0.35796258D-09, -0.35714156D-09, -0.35635212D-09, -0.35556269D-09, -0.35477325D-09, &
-0.35398381D-09, -0.35322595D-09, -0.35243651D-09, -0.35164707D-09, -0.35085764D-09, &
-0.35009978D-09, -0.34931034D-09, -0.34855248D-09, -0.34776304D-09, -0.34700518D-09, &
-0.34624732D-09, -0.34548946D-09, -0.34470002D-09, -0.34394216D-09, -0.34318430D-09, &
-0.34242644D-09, -0.34166858D-09, -0.34094230D-09, -0.34018444D-09, -0.33942658D-09, &
-0.33870030D-09, -0.33794244D-09, -0.33718458D-09, -0.33645830D-09, -0.33570044D-09, &
-0.33497415D-09, -0.33424787D-09, -0.33352159D-09, -0.33276373D-09, -0.33203745D-09, &
-0.33131116D-09, -0.33058488D-09, -0.32985860D-09, -0.32913232D-09, -0.32843761D-09, &
-0.32771133D-09, -0.32698504D-09, -0.32629034D-09, -0.32556406D-09, -0.32483777D-09, &
-0.32414307D-09, -0.32344836D-09, -0.32272208D-09, -0.32202738D-09, -0.32133267D-09, &
-0.32060639D-09, -0.31991168D-09, -0.31921698D-09, -0.31852227D-09, -0.31782757D-09, &
-0.31713286D-09, -0.31643816D-09, -0.31576240D-09, -0.31507717D-09, -0.31439194D-09, &
-0.31370986D-09, -0.31302779D-09, -0.31234887D-09, -0.31166996D-09, -0.31099420D-09, &
-0.31031844D-09, -0.30964584D-09, -0.30897639D-09, -0.30830695D-09, -0.30763751D-09, &
-0.30697122D-09, -0.30630810D-09, -0.30564497D-09, -0.30498500D-09, -0.30432503D-09, &
-0.30366822D-09, -0.30301140D-09, -0.30235775D-09, -0.30170410D-09, -0.30105360D-09, &
-0.30040310D-09, -0.29975576D-09, -0.29911158D-09, -0.29846424D-09, -0.29782322D-09, &
-0.29718220D-09, -0.29654117D-09, -0.29590331D-09, -0.29526860D-09, -0.29463074D-09, &
-0.29399919D-09, -0.29336764D-09, -0.29273608D-09, -0.29210769D-09, -0.29148246D-09, &
-0.29085722D-09, -0.29023199D-09, -0.28960991D-09, -0.28898784D-09, -0.28836892D-09, &
-0.28775315D-09, -0.28713739D-09, -0.28652163D-09, -0.28590903D-09, -0.28529643D-09, &
-0.28468698D-09, -0.28408069D-09, -0.28347125D-09, -0.28286812D-09, -0.28226183D-09, &
-0.28166185D-09, -0.28105872D-09, -0.28045875D-09, -0.27986194D-09, -0.27926512D-09, &
-0.27867147D-09, -0.27807781D-09, -0.27748415D-09, -0.27689365D-09, -0.27630631D-09, &
-0.27571897D-09, -0.27513163D-09, -0.27454744D-09, -0.27396326D-09, -0.27338223D-09, &
-0.27280121D-09, -0.27222334D-09, -0.27164547D-09, -0.27107076D-09, -0.27049605D-09, &
-0.26992134D-09, -0.26934979D-09, -0.26878139D-09, -0.26820984D-09, -0.26764460D-09, &
-0.26707621D-09, -0.26651097D-09, -0.26594889D-09, -0.26538681D-09, -0.26482789D-09, &
-0.26426581D-09, -0.26371004D-09, -0.26315428D-09, -0.26259852D-09, -0.26204275D-09, &
-0.26149330D-09, -0.26094070D-09, -0.26039125D-09, -0.25984180D-09, -0.25929551D-09, &
-0.25874922D-09, -0.25820609D-09, -0.25766295D-09, -0.25712298D-09, -0.25657984D-09, &
-0.25604303D-09, -0.25550305D-09, -0.25496939D-09, -0.25443257D-09, -0.25389891D-09, &
-0.25336841D-09, -0.25283475D-09, -0.25230741D-09, -0.25177691D-09, -0.25124956D-09, &
-0.25072537D-09, -0.25020119D-09, -0.24967700D-09, -0.24915597D-09, -0.24863494D-09, &
-0.24811392D-09, -0.24759604D-09, -0.24707817D-09, -0.24656346D-09, -0.24604875D-09, &
-0.24553403D-09, -0.24502248D-09, -0.24451408D-09, -0.24400252D-09, -0.24349413D-09, &
-0.24298889D-09, -0.24248365D-09, -0.24197841D-09, -0.24147317D-09, -0.24097108D-09, &
-0.24047216D-09, -0.23997008D-09, -0.23947431D-09, -0.23897539D-09, -0.23847962D-09, &
-0.23798385D-09, -0.23749124D-09, -0.23699863D-09, -0.23650603D-09, -0.23601657D-09, &
-0.23552712D-09, -0.23504083D-09, -0.23455454D-09, -0.23406824D-09, -0.23358511D-09, &
-0.23310197D-09, -0.23261883D-09, -0.23213886D-09, -0.23165888D-09, -0.23118206D-09, &
-0.23070524D-09, -0.23022842D-09, -0.22975160D-09, -0.22927793D-09, -0.22880743D-09, &
-0.22833377D-09, -0.22786326D-09, -0.22739592D-09, -0.22692541D-09, -0.22646122D-09, &
-0.22599387D-09, -0.22552969D-09, -0.22506550D-09, -0.22460446D-09, -0.22414027D-09, &
-0.22368240D-09, -0.22322137D-09, -0.22276350D-09, -0.22230562D-09, -0.22185091D-09, &
-0.22139619D-09, -0.22094147D-09, -0.22048992D-09, -0.22003836D-09, -0.21958680D-09, &
-0.21913840D-09, -0.21869000D-09, -0.21824160D-09, -0.21779635D-09, -0.21735111D-09, &
-0.21690587D-09, -0.21646378D-09, -0.21602170D-09, -0.21557961D-09, -0.21514069D-09, &
-0.21470176D-09, -0.21426283D-09, -0.21382706D-09, -0.21339129D-09, -0.21295552D-09, &
-0.21252291D-09, -0.21209030D-09, -0.21165769D-09, -0.21122823D-09, -0.21079878D-09, &
-0.21036933D-09, -0.20994303D-09, -0.20951673D-09, -0.20909044D-09, -0.20866730D-09, &
-0.20824416D-09, -0.20782102D-09, -0.20739788D-09, -0.20697790D-09, -0.20655792D-09, &
-0.20614110D-09, -0.20572428D-09, -0.20530745D-09, -0.20489063D-09, -0.20447696D-09, &
-0.20406330D-09, -0.20364963D-09, -0.20323913D-09, -0.20282862D-09, -0.20241811D-09, &
-0.20201076D-09, -0.20160025D-09, -0.20119606D-09, -0.20078871D-09, -0.20038452D-09, &
-0.19998033D-09, -0.19957614D-09, -0.19917510D-09, -0.19877407D-09, -0.19837303D-09, &
-0.19797516D-09, -0.19757728D-09, -0.19717940D-09, -0.19678153D-09, -0.19638681D-09, &
-0.19599209D-09, -0.19560053D-09, -0.19520581D-09, -0.19481425D-09, -0.19442269D-09, &
-0.19403428D-09, -0.19364588D-09, -0.19325748D-09, -0.19286907D-09, -0.19248383D-09, &
-0.19209858D-09, -0.19171334D-09, -0.19132809D-09, -0.19094600D-09, -0.19056392D-09, &
-0.19018499D-09, -0.18980290D-09, -0.18942397D-09, -0.18904504D-09, -0.18866927D-09, &
-0.18829349D-09, -0.18791772D-09, -0.18754195D-09, -0.18716933D-09, -0.18679356D-09, &
-0.18642411D-09, -0.18605149D-09, -0.18568203D-09, -0.18531258D-09, -0.18494312D-09, &
-0.18457366D-09, -0.18420736D-09, -0.18384107D-09, -0.18347792D-09, -0.18311163D-09, &
-0.18274848D-09, -0.18238534D-09, -0.18202220D-09, -0.18166222D-09, -0.18130223D-09, &
-0.18094225D-09, -0.18058542D-09, -0.18022544D-09, -0.17986862D-09, -0.17951495D-09, &
-0.17915812D-09, -0.17880445D-09, -0.17845079D-09, -0.17809712D-09, -0.17774661D-09, &
-0.17739294D-09, -0.17704243D-09, -0.17669508D-09, -0.17634457D-09, -0.17599721D-09, &
-0.17564986D-09, -0.17530251D-09, -0.17495831D-09, -0.17461412D-09, -0.17426992D-09, &
-0.17392573D-09, -0.17358469D-09, -0.17324050D-09, -0.17289946D-09, -0.17256158D-09, &
-0.17222054D-09, -0.17188267D-09, -0.17154479D-09, -0.17120691D-09, -0.17087218D-09, &
-0.17053746D-09, -0.17020274D-09, -0.16986802D-09, -0.16953330D-09, -0.16920173D-09, &
-0.16887017D-09, -0.16853861D-09, -0.16821020D-09, -0.16788180D-09, -0.16755339D-09, &
-0.16722498D-09, -0.16689658D-09, -0.16657133D-09, -0.16624608D-09, -0.16592083D-09, &
-0.16559558D-09, -0.16527349D-09, -0.16495140D-09, -0.16462931D-09, -0.16430722D-09, &
-0.16398513D-09, -0.16366620D-09, -0.16334727D-09, -0.16302833D-09, -0.16271256D-09, &
-0.16239363D-09, -0.16207785D-09, -0.16176523D-09, -0.16144946D-09, -0.16113368D-09, &
-0.16082107D-09, -0.16050845D-09, -0.16019899D-09, -0.15988637D-09, -0.15957691D-09, &
-0.15926745D-09, -0.15895799D-09, -0.15864853D-09, -0.15834223D-09, -0.15803277D-09, &
-0.15772963D-09, -0.15742333D-09, -0.15711702D-09, -0.15681388D-09, -0.15651074D-09, &
-0.15620759D-09, -0.15590445D-09, -0.15560446D-09, -0.15530448D-09, -0.15500449D-09, &
-0.15470450D-09, -0.15440452D-09, -0.15410769D-09, -0.15381086D-09, -0.15351403D-09, &
-0.15321720D-09, -0.15292037D-09, -0.15262670D-09, -0.15233303D-09, -0.15203936D-09, &
-0.15174569D-09, -0.15145518D-09, -0.15116467D-09, -0.15087099D-09, -0.15058364D-09, &
-0.15029313D-09, -0.15000261D-09, -0.14971526D-09, -0.14942790D-09, -0.14914055D-09, &
-0.14885635D-09, -0.14856899D-09, -0.14828480D-09, -0.14800060D-09, -0.14771640D-09, &
-0.14743536D-09, -0.14715116D-09, -0.14687012D-09, -0.14658909D-09, -0.14630805D-09, &
-0.14602701D-09, -0.14574912D-09, -0.14547124D-09, -0.14519336D-09, -0.14491548D-09, &
-0.14463760D-09, -0.14436287D-09, -0.14408499D-09, -0.14381026D-09, -0.14353554D-09, &
-0.14326397D-09, -0.14298925D-09, -0.14271768D-09, -0.14244612D-09, -0.14217455D-09, &
-0.14190298D-09, -0.14163457D-09, -0.14136301D-09, -0.14109460D-09, -0.14082619D-09, &
-0.14056094D-09, -0.14029253D-09, -0.14002728D-09, -0.13975887D-09, -0.13949362D-09, &
-0.13922837D-09, -0.13896628D-09, -0.13870102D-09, -0.13843893D-09, -0.13817684D-09, &
-0.13791474D-09, -0.13765265D-09, -0.13739372D-09, -0.13713162D-09, -0.13687269D-09, &
-0.13661375D-09, -0.13635482D-09, -0.13609904D-09, -0.13584010D-09, -0.13558433D-09, &
-0.13532855D-09, -0.13507277D-09, -0.13481699D-09, -0.13456437D-09, -0.13430859D-09, &
-0.13405597D-09, -0.13380335D-09, -0.13355073D-09, -0.13329811D-09, -0.13304865D-09, &
-0.13279919D-09, -0.13254657D-09, -0.13229711D-09, -0.13205080D-09, -0.13180134D-09, &
-0.13155504D-09, -0.13130557D-09, -0.13105927D-09, -0.13081296D-09, -0.13056666D-09, &
-0.13032351D-09, -0.13007721D-09, -0.12983406D-09, -0.12959092D-09, -0.12934777D-09, &
-0.12910462D-09, -0.12886463D-09, -0.12862149D-09, -0.12838150D-09, -0.12814151D-09, &
-0.12790152D-09, -0.12766153D-09, -0.12742154D-09, -0.12718471D-09, -0.12694788D-09, &
-0.12671105D-09, -0.12647422D-09, -0.12623738D-09, -0.12600055D-09, -0.12576688D-09, &
-0.12553321D-09, -0.12529637D-09, -0.12506270D-09, -0.12483219D-09, -0.12459851D-09, &
-0.12436800D-09, -0.12413432D-09, -0.12390381D-09, -0.12367329D-09, -0.12344278D-09, &
-0.12321226D-09, -0.12298490D-09, -0.12275754D-09, -0.12252703D-09, -0.12229967D-09, &
-0.12207231D-09, -0.12184811D-09, -0.12162075D-09, -0.12139655D-09, -0.12116919D-09, &
-0.12094499D-09, -0.12072079D-09, -0.12049659D-09, -0.12027555D-09, -0.12005135D-09, &
-0.11983031D-09, -0.11960611D-09, -0.11938507D-09, -0.11916402D-09, -0.11894614D-09, &
-0.11872510D-09, -0.11850721D-09, -0.11828617D-09, -0.11806828D-09, -0.11785040D-09, &
-0.11763251D-09, -0.11741779D-09, -0.11719990D-09, -0.11698518D-09, -0.11676729D-09, &
-0.11655256D-09, -0.11633784D-09, -0.11612311D-09, -0.11591154D-09, -0.11569681D-09, &
-0.11548524D-09, -0.11527052D-09, -0.11505895D-09, -0.11484738D-09, -0.11463581D-09, &
-0.11442740D-09, -0.11421583D-09, -0.11400742D-09, -0.11379901D-09, -0.11358744D-09, &
-0.11338218D-09, -0.11317377D-09, -0.11296536D-09, -0.11275695D-09, -0.11255169D-09, &
-0.11234644D-09, -0.11214119D-09, -0.11193593D-09, -0.11173068D-09, -0.11152543D-09, &
-0.11132017D-09, -0.11111808D-09, -0.11091598D-09, -0.11071388D-09, -0.11051179D-09, &
-0.11030969D-09, -0.11010760D-09, -0.10990550D-09, -0.10970656D-09, -0.10950447D-09, &
-0.10930553D-09, -0.10910659D-09, -0.10890765D-09, -0.10870871D-09, -0.10851293D-09, &
-0.10831399D-09, -0.10811821D-09, -0.10791927D-09, -0.10772349D-09, -0.10752771D-09, &
-0.10733193D-09, -0.10713931D-09, -0.10694353D-09, -0.10675091D-09, -0.10655513D-09, &
-0.10636250D-09, -0.10616988D-09, -0.10597726D-09, -0.10578464D-09, -0.10559201D-09, &
-0.10540255D-09, -0.10520992D-09, -0.10502046D-09, -0.10483099D-09, -0.10464153D-09, &
-0.10445206D-09, -0.10426260D-09, -0.10407313D-09, -0.10388683D-09, -0.10369736D-09, &
-0.10351106D-09, -0.10332475D-09, -0.10313844D-09, -0.10295213D-09, -0.10276583D-09, &
-0.10257952D-09, -0.10239321D-09, -0.10221006D-09, -0.10202691D-09, -0.10184061D-09, &
-0.10165746D-09, -0.10147431D-09, -0.10129431D-09, -0.10111117D-09, -0.10092802D-09, &
-0.10074802D-09, -0.10056487D-09, -0.10038488D-09, -0.10020489D-09, -0.10002490D-09, &
-0.99844907D-10, -0.99664915D-10, -0.99488081D-10, -0.99308090D-10, -0.99131256D-10, &
-0.98954422D-10, -0.98774430D-10, -0.98597596D-10, -0.98420762D-10, -0.98247085D-10, &
-0.98070251D-10, -0.97893417D-10, -0.97719741D-10, -0.97546065D-10, -0.97369231D-10, &
-0.97195555D-10, -0.97021878D-10, -0.96848202D-10, -0.96674526D-10, -0.96504007D-10, &
-0.96330331D-10, -0.96159813D-10, -0.95989294D-10, -0.95815618D-10, -0.95645099D-10, &
-0.95474581D-10, -0.95304062D-10, -0.95136701D-10, -0.94966183D-10, -0.94795664D-10, &
-0.94628304D-10, -0.94460943D-10, -0.94290424D-10, -0.94123064D-10, -0.93955703D-10, &
-0.93791500D-10, -0.93624139D-10, -0.93456778D-10, -0.93292575D-10, -0.93125214D-10, &
-0.92961011D-10, -0.92796808D-10, -0.92629448D-10, -0.92465245D-10, -0.92304199D-10, &
-0.92139996D-10, -0.91975793D-10, -0.91814748D-10, -0.91650545D-10, -0.91489500D-10, &
-0.91325297D-10, -0.91164252D-10, -0.91003206D-10, -0.90842161D-10, -0.90681116D-10, &
-0.90523228D-10, -0.90362183D-10, -0.90204295D-10, -0.90043250D-10, -0.89885363D-10, &
-0.89727475D-10, -0.89569588D-10, -0.89411700D-10, -0.89253813D-10, -0.89095925D-10, &
-0.88938038D-10, -0.88783308D-10, -0.88625420D-10, -0.88470691D-10, -0.88312803D-10, &
-0.88158073D-10, -0.88003343D-10, -0.87848614D-10, -0.87693884D-10, -0.87539154D-10, &
-0.87387582D-10, -0.87232852D-10, -0.87081280D-10, -0.86926551D-10, -0.86774979D-10, &
-0.86623407D-10, -0.86471835D-10, -0.86320263D-10, -0.86168691D-10, -0.86017119D-10, &
-0.85865546D-10, -0.85717132D-10, -0.85565560D-10, -0.85417146D-10, -0.85265574D-10, &
-0.85117160D-10, -0.84968745D-10, -0.84820331D-10, -0.84671917D-10, -0.84523503D-10, &
-0.84378246D-10, -0.84229832D-10, -0.84081418D-10, -0.83936161D-10, -0.83790905D-10, &
-0.83642490D-10, -0.83497234D-10, -0.83351977D-10, -0.83206721D-10, -0.83061464D-10, &
-0.82916208D-10, -0.82774109D-10, -0.82628852D-10, -0.82486754D-10, -0.82341497D-10, &
-0.82199398D-10, -0.82057300D-10, -0.81915201D-10, -0.81769944D-10, -0.81627846D-10, &
-0.81488905D-10, -0.81346806D-10, -0.81204707D-10, -0.81065766D-10, -0.80923667D-10, &
-0.80784726D-10, -0.80642627D-10, -0.80503686D-10, -0.80364745D-10, -0.80225804D-10, &
-0.80086863D-10, -0.79947922D-10, -0.79808981D-10, -0.79670040D-10, -0.79534257D-10, &
-0.79395316D-10, -0.79259533D-10, -0.79120592D-10, -0.78984809D-10, -0.78849025D-10, &
-0.78713242D-10, -0.78577459D-10, -0.78441676D-10, -0.78305892D-10, -0.78170109D-10, &
-0.78037483D-10, -0.77901700D-10, -0.77769075D-10, -0.77633291D-10, -0.77500666D-10, &
-0.77368040D-10, -0.77235415D-10, -0.77102789D-10, -0.76970164D-10, -0.76837538D-10, &
-0.76704913D-10, -0.76572287D-10, -0.76439662D-10, -0.76310194D-10, -0.76177569D-10, &
-0.76048101D-10, -0.75918633D-10, -0.75786008D-10, -0.75656540D-10, -0.75527072D-10, &
-0.75397604D-10, -0.75268136D-10, -0.75141826D-10, -0.75012359D-10, -0.74882891D-10, &
-0.74756581D-10, -0.74627113D-10, -0.74500803D-10, -0.74371335D-10, -0.74245025D-10, &
-0.74118715D-10, -0.73992405D-10, -0.73866095D-10, -0.73739785D-10, -0.73613475D-10, &
-0.73487165D-10, -0.73360855D-10, -0.73237703D-10, -0.73111393D-10, -0.72988241D-10, &
-0.72861931D-10, -0.72738778D-10, -0.72615626D-10, -0.72492474D-10, -0.72369322D-10, &
-0.72246169D-10, -0.72123017D-10, -0.71999865D-10, -0.71876713D-10, -0.71756718D-10, &
-0.71633566D-10, -0.71510414D-10, -0.71390419D-10, -0.71270425D-10, -0.71147272D-10, &
-0.71027278D-10, -0.70907283D-10, -0.70787289D-10, -0.70667294D-10, -0.70547300D-10, &
-0.70427305D-10, -0.70307311D-10, -0.70190474D-10, -0.70070479D-10, -0.69950485D-10, &
-0.69833648D-10, -0.69713654D-10, -0.69596817D-10, -0.69479980D-10, -0.69363143D-10, &
-0.69246307D-10, -0.69129470D-10, -0.69012633D-10, -0.68895796D-10, -0.68778960D-10, &
-0.68662123D-10, -0.68545286D-10, -0.68431607D-10, -0.68314770D-10, -0.68201091D-10, &
-0.68087412D-10, -0.67970575D-10, -0.67856896D-10, -0.67743217D-10, -0.67629538D-10, &
-0.67515859D-10, -0.67402180D-10, -0.67288501D-10, -0.67174822D-10, -0.67061143D-10, &
-0.66950622D-10, -0.66836943D-10, -0.66723264D-10, -0.66612743D-10, -0.66502222D-10, &
-0.66388543D-10, -0.66278021D-10, -0.66167500D-10, -0.66056979D-10, -0.65946458D-10, &
-0.65835936D-10, -0.65725415D-10, -0.65614894D-10, -0.65504372D-10, -0.65393851D-10, &
-0.65286488D-10, -0.65175966D-10, -0.65068603D-10, -0.64958082D-10, -0.64850718D-10, &
-0.64740197D-10, -0.64632833D-10, -0.64525470D-10, -0.64418106D-10, -0.64310743D-10, &
-0.64203379D-10, -0.64096016D-10, -0.63988652D-10, -0.63881289D-10, -0.63777083D-10, &
-0.63669720D-10, -0.63562356D-10, -0.63458150D-10, -0.63350787D-10, -0.63246581D-10, &
-0.63142375D-10, -0.63035012D-10, -0.62930806D-10, -0.62826600D-10, -0.62722394D-10, &
-0.62618189D-10, -0.62513983D-10, -0.62409777D-10, -0.62305571D-10, -0.62204523D-10, &
-0.62100318D-10, -0.61996112D-10, -0.61895064D-10, -0.61790858D-10, -0.61689810D-10, &
-0.61588762D-10, -0.61484556D-10, -0.61383508D-10, -0.61282460D-10, -0.61181412D-10, &
-0.61080364D-10, -0.60979316D-10, -0.60878268D-10, -0.60777220D-10, -0.60676172D-10, &
-0.60575124D-10, -0.60477234D-10, -0.60376186D-10, -0.60275138D-10, -0.60177248D-10, &
-0.60079357D-10, -0.59978309D-10, -0.59880419D-10, -0.59782529D-10, -0.59681481D-10, &
-0.59583591D-10, -0.59485700D-10, -0.59387810D-10, -0.59289920D-10, -0.59192030D-10, &
-0.59094139D-10, -0.58999407D-10, -0.58901517D-10, -0.58803626D-10, -0.58708894D-10, &
-0.58611004D-10, -0.58516271D-10, -0.58418381D-10, -0.58323648D-10, -0.58225758D-10, &
-0.58131026D-10, -0.58036293D-10, -0.57941560D-10, -0.57846828D-10, -0.57752095D-10, &
-0.57657363D-10, -0.57562630D-10, -0.57467898D-10, -0.57373165D-10, -0.57278433D-10, &
-0.57186858D-10, -0.57092126D-10, -0.57000551D-10, -0.56905818D-10, -0.56814244D-10, &
-0.56719511D-10, -0.56627936D-10, -0.56536362D-10, -0.56441629D-10, -0.56350054D-10, &
-0.56258480D-10, -0.56166905D-10, -0.56075330D-10, -0.55983755D-10, -0.55892181D-10, &
-0.55800606D-10, -0.55709031D-10, -0.55620614D-10, -0.55529039D-10, -0.55437464D-10, &
-0.55349047D-10, -0.55257473D-10, -0.55169056D-10, -0.55077481D-10, -0.54989064D-10, &
-0.54900647D-10, -0.54812230D-10, -0.54720655D-10, -0.54632238D-10, -0.54543821D-10, &
-0.54455404D-10, -0.54366987D-10, -0.54278570D-10, -0.54190153D-10, -0.54101736D-10, &
-0.54016477D-10, -0.53928060D-10, -0.53839643D-10, -0.53754384D-10, -0.53665967D-10, &
-0.53580707D-10, -0.53492290D-10, -0.53407031D-10, -0.53318614D-10, -0.53233355D-10, &
-0.53148096D-10, -0.53062836D-10, -0.52974419D-10, -0.52889160D-10, -0.52803901D-10, &
-0.52718641D-10, -0.52633382D-10, -0.52548123D-10, -0.52466021D-10, -0.52380762D-10, &
-0.52295503D-10, -0.52210244D-10, -0.52128142D-10, -0.52042883D-10, -0.51960781D-10, &
-0.51875522D-10, -0.51793421D-10, -0.51708161D-10, -0.51626060D-10, -0.51543958D-10, &
-0.51458699D-10, -0.51376598D-10, -0.51294496D-10, -0.51212395D-10, -0.51130293D-10, &
-0.51048192D-10, -0.50966090D-10, -0.50883989D-10, -0.50801887D-10, -0.50719786D-10, &
-0.50640842D-10, -0.50558740D-10, -0.50476639D-10, -0.50397695D-10, -0.50315593D-10, &
-0.50236650D-10, -0.50154548D-10, -0.50075604D-10, -0.49993503D-10, -0.49914559D-10, &
-0.49835615D-10, -0.49753514D-10, -0.49674570D-10, -0.49595626D-10, -0.49516683D-10, &
-0.49437739D-10, -0.49358795D-10, -0.49279851D-10, -0.49200908D-10, -0.49121964D-10, &
-0.49043020D-10, -0.48967234D-10, -0.48888290D-10, -0.48809347D-10, -0.48733561D-10, &
-0.48654617D-10, -0.48575673D-10, -0.48499887D-10, -0.48424101D-10, -0.48345157D-10, &
-0.48269371D-10, -0.48190428D-10, -0.48114642D-10, -0.48038856D-10, -0.47963069D-10, &
-0.47887283D-10, -0.47811497D-10, -0.47735711D-10, -0.47659925D-10, -0.47584139D-10, &
-0.47508353D-10, -0.47432567D-10, -0.47356781D-10, -0.47280995D-10, -0.47208367D-10, &
-0.47132581D-10, -0.47056795D-10, -0.46984167D-10, -0.46908381D-10, -0.46832595D-10, &
-0.46759967D-10, -0.46687338D-10, -0.46611552D-10, -0.46538924D-10, -0.46466296D-10, &
-0.46390510D-10, -0.46317882D-10, -0.46245253D-10, -0.46172625D-10, -0.46099997D-10, &
-0.46027369D-10, -0.45954740D-10, -0.45882112D-10, -0.45809484D-10, -0.45736856D-10, &
-0.45664227D-10, -0.45591599D-10, -0.45522129D-10, -0.45449500D-10, -0.45376872D-10, &
-0.45307401D-10, -0.45234773D-10, -0.45162145D-10, -0.45092674D-10, -0.45023204D-10, &
-0.44950576D-10, -0.44881105D-10, -0.44808477D-10, -0.44739006D-10, -0.44669536D-10, &
-0.44600065D-10, -0.44527437D-10, -0.44457967D-10, -0.44388496D-10, -0.44319026D-10, &
-0.44249555D-10, -0.44180085D-10, -0.44110614D-10, -0.44041144D-10, -0.43974831D-10, &
-0.43905360D-10, -0.43835890D-10, -0.43766419D-10, -0.43700107D-10, -0.43630636D-10, &
-0.43561166D-10, -0.43494853D-10, -0.43425382D-10, -0.43359070D-10, -0.43289599D-10, &
-0.43223286D-10, -0.43153816D-10, -0.43087503D-10, -0.43021190D-10, -0.42954878D-10, &
-0.42885407D-10, -0.42819094D-10, -0.42752781D-10, -0.42686469D-10, -0.42620156D-10, &
-0.42553843D-10, -0.42487530D-10, -0.42421218D-10, -0.42354905D-10, -0.42288592D-10, &
-0.42222279D-10, -0.42159124D-10, -0.42092812D-10, -0.42026499D-10, -0.41960186D-10, &
-0.41897031D-10, -0.41830718D-10, -0.41767563D-10, -0.41701251D-10, -0.41634938D-10, &
-0.41571783D-10, -0.41508628D-10, -0.41442315D-10, -0.41379160D-10, -0.41316005D-10, &
-0.41249692D-10, -0.41186537D-10, -0.41123382D-10, -0.41060227D-10, -0.40997072D-10, &
-0.40930760D-10, -0.40867605D-10, -0.40804450D-10, -0.40741295D-10, -0.40678140D-10, &
-0.40614985D-10, -0.40554987D-10, -0.40491832D-10, -0.40428677D-10, -0.40365522D-10, &
-0.40302367D-10, -0.40242370D-10, -0.40179215D-10, -0.40116060D-10, -0.40056063D-10, &
-0.39992908D-10, -0.39932910D-10, -0.39869755D-10, -0.39809758D-10, -0.39746603D-10, &
-0.39686606D-10, -0.39626609D-10, -0.39563454D-10, -0.39503456D-10, -0.39443459D-10, &
-0.39383462D-10, -0.39320307D-10, -0.39260310D-10, -0.39200312D-10, -0.39140315D-10, &
-0.39080318D-10, -0.39020321D-10, -0.38960323D-10, -0.38900326D-10, -0.38840329D-10, &
-0.38780332D-10, -0.38723492D-10, -0.38663495D-10, -0.38603498D-10, -0.38543500D-10, &
-0.38486661D-10, -0.38426664D-10, -0.38366666D-10, -0.38309827D-10, -0.38249830D-10, &
-0.38192990D-10, -0.38132993D-10, -0.38076153D-10, -0.38016156D-10, -0.37959317D-10, &
-0.37899319D-10, -0.37842480D-10, -0.37785640D-10, -0.37728801D-10, -0.37668803D-10, &
-0.37611964D-10, -0.37555124D-10, -0.37498285D-10, -0.37441445D-10, -0.37384606D-10, &
-0.37327766D-10, -0.37270927D-10, -0.37214087D-10, -0.37157248D-10, -0.37100408D-10, &
-0.37043569D-10, -0.36986729D-10, -0.36929890D-10, -0.36873050D-10, -0.36819369D-10, &
-0.36762529D-10, -0.36705690D-10, -0.36648850D-10, -0.36595168D-10, -0.36538329D-10, &
-0.36484647D-10, -0.36427808D-10, -0.36374126D-10, -0.36317286D-10, -0.36263605D-10, &
-0.36206765D-10, -0.36153083D-10, -0.36099402D-10, -0.36042562D-10, -0.35988880D-10, &
-0.35935199D-10, -0.35878359D-10, -0.35824677D-10, -0.35770996D-10, -0.35717314D-10, &
-0.35663632D-10, -0.35609950D-10, -0.35556269D-10, -0.35502587D-10, -0.35448905D-10, &
-0.35395223D-10, -0.35341541D-10, -0.35287860D-10, -0.35234178D-10, -0.35180496D-10, &
-0.35126814D-10, -0.35076290D-10, -0.35022609D-10, -0.34968927D-10, -0.34915245D-10, &
-0.34864721D-10, -0.34811039D-10, -0.34760515D-10, -0.34706834D-10, -0.34653152D-10, &
-0.34602628D-10, -0.34548946D-10, -0.34498422D-10, -0.34447898D-10, -0.34394216D-10, &
-0.34343692D-10, -0.34290011D-10, -0.34239487D-10, -0.34188963D-10, -0.34138439D-10, &
-0.34084757D-10, -0.34034233D-10, -0.33983709D-10, -0.33933185D-10, -0.33882661D-10, &
-0.33832137D-10, -0.33781613D-10, -0.33727931D-10, -0.33677407D-10, -0.33630041D-10, &
-0.33579517D-10, -0.33528993D-10, -0.33478469D-10, -0.33427945D-10, -0.33377421D-10, &
-0.33326897D-10, -0.33276373D-10, -0.33229007D-10, -0.33178483D-10, -0.33127959D-10, &
-0.33080592D-10, -0.33030068D-10, -0.32979544D-10, -0.32932178D-10, -0.32881654D-10, &
-0.32834288D-10, -0.32783764D-10, -0.32736397D-10, -0.32685873D-10, -0.32638507D-10, &
-0.32587983D-10, -0.32540617D-10, -0.32493251D-10, -0.32442727D-10, -0.32395360D-10 / 

Data roCs/ 3.3/
Data ECsmax/ 3E3/  !!!valores tomados arbitrariamente
Data rfCs/12.0/
End

Double Precision Function V_Cs_Cs(xx)
Parameter(N=4350)
Implicit Real*8(A-H,O-Z)  !! When Fortran was originally developed memory was at a premium. Variables and procedure names could have a maximum of 6 characters, and variables were often implicitly typed. This means that the first letter of the variable name determines its type.
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_Cs/RCs(N),ECs(N),roCs,ECsmax,rfCs

If(xx.le.roCs)Then
  V_Cs_CS=ECsmax
  Return
Endif
If(xx.gt.rfCs)Then
  V_Cs_CS=0
  Return
Endif
Call spline2(RCs,ECs,a,b,c,d,0,0,N)
Call Findi2(RCs,xx,N,ix)
V_Cs_CS= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End


Block Data Inicio_Cs_Cs_tsu
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_Cs_tsu/RCstsu(N),ECstsu(N),roCstsu,ECsmaxtsu,rfCstsu
Data RCstsu/                                                                         &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /

Data ECstsu/                                                                         &
 0.20213707D+05,  0.10583516D+05,  0.52150246D+04,  0.22587388D+04,  0.13455174D+04, &
 0.70512564D+03,  0.27125075D+03, -0.10736351D+02, -0.22293717D+03, -0.30124938D+03, &
-0.33914238D+03, -0.34356323D+03, -0.32777448D+03, -0.27977668D+03, -0.22198985D+03, &
-0.16578189D+03, -0.89048559D+02, -0.48945130D+02, -0.28103978D+02, -0.19893827D+02, &
-0.15157202D+02, -0.11683676D+02 /

Data roCstsu/ 3.3/
Data ECsmaxtsu/ 10E3/  !!!valores tomados arbitrariamente
Data rfCstsu/30.0/

End

Double Precision Function V_Cs_Cs_tsu(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z)  !! When Fortran was originally developed memory was at a premium. Variables and procedure names could have a maximum of 6 characters, and variables were often implicitly typed. This means that the first letter of the variable name determines its type.
Double Precision:: a(N),b(N),c(N),d(N),C_6=-253924406.172384,C_8= 28604481648.8368  
Common/Param_Spline2_Cs_Cs_tsu/RCstsu(N),ECstsu(N),roCstsu,ECsmaxtsu,rfCstsu

If(xx.le.roCstsu)Then
  V_Cs_CS_tsu=ECsmaxtsu
  Return
Endif
If(xx.gt.rfCstsu)Then
  V_Cs_CS_tsu=0
  Return
Endif
If(xx.gt.RCstsu(22))then
V_Cs_CS_tsu=C_6/xx**6+C_8/xx**8
Return
Endif
Call spline2(RCstsu,ECstsu,a,b,c,d,0,0,N)
Call Findi2(RCstsu,xx,N,ix)
V_Cs_CS_tsu= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Rb_Rb
Parameter(N=298)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Rb_Rb/Rrb(N),ECrb(N),rorb,ECrbmax,rfrb
Data Rrb/                                                                            &
  0.10000000D+00,  0.20000000D+00,  0.30000000D+00,  0.40000001D+00,  0.50000001D+00, &
 0.60000001D+00,  0.70000001D+00,  0.80000001D+00,  0.90000001D+00,  0.10000000D+01, &
 0.11000000D+01,  0.12000000D+01,  0.13000000D+01,  0.14000000D+01,  0.15000000D+01, &
 0.16000000D+01,  0.17000000D+01,  0.18000000D+01,  0.19000000D+01,  0.20000000D+01, &
 0.21000000D+01,  0.22000000D+01,  0.23000000D+01,  0.24000000D+01,  0.25000000D+01, &
 0.26000000D+01,  0.27000000D+01,  0.28000000D+01,  0.29000000D+01,  0.30000000D+01, &
 0.31000000D+01,  0.32000000D+01,  0.33000000D+01,  0.34000001D+01,  0.35000001D+01, &
 0.36000001D+01,  0.37000001D+01,  0.38000001D+01,  0.39000001D+01,  0.40000001D+01, &
 0.41000001D+01,  0.42000001D+01,  0.43000001D+01,  0.44000001D+01,  0.45000001D+01, &
 0.46000001D+01,  0.47000001D+01,  0.48000001D+01,  0.49000001D+01,  0.50000001D+01, &
 0.51000001D+01,  0.52000001D+01,  0.53000001D+01,  0.54000001D+01,  0.55000001D+01, &
 0.56000001D+01,  0.57000001D+01,  0.58000001D+01,  0.59000001D+01,  0.60000001D+01, &
 0.61000001D+01,  0.62000001D+01,  0.63000001D+01,  0.64000001D+01,  0.65000001D+01, &
 0.66000001D+01,  0.67000001D+01,  0.68000001D+01,  0.69000001D+01,  0.70000001D+01, &
 0.71000001D+01,  0.72000001D+01,  0.73000001D+01,  0.74000001D+01,  0.75000001D+01, &
 0.76000001D+01,  0.77000001D+01,  0.78000001D+01,  0.79000001D+01,  0.80000001D+01, &
 0.81000001D+01,  0.82000001D+01,  0.83000001D+01,  0.84000001D+01,  0.85000001D+01, &
 0.86000001D+01,  0.87000001D+01,  0.88000001D+01,  0.89000001D+01,  0.90000001D+01, &
 0.91000001D+01,  0.92000001D+01,  0.93000001D+01,  0.94000001D+01,  0.95000001D+01, &
 0.96000001D+01,  0.97000001D+01,  0.98000001D+01,  0.99000001D+01,  0.10000000D+02, &
 0.10100000D+02,  0.10200000D+02,  0.10300000D+02,  0.10400000D+02,  0.10500000D+02, &
 0.10600000D+02,  0.10700000D+02,  0.10800000D+02,  0.10900000D+02,  0.11000000D+02, &
 0.11100000D+02,  0.11200000D+02,  0.11300000D+02,  0.11400000D+02,  0.11500000D+02, &
 0.11600000D+02,  0.11700000D+02,  0.11800000D+02,  0.11900000D+02,  0.12000000D+02, &
 0.12100000D+02,  0.12200000D+02,  0.12300000D+02,  0.12400000D+02,  0.12500000D+02, &
 0.12600000D+02,  0.12700000D+02,  0.12800000D+02,  0.12900000D+02,  0.13000000D+02, &
 0.13100000D+02,  0.13200000D+02,  0.13300000D+02,  0.13400000D+02,  0.13500000D+02, &
 0.13600000D+02,  0.13700000D+02,  0.13800000D+02,  0.13900000D+02,  0.14000000D+02, &
 0.14100000D+02,  0.14200000D+02,  0.14300000D+02,  0.14400000D+02,  0.14500000D+02, &
 0.14600000D+02,  0.14700000D+02,  0.14800000D+02,  0.14900000D+02,  0.15000000D+02, &
 0.15100000D+02,  0.15200000D+02,  0.15300000D+02,  0.15400000D+02,  0.15500000D+02, &
 0.15600000D+02,  0.15700000D+02,  0.15800000D+02,  0.15900000D+02,  0.16000000D+02, &
 0.16100000D+02,  0.16200000D+02,  0.16300000D+02,  0.16400000D+02,  0.16500000D+02, &
 0.16600000D+02,  0.16700000D+02,  0.16800000D+02,  0.16900000D+02,  0.17000000D+02, &
 0.17100000D+02,  0.17200000D+02,  0.17300000D+02,  0.17400000D+02,  0.17500000D+02, &
 0.17600000D+02,  0.17700000D+02,  0.17800000D+02,  0.17900000D+02,  0.18000000D+02, &
 0.18100000D+02,  0.18200000D+02,  0.18300000D+02,  0.18400000D+02,  0.18500000D+02, &
 0.18600000D+02,  0.18700000D+02,  0.18800000D+02,  0.18900000D+02,  0.19000000D+02, &
 0.19100000D+02,  0.19200000D+02,  0.19300000D+02,  0.19400000D+02,  0.19500000D+02, &
 0.19600000D+02,  0.19700000D+02,  0.19800000D+02,  0.19900000D+02,  0.20000000D+02, &
 0.20100000D+02,  0.20200000D+02,  0.20300000D+02,  0.20400000D+02,  0.20500000D+02, &
 0.20600000D+02,  0.20700000D+02,  0.20800000D+02,  0.20900000D+02,  0.21000000D+02, &
 0.21100000D+02,  0.21200000D+02,  0.21300000D+02,  0.21400000D+02,  0.21500000D+02, &
 0.21600000D+02,  0.21700000D+02,  0.21800000D+02,  0.21900000D+02,  0.22000000D+02, &
 0.22100000D+02,  0.22200000D+02,  0.22300000D+02,  0.22400000D+02,  0.22500000D+02, &
 0.22600000D+02,  0.22700000D+02,  0.22800000D+02,  0.22900000D+02,  0.23000000D+02, &
 0.23100000D+02,  0.23200000D+02,  0.23300000D+02,  0.23400000D+02,  0.23500000D+02, &
 0.23600000D+02,  0.23700000D+02,  0.23800000D+02,  0.23900000D+02,  0.24000000D+02, &
 0.24100000D+02,  0.24200000D+02,  0.24300000D+02,  0.24400000D+02,  0.24500000D+02, &
 0.24600000D+02,  0.24700000D+02,  0.24800000D+02,  0.24900000D+02,  0.25000000D+02, &
 0.25100000D+02,  0.25200000D+02,  0.25300000D+02,  0.25400000D+02,  0.25500000D+02, &
 0.25600000D+02,  0.25700000D+02,  0.25800000D+02,  0.25900000D+02,  0.26000000D+02, &
 0.26100000D+02,  0.26200000D+02,  0.26300000D+02,  0.26400000D+02,  0.26500000D+02, &
 0.26600000D+02,  0.26700000D+02,  0.26800000D+02,  0.26900000D+02,  0.27000000D+02, &
 0.27100000D+02,  0.27200000D+02,  0.27300000D+02,  0.27400000D+02,  0.27500000D+02, &
 0.27600000D+02,  0.27700000D+02,  0.27800000D+02,  0.27900000D+02,  0.28000000D+02, &
 0.28100000D+02,  0.28200000D+02,  0.28300000D+02,  0.28400000D+02,  0.28500000D+02, &
 0.28600000D+02,  0.28700000D+02,  0.28800000D+02,  0.28900000D+02,  0.29000000D+02, &
 0.29100000D+02,  0.29200000D+02,  0.29300000D+02,  0.29400000D+02,  0.29500000D+02, &
 0.29600000D+02,  0.29700000D+02,  0.29800000D+02  /

Data ECrb/                                                                           &
  0.46054087D+06,  0.45680239D+06,  0.45060884D+06,  0.44201607D+06,  0.43110232D+06,&
 0.41796813D+06,  0.40273642D+06,  0.38555244D+06,  0.36658380D+06,  0.34602045D+06, &
 0.32407468D+06,  0.30098114D+06,  0.27699682D+06,  0.25240106D+06,  0.22749553D+06, &
 0.20260429D+06,  0.17807370D+06,  0.15427249D+06,  0.13159174D+06,  0.11044486D+06, &
 0.91263465D+05,  0.74413806D+05,  0.60044868D+05,  0.48087909D+05,  0.38335354D+05, &
 0.30513279D+05,  0.24330053D+05,  0.19504170D+05,  0.15777870D+05,  0.12922141D+05, &
 0.10737041D+05,  0.90524715D+04,  0.77316783D+04,  0.66703246D+04,  0.57918339D+04, &
 0.50419027D+04,  0.43834334D+04,  0.37923020D+04,  0.32540017D+04,  0.27610655D+04, &
 0.23111010D+04,  0.19047749D+04,  0.15433128D+04,  0.12268988D+04,  0.95407674D+03, &
 0.72168897D+03,  0.52528745D+03,  0.36001282D+03,  0.22127338D+03,  0.10508502D+03, &
 0.82022293D+01, -0.71929923D+02, -0.13739448D+03, -0.18999913D+03, -0.23139700D+03, &
-0.26315583D+03, -0.28676228D+03, -0.30354604D+03, -0.31462335D+03, -0.32091948D+03, &
-0.32323099D+03, -0.32226301D+03, -0.31863204D+03, -0.31286433D+03, -0.30540469D+03, &
-0.29663077D+03, -0.28686447D+03, -0.27637891D+03, -0.26540408D+03, -0.25413218D+03, &
-0.24272259D+03, -0.23130635D+03, -0.21999006D+03, -0.20885939D+03, -0.19798215D+03, &
-0.18741102D+03, -0.17718589D+03, -0.16733590D+03, -0.15788108D+03, -0.14883393D+03, &
-0.14020065D+03, -0.13198230D+03, -0.12417575D+03, -0.11677457D+03, -0.10976974D+03, &
-0.10315019D+03, -0.96903351D+02, -0.91015511D+02, -0.85472188D+02, -0.80258387D+02, &
-0.75358826D+02, -0.70758130D+02, -0.66441010D+02, -0.62392402D+02, -0.58597587D+02, &
-0.55042261D+02, -0.51712601D+02, -0.48595298D+02, -0.45677586D+02, -0.42947248D+02, &
-0.40392625D+02, -0.38002614D+02, -0.35766675D+02, -0.33674833D+02, -0.31717671D+02, &
-0.29886316D+02, -0.28172426D+02, -0.26568167D+02, -0.25066194D+02, -0.23659624D+02, &
-0.22342015D+02, -0.21107346D+02, -0.19949989D+02, -0.18864696D+02, -0.17846578D+02, &
-0.16891082D+02, -0.15993971D+02, -0.15151311D+02, -0.14359445D+02, -0.13614976D+02, &
-0.12914754D+02, -0.12255857D+02, -0.11635571D+02, -0.11051378D+02, -0.10500938D+02, &
-0.99820759D+01, -0.94927698D+01, -0.90311379D+01, -0.85954284D+01, -0.81840102D+01, &
-0.77953630D+01, -0.74280698D+01, -0.70808092D+01, -0.67523489D+01, -0.64415391D+01, &
-0.61473074D+01, -0.58686528D+01, -0.56046408D+01, -0.53543992D+01, -0.51171134D+01, &
-0.48920224D+01, -0.46784152D+01, -0.44756273D+01, -0.42830371D+01, -0.41000628D+01, &
-0.39261601D+01, -0.37608187D+01, -0.36035608D+01, -0.34539381D+01, -0.33115303D+01, &
-0.31759427D+01, -0.30468050D+01, -0.29237691D+01, -0.28065081D+01, -0.26947144D+01, &
-0.25880989D+01, -0.24863895D+01, -0.23893298D+01, -0.22966787D+01, -0.22082087D+01, &
-0.21237057D+01, -0.20429675D+01, -0.19658035D+01, -0.18920338D+01, -0.18214886D+01, &
-0.17540077D+01, -0.16894396D+01, -0.16276411D+01, -0.15684770D+01, -0.15118194D+01, &
-0.14575474D+01, -0.14055465D+01, -0.13557084D+01, -0.13079305D+01, -0.12621158D+01, &
-0.12181721D+01, -0.11760124D+01, -0.11355540D+01, -0.10967185D+01, -0.10594317D+01, &
-0.10236231D+01, -0.98922565D+00, -0.95617595D+00, -0.92441367D+00, -0.89388151D+00, &
-0.86452507D+00, -0.83629261D+00, -0.80913497D+00, -0.78300538D+00, -0.75785939D+00, &
-0.73365466D+00, -0.71035094D+00, -0.68790988D+00, -0.66629500D+00, -0.64547153D+00, &
-0.62540635D+00, -0.60606789D+00, -0.58742608D+00, -0.56945224D+00, -0.55211902D+00, &
-0.53540033D+00, -0.51927128D+00, -0.50370811D+00, -0.48868816D+00, -0.47418976D+00, &
-0.46019224D+00, -0.44667585D+00, -0.43362170D+00, -0.42101176D+00, -0.40882876D+00, &
-0.39705623D+00, -0.38567837D+00, -0.37468009D+00, -0.36404695D+00, -0.35376512D+00, &
-0.34382136D+00, -0.33420300D+00, -0.32489789D+00, -0.31589441D+00, -0.30718141D+00, &
-0.29874820D+00, -0.29058455D+00, -0.28268063D+00, -0.27502703D+00, -0.26761472D+00, &
-0.26043502D+00, -0.25347962D+00, -0.24674054D+00, -0.24021012D+00, -0.23388099D+00, &
-0.22774608D+00, -0.22179862D+00, -0.21603207D+00, -0.21044017D+00, -0.20501690D+00, &
-0.19975645D+00, -0.19465328D+00, -0.18970201D+00, -0.18489751D+00, -0.18023482D+00, &
-0.17570916D+00, -0.17131596D+00, -0.16705080D+00, -0.16290942D+00, -0.15888774D+00, &
-0.15498181D+00, -0.15118785D+00, -0.14750218D+00, -0.14392130D+00, -0.14044179D+00, &
-0.13706040D+00, -0.13377397D+00, -0.13057946D+00, -0.12747393D+00, -0.12445456D+00, &
-0.12151863D+00, -0.11866349D+00, -0.11588662D+00, -0.11318556D+00, -0.11055796D+00, &
-0.10800153D+00, -0.10551406D+00, -0.10309345D+00, -0.10073763D+00, -0.98444621D-01, &
-0.96212512D-01, -0.94039456D-01, -0.91923666D-01, -0.89863417D-01, -0.87857043D-01, &
-0.85902934D-01, -0.83999533D-01, -0.82145335D-01, -0.80338886D-01, -0.78578780D-01, &
-0.76863656D-01, -0.75192200D-01, -0.73563140D-01, -0.71975245D-01, -0.70427325D-01, &
-0.68918230D-01, -0.67446843D-01, -0.66012088D-01, -0.64612919D-01, -0.63248328D-01, &
-0.61917335D-01, -0.60618995D-01, -0.59352389D-01, -0.58116631D-01, -0.56910860D-01, &
-0.55734243D-01, -0.54585973D-01, -0.53465268D-01, -0.52371372D-01, -0.51303549D-01, &
-0.50261089D-01, -0.49243303D-01, -0.48249521D-01  /

Data rorb/ 5/
Data ECrbmax/ 10E3/  !!!valores tomados arbitrariamente
Data rfrb/30.0/

End

Double Precision Function V_Rb_Rb(xx)
Parameter(N=298)
Implicit Real*8(A-H,O-Z)  !! When Fortran was originally developed memory was at a premium. Variables and procedure names could have a maximum of 6 characters, and variables were often implicitly typed. This means that the first letter of the variable name determines its type.
Double Precision:: a(N),b(N),c(N),d(N),C_6=-253924406.172384,C_8= 28604481648.8368  
Common/Param_Spline2_Rb_Rb/Rrb(N),ECrb(N),rorb,ECrbmax,rfrb

If(xx.le.rorb)Then
  V_Rb_Rb=ECrbmax
  Return
Endif
If(xx.gt.rfrb)Then
  V_Rb_Rb=0
  Return
Endif
Call spline2(Rrb,ECrb,a,b,c,d,0,0,N)
Call Findi2(Rrb,xx,N,ix)
V_Rb_Rb= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Cs_ssg
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_ssg/RCs(N),ECcs(N),rocs,ECcsmax,rfcs
Data RCs/                                                                         &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02  /

Data ECcs/                                                                        &
 0.10533623D+05,  0.18359160D+04, -0.31340672D+04, -0.50678735D+04, -0.52503914D+04, &
-0.50858727D+04, -0.46924170D+04, -0.41666515D+04, -0.33412156D+04, -0.27734521D+04, &
-0.22577915D+04, -0.18027597D+04, -0.14175141D+04, -0.86017119D+03, -0.51660795D+03, &
-0.30061783D+03, -0.12125761D+03, -0.55576405D+02, -0.29367078D+02, -0.20209602D+02, &
-0.15157202D+02, -0.11683676D+02  /

Data rocs/ 1/
Data ECcsmax/ 10E3/  !!!valores tomados arbitrariamente
Data rfcs/30.0/
End

Double Precision Function V_Cs_ssg(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z)  !! When Fortran was originally developed memory was at a premium. Variables and procedure names could have a maximum of 6 characters, and variables were often implicitly typed. This means that the first letter of the variable name determines its type.
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_ssg/RCs(N),ECcs(N),rocs,ECcsmax,rfcs

If(xx.le.rocs)Then
  V_Cs_ssg=ECcsmax
  Return
Endif
If(xx.gt.rfcs)Then
  V_Cs_ssg=0
  Return
Endif
Call spline2(RCs,ECcs,a,b,c,d,0,0,N)
Call Findi2(RCs,xx,N,ix)
V_Cs_ssg= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Cs_tsu
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_tsu/RCstsu(N),ECcstsu(N),rocstsu,ECcsmaxtsu,rfcstsu
Data RCstsu/                                                                         &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02   /

Data ECcstsu/                                                                        &
 0.20213707D+05,  0.10583516D+05,  0.52150246D+04,  0.22587388D+04,  0.13455174D+04, &
 0.70512564D+03,  0.27125075D+03, -0.10736351D+02, -0.22293717D+03, -0.30124938D+03, &
-0.33914238D+03, -0.34356323D+03, -0.32777448D+03, -0.27977668D+03, -0.22198985D+03, &
-0.16578189D+03, -0.89048559D+02, -0.48945130D+02, -0.28103978D+02, -0.19893827D+02, &
-0.15157202D+02, -0.11683676D+02 /

Data rocstsu/ 1/
Data ECcsmaxtsu/ 10E3/  !!!valores tomados arbitrariamente
Data rfcstsu/30.0/
End

Double Precision Function V_Cs_tsu(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_tsu/RCstsu(N),ECcstsu(N),rocstsu,ECcsmaxtsu,rfcstsu

If(xx.le.rocstsu)Then
  V_Cs_tsu=ECcsmaxtsu
  Return
Endif
If(xx.gt.rfcstsu)Then
  V_Cs_tsu=0
  Return
Endif
Call spline2(RCstsu,ECcstsu,a,b,c,d,0,0,N)
Call Findi2(RCstsu,xx,N,ix)
V_Cs_tsu= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End


Block Data Inicio_Cs_ssg_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_ssg_ps/RCssgps(N),ECssgps(N),rocsssgps,ECcsmaxssgps,rfcsssgps
Data RCssgps/                                                                         &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /
Data ECssgps/                                                                        &
 0.27053394D+05,  0.19978139D+05,  0.16649554D+05,  0.14393973D+05,  0.13569800D+05, &
 0.12966354D+05,  0.12553952D+05,  0.12299753D+05,  0.12166496D+05,  0.12193653D+05, &
 0.12314595D+05,  0.12527743D+05,  0.12817309D+05,  0.13547696D+05,  0.14378184D+05, &
 0.15184358D+05,  0.16359357D+05,  0.16848493D+05,  0.16958067D+05,  0.16939436D+05, &
 0.16892701D+05,  0.16845966D+05 /

Data rocsssgps/ 1/
Data ECcsmaxssgps/ 10E7/  !!!valores tomados arbitrariamente
Data rfcsssgps/30.0/
End

Double Precision Function V_Cs_ssg_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_ssg_ps/RCssgps(N),ECssgps(N),rocsssgps,ECcsmaxssgps,rfcsssgps

If(xx.le.rocsssgps)Then
  V_Cs_ssg_ps=ECcsmaxssgps
  Return
Endif
If(xx.gt.rfcsssgps)Then
  V_Cs_ssg_ps=0
  Return
Endif
Call spline2(RCssgps,ECssgps,a,b,c,d,0,0,N)
Call Findi2(RCssgps,xx,N,ix)
V_Cs_ssg_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Cs_ssu_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_ssu_ps/RCssups(N),ECssups(N),rocsssups,ECcsmaxssups,rfcsssups
Data RCssups/                                                                         &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /
Data ECssups/                                                                        &
 0.27458849D+05,  0.19688573D+05,  0.14042516D+05,  0.10459732D+05,  0.94479889D+04, &
 0.88562265D+04,  0.85928702D+04,  0.85786603D+04,  0.88631736D+04,  0.92032633D+04, &
 0.96194548D+04,  0.10093117D+05,  0.10599620D+05,  0.11627784D+05,  0.12610476D+05, &
 0.13485804D+05,  0.14772272D+05,  0.15491923D+05,  0.15866748D+05,  0.16072318D+05, &
 0.16199575D+05,  0.16285781D+05 /

Data rocsssups/ 1/
Data ECcsmaxssups/ 10E7/  !!!valores tomados arbitrariamente
Data rfcsssups/30.0/
End

Double Precision Function V_Cs_ssu_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_ssu_ps/RCssups(N),ECssups(N),rocsssups,ECcsmaxssups,rfcsssups

If(xx.le.rocsssups)Then
  V_Cs_ssu_ps=ECcsmaxssgps
  Return
Endif
If(xx.gt.rfcsssups)Then
  V_Cs_ssu_ps=0
  Return
Endif
Call spline2(RCssups,ECssups,a,b,c,d,0,0,N)
Call Findi2(RCssups,xx,N,ix)
V_Cs_ssu_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Cs_spg_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_spg_ps/RCspgps(N),ECspgps(N),rocspgps,ECcsmaxspgps,rfcspgps
Data RCspgps/                                                                         &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /
Data ECspgps/                                                                        &
 0.33044278D+05,  0.24908334D+05,  0.19633313D+05,  0.16773654D+05,  0.15946955D+05, &
 0.15417716D+05,  0.15108572D+05,  0.14955105D+05,  0.14911844D+05,  0.14941211D+05, &
 0.14998051D+05,  0.15077626D+05,  0.15168569D+05,  0.15349824D+05,  0.15526974D+05, &
 0.15693072D+05,  0.15957375D+05,  0.16139578D+05,  0.16262098D+05,  0.16341989D+05, &
 0.16398197D+05,  0.16439248D+05 /

Data rocspgps/ 1/
Data ECcsmaxspgps/ 10E7/  !!!valores tomados arbitrariamente
Data rfcspgps/30.0/
End

Double Precision Function V_Cs_spg_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_spg_ps/RCspgps(N),ECspgps(N),rocspgps,ECcsmaxspgps,rfcspgps

If(xx.le.rocspgps)Then
  V_Cs_spg_ps=ECcsmaxspgps
  Return
Endif
If(xx.gt.rfcspgps)Then
  V_Cs_spg_ps=0
  Return
Endif
Call spline2(RCspgps,ECspgps,a,b,c,d,0,0,N)
Call Findi2(RCspgps,xx,N,ix)
V_Cs_spg_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End


Block Data Inicio_Cs_spu_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_spu_ps/RCspups(N),ECspups(N),rocspups,ECcsmaxspups,rfcspups
Data RCspups/                                                                        &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /
Data ECspups/                                                                        &
 0.23340511D+05,  0.18208220D+05,  0.15642863D+05,  0.13921258D+05,  0.13552749D+05, &
 0.13481068D+05,  0.13628219D+05,  0.13931047D+05,  0.14491548D+05,  0.14912792D+05, &
 0.15314773D+05,  0.15683283D+05,  0.16002531D+05,  0.16462615D+05,  0.16728498D+05, &
 0.16864281D+05,  0.16927121D+05,  0.16891438D+05,  0.16840282D+05,  0.16791969D+05, &
 0.16753444D+05,  0.16724077D+05  /

Data rocspups/ 1/
Data ECcsmaxspups/ 10E7/  !!!valores tomados arbitrariamente
Data rfcspups/30.0/
End


Double Precision Function V_Cs_spu_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_spu_ps/RCspups(N),ECspups(N),rocspups,ECcsmaxspups,rfcspups

If(xx.le.rocspups)Then
  V_Cs_spu_ps=ECcsmaxspups
  Return
Endif
If(xx.gt.rfcspups)Then
  V_Cs_spu_ps=0
  Return
Endif
Call spline2(RCspups,ECspups,a,b,c,d,0,0,N)
Call Findi2(RCspups,xx,N,ix)
V_Cs_spu_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Cs_tpg_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_tpg_ps/RCtpgps(N),ECtpgps(N),roctpgps,ECcsmaxtpgps,rfctpgps
Data RCtpgps/                                                                        &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02  /
Data ECtpgps/                                                                        &
 0.31817176D+05,  0.24118581D+05,  0.19863513D+05,  0.17369837D+05,  0.16676079D+05, &
 0.16285781D+05,  0.16128525D+05,  0.16136736D+05,  0.16321464D+05,  0.16496403D+05, &
 0.16668185D+05,  0.16819757D+05,  0.16936278D+05,  0.17049957D+05,  0.17069535D+05, &
 0.17051852D+05,  0.16978276D+05,  0.16904385D+05,  0.16843440D+05,  0.16792600D+05, &
 0.16753444D+05,  0.16724077D+05   /

Data roctpgps/ 1/
Data ECcsmaxtpgps/ 10E7/  !!!valores tomados arbitrariamente
Data rfctpgps/30.0/
End

Double Precision Function V_Cs_tpg_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_tpg_ps/RCtpgps(N),ECtpgps(N),roctpgps,ECcsmaxtpgps,rfctpgps

If(xx.le.roctpgps)Then
  V_Cs_tpg_ps=ECcsmaxtpgps
  Return
Endif
If(xx.gt.rfctpgps)Then
  V_Cs_tpg_ps=0
  Return
Endif
Call spline2(RCtpgps,ECtpgps,a,b,c,d,0,0,N)
Call Findi2(RCtpgps,xx,N,ix)
V_Cs_tpg_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Cs_tpu_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_tpu_ps/RCtpups(N),ECtpups(N),roctpups,ECcsmaxtpups,rfctpups
Data RCtpups/                                                                        &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /
Data ECtpups/                                                                        &
 0.20371279D+05,  0.11877878D+05,  0.78078534D+04,  0.65352800D+04,  0.65633840D+04, &
 0.68987371D+04,  0.74538696D+04,  0.81555217D+04,  0.92806282D+04,  0.10109853D+05, &
 0.10924553D+05,  0.11706728D+05,  0.12433958D+05,  0.13655375D+05,  0.14552492D+05, &
 0.15164780D+05,  0.15817803D+05,  0.16105474D+05,  0.16253888D+05,  0.16340095D+05, &
 0.16397882D+05,  0.16438932D+05  /

Data roctpups/ 1/
Data ECcsmaxtpups/ 10E7/  !!!valores tomados arbitrariamente
Data rfctpups/30.0/
End

Double Precision Function V_Cs_tpu_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_tpu_ps/RCtpups(N),ECtpups(N),roctpups,ECcsmaxtpups,rfctpups

If(xx.le.roctpups)Then
  V_Cs_tpu_ps=ECcsmaxtpups
  Return
Endif
If(xx.gt.rfctpups)Then
  V_Cs_tpu_ps=0
  Return
Endif
Call spline2(RCtpups,ECtpups,a,b,c,d,0,0,N)
Call Findi2(RCtpups,xx,N,ix)
V_Cs_tpu_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End


Block Data Inicio_Cs_tsg_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_tsg_ps/RCtsgps(N),ECtsgps(N),roctsgps,ECcsmaxtsgps,rfctsgps
Data RCtsgps/                                                                        &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /
Data ECtsgps/                                                                        &
 0.26690253D+05,  0.19898564D+05,  0.16728498D+05,  0.14695223D+05,  0.13457385D+05, &
 0.12646474D+05,  0.12184811D+05,  0.11979873D+05,  0.11992188D+05,  0.12128603D+05, &
 0.12332910D+05,  0.12586161D+05,  0.12866569D+05,  0.13436859D+05,  0.13977150D+05, &
 0.14455234D+05,  0.15172043D+05,  0.15624549D+05,  0.15905904D+05,  0.16083054D+05, &
 0.16202417D+05,  0.16286729D+05   /

Data roctsgps/ 1/
Data ECcsmaxtsgps/ 10E7/  !!!valores tomados arbitrariamente
Data rfctsgps/30.0/
End


Double Precision Function V_Cs_tsg_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_tsg_ps/RCtsgps(N),ECtsgps(N),roctsgps,ECcsmaxtsgps,rfctsgps

If(xx.le.roctsgps)Then
  V_Cs_tsg_ps=ECcsmaxtsgps
  Return
Endif
If(xx.gt.rfctsgps)Then
  V_Cs_tsg_ps=0
  Return
Endif
Call spline2(RCtsgps,ECtsgps,a,b,c,d,0,0,N)
Call Findi2(RCtsgps,xx,N,ix)
V_Cs_tsg_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

Block Data Inicio_Cs_tsu_ps
Parameter(N=22)
Implicit Real*8(A-H,O-Z)
!!!Energy in K and R in \AA
Common/Param_Spline2_Cs_tsu_ps/RCtsups(N),ECtsups(N),roctsups,ECcsmaxtsups,rfctsups
Data RCtsups/                                                                        &
 0.28046391D+01,  0.33338163D+01,  0.38629935D+01,  0.43921707D+01,  0.46567593D+01, &
 0.49213479D+01,  0.51859365D+01,  0.54505250D+01,  0.58209491D+01,  0.60855377D+01, &
 0.63501263D+01,  0.66147149D+01,  0.68793035D+01,  0.74084806D+01,  0.79376578D+01, &
 0.84668350D+01,  0.95251894D+01,  0.10583544D+02,  0.11641898D+02,  0.12700253D+02, &
 0.13758607D+02,  0.14816961D+02 /
Data ECtsups/                                                                        &
 0.38036366D+05,  0.28540695D+05,  0.23256831D+05,  0.20252547D+05,  0.19363009D+05, &
 0.18789877D+05,  0.18459577D+05,  0.18305479D+05,  0.18284953D+05,  0.18327899D+05, &
 0.18366107D+05,  0.18354424D+05,  0.18261270D+05,  0.17938864D+05,  0.17652456D+05, &
 0.17455412D+05,  0.17232791D+05,  0.17115322D+05,  0.17031642D+05,  0.16959014D+05, &
 0.16897753D+05,  0.16847229D+05  /

Data roctsups/ 1/
Data ECcsmaxtsups/ 10E7/  !!!valores tomados arbitrariamente
Data rfctsups/30.0/
End

Double Precision Function V_Cs_tsu_ps(xx)
Parameter(N=22)
Implicit Real*8(A-H,O-Z) 
Double Precision:: a(N),b(N),c(N),d(N)
Common/Param_Spline2_Cs_tsu_ps/RCtsups(N),ECtsups(N),roctsups,ECcsmaxtsups,rfctsups

If(xx.le.roctsups)Then
  V_Cs_tsu_ps=ECcsmaxtsups
  Return
Endif
If(xx.gt.rfctsups)Then
  V_Cs_tsu_ps=0
  Return
Endif
Call spline2(RCtsups,ECtsups,a,b,c,d,0,0,N)
Call Findi2(RCtsups,xx,N,ix)
V_Cs_tsu_ps= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return 
End

