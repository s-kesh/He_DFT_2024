 &INPUT
 NX      =         320,
 NY      =         320,
 NZ      =         320,
 XMAX    =   64.0000000000000     ,
 YMAX    =   64.0000000000000     ,
 ZMAX    =   64.0000000000000     ,
 N_IMP   =           2,
 FILE4   = imp.input                                                                       ,
 FILE5   = potential_multi_impurity.input                                                  ,
 FILE6   = potential_multi_impurity.res                                                    ,
 FILE7   = potential_multi_impurity.out                                                    ,
 SELEC_AR_AR     = Na2_sspg                                                                        
 /
 &IMP
 RIMP    = 2*0.000000000000000E+000  ,   1.53600000000000     ,  -1.53600000000000     , 2*42.0000000000000        ,
 SELEC_GS_K      = NA                                                                              NA                                
                                               ,
 UMAX_GS_K       = 2*10000.0000000000        ,
 R_CUTOFF_GS_K   = 2*3.00000000000000        ,
 FAC_X_Y_Z       =   1.00000000000000     
 /
#
# Multi impurity potential, the following parameters had been used:
# N_imp ...........:   2
# rimp(1).............:   0.000000E+00   1.536000E+00   4.200000E+01
# selec_gs_k(1).......:NA                                                                              
# umax_gs_k(1),........:   1.000000E+04
# r_cutoff_gs_k(1)....:   3.000000E+00
# rimp(2).............:   0.000000E+00  -1.536000E+00   4.200000E+01
# selec_gs_k(2).......:NA                                                                              
# umax_gs_k(2),........:   1.000000E+04
# r_cutoff_gs_k(2)....:   3.000000E+00
#
         320         320         320  0.400000000000000     
  0.400000000000000       0.400000000000000     
Na2_sspg: aquest potencial no està definit, de part de: Selec_Pot
