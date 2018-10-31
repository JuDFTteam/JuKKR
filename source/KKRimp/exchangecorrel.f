!-------------------------------------------------------------------------------
!> Summary: Module for the exchange correlation potentials
!> Author: 
!> This module gathers the different exchange correlation potentials
!-------------------------------------------------------------------------------
        MODULE MOD_EXCHANGECORRELATION
        CONTAINS    
!           include 'XC/vxcdrvnew.f90'

          include 'XC/corlsd.f'
          include 'XC/cpw91.f'
          include 'XC/cylm02.f'
          include 'XC/exch91.f'
          include 'XC/gcor91.f'
          include 'XC/gradr.f'
          include 'XC/gradrl.f'
          include 'XC/gxcpt.f'
          include 'XC/lebedev.f'
          include 'XC/mkxcpe.f' 
          include 'XC/sphere_gga.f'
          include 'XC/sphere_nogga.f'
          include 'XC/spher.f'
          include 'XC/trarea.f'
          include 'XC/vosko.f'
          include 'XC/vxcgga.f'
          include 'XC/vxclm.f'    
          include 'XC/vxcspo.f'  
          include 'XC/rinit.f'
        END MODULE MOD_EXCHANGECORRELATION 


