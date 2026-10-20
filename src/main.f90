program main_drv
  use types      
  use global_variables
  use input_values
  use math_subroutines
  use leg_integrator
  use dmc_cal
  implicit none

   call reading_data 
   if (.not.fixran) call seed_cal(idum)

   sim_type: select case (simtyp)
            case (1:2)
             call dmc_drv
            case (3)
             write(*,*) 'DVR not implemented yet'
            case (4)
             write(*,*) 'Adiabatic diagonalization not implemented yet'
            case default
             write(*,*) 'Simulation type not recognized'
            end select   sim_type






end program main_drv



