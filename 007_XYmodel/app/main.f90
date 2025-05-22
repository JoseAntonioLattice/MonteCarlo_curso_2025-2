program xymodel
  
  use parameters
  use pbc
  use spin_field
  use dynamics
  implicit none

  integer(i4), parameter :: nb = 20
  real(dp) :: beta(nb), bi, bf, db
  integer(i4) :: i_b

  bf = 10.0_dp
  bi = 0.1_dp
  db = (bf-bi)/nb
  
  ! Read input parameters file
  call read_input()

  ! Allocate arrays
  allocate(spin(2,L))
  allocate(energy_density_array(N_measurements))
  allocate(action_array(N_measurements))
  allocate(topological_charge_array(N_measurements))
  allocate(magnetization_array(N_measurements))
  allocate(correlation_array(0:L,N_measurements))

  ! Set periodic boundary conditions
  call set_pbc(L)

  beta = [(bi + (i_b-1)*(bf-bi)/(nb-1), i_b = 1, nb)]
  
  ! Start simulation
  call mc_simulation(spin,trim(start),trim(algorithm),beta,N_thermalization,N_measurements,N_skip)
  
     
end program xymodel
