module spin_field

  use precision
  implicit none

  ! configuration
  real(dp), allocatable, dimension(:,:) :: spin
  ! observables
  real(dp), allocatable, dimension(:) :: action_array, energy_density_array, &
       magnetization_array, topological_charge_array, &
       correlation_array(:,:)
  
end module spin_field
