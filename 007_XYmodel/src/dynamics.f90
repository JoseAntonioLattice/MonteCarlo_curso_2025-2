module dynamics

  use precision
  use pbc, only : ip, im
  use random
  use constants
  use measurements
  use update_algorithms
  implicit none
  
contains
  
  subroutine hot_start(spin)
    real(dp), dimension(:,:), intent(out) :: spin
    real(dp) :: phi(size(spin(1,:)))
    integer(i4) :: i
    call random_number(phi)

    phi = 2*pi*(2*phi-1.0_dp)
    do i = 1, size(spin(1,:))
       spin(:,i) = [cos(phi(i)), sin(phi(i))] 
    end do
    
  end subroutine hot_start

  subroutine cold_start(spin)
    real(dp), dimension(:,:), intent(out) :: spin

    spin(1,:) = 1.0_dp
    spin(2,:) = 0.0_dp
  end subroutine cold_start

  subroutine start_configuration(spin, start)
    real(dp), dimension(:,:), intent(out) :: spin
    character(*), intent(in) :: start

    select case(start)
    case("hot")
       call hot_start(spin)
    case("cold")
       call cold_start(spin)
    end select

  end subroutine start_configuration
  
  subroutine thermalization(spin,N_thermalization,algorithm, beta)
    real(dp), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_thermalization
    integer(i4) :: i_th
    real(dp) :: h

    do i_th = 1, N_thermalization
       call sweeps(spin,beta,algorithm)
       !if( mod(i_th, 100) == 0) call normalize(spin) 
    end do
    
  end subroutine thermalization

  subroutine take_measurements(spin,N_measurements,N_skip,algorithm,beta)
    use spin_field, only : topological_charge_array, energy_density_array, correlation_array
    real(dp), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_measurements, N_skip
    integer(i4) :: i_m, i_skip
   
    do i_m = 1, N_measurements
       do i_skip = 1, N_skip
          call sweeps(spin,beta,algorithm)
       end do
       !call normalize(spin) 
       topological_charge_array(i_m) = topological_charge(spin)
       energy_density_array(i_m) = hamiltonian(spin)/size(spin(1,:))
       correlation_array(0:,i_m) = correlation(spin)
    end do
    
    
  end subroutine take_measurements
 
  
  subroutine sweeps(spin,beta,algorithm)
    real(dp), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer :: i

    select case(algorithm)
    case("metropolis")
       do i = 1, size(spin(1,:))
          call metropolis(spin,beta,i)
       end do
    case("wolff")
       call wolff(spin,beta)
    case default
       stop "Not a valid algorithm. Try 'metropolis' or 'wolff'."
    end select
    
  end subroutine 


  subroutine mc_simulation(spin,start,algorithm,beta,N_thermalization,N_measurements,N_skip)

    use spin_field, only : topological_charge_array, energy_density_array, correlation_array
    use statistics
    
    real(dp), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: beta(:)
    character(*), intent(in) :: start, algorithm
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    integer(i4) :: i_b, i, outcorrelation 
    
    call start_configuration(spin,trim(start))

    open(newunit = outcorrelation, file = "data/correlation.dat", status = "unknown", action = "write")
    do i_b = 1, size(beta)
       call thermalization(spin,N_thermalization,trim(algorithm),beta(i_b))
       call take_measurements(spin,N_measurements,N_skip,algorithm,beta(i_b))

       write(outcorrelation,*) "#",beta(i_b)
       do i = 0, size(spin(1,:))
          write(outcorrelation,*) avr(correlation_array(i,:)), std_err(correlation_array(i,:))
       end do
       write(outcorrelation,*) " "
       write(outcorrelation,*) " "
       write(outcorrelation,*) " "
       
       print*, beta(i_b), avr(energy_density_array), std_err(energy_density_array), &
            avr(topological_charge_array), std_err(topological_charge_array), &
            avr(topological_charge_array**2)/size(spin(1,:)), std_err(topological_charge_array**2)/size(spin(1,:))
    end do
  end subroutine mc_simulation

  subroutine normalize(spin)
    real(dp), dimension(:,:), intent(inout) :: spin
    integer(i4) :: i

    do i = 1, size(spin(1,:))
       spin(:,i) = spin(:,i)/sqrt(norm2(spin(:,i)))
    end do
  end subroutine normalize
  
end module dynamics
