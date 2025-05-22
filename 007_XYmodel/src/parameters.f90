module parameters

  use precision
  implicit none

  integer :: L                ! Lattice size
  integer :: N_measurements   ! # of measurements
  integer :: N_skip           ! sweeps between measurements
  integer :: N_thermalization ! # of sweeps of thermalization
  character(100) :: algorithm ! Algorithm: Metropolis
  character(100) :: start     ! Hot or cold start

  namelist /input_parameters/ L, N_measurements, N_skip, N_thermalization, &
       algorithm, start
  
contains

  subroutine read_input()

    character(100):: input_file
    integer(i4) :: inunit
    
    write(*,*) "Please enter the parameters file: "
    read(*,'(a)') input_file
    write(*,'(a)') trim(input_file)

    open( newunit = inunit, file = trim(input_file), status = 'old')
    read(inunit, nml = input_parameters)
    close(inunit)
    write(*,nml=input_parameters)
    
  end subroutine read_input
 

end module parameters
