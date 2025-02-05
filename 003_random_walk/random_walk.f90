program random_walk_program

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  integer(i4), parameter :: Max_steps = 1000, N_walkers = 1000, N_experiments = 1000
  integer(i4) :: Max_steps_array(500)
  real(dp), parameter :: p = 0.5_dp

  integer(i4) :: walker(N_walkers)
  real(dp) :: probability_of_return(N_experiments)
  integer(i4) :: i_steps, i_walker, i_experiments,counter,k,i, ounit
  
  real(dp) :: r

  max_steps_array = [(i*10,i=1,size(Max_steps_array))]
  !max_steps_array(31:50) = [(i*100,i=31,50)]
  open(newunit = ounit, file = 'probabilities.dat')

  max_step: do k = 1, size(max_steps_array)
     experiments: do i_experiments = 1, N_experiments
        counter = 0
        walker = 0
        walkers: do i_walker = 1, N_walkers
           steps: do i_steps = 1, Max_steps_array(k)
              call random_number(r)
              if( r <= p )then
                 walker(i_walker) = walker(i_walker) + 1
              else
                 walker(i_walker) = walker(i_walker) - 1
              end if

              if(walker(i_walker) == 0)then
                 counter = counter + 1
                 cycle walkers
              end if
           end do steps
        end do walkers
        probability_of_return(i_experiments) = counter/real(N_walkers,dp)
        
     end do experiments
     write(ounit,*) Max_steps_array(k),avr(probability_of_return),&
          jackknife(probability_of_return,10), stderr(probability_of_return)
     flush(ounit)
  end do max_step
contains

  function avr(x)
    real(dp) :: avr, x(:)
    avr = sum(x)/size(x)
  end function avr


  function var(x)
    real(dp) :: var, x(:), xbar
    xbar = avr(x)
    var = 1.0_dp/(size(x) - 1) * sum( (x - xbar)**2 )
  end function var

  function stderr(x)
    real(dp) :: stderr, x(:)
    stderr = sqrt(var(x)/size(x))
  end function stderr


  function jackknife(x,bins)
    real(dp) :: jackknife, x(:)
    integer(i4), intent(in) :: bins
    integer(i4) :: MM, NN, i
    real(dp) :: xbar, sum_x
    real(dp) :: x_m(bins)


    NN = size(x)
    MM = NN/bins

    xbar = avr(x)
    sum_x = sum(x)
    x_m = 1.0_dp/(NN-MM) * [(sum_x - sum(x(MM*(i-1)+1:MM*i)),i=1,bins)]

    jackknife = sqrt( real(bins - 1,dp)/bins * sum( (x_m - xbar)**2) )
  end function jackknife
  
end program random_walk_program
