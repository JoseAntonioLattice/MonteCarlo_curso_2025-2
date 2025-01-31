program buffon

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  real(dp) :: y, theta
  real(dp), parameter :: D = 1.0_dp, L = 1.0_dp 
  
  integer(i4) :: i,j,counter
  integer(i4), parameter :: N = 10**5, M = 10**3

  real(dp), parameter :: pi = acos(-1.0_dp)

  real(dp) :: pi_approximation(M)

  ! Do M experiments
  stats: do j = 1, M
     counter = 0 
     MC: do i = 1, N
        !Throw needle
        call random_number(y)
        call random_number(theta)
        theta = pi*theta
        ! Check if needle crosses horizontal lines
        if (y + 0.5*L*sin(theta) > D .or. y - 0.5*L*sin(theta) < 0.0_dp) counter = counter + 1
     end do MC
     
     pi_approximation(j) =  2*real(N,dp)/counter
     !print*, pi_approximation(j)
  end do stats

  print*, avr(pi_approximation), stderr(pi_approximation), &
       jackknife(pi_approximation,5),jackknife(pi_approximation,10),jackknife(pi_approximation,20)

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
  
end program buffon
