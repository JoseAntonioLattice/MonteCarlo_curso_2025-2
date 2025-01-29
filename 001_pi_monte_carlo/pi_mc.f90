program pi_mc

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  integer(i4), parameter :: n = 10**2, m = 10**5
  
  real(dp) :: x, y, pi(n)
  integer(i4) :: i,j, counter


  main: do i = 1, n

     counter = 0
     pi_loop : do j = 1, m

        call random_number(x)
        call random_number(y)

        if (x**2 + y**2 <= 1.0_dp) counter = counter + 1
        
     end do pi_loop
     pi(i) = 4*real(counter,dp)/m
     !print*, pi(i)
  end do main


  print'("pi = ",f7.5," +/- ",f7.5)', avr(pi), stderr(pi)

contains

  function avr(x)
    real(dp) :: avr, x(:)

    avr = sum(x)/size(x)
  end function avr

  function var(x)
    real(dp) :: var, x(:)

    var = sum((x - avr(x))**2)/(size(x) - 1)
  end function var

  function stderr(x)
    real(dp) :: stderr, x(:)

    stderr = sqrt( var(x)/size(x) )
  end function stderr
  
end program pi_mc
