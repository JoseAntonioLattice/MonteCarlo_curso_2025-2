program tarea

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  integer(i4) :: inunit, i, k
  real(dp) :: r(100)


  open(newunit = inunit,file = "datos.dat")

  do i = 1, 100
     read(inunit,*) k, r(i)
  end do

   print'("<r> = ",f10.5," +/- ",f7.5)', avr(r), stderr(r)

  
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
  
end program tarea
