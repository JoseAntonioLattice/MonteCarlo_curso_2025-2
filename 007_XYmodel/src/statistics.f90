module statistics

  use precision
  implicit none

contains

  function avr(x)
    real(dp) :: avr
    real(dp), intent(in) :: x(:)

    avr = sum(x)/size(x)
    
  end function avr

  function var(x)
    real(dp) :: var
    real(dp), intent(in) :: x(:)

    var = sum((x-avr(x))**2)/(size(x)-1)
    
  end function var

  function std_err(x)
    real(dp) :: std_err
    real(dp), intent(in) :: x(:)

    std_err = sqrt(var(x)/size(x))

  end function std_err

end module statistics
