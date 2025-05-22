module random

  use precision
  implicit none

contains

  function uran(a,b)
    real(dp) :: uran
    real(dp), intent(in) :: a, b
    real(dp) :: r

    call random_number(r)

    uran = r*(b - a) + a 
   
  end function uran
  
end module random
