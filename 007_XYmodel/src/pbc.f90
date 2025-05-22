module pbc

  use precision
  implicit none

  integer(i4), allocatable, dimension(:) :: ip, im

  
contains
  
  subroutine set_pbc(L)
    integer(i4), intent(in) :: L
    integer(i4) :: i

    allocate(ip(L),im(L))

    do i = 1, L
       ip(i) = i + 1
       im(i) = i - 1
    end do
    ip(L) = 1
    im(1) = L
  end subroutine set_pbc
    
end module pbc
