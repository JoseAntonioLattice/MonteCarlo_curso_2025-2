module update_algorithms

  use precision
  use constants
  use random
  use measurements, only : DH
  use id_cluster
  implicit none

contains

  subroutine metropolis(spin,beta,i)
    real(dp), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: beta
    integer(i4), intent(in) :: i
    real(dp) :: r, u, p, spin_p(2)
    
    call random_number(r)
   
    u = uran(0.0_dp,2*pi)
    spin_p = [cos(u),sin(u)] 
    p = min(1.0_dp,exp(-beta*DH(spin,spin_p,i)))
    
    if(r < p) spin(:,i) = spin_p
    
  end subroutine metropolis

  subroutine wolff(spin,beta)
    real(dp), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: beta
    integer(i4) :: clusters(size(spin(1,:))), bond(size(spin(1,:)))

    real(dp) :: r(2), phi

    ! Select Wolff direction
    phi = uran(0.0_dp,2*pi)
    r = [cos(phi),sin(phi)]

    
    call put_bonds(bond,spin,r,beta) 
    call identify_clusters(bond,clusters)
    call flip_clusters(spin,clusters,r)
    
  end subroutine wolff

  subroutine flip_clusters(spin,clusters,r)
    real(dp), dimension(:,:), intent(inout) :: spin
    integer(i4), intent(in) :: clusters(:)
    real(dp), intent(in) :: r(2)
    integer :: i
    real(dp) :: u
    real(dp), allocatable :: p(:)


    allocate(p(maxval(clusters)))
    call random_number(p)  
    do i = 1, size(spin(1,:))
       if(p(clusters(i)) < 0.5_dp) spin(:,i) = flip(spin(:,i),r)
    end do
    deallocate(p)
  end subroutine flip_clusters

  function flip(spin,r)
    real(dp), dimension(2) :: flip
    real(dp), dimension(2), intent(in) :: spin, r

    flip = spin - 2*r*dot_product(spin,r)
    
  end function flip

  subroutine put_bonds(bond,spin,r,beta)
    real(dp), dimension(:,:), intent(in) :: spin
    integer(i4), intent(out) :: bond(size(spin(2,:)))
    real(dp), intent(in) :: r(2), beta
    integer(i4) :: i
    real(dp) :: DHx, u
    
    bond = 0
    do i = 1, size(bond)
       DHx = dot_product(-flip(spin(:,i),r)+spin(:,i),spin(:,ip(i))) !+ dot_product(flip(spin(:,i),r),spin(:,ip(i)))
       if( DHx > 0.0_dp )then
          call random_number(u)
          if(u <= 1.0_dp - exp(-beta*DHx)) bond(i) = 1
       end if
    end do
    
  end subroutine put_bonds
end module update_algorithms
