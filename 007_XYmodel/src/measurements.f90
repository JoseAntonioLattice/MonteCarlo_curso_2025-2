module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use constants
  use pbc, only : im, ip
  implicit none
  

contains

  function topological_charge(spin)
    real(dp) :: topological_charge
    real(dp), dimension(:,:), intent(in) :: spin
    real(dp), dimension(size(spin(1,:))) :: phi, dphi
    integer(i4) :: L, i, n

    L = size(spin(1,:))
    
    do i = 1, L
       phi(i) = atan2(spin(2,i),spin(1,i))
    end do

    do i = 1, L
       dphi(i) = phi(ip(i)) - phi(i)
       if(dphi(i) >=  0.0_dp) then
          dphi(i) = mod(dphi(i),2*pi)
          if(dphi(i) > pi) dphi(i) = dphi(i) - 2*pi
       elseif(dphi(i) < 0.0_dp) then
          n = abs(dphi(i)/(2*pi))
          dphi(i) = dphi(i) + n*2*pi
          if(dphi(i) < -pi) dphi(i) = dphi(i) + 2*pi
       end if
    end do

    topological_charge = 0.0_dp
    do i = 1, L
       topological_charge = topological_charge + dphi(i)
    end do
    topological_charge = topological_charge/(2*pi)
    
  end function topological_charge

  function action(spin,beta)
    real(dp) :: action
    real(dp), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: beta

    action = beta*hamiltonian(spin)
    
  end function action

  function hamiltonian(spin)
    real(dp) :: hamiltonian
    real(dp), dimension(:,:), intent(inout) :: spin
    integer(i4) :: i

    hamiltonian = 0.0_dp
    do i = 1, size(spin(1,:))
       hamiltonian = hamiltonian - dot_product(spin(:,i),spin(:,ip(i)))
    end do  
    
  end function hamiltonian


  function DH(spin,spin_p,i)
    real(dp) :: DH
    real(dp), dimension(:,:), intent(in) :: spin
    integer(i4), intent(in) :: i
    real(dp), dimension(2), intent(in) :: spin_p
    real(dp) :: r
        
    DH = dot_product(spin(:,ip(i))+spin(:,im(i)), spin(:,i) - spin_p)
    
  end function DH
  
  function correlation(spin)
    real(dp), dimension(:,:), intent(in) :: spin
    real(dp), dimension(0:size(spin(1,:))) :: correlation
    integer(i4) :: i
    
    do i = 0, size(spin(1,:))-1
       correlation(i) = dot_product(spin(:,1),spin(:,i+1))
    end do
    correlation(size(spin(1,:))) = dot_product(spin(:,1),spin(:,1))

  end function correlation
  
end module measurements
