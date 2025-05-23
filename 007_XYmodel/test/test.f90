program test

  use id_cluster
  use pbc
  implicit none

  integer, parameter :: L = 50
  integer, dimension(L) :: bonds, cluster
  real :: r(L)

  call random_number(r)
  
  call set_pbc(L)
  
  bonds = nint(r)![0,1,0,1,1,0,0,0,1,1]
  call identify_clusters(bonds,cluster)
  print"(i4,*(i3))", bonds
  print"(*(i3))", cluster

  
  !bonds = 1
  !call identify_clusters(bonds,cluster)
  !print*, cluster
  
end program test
