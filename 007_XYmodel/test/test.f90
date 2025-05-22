program test

  use id_cluster
  use pbc
  implicit none

  integer, parameter :: L = 10
  integer, dimension(L) :: bonds, cluster

  call set_pbc(L)
  
  bonds = [0,1,0,1,1,0,0,0,1,1]
  call identify_clusters(bonds,cluster)
  print*, cluster

  
  bonds = 1
  call identify_clusters(bonds,cluster)
  print*, cluster
  
end program test
