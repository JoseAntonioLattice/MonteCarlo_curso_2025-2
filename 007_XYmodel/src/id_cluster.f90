module id_cluster

  
  use precision
  use pbc, only : ip, im
  implicit none

contains
  
  subroutine identify_clusters(bond,clusters)
    integer(i4), intent(in) :: bond(:)
    integer(i4), intent(out) :: clusters(:)
    integer(i4) ::  i, label

    label = 1
    clusters = 0
    clusters(1) = label 
    do i = size(bond), 1, -1
       if(bond(i) == 1) then
          clusters(i) = label
       else
          exit
       end if
    end do

    do i = 1, size(bond) - 1
       if(clusters(ip(i)) == 0) then
          if(bond(i) == 1) then
             clusters(ip(i)) = label
          else
             label = label + 1
             clusters(ip(i)) = label
          end if
       else
          exit
       end if
    end do
    
  end subroutine identify_clusters

  
end module id_cluster
