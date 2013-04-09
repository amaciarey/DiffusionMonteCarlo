module pbc_mod

use global_mod, only: dim,Np,Lbox,LboxHalf

implicit none

contains

!-----------------------------------------------------------------------

  subroutine MinimumImage(xij,rij2)

    implicit none

    real (kind=8)    :: rij2
    integer (kind=4) :: k

    real (kind=8),dimension (dim) :: xij

    rij2 = 0.d0
    
    do k=1,dim
       
       if (xij(k)> LboxHalf(k)) xij(k) = xij(k)-Lbox(k)
       if (xij(k)<-LboxHalf(k)) xij(k) = xij(k)+Lbox(k)
       
       rij2 = rij2+xij(k)*xij(k)
       
    end do
    
    return
  end subroutine MinimumImage

!-----------------------------------------------------------------------

  subroutine BoundaryConditions(k,x)

    implicit none

    real (kind=8)    :: x
    integer (kind=4) :: k

    if (x<0.d0   ) x = x+Lbox(k)
    if (x>Lbox(k)) x = x-Lbox(k)

    return
  end subroutine BoundaryConditions

!-----------------------------------------------------------------------

end module pbc_mod
