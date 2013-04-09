module global_mod

  logical          :: table=.True.
  real (kind=8)    :: V0
  real (kind=8)    :: pi
  real (kind=8)    :: Rm,Amp
  real (kind=8)    :: rbin,dr
  real (kind=8)    :: rcut,rcut2
  integer (kind=4) :: dim,Np,Nbin,Nmax

  real (kind=8),dimension (:),allocatable   :: Lbox,LboxHalf,qbin
  real (kind=8),dimension (:,:),allocatable :: Req
  

end module global_mod
