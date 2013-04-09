module dmc_mod

use global_mod
use random_mod
use sample_mod
use system_mod
use pbc_mod

implicit none

contains

!-----------------------------------------------------------------------
!
! This subroutine generates a table of values of the two body Jastrow
! function in the form WF(i) = Log(f(r_i)) in order to speed up the
! calculation.
!-----------------------------------------------------------------------

  subroutine Jastrow_Table(rmax,Rm,WF)

    implicit none

    real (kind=8)    :: Rm,rmax
    real (kind=8)    :: r
    integer (kind=4) :: i

    real (kind=8),dimension (0:Nmax+1) :: WF

    dr = rmax/real(Nmax-1)

    open (unit=1,file='jastrow.out')

    do i=1,Nmax

       r     = (i-1)*dr
       WF(i) = LogPsi(0,Rm,r)
       write (1,'(20g20.10e3)') r,exp(WF(i))
 
    end do

    close (unit=1)

    WF(0)      = WF(2)
    WF(Nmax+1) = WF(Nmax)

    return
  end subroutine Jastrow_Table
	
!-----------------------------------------------------------------------  
!
! This subroutine defines the initial set of walkers of the simulation.
! The position, drift force and local energy for each walker is stored
! in variables Walker, DriftForce and Eloc respectively. The initial 
! set of walkers can be set randomly (opt=0) or alternatively read from
! a previous calculation (opt=1).
! The energy shift E0 is initially defined as the mean energy per walker
! in order to avoid problems with the population growth.
!
!-----------------------------------------------------------------------
  
  subroutine init(LogWF,NwStart,NwMax,opt,Walker,DriftForce,Eloc,E0)

    implicit none

    real (kind=8)    :: E,Kin,Kf,Pot,E0
    integer (kind=4) :: ip
    integer (kind=4) :: NwStart,NwMax,iw,opt
    integer (kind=4) :: k
    
    real (kind=8),dimension(0:Nmax+1)       :: LogWF
    real (kind=8),dimension(dim,Np,NwMax,2) :: Walker,DriftForce
    real (kind=8),dimension(NwMax,2)        :: Eloc
    real (kind=8),dimension(dim,Np)         :: R,F
    real (kind=8),dimension(Nbin)           :: gr,rho,nrho

    if (opt==0) then

       !Random generation

       E0 = 0.d0

       do iw=1,NwStart
   
          do ip=1,Np
             do k=1,dim
                R(k,ip) = Lbox(k)*grnd()
             end do
          end do
          
          call Forces(.true.,LogWf,R,E,Kin,Kf,Pot,F,gr,rho,nrho)
   
          do ip=1,Np
             do k=1,dim
                Walker(k,ip,iw,1)     = R(k,ip)
                DriftForce(k,ip,iw,1) = F(k,ip)
             end do
          end do

          Eloc(iw,1) = E
          E0         = E0+E
          
       end do

       E0 = E0/real(NwStart)

    else

       !Reading a configuration from a VMC or another many body 
       !calculation
       
       open (unit=2,file='config_vmc.in',status='old')
       
       E0 = 0.d0
       
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*) NwStart

       do iw=1,NwStart
      
          read (2,*) E
     
          Eloc(iw,1) = E
          E0         = E0+E
   
          do ip=1,Np
         
             read (2,'(30g20.10e3)') (R(k,ip),k=1,dim)
                     
          end do
          
          call Forces(.false.,LogWf,R,E,Kin,Kf,Pot,F,gr,rho,nrho)
          
          do ip=1,Np
             do k=1,dim
                Walker(k,ip,iw,1)     = R(k,ip)
                DriftForce(k,ip,iw,1) = F(k,ip)
             end do
          end do
          
       end do

       E0 = E0/real(NwStart)
       
       close (unit=2)

    end if

    return
  end subroutine init

!-----------------------------------------------------------------------
!
! This subroutine performs the gaussian diffussion for a single walker.
!
!-----------------------------------------------------------------------

  subroutine GaussianDiffusion(sigma,mu,R)

    implicit none

    real (kind=8)    :: sigma,mu
    real (kind=8)    :: gauss1,gauss2
    integer (kind=4) :: k,ip

    real (kind=8),dimension (dim,Np) :: R

    do ip=1,Np/2

       do k=1,dim
          
          call rangauss(sigma,mu,gauss1,gauss2)

          R(k,2*ip-1) = R(k,2*ip-1)+gauss1 
          R(k,2*ip)   = R(k,2*ip)+gauss2

          call BoundaryConditions(k,R(k,2*ip-1))
          call BoundaryConditions(k,R(k,2*ip))
               
       end do

    end do

    return
  end subroutine GaussianDiffusion

!-----------------------------------------------------------------------
!
! This soubroutine performs the drift movement for a single walker.
!
!-----------------------------------------------------------------------

  subroutine DriftIntegration(LogWF,R,F,dt)

    implicit none

    real (kind=8)    :: dt
    real (kind=8)    :: E,Kin,Kf,Pot
    integer (kind=4) :: ip
    integer (kind=4) :: k
        
    real (kind=8),dimension (0:Nmax+1) :: LogWF
    real (kind=8),dimension (dim,Np)   :: R,Rp,F,Fp
    real (kind=8),dimension (Nbin)     :: gr,rho,nrho
    
    !Predictor step

    do ip=1,Np
            
       do k=1,dim

          Rp(k,ip) = R(k,ip)+dt*F(k,ip)

          call BoundaryConditions(k,Rp(k,ip))
               
       end do

    end do

    !Force evaluation
    
    call Forces(.false.,LogWf,R,E,Kin,Kf,Pot,Fp,gr,rho,nrho)

    !Corrector step

    do ip=1,Np

       do k=1,dim

          R(k,ip) = Rp(k,ip)+0.5d0*dt*(Fp(k,ip)-F(k,ip))

          call BoundaryConditions(k,R(k,ip))

       end do

    end do

    return
  end subroutine DriftIntegration

!-----------------------------------------------------------------------

end module dmc_mod
