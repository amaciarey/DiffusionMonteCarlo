module sample_mod

use global_mod
use random_mod
use system_mod
use pbc_mod

implicit none

contains

!-----------------------------------------------------------------------
!
! This subroutine performs the evaluation of the drift force for a 
! single walker that is needed for the drift movement and the energy
! evaluation.
! If the control variable samp is set to .TRUE. then the routine 
! evaluates also the Energies, pair correlation function and one body
! density matrix.
!
!-----------------------------------------------------------------------

  subroutine Forces(samp,LogWf,R,E,Kin,Kf,Pot,F,gr,rho,nrho)

    implicit none

    logical          :: samp
    real (kind=8)    :: E,Kin,Kf,Pot
    real (kind=8)    :: Interpolate
    real (kind=8)    :: rij,rij2,fij
    real (kind=8)    :: ur,dudr,d2udr2
    real (kind=8)    :: LapLogPsi
    real (kind=8)    :: sum_ur
        
    integer (kind=4) :: ip,jp
    integer (kind=4) :: k
    integer (kind=4) :: ibin
    
    real (kind=8),dimension (0:Nmax+1) :: LogWF
    real (kind=8),dimension (dim,Np)   :: R,F
    real (kind=8),dimension (Nbin)     :: gr,rho,nrho
    real (kind=8),dimension (dim)      :: xij,xmc
    real (kind=8),dimension (Np)       :: aux1rho,aux2rho

    Pot       = 0.d0
    LapLogPsi = 0.d0

    do ip=1,Np
       do k=1,dim
          F(k,ip) = 0.d0
       end do
       aux1rho(ip) = 0.d0
    end do

    do ip=1,Np-1
       do jp=ip+1,Np
      
          !Calculation of the distance between each pair of dipoles
          !choosing only the nearest image of each particle
          
          do k=1,dim

             xij(k) = R(k,ip)-R(k,jp)

          end do

          call MinimumImage(xij,rij2)
          
          !Calculation of the local kinetic energy
          !We use in this routine the following prescription:
          
          !    K = -0.5*Laplacian(Psi)/Psi
          !    K = -0.5*(Laplacian(Log(Psi))+F*F)
          
          !where F is the so called drift force given by:
          
          !    F = Grad(Log(Psi))

          if (rij2<=rcut2) then
         
             rij = sqrt(rij2)

             if (table) then
                
                dudr = Interpolate(1,Nmax,dr,LogWF,rij)
              
             else

                dudr = LogPsi(1,Rm,rij)

             end if

             do k=1,dim
                fij     = dudr*xij(k)/rij
                F(k,ip) = F(k,ip)+fij
                F(k,jp) = F(k,jp)-fij
             end do      

             if (samp) then
             
                !Density matrix
                
                if (table) then
                
                   ur     = Interpolate(0,Nmax,dr,LogWF,rij)
                   d2udr2 = Interpolate(2,Nmax,dr,LogWF,rij) 
                   
                else
                
                   ur     = LogPsi(0,Rm,rij)
                   d2udr2 = LogPsi(2,Rm,rij)
                
                end if
             
                aux1rho(ip) = aux1rho(ip)+ur
                aux1rho(jp) = aux1rho(jp)+ur
             
                !Kinetic energy
             
                LapLogPsi = LapLogPsi+(real(dim-1)*dudr/rij+d2udr2)   
             
                !Calculation of the local potential energy 
             
                Pot = Pot+Potential(rij)
             
                !Construction of the histogram for g(x,y)
             
                ibin     = int(rij/rbin)+1
                gr(ibin) = gr(ibin)+2.d0
                
             end if

          end if

       end do
    end do

    !Kinetic energy and drift force calculation

    !Density matrix evaluation

    if (samp) then
       
       sum_ur = 0.d0

       do k=1,dim
          xmc(k) = Lbox(k)*grnd()
       end do
       
       do ip=1,Np
          
          aux2rho(ip) = 0.d0
          
          do k=1,dim
             
             xij(k) = xmc(k)-R(k,ip)
             
          end do
          
          call MinimumImage(xij,rij2)
          
          if (rij2<=rcut2) then
             
             rij = sqrt(rij2)
             
             if (table) then
                
                ur = Interpolate(0,Nmax,dr,LogWF,rij)
                
             else
                
                ur = LogPsi(0,Rm,rij)
                
             end if
             
             aux2rho(ip) = aux2rho(ip)+ur
             sum_ur      = sum_ur+ur
             
          end if
          
       end do
       
       do ip=1,Np
          
          do k=1,dim
             
             xij(k) = xmc(k)-R(k,ip)
             
          end do
          
          call MinimumImage(xij,rij2)
          
          if (rij2<=rcut2) then
             
             rij  = sqrt(rij2)
             ibin = int(rij/rbin)+1
             
             rho(ibin)  = rho(ibin)+exp(sum_ur-aux2rho(ip)-aux1rho(ip))
             nrho(ibin) = nrho(ibin)+1.d0
             
          end if
          
       end do
       
       Kin = 2.d0*LapLogPsi
       Kf  = 0.d0
       
       do ip=1,Np
          do k=1,dim
             Kf = Kf+F(k,ip)*F(k,ip)
          end do
       end do
       
       !Computing energy
       
       Kin = -0.5d0*(Kin+Kf)
       Kf  = 0.5d0*Kf
       E   = Kin+Pot
       
    end if

    return
  end subroutine Forces

!-----------------------------------------------------------------------
!
! This subroutine performs the evaluation of the static structure factor
! by using the identity:
!
! S(k) = \sum_{i}^{Np} exp(ii*k*r_i) 
!
! where {k} is a wave vector compatible with the periodic boundary 
! conditions and r_i is the position vector of a particle in a single 
! walker.
!
!-----------------------------------------------------------------------

  subroutine StructureFactor(Nk,R,Sk)

    implicit none

    real (kind=8)    :: qr
    real (kind=8)    :: SumCos,SumSin
    integer (kind=4) :: iq,Nk
    integer (kind=4) :: k
    integer (kind=4) :: ip

    real (kind=8),dimension(dim)    :: x
    real (kind=8),dimension(dim,Np) :: R
    real (kind=8),dimension(dim,Nk) :: Sk

    SumCos = 0.d0
    SumSin = 0.d0

    do iq=1,Nk

       do k=1,dim
   
          SumCos = 0.d0
          SumSin = 0.d0

          do ip=1,Np
        
             x(k) = R(k,ip)
             qr   = real(iq)*qbin(k)*x(k)
             
             SumCos = SumCos+cos(qr)
             SumSin = SumSin+sin(qr)
             
          end do
          
          Sk(k,iq) = Sk(k,iq)+(SumCos*SumCos+SumSin*SumSin)
          
       end do

    end do

    return
  end subroutine StructureFactor

!-----------------------------------------------------------------------
!
! This subroutine performs the normalization of the structural functions
! sampled during the MC simulation: the pair correlation function g(r),
! the static structure factor S(k) and the one body density matrix n(r).
!
!-----------------------------------------------------------------------

  subroutine Normalize(pur_est,density,Nk,ngr,gr,Sk,rho,nrho)

    implicit none

    real (kind=8)    :: r8_gamma
    real (kind=8)    :: density
    real (kind=8)    :: nid,r,norm
    real (kind=8)    :: k_n
    integer (kind=4) :: ibin,k
    integer (kind=4) :: ngr,Nk
    logical          :: pur_est
        
    real (kind=8),dimension (Nbin)   :: gr
    real (kind=8),dimension (dim,Nk) :: Sk
    real (kind=8),dimension (Nbin)   :: rho
    real (kind=8),dimension (Nbin)   :: nrho

    k_n  = pi**(0.5d0*dim)/r8_gamma(0.5d0*dim+1.d0)
    norm = real(Np)*real(ngr)
    
    if (pur_est) then
       open (unit=111,file='gr_pur.out',status='unknown')
    else
       open (unit=111,file='gr_dmc.out',status='unknown')
       open (unit=222,file='nr_dmc.out',status='unknown')
    end if
       
    do ibin=1,Nbin
       r   = (real(ibin)-0.5d0)*rbin
       nid = density*k_n*((r+0.5d0*rbin)**dim-(r-0.5d0*rbin)**dim)
       if (pur_est) then
          write (111,'(20g20.10e3)') r,gr(ibin)/(nid*norm)
       else
          write (111,'(20g20.10e3)') r,gr(ibin)/(nid*norm)
          write (222,'(20g20.10e3)') r,rho(ibin)/nrho(ibin)
          gr(ibin)  = gr(ibin)/(nid*norm)
          rho(ibin) = rho(ibin)/nrho(ibin) 
       end if
    end do

    close (unit=111)
    close (unit=222)

    if (pur_est) then
       open (unit=333,file='sk_pur.out',status='unknown')
    else
       open (unit=333,file='sk_dmc.out',status='unknown')
    end if

    do ibin=1,Nk
       if (pur_est) then
          write (333,'(20g20.10e3)') (ibin*qbin(k),Sk(k,ibin)/norm,k=1,dim)
       else
          write (333,'(20g20.10e3)') (ibin*qbin(k),Sk(k,ibin)/norm,k=1,dim)
          do k=1,dim
             Sk(k,ibin) = Sk(k,ibin)/norm
          end do
       end if
    end do
    
    close (unit=333)
    
    return
  end subroutine Normalize

!-----------------------------------------------------------------------

  subroutine AccumGr(gr,AvGr,AvGr2)

    implicit none
    
    integer (kind=4) :: j

    real (kind=8),dimension(Nbin) :: gr,AvGr,AvGr2

    do j=1,Nbin
       AvGr(j)  = AvGr(j)+gr(j)
       AvGr2(j) = AvGr2(j)+gr(j)*gr(j)
    end do

    return
  end subroutine AccumGr

!-----------------------------------------------------------------------

  subroutine AccumSk(Nk,Sk,AvSk,AvSk2)

    implicit none
    
    integer (kind=4) :: j,k,Nk

    real (kind=8),dimension(dim,Nk) :: Sk,AvSk,AvSk2

    do j=1,Nk
       do k=1,dim
          AvSk(k,j)  = AvSk(k,j)+Sk(k,j)
          AvSk2(k,j) = AvSk2(k,j)+Sk(k,j)*Sk(k,j)
       end do
    end do

    return
  end subroutine AccumSk

!-----------------------------------------------------------------------

  subroutine AccumNr(rho,AvRho,AvRho2)

    implicit none
    
    integer (kind=4) :: j

    real (kind=8),dimension(Nbin) :: rho,AvRho,AvRho2
    
    do j=1,Nbin
       AvRho(j)   = AvRho(j)+rho(j)
       AvRho2(j)  = AvRho2(j)+rho(j)*rho(j)
    end do

    return
  end subroutine AccumNr

!-----------------------------------------------------------------------

  subroutine NormAvGr(Nitem,AvGr,AvGr2,VarGr)

    implicit none

    real (kind=8)    :: r
    integer (kind=4) :: Nitem,j

    real (kind=8),dimension(Nbin) :: AvGr,AvGr2,VarGr
    
    open (unit=11,file='gr_dmc.out')

    do j=1,Nbin
       r        = (real(j)-0.5d0)*rbin
       AvGr(j)  = AvGr(j)/real(Nitem)
       AvGr2(j) = AvGr2(j)/real(Nitem)
       VarGr(j) = Var(Nitem,AvGr(j),AvGr2(j))
       write (11,'(20g20.10e3)') r,AvGr(j),VarGr(j)
    end do
    
    close (unit=11)

    return
  end subroutine NormAvGr

!-----------------------------------------------------------------------

  subroutine NormAvSk(Nitem,Nk,AvSk,AvSk2,VarSk)

    implicit none

    integer (kind=4) :: Nitem,Nk,j,k

    real (kind=8),dimension(dim,Nk) :: AvSk,AvSk2,VarSk
    
    open (unit=11,file='sk_dmc.out')

    do j=1,Nk
       do k=1,dim
          AvSk(k,j)  = AvSk(k,j)/real(Nitem)
          AvSk2(k,j) = AvSk2(k,j)/real(Nitem)
          VarSk(k,j) = Var(Nitem,AvSk(k,j),AvSk2(k,j))
       end do
       write (11,'(20g20.10e3)') (j*qbin(k),AvSk(k,j),VarSk(k,j),k=1,dim)
    end do

    close (unit=11)

    return
  end subroutine NormAvSk

!-----------------------------------------------------------------------

  subroutine NormAvNr(Nitem,AvRho,AvRho2,VarRho)

    implicit none

    real (kind=8)    :: r
    integer (kind=4) :: Nitem,j

    real (kind=8),dimension(Nbin) :: AvRho,AvRho2,VarRho
    
    open (unit=11,file='nr_dmc.out')

    do j=1,Nbin
       r         = (real(j)-0.5d0)*rbin
       AvRho(j)  = AvRho(j)/real(Nitem)
       AvRho2(j) = AvRho2(j)/real(Nitem)
       VarRho(j) = Var(Nitem,AvRho(j),AvRho2(j))
       write (11,'(20g20.10e3)') r,AvRho(j),VarRho(j)
    end do

    close (unit=11)

    return
  end subroutine NormAvNr

!-----------------------------------------------------------------------
!
! This subroutine simply accumulates the value of the energies at each 
! step/block in order to obtain average values.
!
!-----------------------------------------------------------------------

  subroutine Accumulate(E,Ec,Kf,Ep,SumE,SumEc,SumKf,SumEp)

    implicit none

    real (kind=8) :: E,SumE
    real (kind=8) :: Ec,SumEc
    real (kind=8) :: Kf,SumKf
    real (kind=8) :: Ep,SumEp
    
    SumE  = SumE+E
    SumEc = SumEc+Ec
    SumKf = SumKf+Kf
    SumEp = SumEp+Ep
    
    return
  end subroutine Accumulate
  
!-----------------------------------------------------------------------
!
! This subroutine evaluates the average of a magnitude along a step, 
! block or the whole simulation.
!
!-----------------------------------------------------------------------

  subroutine NormalizeAv(Nitem,SumE,SumEc,SumKf,SumEp)

    implicit none

    real (kind=8)    :: SumE,SumEc,SumKf,SumEp
    integer (kind=4) :: Nitem

    SumE  = SumE/real(Nitem)
    SumEc = SumEc/real(Nitem)
    SumKf = SumKf/real(Nitem)
    SumEp = SumEp/real(Nitem)

    return
  end subroutine NormalizeAv

!-----------------------------------------------------------------------
!
! This function simply evaluates the variance of a variable.
!
!-----------------------------------------------------------------------

  function Var(Nitem,Sum,Sum2)

    implicit none
    
    real (kind=8)    :: Var
    real (kind=8)    :: Sum,Sum2
    integer (kind=4) :: Nitem

    Var = sqrt((Sum2-Sum*Sum)/real(Nitem))

    return
  end function Var

!-----------------------------------------------------------------------

end module sample_mod
