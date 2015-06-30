!PROGRAM DMC

program dmc

use global_mod
use random_mod
use sample_mod
use dmc_mod

implicit none

real (kind=8)    :: density,E0
real (kind=8)    :: gamma
real (kind=8)    :: dt,sigma,mu
integer (kind=4) :: seed
integer (kind=4) :: opt
integer (kind=4) :: istep,Nstep
integer (kind=4) :: iblock,Nblock
integer (kind=4) :: iw,NwStart,NwEnd,NwMax,NwCrit
integer (kind=4) :: ip
integer (kind=4) :: k
integer (kind=4) :: isons,Nsons
integer (kind=4) :: ic,ngr,ngr_pur
integer (kind=4) :: ik,Nk

real (kind=8)    :: Ep_pur,BlockEp_pur
real (kind=8)    :: Kin,Kf,Pot,E
real (kind=8)    :: AvE,AvE2,VarE
real (kind=8)    :: AvEc,AvEc2,VarEc
real (kind=8)    :: AvKf,AvKf2,VarKf
real (kind=8)    :: AvEp,AvEp2,VarEp
real (kind=8)    :: BlockAvE,BlockAvE2,BlockVarE
real (kind=8)    :: BlockAvEc,BlockAvEc2,BlockVarEc
real (kind=8)    :: BlockAvKf,BlockAvKf2,BlockVarKf
real (kind=8)    :: BlockAvEp,BlockAvEp2,BlockVarEp
real (kind=8)    :: WalkerE,WalkerEc,WalkerKf,WalkerEp
real (kind=8)    :: end,begin
integer (kind=4) :: in,io                                               

real (kind=8),dimension (:),allocatable       :: LogWF
real (kind=8),dimension (:,:,:,:),allocatable :: Walker,DriftForce
real (kind=8),dimension (:,:),allocatable     :: R,F,Eloc
real (kind=8),dimension (:),allocatable       :: gr,gr_i,gr_pur
real (kind=8),dimension (:),allocatable       :: rho,rho_i
real (kind=8),dimension (:),allocatable       :: nrho,nrho_i
real (kind=8),dimension (:,:),allocatable     :: Sk,Sk_i,Sk_pur
real (kind=8),dimension (:,:),allocatable     :: AvSk,AvSk2,VarSk
real (kind=8),dimension (:),allocatable       :: AvNr,AvNr2,VarNr
real (kind=8),dimension (:),allocatable       :: AvGr,AvGr2,VarGr

real (kind=8),dimension (:,:),allocatable     :: Ep_w,Ep_rw
real (kind=8),dimension (:,:,:),allocatable   :: gr_w,gr_rw
real (kind=8),dimension (:,:,:,:),allocatable :: Sk_w,Sk_rw

!Reading input parameters

open (unit=1,file='dmc.in')

read (1,*) dim
read (1,*) Np
read (1,*) density
read (1,*) V0
read (1,*) Nmax
read (1,*) Nblock
read (1,*) Nstep
read (1,*) dt
read (1,*) gamma
read (1,*) NwStart
read (1,*) NwMax
read (1,*) NwCrit
read (1,*) Rm
read (1,*) Nbin
read (1,*) Nk
read (1,*) opt
read (1,*) seed

close (unit=1)

!Definition of useful parameters

allocate (Lbox(dim),LboxHalf(dim),qbin(dim))

pi = acos(-1.d0)

if (opt==1) then
   
   open (unit=11,file='config_vmc.in',status='old')
   read (11,*) Np
   read (11,*) (Lbox(k),k=1,dim)
   read (11,*) density
   close (unit=11)

else

   Amp = 0.d0

   do k=1,dim
   
      Lbox(k) = (real(Np)/density)**(1.d0/dim)
            
   end do
   
end if

do k=1,dim

   LboxHalf(k) = 0.5d0*Lbox(k)
   qbin(k)     = 2.d0*pi/Lbox(k)

end do

rcut  = minval(LboxHalf)
rcut2 = rcut*rcut
rbin  = rcut/real(Nbin)
sigma = sqrt(dt)
mu    = 0.d0
in    = 1
io    = 2

allocate (Walker(dim,Np,NwMax,2),DriftForce(dim,Np,NwMax,2),Eloc(NwMax,2))

allocate (LogWf(0:Nmax+1))

!Initial tasks:
! 
! - Build table for Jastrow Wave Function
! - Initialize random number generator
! - Read/Build the initial configuration of the walkers

call Jastrow_Table(rcut,Rm,LogWF)  

call sgrnd(seed)

call init(LogWF,NwStart,NwMax,opt,Walker,DriftForce,Eloc,E0)

!Priting simulation parameters

print *, '============================================================='
print *, '# DIFFUSION MONTE CARLO SIMULATION:                          '
print *, ''
print *, '   >> Simulation parameters: '
print *, ''
print *, '      -> Density                       :',density
print *, '      -> Number of particles           :',Np
print *, '      -> Variational parameter         :',Rm
print *, '      -> Time step                     :',dt
print *, '      -> Size of the simulation cell   :',Lbox
print *, '      -> Number of calculation blocks  :',Nblock
print *, '      -> Number of iteration per block :',Nstep
print *, '      -> Initial number of walkers     :',NwStart
print *, '      -> Critical number of walkers    :',NwCrit
print *, '      -> Maximum number of walkers     :',NwMax
print *, '=============================================================' 

allocate (gr(Nbin),gr_i(Nbin),gr_pur(Nbin))
allocate (Sk(dim,Nk),Sk_i(dim,Nk),Sk_pur(dim,Nk))
allocate (rho(Nbin),rho_i(Nbin))
allocate (nrho(Nbin),nrho_i(Nbin))
allocate (Ep_w(NwMax,2),Ep_rw(NwMax,2))
allocate (gr_w(Nbin,NwMax,2),gr_rw(Nbin,NwMax,2))
allocate (Sk_w(dim,Nk,NwMax,2),Sk_rw(dim,Nk,NwMax,2))
allocate (AvGr(Nbin),AvGr2(Nbin),VarGr(Nbin))
allocate (AvNr(Nbin),AvNr2(Nbin),VarNr(Nbin))
allocate (AvSk(dim,Nk),AvSk2(dim,Nk),VarSk(dim,Nk))
       
!Initialize the average observables and its variances

AvE   = 0.d0
AvEc  = 0.d0
AvKf  = 0.d0
AvEp  = 0.d0

AvE2  = 0.d0
AvEc2 = 0.d0
AvKf2 = 0.d0
AvEp2 = 0.d0

AvGr  = 0.d0
AvSk  = 0.d0
AvNr  = 0.d0

AvGr2 = 0.d0
AvSk2 = 0.d0
AvNr2 = 0.d0

VarGr = 0.d0
VarSk = 0.d0
VarNr = 0.d0

Ep_pur  = 0.d0
gr_pur  = 0.d0
Sk_pur  = 0.d0
ngr_pur = 0

Ep_rw = 0.d0
gr_rw = 0.d0
Sk_rw = 0.d0

!Beginning of the Montecarlo calculation 

open (unit=3,file='e_dmc.out') 

write (3,*) '# Block -',' <E>/Np -',' <E0/Np> -',' <K>/Np -',&
     & ' <Kf>/Np -',' <V>/Np -',' <Vpure>/Np'

allocate (R(dim,Np),F(dim,Np))

do iblock=1,Nblock

   call cpu_time(begin)

   !Initialize variables of each block

   BlockAvE   = 0.d0
   BlockAvEc  = 0.d0
   BlockAvKf  = 0.d0
   BlockAvEp  = 0.d0
   
   BlockAvE2  = 0.d0
   BlockAvEc2 = 0.d0
   BlockAvKf2 = 0.d0
   BlockAvEp2 = 0.d0
   
   Ep_w = 0.d0
   gr_w = 0.d0
   Sk_w = 0.d0

   gr   = 0.d0
   Sk   = 0.d0
   rho  = 0.d0
   nrho = 0.d0
   ngr  = 0  
 
   do istep=1,Nstep
      
      NwEnd = 0

      WalkerE  = 0.d0
      WalkerEc = 0.d0
      WalkerKf = 0.d0
      WalkerEp = 0.d0
   
      do iw=1,NwStart

         do ip=1,Np
            do k=1,dim
               
               R(k,ip) = Walker(k,ip,iw,in)
               F(k,ip) = DriftForce(k,ip,iw,in)

            end do
         end do     

         !Gaussian diffusion

         call GaussianDiffusion(sigma,mu,R)
         
         !Drift integration 

         call DriftIntegration(LogWF,R,F,dt)
        
         !Evaluating properties of the system

         gr_i   = 0.d0
         Sk_i   = 0.d0
         rho_i  = 0.d0
         nrho_i = 0.d0

         call Forces(.true.,LogWF,R,E,Kin,Kf,Pot,F,gr_i,rho_i,nrho_i)
         
         call StructureFactor(Nk,R,Sk_i)
         
         !Branching
        
         Nsons = int(exp(dt*(E0-0.5d0*(E+Eloc(iw,in))))+grnd())

         do isons=1,Nsons

            NwEnd = NwEnd+1
            
            if (NwEnd>NwMax) then
               print *, 'THE POPULATION HAS GROWN BEYOND THE MAXIMUM!!'
               stop
            end if
  
            do ip=1,Np
               do k=1,dim
                  
                  Walker(k,ip,NwEnd,io)     = R(k,ip)
                  DriftForce(k,ip,NwEnd,io) = F(k,ip)              

               end do
            end do

            Eloc(NwEnd,io) = E

            !Branching of the evaluated magnitudes for calculate the 
            !pure estimators

            Ep_rw(NwEnd,io) = Ep_rw(iw,in)
            Ep_w(NwEnd,io)  = Ep_w(iw,in)+Pot

            do ic=1,Nbin

               gr_rw(ic,NwEnd,io) = gr_rw(ic,iw,in) 
               gr_w(ic,NwEnd,io)  = gr_w(ic,iw,in)+gr_i(ic)

            end do
             
            do ik=1,Nk
               do k=1,dim
                  Sk_rw(k,ik,NwEnd,io) = Sk_rw(k,ik,iw,in)
                  Sk_w(k,ik,NwEnd,io)  = Sk_w(k,ik,iw,in)+Sk_i(k,ik)
               end do
            end do

         end do

         !Evaluate averages over the ensemble of walkers
         
         !Energies

         call Accumulate(Nsons*E,Nsons*Kin,Nsons*Kf,Nsons*Pot,&
                        &WalkerE,WalkerEc,WalkerKf,WalkerEp)

         !Structural properties
         
         ngr = ngr+Nsons

         do ic=1,Nbin
            gr(ic)   = gr(ic)+Nsons*gr_i(ic)
            rho(ic)  = rho(ic)+Nsons*rho_i(ic)
            nrho(ic) = nrho(ic)+Nsons*nrho_i(ic)
         end do

         do ik=1,Nk
            do k=1,dim
               Sk(k,ik) = Sk(k,ik)+Nsons*Sk_i(k,ik)
            end do
         end do

      end do

      !Normalization of the averages in each iteration

      !Energies

      call NormalizeAv(NwEnd,WalkerE,WalkerEc,WalkerKf,WalkerEp)
 
      !Adjusting the reference energy to control the population

      if (NwEnd >= 1.1d0*NwCrit) then
         E0 = WalkerE+gamma*(1.d0-real(NwEnd)/real(NwCrit))/dt
      else if (NwEnd <= 0.9d0*NwCrit) then
         E0 = WalkerE+gamma*(1.d0-real(NwEnd)/real(NwCrit))/dt
      else
         E0 = WalkerE
      end if

      !Swap old and new generations of walkers and update the number
      !of walkers and the reference energy

      in      = 3-in
      io      = 3-io
      NwStart = NwEnd
      
      if (iblock>1) ngr_pur = ngr_pur+NwEnd

      !Accumulating the averages to evaluate the block average

      call Accumulate(WalkerE,WalkerEc,WalkerKf,WalkerEp,&
                     &BlockAvE,BlockAvEc,BlockAvKf,BlockAvEp)

      call Accumulate(WalkerE**2,WalkerEc**2,WalkerKf**2,WalkerEp**2,&
                     &BlockAvE2,BlockAvEc2,BlockAvKf2,BlockAvEp2)

   end do

   !Normalization of block averages and calculation of variances

   call NormalizeAv(Nstep,BlockAvE,BlockAvEc,BlockAvKf,BlockAvEp)
   
   call NormalizeAv(Nstep,BlockAvE2,BlockAvEc2,BlockAvKf2,BlockAvEp2)

   BlockVarE  = Var(Nstep,BlockAvE,BlockAvE2) 
   BlockVarEc = Var(Nstep,BlockAvEc,BlockAvEc2)
   BlockVarKf = Var(Nstep,BlockAvKf,BlockAvKf2)
   BlockVarEp = Var(Nstep,BlockAvEp,BlockAvEp2)

   !Accumulation of total averages
    
   call Accumulate(BlockAvE,BlockAvEc,BlockAvKf,BlockAvEp,&
                  &AvE,AvEc,AvKf,AvEp)

   call Accumulate(BlockAvE**2,BlockAvEc**2,BlockAvKf**2,BlockAvEp**2,&
                  &AvE2,AvEc2,AvKf2,AvEp2)

   !Pure estimators

   do iw=1,NwStart
      
      Ep_pur = Ep_pur+Ep_rw(iw,in)

      do ic=1,Nbin
         gr_pur(ic) = gr_pur(ic)+gr_rw(ic,iw,in)
      end do

      do ik=1,Nk
         do k=1,dim
            Sk_pur(k,ik) = Sk_pur(k,ik)+Sk_rw(k,ik,iw,in)
         end do      
      end do

   end do
   
   !Output partial results

   call Normalize(.false.,density,Nk,ngr,gr,Sk,rho,nrho)

   call AccumGr(gr,AvGr,AvGr2)
   call AccumSk(Nk,Sk,AvSk,AvSk2)
   call AccumNr(rho,AvNr,AvNr2)

   if (mod(iblock,5)==0) then
      call Normalize(.true.,density,Nk,ngr_pur,gr_pur,Sk_pur,rho,nrho)
   end if

   !Pure estimators

   do iw=1,NwStart

      Ep_rw(iw,in) = Ep_w(iw,in)
      
      do ic=1,Nbin
         gr_rw(ic,iw,in) = gr_w(ic,iw,in)
      end do
      
      do ik=1,Nk
         do k=1,dim
            Sk_rw(k,ik,iw,in) = Sk_w(k,ik,iw,in)
         end do
      end do

   end do

   BlockEp_pur = Ep_pur/real(ngr_pur)
   
   if (iblock>1) then

      write (3,'(10g20.10e3)') real(iblock),BlockAvE/Np,E0/Np,BlockAvEc/Np,&
         & BlockAvKf/Np,BlockAvEp/Np,BlockEp_pur/Np
      
   end if

   call cpu_time(end)

   !Printing results of the block

   print *, '----------------------------------------------------------'
   print *, '# Block number :',iblock
   print *, ''
   print *, '  > <E>        =',BlockAvE/Np,'+/-',BlockVarE/Np
   print *, '  > <K>        =',BlockAvEc/Np,'+/-',BlockVarEc/Np
   print *, '  > <Kf>       =',BlockAvKf/Np,'+/-',BlockVarKf/Np
   print *, '  > <V>        =',BlockAvEp/Np,'+/-',BlockVarEp/Np
   print *, '  > <Vp>       =',BlockEp_pur/Np
   print *, '  > Walkers    =',NwEnd
   print *, '  > Time/block =',end-begin
  
end do

close (unit=12)

!Normalizing the total averages of the simulation

call NormalizeAv(Nblock,AvE,AvEc,AvKf,AvEp)
call NormalizeAv(Nblock,AvE2,AvEc2,AvKf2,AvEp2)

VarE  = Var(Nblock,AvE,AvE2)
VarEc = Var(Nblock,AvEc,AvEc2)
VarKf = Var(Nblock,AvKf,AvKf2)
VarEp = Var(Nblock,AvEp,AvEp2)

!Normalize the histogram for the pair correlation function

call NormAvGr(Nblock,AvGr,AvGr2,VarGr)
call NormAvSk(Nblock,Nk,AvSk,AvSk2,VarSk)
call NormAvNr(Nblock,AvNr,AvNr2,VarNr)

call Normalize(.true.,density,Nk,ngr_pur,gr_pur,Sk_pur,rho,nrho)

!Print final results

print *, '============================================================='
print *, '# FINAL RESULTS OF THE SIMULATION:'
print *, ''
print *, '  > <E>  =',AvE/Np,'+/-',VarE/Np
print *, '  > <K>  =',AvEc/Np,'+/-',VarEc/Np
print *, '  > <Kf> =',AvKf/Np,'+/-',VarKf/Np
print *, '  > <V>  =',AvEp/Np,'+/-',VarEp/Np
print *, '============================================================='

deallocate (Walker,Eloc,DriftForce,F,R)
deallocate (rho,rho_i)
deallocate (nrho,nrho_i)
deallocate (gr,gr_i,gr_pur)
deallocate (Sk,Sk_i,Sk_pur)
deallocate (gr_w,gr_rw)
deallocate (Sk_w,Sk_rw)
deallocate (AvGr,AvGr2,VarGr)
deallocate (AvSk,AvSk2,VarSk)
deallocate (AvNr,AvNr2,VarNr)

end program dmc
