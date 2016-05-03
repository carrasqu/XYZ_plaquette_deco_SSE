
! Kagome lattice XYZ model SSE program     
!-------------------------------------------------

!--------------------
module configuration
!--------------------

save

!lattice has L=3*lx*lysites


integer :: lx  
integer :: ly 
integer :: lz ! not used in this code
integer :: L  ! total number of sites or spins
integer :: nb ! number of plquettes
integer :: nver ! number of vertices
integer(4) :: nhl ! number of sites or spins again
integer(4) :: ndim ! number of sites in the bravais lattice of the kagome
integer(4) :: nbase ! number of lattice sites in the unit cell
integer(4) :: nindm,nindms,maxmulti,maxmultis ! dimensions of tables of neighbours
integer(4), allocatable :: imulti(:),ivicn(:,:,:),ivicns(:,:,:),imultis(:)

integer :: nh ! order of the SSE expansion 
integer :: mm ! max order considered in the SSE expansion 
integer :: dd ! dimensionality
integer :: disseed ! seed for the disorder  configuration
integer :: nl ! number of loops per update
integer(8) :: maxloop ! maximum length of the loop operators
integer :: mloopc ! coefficient for maximum loop length maxloop=mloopc*<n>
integer :: termal ! termalization call sign
integer(8) :: totvisit
integer(8) :: looplength
integer :: spino
integer :: nh0
integer :: densiono ! measure density?
integer :: info ! lapack info variable
integer :: inittype ! type of initialization (mf based=1 , random=0)
integer :: thermalization ! whether the simulation is thermalizing or not
integer :: full ! what part of nk to measure (0 only along x. 1 full nk)
integer :: directed
integer :: restart ! restart=1(restart) 0(start from scratch)
integer :: kstart ! number of times writeresults have been executed 

real(8) :: beta  ! inverse temperature 
real(8) :: aprob ! nb*beta
real(8) :: mu    ! chemical potential
real(8) :: cadd  ! sse constant makes all config postive definite 
real(8) :: delta ! bound of the disordered potential
real(8) :: t     ! "hopping" matrix element for S^+S^- +  S^-S^+ (down triangle)
real(8) :: t2     ! "hopping" matrix element for S^+S^- +  S^-S^+(up   triangle) 

real(8) :: tp    ! "hopping" matrix element for S^+S^+ +  S^-S^- (down triangle)
real(8) :: tp2    ! "hopping" matrix element for S^+S^+ +  S^-S^- (up triangle)
real(8) :: Jz    ! SzSz term (or NN  interaction V for hardcore bosons )(down triangle)
real(8) :: Jz2    ! SzSz term (or NN  interaction V for hardcore bosons )(up triangle)
real(8) :: sup   ! amplitude of the superlattice potential is present
real(8) :: vo    ! strength of the optical trap potential vo*(r)**2, r distance from the center of the trap
real(8) :: timeli,timelo,timeme,timew

real(8), allocatable :: epsdis(:) ! espdis(i)=mu+ epsilon_i : disorder configuration plus overall chemical potential
real(8), allocatable :: prob(:,:,:,:) ! probabilities for the loop updates
                       !prob(bond,vertex type(st bond p), entrance leg, exit leg) en principio.  

real(8), allocatable :: tvector(:),tpvector(:),jzvector(:) ! vectors of the hamiltonian parameters

!real(8), allocatable :: wdatgf(:,:),wdat2gf(:,:)                        

real(8), allocatable :: vtex(:,:) ! weight of the vertices
real(8) :: NNN
real(8) :: z
real(8), allocatable :: gf(:) ! one body density matrix
!real(8), allocatable :: eigenval(:) ! eigenvalues of gf
!real(8), allocatable :: workgf(:)

real(kind=8), allocatable :: loc(:) ! local density
real(kind=8), allocatable :: loc0(:) ! local density reduced to the zeroth node


integer, allocatable :: tryloop(:,:)

real(8) :: densgf

integer, allocatable :: spin(:),tspin(:) ! spin (-1 or 1 ) or (boson 0 or 1 respectively) configuration
integer, allocatable :: bsites(:,:) ! bond sites which is the lattice information
integer, allocatable :: bbond(:,:) 
integer, allocatable :: sm(:) ! operator string

integer, allocatable :: frstspin(:) ! temporary; helps constructing vtx
integer, allocatable :: lastspin(:) ! temporary; helps constructing vtx
integer, allocatable :: vtx(:) ! linked vertex storage
!integer, allocatable :: A(:,:,:)

integer, allocatable :: vtype(:) ! type of vertex
integer, allocatable :: vtypel(:) ! temporary type of vertex during loop construction
 
integer, allocatable :: newvtx(:,:,:) ! newvtx(lentrance,lexit,vtype before the change)
integer, allocatable :: optype(:) ! 0 for diagonal 1 for off-diagonal operators,it is function of the vertex type vtype(:)
integer, allocatable :: legspin(:,:) ! value of the spin or boson density for a given leg and type of vertex


end module configuration

!----------------------!
 module measurementdata
!----------------------!
 save

 real(8) :: enrg1=0.d0
 real(8) :: enrg2=0.d0
 real(8) :: kinetic=0.0d0
 real(8) :: stiff=0.d0
 real(8) :: den=0.0d0
 real(8) :: comp2=0.0d0
 real(8) :: comp=0.0d0
 real(8) :: compress=0.0
 real(8) :: ecomp=0.0d0
 real(8) :: gfden=0.0d0
 real(8) :: token(8)
 real(8) :: token0(8)
 real(8) :: data1(8)=0.d0
 real(8) :: data2(8)=0.d0
 real(kind=8) pi

 ! for measurements of the full greens function g(i,j) (memory intensive)
 !! temporary gij
 !real(8), allocatable :: gfm(:,:)
 !real(8), allocatable :: gfm0(:,:)
 !real(8), allocatable :: datagfm(:,:)
 !real(8), allocatable :: data2gfm(:,:) 
 !real(8), allocatable :: databoot(:,:)   

 ! one body density matrix 
 real(8), allocatable :: gf0(:)
 real(8), allocatable :: datagf(:)
 real(8), allocatable :: data2gf(:)
 real(8), allocatable :: databootgf(:)

 ! local density variables
 real(kind=8),allocatable :: dataloc(:)
 real(kind=8),allocatable :: dataloc2(:)

 ! estimate error in condensate fraction
 integer :: boots

 ! momentum distribution  

 real(kind=8),allocatable :: rnk(:,:,:)    ! real part of nk
 real(kind=8),allocatable :: ink(:,:,:)    ! imag part of nk 
 real(kind=8),allocatable :: ernk(:,:,:)   ! error real part
 real(kind=8),allocatable :: arnk(:,:,:)    ! aver real part of nk
 real(kind=8),allocatable :: aink(:,:,:)    ! aver imag part of nk 
 real(kind=8),allocatable :: rnk2(:,:,:)   ! error real part  
 real(kind=8),allocatable :: ink2(:,:,:)   ! error imag part 
  
 real(kind=8),allocatable :: kkx(:,:)     ! k values
 real(kind=8),allocatable :: ord(:,:)  ! coordinates
 real(kind=8),allocatable :: distances(:,:)  ! coordinates
 integer(4),allocatable :: cdistances(:)
 integer(4) indist                           ! number of distances considered
 integer(4) nqv                              ! number of momenta considered   
 integer(4),allocatable :: site_to_d(:,:) 
 real(kind=8),allocatable :: q(:,:)    !
 real(kind=8)ba1(2),ba2(2)

! four point fourier transform
 integer(4) fqmeas ! how often to measure sq. should be less than msteps
 integer(4) fqmsteps
! real(kind=8),allocatable :: fourq(:,:,:)    ! real part of Bq
! real(kind=8),allocatable :: fourq0(:,:,:)
! real(kind=8),allocatable :: datafourq(:,:,:)      
! real(kind=8),allocatable :: data2fourq(:,:,:)  

! real(kind=8) fourb    ! real part of Bq
! real(kind=8) fourb0
! real(kind=8)datafourb
! real(kind=8) data2fourb

real(kind=8),allocatable :: real_bond(:,:,:)
real(kind=8),allocatable :: real_bond0(:,:,:) 
real(kind=8),allocatable :: datareal_bond(:,:,:)
real(kind=8),allocatable :: data2real_bond(:,:,:)

! momentum distribution directly in the loop update (expensive)
! real(kind=8),allocatable :: rN(:,:,:) ! momentum distribution
! real(kind=8),allocatable :: iN(:,:,:) ! momentum distribution
! complex(kind=16),allocatable :: cN(:,:,:) ! momentum distribution 
! real(kind=8),allocatable :: rN0(:,:,:) ! momentum distribution
! real(kind=8),allocatable :: iN0(:,:,:) ! momentum distribution
! real(kind=8),allocatable :: datarN(:,:,:) ! momentum distribution
! real(kind=8),allocatable :: dataiN(:,:,:) ! momentum distribution
! real(kind=8),allocatable :: data2rN(:,:,:) ! momentum distribution
! real(kind=8),allocatable :: data2iN(:,:,:) ! momentum distribution 
  

! Structure factor Sq
integer(4) sqmeas ! how often to measure sq. should be less than msteps
integer(4) sqmsteps  
real(kind=8),allocatable :: rnq(:,:)    ! for accumulating the measurements (  \sum_i S^z_i cos(q . r_i) )
real(kind=8),allocatable :: inq(:,:)    ! for accumulating the measurements (  \sum_i S^z_i sin(q . r_i) )
real(kind=8),allocatable :: cost(:,:) ! cosine table 
real(kind=8),allocatable :: sint(:,:) ! sin table
real(kind=8),allocatable :: sq(:,:,:)     ! sq(x) structure factor at configuration x 
real(kind=8),allocatable :: sq0(:,:,:)    ! after reducing over different nodes
real(kind=8),allocatable :: sqt(:,:,:)    ! sq temp during subroutine check
real(kind=8),allocatable :: datasq(:,:,:) ! accumulating measurements over different bins
real(kind=8),allocatable :: data2sq(:,:,:) ! accumulating measurement^2 to compute errors
real(kind=8),allocatable :: rhoi(:)
real(kind=8),allocatable :: cnq(:,:)
real(kind=8),allocatable :: snq(:,:)
real(kind=8),allocatable :: rdsq(:,:)  
real(kind=8),allocatable :: idsq(:,:)  
real(kind=8),allocatable :: rdsq0(:,:) 
real(kind=8),allocatable :: idsq0(:,:)
real(kind=8),allocatable :: datardsq(:,:) 
real(kind=8),allocatable :: dataidsq(:,:) 
real(kind=8),allocatable :: data2rdsq(:,:) 
real(kind=8),allocatable :: data2idsq(:,:)  
real(kind=8),allocatable :: bootloc(:)
real(kind=8) on

! local triangle compressibility
real(kind=8)tri(2),tri0(2),trir(2),datatri(2),data2tri(2)

! Bq= 1/L \sum_{alpha,beta} exp(iq(r_alpha-r_beta)) < (S+S+)_alpha  (S-S-)_beta > 

 ! <n_i n_j>
 real(8), allocatable :: cf(:) 
 real(8), allocatable :: cf0(:)
 real(8), allocatable :: datacf(:)
 real(8), allocatable :: data2cf(:)
 real(8), allocatable :: databootcf(:)
 integer(4),allocatable :: site_to_bondt(:,:,:)
  
 integer :: lk
 real(kind=8) :: kmax
 complex(8) ima

 character(len=1024) :: filename
 character(len=1024) :: format_string
 end module measurementdata
!--------------------------!



!============================!
 program hcb_sse
!============================!
 use configuration
 use measurementdata
 use mpi
 implicit none
 
 !include 'mpif.h'
 integer(4) stats(MPI_STATUS_SIZE)

 integer :: nbins,msteps,isteps,iseedd,iseeddcopy,k,m,irstart,mstepsq,kkk

 !parallelization variables
 integer :: size,rank,ierr

 real(8) :: drand1,avt,avn, time1,time,avn0,avt0,con
 integer :: countones,countonesD,front ! counts how many times the increase in nl is one increasing and decreasing 

 !parallel initialization


 integer :: values(1:8)
 integer, dimension(:), allocatable :: seed
 real(8) :: r

  call mpi_init(ierr)

  if (ierr.ne.mpi_success) write(6,*)'some problem initializing parallel'

  call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

  

 !input file reading. All processors read the file.
 open(unit=10,file='input.in',form='formatted',status='unknown')

 !read(10,*) dd
 dd=2 ! works only for the kagome
 
 read(10,*) lx 
 read(10,*) ly
 !read(10,*) lz 

 read(10,*) beta
 read(10,*) mu
 read(10,*) t
 read(10,*) t2 
 read(10,*) tp
 read(10,*) tp2
 read(10,*) Jz
 read(10,*) Jz2
 read(10,*) nbins
 read(10,*) msteps
 read(10,*) sqmeas 
 read(10,*) isteps
 read(10,*) iseedd
 read(10,*) disseed
 read(10,*) delta
 read(10,*) cadd 
 read(10,*) sup
 read(10,*) vo
 read(10,*) densiono
 read(10,*) boots
 read(10,*) inittype
 read(10,*) nqv
 read(10,*) kmax
 read(10,*) full
 read(10,*)directed
 read(10,*) mloopc
 read(10,*) restart 

  if(abs(Jz)>=abs(Jz2))then 
   cadd=4.0+3.0*abs(mu)/4.0+3.0*abs(Jz)/4.0
  elseif(abs(Jz2)>abs(Jz))then
   cadd=4.0+3.0*abs(mu)/4.0+3.0*abs(Jz2)/4.0
  end if

 fqmeas=1 
 !sqmeas
 ! nbins: number of bins, msteps: monte carlo steps, isteps: thermalization steps,iseed: seed for the random generator,
 ! mu:chemical potential,disseed: seed for the disorder, sup=amplitude staggered potential
 ! vo: the strength of the parabolic potential. 
 ! densiono =1 measure loc density, 0 to avoid writing it
 ! inittype= 1 initialization based on the  meanfield density. 0 for random initialization ! lk: number of k points in the momentum distribution along x (lx/2 should be used if pbc are used and one wants the right quasimomenta)
 ! kmax : maximum kx value considered. (kmax should bi \pi if pbc are used to
 ! get the right quasimomenta)
 ! directed=0: directed loop updates with no bounces when possible. If not possible, minimal bounces. directed>0 : heat bath solution 
 ! mloopc coefficient for maximum loop length maxloop=mloopc*<n> 

 close(10)


! if(rank==0)then ! rank 0 writes the energy each msteps
! open(unit=13,file='fort.12',form='unformatted',status='unknown') ! energy at  each mc step to check binning
! end if

!initialization ------------------------------------------------------- 

 if(dd==1)then
  ly=1 
  lz=1
 elseif(dd==2)then
  lz=1   
 endif         

 z=2 

 pi=2.0d0*asin(1.0d0)
 ima=cmplx(0.0d0,1.0d0)

 call makelattice()  !construct the lattice
 
 call newvertex()  

 aprob=beta*nb

 if(directed==0)then
  call probupdates() ! construct the probabilities used in the loop updates. Directed loop updates
 elseif(directed>0)then
  call probupdates_heat_bath() ! heat bath solution for the updates
 end if

! write(6,*)'iseedd before',iseedd
! call restarts(iseedd)
! iseeddcopy=iseedd
 ! this bit makes a different seed for each processor
!  irstart=iseedd
!  do k=1,rank
!   irstart=ran(iseedd)*2**29
!  enddo
!  iseedd=2*irstart+1

 ! this bit makes a different seed for each processor

 call date_and_time(values=values)

 call random_seed(size=k)
 allocate(seed(1:k))
 seed(:) = values(8)
 call random_seed(put=seed)
 call random_number(r)

 irstart=0
  do k=1,rank+2
   call random_number(r)
   irstart=irstart+10000*r
  enddo
  iseedd=2*irstart+1


  write(6,*)'seed, processor', iseedd, rank




 ! write(6,*) 'seed',  iseedd,rank

 call rand_init(iseedd) ! initializaation of random number
  
 if(restart==0)then

  kstart=1

  call initconf(rank) ! initial random configuration


! write(6,*) 'configuration', spin, rank


! write(6,*)'epsdis',epsdis,rank
!----------------------------------------------------------------------

  nl=200
  maxloop=30000
  avn=0.0d0
  avt=0.0d0
  timeli=0.0d0
  timelo=0.0d0  
  timeme=0.0d0
  timew=0.0d0 
  con=2.0d0
  front=10
  countones=0
  countonesD=0
  allocate(tryloop(nl,4))

  ! Thermalization
  thermalization=1
  if(rank==0)kkk=0
  do  !k=1,isteps ( now it is automatic)
   do m=1,msteps

    call diagonalupdate()
    call loopupdate()
    ! collects the maximum nh in order to enlarge the cutoff
    call mpi_allreduce(nh,nh0,1,mpi_integer4,mpi_max,mpi_comm_world,ierr)
    call adjustcutoff(k)
    avn=avn+dble(nh)
    avt=avt+dble(totvisit)

   end do
  
    avn=avn/dble(msteps)
    avt=avt/dble(msteps)    

    call mpi_reduce(avn,avn0,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(avt,avt0,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

    !if(rank==0)then
    ! avn0=avn0/dble(size)
    ! avt0=avt0/dble(size)
    ! write(6,*)'avn,avt', avn0,avt0,nl  
    ! if(2.0d0*avn0-avt0>0)then
    !  nl=nl+1
    ! else
    !  nl=nl-1      
    !  if(nl==0)nl=1  
    ! end if
     if(rank==0)then

      avn0=avn0/dble(size)
      avt0=avt0/dble(size)
      write(6,*)'avn,avt', avn0,avt0,nl,kkk

      if(con*avn0-avt0>0)then
       nl=nl+ceiling(front*abs(con*avn0-avt0)/abs(con*avn0+avt0))
       write(6,*)'add', front*abs(con*avn0-avt0)/abs(con*avn0+avt0),countones
       k=ceiling(front*abs(con*avn0-avt0)/abs(con*avn0+avt0))
       if(k==1.or.k==2.or.k==3)countones=countones+1
      else
       nl=nl-1*ceiling(front*abs(con*avn0-avt0)/abs(con*avn0+avt0))
       write(6,*)'subs', front*abs(con*avn0-avt0)/abs(con*avn0+avt0),countonesD
            k=ceiling(front*abs(con*avn0-avt0)/abs(con*avn0+avt0))
          if(k==1.or.k==2.or.k==3)countonesD=countonesD+1
       if(nl<=0)nl=1
      end if
     
     kkk=kkk+1 
     maxloop=int(mloopc*avn0)
     avn0=0.0d0
     avt0=0.0d0
    end if

    call mpi_bcast(nl,1,mpi_integer4,0,mpi_comm_world,ierr)
    call mpi_bcast(maxloop,1,mpi_integer8,0,mpi_comm_world,ierr)
    call mpi_bcast(countonesD,1,mpi_integer4,0,mpi_comm_world,ierr)
    call mpi_bcast(countones,1,mpi_integer4,0,mpi_comm_world,ierr) 
    call mpi_bcast(kkk,1,mpi_integer4,0,mpi_comm_world,ierr)

    avn=0.0d0
    avt=0.0d0

   deallocate(tryloop)
   allocate(tryloop(nl,4))
   if(countones>5.and.countonesD>5)then
    write(6,*)'Thermalization seems completed',rank
    exit
   end if
  
  
   if(kkk>isteps)then
    write(6,*)'Thermalization exceeded istep.. kkk, isteps',kkk,isteps
    exit
   end if

  end do

  if(rank==0)then
   write(6,*)'Finished equilibration, M = ',mm, 'maxloop=',maxloop, 'number of loops',nl
  end if


 elseif(restart==1)then

  call restarts(iseedd,rank)

 end if

 thermalization=0
 
 !rN=0.0d0
 !iN=0.0d0
 !rN0=0.0d0
 !iN0=0.0d0 
 !cN=(0.0d0,0.0d0) 
 gf=0.0d0
 cf=0.0d0

 
! gfm=0.0d0
 
 !------- mc run ------------
 do k=kstart,kstart-1+nbins
  
  sqmsteps=0
  fqmsteps=0
  do m=1,msteps
  
   call diagonalupdate()
   
   call loopupdate()
 

   if(looplength>=0)then
    call check(m)
   end if
 
   !!!!call mpi_reduce(nh,nh0,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr) 
 
   !!!!if(rank==0)then 
   !!!! write(13) 1.0d0,1.0d0, -(1.0d0/beta)*dble(nh0)/dble(size)+cadd*dble(nb)
   !!!!end if 

  end do
  
  
  call writeresults(msteps,k,size,rank,nbins,iseeddcopy)

 
 end do



 !---------------------------
 if(rank==0)then
  close(12)
  close(13)
 end if
 close(37)
 
 call mpi_finalize(ierr)
 stop

 end program hcb_sse
!================================!

! MEasures momentum distribution

subroutine nq(li,lo,factor,typ)
use configuration
use measurementdata
implicit none
integer(4)li,lo,ki,ko,factor,i,typ
real(8)ct,aa


ct=dble(factor*typ)*1.0d0/dble(ndim)

ki=mod(li,nbase)
ko=mod(lo,nbase)
if(ki==0)ki=nbase
if(ko==0)ko=nbase



do i=1,nqv
 !!rN(i,ki,ko)=rN(i,ki,ko)+ct*cos(2.0d0*pi*(q(i,1)*(ord(li,1)-ord(lo,1)) + q(i,2)*(ord(li,2)-ord(lo,2))) )
 !!iN(i,ki,ko)=iN(i,ki,ko)+ct*sin(2.0d0*pi*(q(i,1)*(ord(li,1)-ord(lo,1)) + q(i,2)*(ord(li,2)-ord(lo,2)) ) )
  !!cN(i,ki,ko)=cN(i,ki,ko)+ct*exp(2.0d0*pi*ima*(q(i,1)*(ord(li,1)-ord(lo,1)) + q(i,2)*(ord(li,2)-ord(lo,2))) )
end do



! Checks: kinetic energy from worm 
!aa=sqrt(abs( ord(li,1)-ord(lo,1) )**2+abs(ord(li,2)-ord(lo,2))**2)
!if(li==1.and.lo==2)then
!cN(1,1,1)=cN(1,1,1)+1.0d0
!end if


end subroutine nq


subroutine makelattice() 
!------------------------------------------------------------------------------
! Constructs the list of sites bsites(1,b) and bsites(2,b)  bsites(3,b)  in plaquette b
!------------------------------------------------------------------------------

use configuration
use measurementdata
implicit none

integer is,x1,x2,y1,y2,z1,z2,i,bond,j,qv,typeb,wh,sub,k,kk,nbonds
real(8) drand1,qx,qy,qz
complex(16) valuec




is=0 

! ndim number of lattice sites in the bravais lattice
! nbase number of sites in the unit cell.
! total number of sites L=ndim*nbase


read(9) ndim,nbase,nindm,maxmulti

L=nbase*ndim
nhl=L
nb=2*L/nbase ! number of triangular plaquettes ( works for the kagome lattice)

ALLOCATE(ivicn(nhl,maxmulti,nindm))
allocate(ord(L,dd))

read(9) ivicn,ord




! momentum

if(full==1)nqv=ndim

if(nqv>ndim)then
 write(6,*) 'nqv should be less or equal than ndim'
 stop
end if
open(unit=27,file='fort.27',form='formatted',status='unknown')

allocate(q(nqv,2),rnq(nqv,nbase),inq(nqv,nbase),cost(nqv,L),sint(nqv,L),sq(nqv,nbase,nbase),&
 sqt(nqv,nbase,nbase),sq0(nqv,nbase,nbase),rdsq(nqv,nbase),idsq(nqv,nbase),rdsq0(nqv,nbase),&
 idsq0(nqv,nbase),datardsq(nqv,nbase),dataidsq(nqv,nbase),data2rdsq(nqv,nbase),&
 data2idsq(nqv,nbase),datasq(nqv,nbase,nbase), data2sq(nqv,nbase,nbase),rhoi(L),cnq(nqv,nbase),&
snq(nqv,nbase) )

do j=1,nqv
 read(27,*)q(j,1),q(j,2)
end do

close(27)

!allocate(iN(nqv,nbase,nbase),rN(nqv,nbase,nbase),iN0(nqv,nbase,nbase),rN0(nqv,nbase,nbase),&
!dataiN(nqv,nbase,nbase),datarN(nqv,nbase,nbase),data2iN(nqv,nbase,nbase),data2rN(nqv,nbase,nbase),cN(nqv,nbase,nbase)



!iN=0.0d0
!rN=0.0d0
!iN0=0.0d0
!rN0=0.0d0
!dataiN=0.0d0
!datarN=0.0d0
!data2iN=0.0d0
!data2rN=0.0d0

!allocate(fourq(6,6,nqv),fourq0(6,6,nqv),datafourq(6,6,nqv),data2fourq(6,6,nqv))

!fourq=0.0d0
!fourq0=0.0d0
!datafourq=0.0d0
!data2fourq=0.0d0

!fourb=0.0d0
!fourb0=0.0d0
!datafourb=0.0d0
!data2fourb=0.0d0

sq=0.0d0
sq0=0.0d0
sqt=0.0d0
datasq=0.0d0
data2sq=0.0d0
rdsq=0.0d0
idsq=0.0d0
rdsq0=0.0d0
idsq0=0.0d0
datardsq=0.0d0
dataidsq=0.0d0
data2rdsq=0.0d0
data2idsq=0.0d0



allocate(bsites(3,nb)) ! 3 from having 3 spins per plaquette

! cosine and sine tables for structure factors
on=1.0d0/dble(ndim)
do i=1,L 
 do j=1,nqv
  cost(j,i)=on*cos(2.0d0*pi*(q(j,1)*ord(i,1)+q(j,2)*ord(i,2)))
  sint(j,i)=on*sin(2.0d0*pi*(q(j,1)*ord(i,1)+q(j,2)*ord(i,2)))
 end do
enddo


! constructing the plquettes. Works for the kagome defined by the provided
! geometry.d but (hopefully later?) can be generalized.

bond=1 ! bond means plaquette in this code
j=1

do i=1,L
 
   
   if(mod(i,nbase)==1)then
    bsites(1,bond)=i
    bsites(2,bond)=ivicn(i,1,1)
    bsites(3,bond)=ivicn(i,2,1)   
    bond=bond+1
   elseif(mod(i,nbase)==2)then
    bsites(1,bond)=i
    bsites(2,bond)=ivicn(i,3,1)
    bsites(3,bond)=ivicn(i,4,1)   
    bond=bond+1 
   end if 

end do


!write(6,*)'number of plquettes', bond-1, nb
!do i=1,nb
!write(6,*) 'plaquette ',i,'sites',bsites(1,i),bsites(2,i),bsites(3,i)
!end do
!stop


! allocate momentum distribution


allocate(rnk(nqv,nbase,nbase),ink(nqv,nbase,nbase),ernk(nqv,nbase,nbase),rnk2(nqv,nbase,nbase)&
     ,ink2(nqv,nbase,nbase),arnk(nqv,nbase,nbase),aink(nqv,nbase,nbase))


call cdist() ! organizes the different distances in the torus (PBC)



allocate(epsdis(L),gf(indist),gf0(indist),datagf(indist),data2gf(indist),databootgf(indist),&
         dataloc(L),dataloc2(L), loc(L),loc0(L),bootloc(L))

allocate(frstspin(L))
allocate(lastspin(L))

allocate(spin(L),tspin(L))



allocate(cf(indist),cf0(indist),datacf(indist),data2cf(indist),databootcf(indist))

allocate(real_bond(6,6,indist),real_bond0(6,6,indist),datareal_bond(6,6,indist),data2real_bond(6,6,indist))

real_bond=0.0d0
real_bond0=0.0d0 
datareal_bond=0.0d0
data2real_bond=0.0d0 



!allocate(gfm(L,L),gfm0(L,L),datagfm(L,L),data2gfm(L,L),databoot(L,L)) ! full g(i,j)

gf=0.0d0
gf0=0.0d0
datagf=0.0d0
data2gf=0.0d0

dataloc=0.0d0
dataloc2=0.0d0
loc=0.0d0
loc0=0.0d0

cf=0.0d0
cf0=0.0d0
datacf=0.0d0
data2cf=0.0d0

tri=0.0d0
tri0=0.0d0
datatri=0.0d0
data2tri=0.0d0





!gfm=0.0d0
!datagfm=0.0d0
!data2gfm=0.0d0



! disorder configuration in the local chemical potential
! call rand_init(disseed)

      
 epsdis=mu


!stop



! organizes the bonds

allocate(bbond(2,2*L),tvector(2*L),tpvector(2*L),jzvector(2*L))
allocate(site_to_bondt(L,L,2))
bond=1
j=1

do i=1,L
   
 if(mod(i,nbase)==1)then

! Bonds
     !             -
     !            - -
     !          2    5
     !        -       -
! -----4---j----1-------
!  -        -
!   6     3
!    -   -
!      -
      wh=1
      k=ivicn(i,wh,1)
      bbond(1,bond)=i
      bbond(2,bond)=k
      bond=bond+1
   
      wh=2
      k=ivicn(i,wh,1)
      bbond(1,bond)=i
      bbond(2,bond)=k
      bond=bond+1 
 
      wh=3
      k=ivicn(i,wh,1)
      bbond(1,bond)=i
      bbond(2,bond)=k
      bond=bond+1 

      wh=4
      k=ivicn(i,wh,1)
      bbond(1,bond)=i
      bbond(2,bond)=k
      bond=bond+1


      kk=ivicn(i,1,1)
      k=ivicn(kk,2,1)
      bbond(1,bond)=kk
      bbond(2,bond)=k
      bond=bond+1 

      kk=ivicn(i,4,1)
      k=ivicn(kk,4,1)
      bbond(1,bond)=kk
      bbond(2,bond)=k
      bond=bond+1

  end if

end do 

write(6,*) 2*L,bond-1
site_to_bondt=-1
bond=bond-1

sub=1
do i=1,bond
 typeb=mod(i,6)
 if(typeb==0)typeb=6
 !write(6,*)'typeb',typeb,sub
 site_to_bondt(bbond(1,i),bbond(2,i),1)=typeb
 site_to_bondt(bbond(2,i),bbond(1,i),1)=site_to_bondt(bbond(1,i),bbond(2,i),1)
 site_to_bondt(bbond(1,i),bbond(2,i),2)=nbase*(sub-1)+1
 site_to_bondt(bbond(2,i),bbond(1,i),2)=site_to_bondt(bbond(1,i),bbond(2,i),2)
 if(typeb==6)sub=sub+1
 
 if(typeb==2.or.typeb==5.or.typeb==1)then ! up triangle \/ 
  tvector(i)=t2
  tpvector(i)=tp2
  jzvector(i)=Jz2
 elseif(typeb==4.or.typeb==3.or.typeb==6)then
  tvector(i)=t
  tpvector(i)=tp
  jzvector(i)=Jz
 end if
  
end do
write(6,*)ndim,sub-1

!write(6,*)'site i j site to bond'
!do i=1,L
! do j=1,L
! if(site_to_bondt(i,j,1)>0)then
! write(6,*) i,j,site_to_bondt(i,j,1)  
! endif
! end do
!end do
!stop
end subroutine makelattice 

!------------------------
!initial configuration
!-----------------------
 subroutine initconf(rank)
 use measurementdata
 use configuration
 implicit none
 integer :: i,rank 
 real(8) :: drand1, deni,maxni


if(inittype==1)then
 open(unit=38,file='mfdensity.dat',form='unformatted',status='unknown')
 read(38)loc
 close(38)
 maxni=maxval(loc)

 do i=1,L
 
   ! initialize according to mean-field density
   deni=drand1()*maxni

   if(deni<=loc(i))then
    spin(i)=1
   else
    spin(i)=-1
   end if
  !write(6,*) 'Sz', spin(i), 'boson occupancy n_i=' , (iabs(spin(i))+spin(i) )/2

 end do
loc=0.0d0

!  open(unit=18,file='localdensity.dat',form='unformatted',status='unknown')
!   loc=(1+spin)/2.0d0
!  write(18)loc
!  close(18)
! stop
!!hola

elseif(inittype==0)then
 
 do i=1,L
 ! random initialization
  spin(i)=2*int(2.*drand1())-1
 end do


end if

 mm=5 !initial guess for mm the max order of the SSE expansion

 if(rank==0)then
  open(unit=12,file='fort.12',form='formatted',status='unknown')
  open(unit=13,file='fort.13',form='formatted',status='unknown')
 end if

  ! file to write configuration
  WRITE(filename, fmt = '(A3,I0,A4)')'conf',rank,'.dat'
  !print *, trim(filename)
  open(unit=37,file=filename,form='unformatted',status='unknown')

 allocate(sm(mm))
 sm(:)=0 !all mm operators in the initial configuration are of the [0,0]  type => nh=0
 nh=0 
 
 
 allocate(vtx(0:6*mm-1))
 !allocate(A(mm,L,2))
 allocate(vtype(mm),vtypel(mm))
 end subroutine initconf
!-----------------------

 subroutine diagonalupdate()
 use configuration
 implicit none
 integer :: i, b,op,cc
 real(8) :: drand1, ratio,nran,Jztt
 integer :: spint(3), typo,legt(0:5),vertex
 tspin=spin

 


 !vtype(:)=0

 do i=1,mm
  op=sm(i)

  if(op==0)then ! [0,0] operator found
   b=int(drand1()*nb+1)     !select a plaquette at random among the possible nb 

   if(mod(b,2)==1)then
    Jztt=Jz2
   elseif(mod(b,2)==0)then
    Jztt=Jz
   end if

   spint(1)=spin(bsites(1,b))
   spint(2)=spin(bsites(2,b))
   spint(3)=spin(bsites(3,b))  

   ratio=cadd+ 0.5d0/z*(epsdis(bsites(1,b))*spin(bsites(1,b))+ &
             epsdis(bsites(2,b))*spin(bsites(2,b))+ epsdis(bsites(3,b))*spin(bsites(3,b))     )-&
             0.25*Jztt*( spin(bsites(1,b))*spin(bsites(2,b))+ spin(bsites(1,b))*spin(bsites(3,b)) +&
              spin(bsites(2,b))*spin(bsites(3,b)))

   ratio=aprob*ratio/dble(mm-nh)   
 
  
   ! metropolis        
   nran=drand1()           
   if(ratio>=1.0d0 .or. ratio>=nran  )then
    
    sm(i)=2*b  !insert diagonal operator      
    nh=nh+1    ! increases the order of the SSE expansion 

   end if  

  elseif(mod(op,2)==0)then ! diagonal operator found
   !op=2*b
   b=op/2

   if(mod(b,2)==1)then
    Jztt=Jz2
   elseif(mod(b,2)==0)then
    Jztt=Jz
   end if 

   spint(1)=spin(bsites(1,b))
   spint(2)=spin(bsites(2,b))
   spint(3)=spin(bsites(3,b))  

   ratio=cadd+ 0.5d0/z*(epsdis(bsites(1,b))*spin(bsites(1,b))+ &
             epsdis(bsites(2,b))*spin(bsites(2,b))+ epsdis(bsites(3,b))*spin(bsites(3,b))     )-&
             0.25*Jztt*(spin(bsites(1,b))*spin(bsites(2,b))+ spin(bsites(1,b))*spin(bsites(3,b)) +&
             spin(bsites(2,b))*spin(bsites(3,b)))

   ratio=dble(mm-nh+1)/(aprob*ratio)
        
       
         

   !metropolis      
   nran=drand1()
   if(ratio>=1.0d0 .or. ratio>=nran  )then

    sm(i)=0  !remove diagonal operator      
    nh=nh-1  ! decreases the order of the SSE expansion

   end if


  else ! off-diagonal operator found
          
   ! op=2*b+1  
   !propagate the state |alpha(p)> 

   b=op/2 ! (its integer part) 

   spint(1)=spin(bsites(1,b))
   spint(2)=spin(bsites(2,b))
   spint(3)=spin(bsites(3,b))
  
  ! write(6,*) 'should be zero', spint(1)-legspin(0,vtype(i))+ &
  !                             spint(2)-legspin(1,vtype(i))+ &
  !                             spint(3)-legspin(2,vtype(i))
  ! write(6,*) 'shoild be >8',vtype(i)

   spin(bsites(1,b))=legspin(3,vtype(i))
   spin(bsites(2,b))=legspin(4,vtype(i))
   spin(bsites(3,b))=legspin(5,vtype(i))

  end if 

  !construction of the types of vertices
  if(sm(i)>0)then
 
   cc=iabs(spint(1)-spin(bsites(1,b)))+iabs(spint(2)-spin(bsites(2,b)))+iabs(spint(3)-spin(bsites(3,b)))


   if(cc==0)then ! diagonal 
  
     legt(0)=spint(1)
     legt(1)=spint(2)
     legt(2)=spint(3)
     legt(3)=spint(1)
     legt(4)=spint(2)
     legt(5)=spint(3)
   
     call searchvtx(vertex,legt) 
     vtype(i)=vertex
  
   end if

  elseif(sm(i)==0)then

  vtype(i)=0

  end if 


 end do

if(thermalization==0)then
!=== handy for debugging=== 
 typo=0
 do i=1,mm
 
  if(vtype(i).ne.0)then
 ! write(6,*)vtype(i) !,i

  typo=typo+optype(vtype(i))
  
  end if
 end do

if(sum(tspin-spin)/=0)then 
 write(6,*) 'diagonal check', sum(tspin-spin), 'typo=',typo
endif
endif
! write(6,*)'order of exp n', nh

!write(6,*)'sm after',sm
! write(6,*)'type after',vtype
!================================= 
 end subroutine diagonalupdate
 !------------------------------------------------------------

 subroutine loopupdate()
 use configuration
 use measurementdata
 implicit none
 integer :: i,n,b,op,s1,s2,s3,v0,v1,v2,v3,p,ir,br,br1,k,aa,kk,cc,pbef,plast,io
 integer :: jo,j,li,le,ex,break,po,bo,igfo,lii,spinprop
 integer :: lf,pf,bf,pff,bff,lff, igff,igf,lo,typ,it,typoo
 integer :: crossed, go,lfi,jf,nobounce,dist
 real(8) drand1, ratio,nran !,time,time1


 
 densgf=0.0d0
 tryloop=-1
 
 do i=1,nl
   br=drand1()*(mm)+1
   br1=drand1()*L+1
   tryloop(i,1)=br1
   tryloop(i,2)=br 

   if(drand1()<0.5)then
    tryloop(i,3)=1 ! going up
   else
    tryloop(i,3)=0 ! going down         
   end if        

 end do
!------ construction of the linked vertex list -------

 !time=mclock()

 !A=-100 ! Given a random point in space-time, A tells you which is the closest leg going up or down in time
 ! if given the random point A is negative(-100), then that site has only trivial
 ! operators acting on it. 
 ! A is a huge memory waste, not used anymore.  

 frstspin(:)=-1
 lastspin(:)=-1

 vtx=-2
 do p=1,mm

  v0=6*(p-1)
  !p=v0/4+1 ! position of the operator
  op=sm(p) ! type of operator

  if(op/=0)then
   b=op/2  ! plaquette on which operator acts
   s1=bsites(1,b) 
   s2=bsites(2,b)
   s3=bsites(3,b) 
   v1=lastspin(s1) 
   v2=lastspin(s2)
   v3=lastspin(s3)
!   write(6,*) 'v1,s1,v2,s2,v0',v1,s1,v2,s2,v0

   if(v1/=-1)then  ! found a link to v0 
    pbef=v1/6+1
    vtx(v1)=v0
    vtx(v0)=v1

    do i=1,nl
           
     if(tryloop(i,1)==s1)then 
      if(tryloop(i,2)>=pbef+1.and.tryloop(i,2)<=p)then
       if(tryloop(i,3)==1)then
        tryloop(i,4)=v0
       elseif(tryloop(i,3)==0)then
        tryloop(i,4)=v1
       end if        
      end if
     end if       
   
    end do      
    

!    A(pbef+1:p,s1,1)=v0 ! v0 is up
!    A(pbef+1:p,s1,2)=v1 ! v1 is down 
   else
    frstspin(s1)=v0  ! no link so it is first visit to that site
   end if     

   if(v2/=-1)then  ! found a link to v0+1 
    pbef=v2/6+1
    vtx(v2)=v0+1
    vtx(v0+1)=v2
    
    do i=1,nl
           
     if(tryloop(i,1)==s2)then 
      if(tryloop(i,2)>=pbef+1.and.tryloop(i,2)<=p)then
       if(tryloop(i,3)==1)then
        tryloop(i,4)=v0+1
       elseif(tryloop(i,3)==0)then
        tryloop(i,4)=v2
       end if        
      end if
     end if       
   
    end do    

!    A(pbef+1:p,s2,1)=v0+1 ! v0+1 is up
!    A(pbef+1:p,s2,2)=v2   ! v2 is down

   else
    frstspin(s2)=v0+1 ! no link so it is first visit to 
   end if
   
   if(v3/=-1)then  ! found a link to v0+1 
    pbef=v3/6+1
    vtx(v3)=v0+2
    vtx(v0+2)=v3
   
    do i=1,nl
   
     if(tryloop(i,1)==s3)then
      if(tryloop(i,2)>=pbef+1.and.tryloop(i,2)<=p)then
       if(tryloop(i,3)==1)then
        tryloop(i,4)=v0+2
       elseif(tryloop(i,3)==0)then
        tryloop(i,4)=v3
       end if
      end if
     end if
  
    end do

!    A(pbef+1:p,s2,1)=v0+1 ! v0+1 is up
!    A(pbef+1:p,s2,2)=v2   ! v2 is down

   else
    frstspin(s3)=v0+2 ! no link so it is first visit to 
   end if
  
   
   ! new last visited spins in vertex p 
   lastspin(s1)=v0+3
   lastspin(s2)=v0+4
   lastspin(s3)=v0+5
   

  else 
  
   vtx(v0:v0+5)=-1 ! -1 indicates no links in operators.


  end if        



 end do

! links across the "imaginary time boundary" p=1 <--> p=mm

 do i=1,L

  it=frstspin(i)
  io=lastspin(i)
  pbef=it/6+1
  plast=io/6+1

  if(it/=-1)then
   vtx(io)=it
   vtx(it)=io

   if(plast==mm)then

    do j=1,nl
     if(tryloop(j,1)==i)then
      if(tryloop(j,2)>=1.and.tryloop(j,2)<=pbef)then
        if(tryloop(j,3)==1)then
         tryloop(j,4)=it
        elseif(tryloop(j,3)==0)then
         tryloop(j,4)=io
        end if
      end if
     end if 
    end do    


   else
   
    do j=1,nl
     if(tryloop(j,1)==i)then
      if(tryloop(j,2)>=1.and.tryloop(j,2)<=pbef)then
       if(tryloop(j,3)==1)then
        tryloop(j,4)=it 
       elseif(tryloop(j,3)==0)then
        tryloop(j,4)=io
       end if        
      elseif(tryloop(j,2)>=plast+1.and.tryloop(j,2)<=mm)then
       if(tryloop(j,3)==1)then
        tryloop(j,4)=it        
       elseif(tryloop(j,3)==0)then
        tryloop(j,4)=io
       end if
      end if 
     end if        
    end do     



   end if        
  
     

   !if(plast==mm)then
   ! A(1:pbef,i,1)=it
   ! A(1:pbef,i,2)=io
   !else
   ! A(1:pbef,i,1)=it
   ! A(1:pbef,i,2)=io

   ! A(plast+1:mm,i,1)=it
   ! A(plast+1:mm,i,2)=io
           
  ! end if        

  end if        

 end do

  !time1=mclock()
  !timeli=timeli+time1-time
  !write(6,*)'linked vertex list', timeli,time1-time
  !time=mclock()

!  write(6,*) 'i, vtype,plaquette'
!  do i=1,mm
!   write(6,*)i,vtype(i),sm(i)/2
!  end do 

!  write(6,*) 'vtx table'
!  do i=0,6*mm-1

!    write(6,*)i, vtx(i)
!  end do

! write(6,*) 'proposals'
! do i=1,nl
!  write(6,*) i,tryloop(i,1),tryloop(i,2),tryloop(i,3),tryloop(i,4)
! end do
! write(6,*) 'spin',spin


!----- loop construction and update------------------------

 totvisit=0

 vtypel=vtype
! write(6,*) '======= start operator loops======'
loop0: do i=1,nl

  !select a random point in space and time variables
 ! br=drand1()*(mm)+1
 ! br1=drand1()*L+1 

   br1=tryloop(i,1) ! space
   br=tryloop(i,2)  ! time 

   !jo=A(br,br1,1)

   jo=tryloop(i,4)

  if(jo<0)then ! no operators in that worldline
   typ=(1+ spin(br1))/2
   densgf=densgf+typ
   !call nq(br1,br1,2,typ) ! measuring momentum dist directly on the loopupdate is expensive
   dist=site_to_d(br1,br1)
   gf(dist)=gf(dist)+typ
   !gfm(br1,br1)= gfm(br1,br1)+typ ! memory expensive 
  
   cycle

  end if

  !if(drand1()<0.5)then ! worm goes up
   
   !jo=A(br,br1,1)         
   po=jo/6+1
   bo=sm(po)/2
   li=mod(jo,6)
   typ=(1+legspin(li,vtypel(po)))/2 
   dist=site_to_d(br1,br1) 
   gf(dist)=gf(dist)+typ 
   !gfm(br1,br1)=gfm(br1,br1)+typ    
   densgf=densgf+typ
   !call nq(br1,br1,2,typ)


  !else !worm goes down

   !jo=A(br,br1,2)       
   !po=jo/4+1
   !bo=sm(po)/2  
   !li=mod(jo,4)
   !typ=(1+legspin(li,vtypel(po)))/2
   !gf(br1,br1)=gf(br1,br1)+typ
   !densgf=densgf+typ 

  !end if        

  !lo=br1 
  j=jo

  break=0
  looplength=0

  if(li<3)then ! find the time level of the selected closest leg to the initial  point                   
    igfo=po
  else
    if(po==mm)then
     igfo=1
    else
     igfo=po+1
    end if 
  end if     

  igff=igfo

 
  if(li==0.or.li==3)then
   lo=bsites(1,bo)
  elseif(li==1.or.li==4)then
   lo=bsites(2,bo)
  elseif(li==2.or.li==5)then
   lo=bsites(3,bo) 
  endif

  lo=br1

 
loop1:  do while(break==0)  ! loop operator construction

  p=j/6+1 !vertex number
    
  b=sm(p)/2 ! plaquette where the operator is located

  li=mod(j,6) !entrance leg index
 

! **** selection of exit leg ******
  nran=drand1()
  le=0
  ex=0 
loopleg:  do while(ex==0)

   if(vtypel(p)<0) write(6,*) 'p,vtypel(p)', p, vtypel(p)


   if(nran<prob(b,vtypel(p),li,le))then
    ex=1
    exit loopleg
   end if 
   le=le+1

  end do loopleg

  ! le is the exit leg

      if(newvtx(li,le,vtypel(p))<0)then
       write(6,*) 'li,le,vtypel(p),nran', li,le,vtypel(p),nran
       
       write(6,*) 'nb,b,vtypel,li',nb,b,vtypel(p),li
       write(6,*) 'probs',prob(b,vtypel(p),li,0),prob(b,vtypel(p),li,1), &
       prob(b,vtypel(p),li,2),prob(b,vtypel(p),li,3)


       stop
      end if


   vtypel(p)=newvtx(li,le,vtypel(p)) ! update vertex   


   j=6*(p-1)+le
   jf=j


   if(le<3)then !finds the time level of the head of the worm
    igf=p
   else
    if(p==mm)then
     igf=1
    else
     igf=p+1
    end if
   end if

   if(le==0.or.le==3)then
     lf=bsites(1,b) 
   elseif(le==1.or.le==4)then
     lf=bsites(2,b)
   elseif(le==2.or.le==5)then
     lf=bsites(3,b)  
   end if 

    
    if(igf==br)then ! measure green's function after crossing the vertex p

      if(lf==lo)then
      !density has been measured already. 
      else 
        dist=site_to_d(lf,lo) 
        gf(dist)=gf(dist)+1
        !gfm(lf,lo)=gfm(lf,lo)+1 
        !call nq(lf,lo,1,1)
      end if

    endif


   if(igf==br.and.lf==br1)then ! head finds tail and closes the loop?
    exit loop1
   end if        

   j=vtx(j)

   pff=j/6+1 !vertex number

   bff=sm(pff)/2 ! bond where the operator is located

   lff=mod(j,6) !entrance leg index


   if(lff<3)then ! find the time level of the head after traveling.
    igff=pff
   else
    if(pff==mm)then
    igff=1
    else
    igff=pff+1
    end if
   end if
   
   if(lff==0.or.lff==3)then
    lfi=bsites(1,bff)
   elseif(lff==1.or.lff==4)then
    lfi=bsites(2,bff)
   elseif(lff==2.or.lff==5)then
    lfi=bsites(3,bff) 
   end if 

   
   ! this bit finds whether the head is going up or down in time and
   ! whether it has crossed the time boundary or not. 

   if(le==3.or.le==4.or.le==5)then
    go=1 ! going up
    if(igff-igf<0 )then
     crossed=1
    else
     crossed=0
    end if  
   elseif(le==0.or.le==1.or.le==2)then
    go=0 ! going down
    if(igff-igf>0)then
     crossed=1
    else
     crossed=0
    end if  
   end if
   
   aa=igff-igf

  ! ========= handy for debugging =============
  ! write(6,*)
  ! write(6,*) 'mm=',mm 
  ! write(6,*) 'jo=',jo,'jf=',jf,'j=',j
  ! write(6,*) 'li=',li,'le=',le,'lff=',lff
  ! write(6,*) 'po=',po, 'pf=',p,'pff=',pff
  ! write(6,*) 'igfo=',igfo,'igf=',igf,'igff=',igff 
  ! write(6,*) 'lo=',lo,'lf=',lf,'lfi=',lfi
  ! write(6,*) 'aa=',aa,'crossed=',crossed,'go=',go
  ! write(6,*) 
  ! ==========================================

  igfo=br  

  ! green's function measurements and checks if the head reached the tail.

  if(crossed==0)then

     if(aa>0)then 

          if(igfo<=igff.and.igfo>igf)then

            if(lfi==lo)then
             break=1    ! loop closes
            else
              dist=site_to_d(lfi,lo) 
              gf(dist)=gf(dist)+1     
              !gfm(lfi,lo)= gfm(lfi,lo)+1
              !call nq(lfi,lo,1,1)
            end if
          end if
  
     elseif(aa<0)then
    
          if(igfo<igf.and.igfo>=igff)then

            if(lfi==lo)then
             break=1 
            else
              dist=site_to_d(lfi,lo) 
              gf(dist)=gf(dist)+1
              !gfm(lfi,lo)= gfm(lfi,lo)+1  
              !call nq(lfi,lo,1,1)  
            end if
          end if


     end if      

  elseif(crossed==1)then
   
     if(aa>0)then

          if(igfo>=igff)then

            if(lfi==lo)then
             break=1 
            else
              dist=site_to_d(lfi,lo) 
              gf(dist)=gf(dist)+1
              !gfm(lfi,lo)= gfm(lfi,lo)+1
              !call nq(lfi,lo,1,1)
            end if 

          elseif(igfo<igf)then

            if(lfi==lo)then
             break=1 
            else
             dist=site_to_d(lfi,lo) 
             gf(dist)=gf(dist)+1
             !gfm(lfi,lo)= gfm(lfi,lo)+1
             !call nq(lfi,lo,1,1)
            end if
          end if       

     elseif(aa<0)then

          if(igfo<=igff)then

            if(lfi==lo)then
             break=1 
            else
              dist=site_to_d(lfi,lo) 
              gf(dist)=gf(dist)+1
              !gfm(lfi,lo)= gfm(lfi,lo)+1 
              !call nq(lfi,lo,1,1)
            end if 
          elseif(igfo>igf)then
            if(lfi==lo)then
             break=1 
            else
             dist=site_to_d(lfi,lo) 
             gf(dist)=gf(dist)+1
             !gfm(lfi,lo)= gfm(lfi,lo)+1 
             !call nq(lfi,lo,1,1) 
            end if    

          end if

     end if


  end if


   if(li/=le)then 
     looplength=looplength+1 
   end if
  
   ! if the loop does not close after a long construction  it should be killed.

   if(looplength>=maxloop)then 
    write(6,*)'broken===========**************************************'
     looplength=-1
    exit loop0
    exit loop1
  end if



  end do loop1


    totvisit=totvisit+looplength 

 end do loop0

 !time1=mclock()
 !timelo=timelo+time1-time
 !write(6,*) 'loop construction', timelo, time1-time

 if(looplength>=0)then
  
   vtype=vtypel
  
  !update of the operator chain after the loops

  do i=1,mm

   if(sm(i)>0)then 

    b=sm(i)/2
    sm(i)=2*b+optype(vtype(i)) 
    
   end if        

  end do 
 
  ! update of spin configuration

  do i=1,L
 
   if(frstspin(i)/=-1)then
  
    p=frstspin(i)/6+1
    le=mod(frstspin(i),6)  
    spin(i)=legspin(le,vtype(p))

   else ! free spin

    !if(thermalization==0)then       
     if (drand1()<0.5) spin(i)=-spin(i)       
    !end if 

   end if        

  end do 
 
 end if        

!write(6,*) 'spin',spin
!write(6,*) 'spin 2', tspin
!write(6,*) 'after all loops'
!do i=1,mm
!write(6,*)  i,'vtp',vtype(i),'sm',sm(i), 'bond=',sm(i)/2
!write(6,*) 
!end do
!write(6,*)  'states on the legs'
!do i=1,mm
!write(6,*)i, legspin(0:3,vtype(i))
!end do


 !----------------------------------------------------------
 end subroutine loopupdate 
                      
                     
 subroutine bouncefree(ma,i,vlist,ok)
 use configuration
 implicit none
 real(8) ma(6,6),w1,w2,w3,w4,w5,w6
 real(8)dw(6),dwr(6),ap(6,6),C,D
 integer(4)indx(6),ii,jj,i,vlist(6)
 integer :: ok,IER


 ! Syljuasen solution
 ! 
 ! 0 1 1 1 1 1 1
 ! 1 0 1 0 0 0 0
 ! 1 1 0 0 0 0 0 
 ! 1 0 0 0 0 1 0
 ! 1 0 0 0 1 0 1 
 ! 1 0 0 0 0 1 0 

 do ii=1,6
  dw(ii)=vtex(i,vlist(ii))
 end do

dwr=dw

indx=-1

! write(6,*)'bfore', dw

 
 call DPSORT(dw,6,indx,-2,IER) 
!  write(6,*)'after', dw
!  write(6,*) 'index',indx
   
  !DSORTX (x, incx, n, indx)
 ! call dsortx(dw,1,4,indx)
  

 ap=0.0d0

 ap(indx(1),indx(6))=dw(6)/2.0
 ap(indx(1),indx(5))=(dw(5)-dw(6))/2.0
 ap(indx(1),indx(4))=dw(4)-dw(5)/2.0
 ap(indx(4),indx(5))=dw(5)/2.0
 ap(indx(5),indx(6))=dw(6)/2.0
 ap(indx(1),indx(2))=(dw(1)+dw(2)-dw(3)-dw(4))/2.0d0
 ap(indx(1),indx(3))=(dw(1)-dw(2)+dw(3)-dw(4))/2.0d0
 ap(indx(2),indx(3))=(-dw(1)+dw(2)+dw(3)+dw(4))/2.0d0

 
 ap(indx(6),indx(1))=ap(indx(1),indx(6))
 ap(indx(5),indx(1))=ap(indx(1),indx(5))
 ap(indx(4),indx(1))=ap(indx(1),indx(4))
 ap(indx(5),indx(4))=ap(indx(4),indx(5))
 ap(indx(6),indx(5))=ap(indx(5),indx(6))
 ap(indx(2),indx(1))=ap(indx(1),indx(2))
 ap(indx(3),indx(1))=ap(indx(1),indx(3)) 
 ap(indx(3),indx(2))=ap(indx(2),indx(3))
 

 
 ma=ap
  
! write(6,*) dwr(1),sum(ap(1,:))
! write(6,*) dwr(2),sum(ap(2,:))
! write(6,*) dwr(3),sum(ap(3,:))
! write(6,*) dwr(4),sum(ap(4,:))
! write(6,*) dwr(5),sum(ap(5,:))
! write(6,*) dwr(6),sum(ap(6,:))
 

! write(6,*)'MATRIX',ap
 
 
 ok=1 

 do ii=1,6
  do jj=1,6
   if(ma(ii,jj)<0.0d0)ok=0
!    write(6,*) ii,jj,ap(ii,jj)
  enddo
 enddo

! write(6,*) 'a23,a13,a12',  ap(indx(2),indx(3)),ap(indx(1),indx(3)),ap(indx(1),indx(2))


 ! Bounce from the largest vertex only solution
 ! 
 ! 1 1 1 1 1 1 1
 ! 1 0 0 0 0 0 0
 ! 1 0 0 0 0 0 0 
 ! 1 0 0 0 0 0 0
 ! 1 0 0 0 0 0 0 
 ! 1 0 0 0 0 0 0 
 
 if(ok==0)then
  !write(6,*)'a23',ap(indx(2),indx(3))
  ap=0.0d0

  ap(indx(1),indx(1))=dw(1)-(dw(2)+dw(3)+dw(4)+dw(5)+dw(6))

  ap(indx(1),indx(2))=dw(2)
  ap(indx(1),indx(3))=dw(3)
  ap(indx(1),indx(4))=dw(4)
  ap(indx(1),indx(5))=dw(5)
  ap(indx(1),indx(6))=dw(6)
  
  ap(indx(2),indx(1))=dw(2)
  ap(indx(3),indx(1))=dw(3)
  ap(indx(4),indx(1))=dw(4)
  ap(indx(5),indx(1))=dw(5)
  ap(indx(6),indx(1))=dw(6) 

  ma=ap

! write(6,*)' W1 is large'
  
! write(6,*) dwr(1),sum(ap(1,:))

! write(6,*) dwr(2),sum(ap(2,:))

! write(6,*) dwr(3),sum(ap(3,:))

! write(6,*) dwr(4),sum(ap(4,:))

! write(6,*) dwr(5),sum(ap(5,:))

! write(6,*) dwr(6),sum(ap(6,:))

!write(6,*) 'MATRIX W1 is large'
! do ii=1,6
!  do jj=1,6
!    write(6,*) ii,jj,ap(ii,jj)
!  enddo
! enddo

 end if


 ! double check everything is OK, but it should already be. 
 ok=1

 do ii=1,6
  do jj=1,6
   if(ma(ii,jj)<0.0d0)ok=0
  enddo
 enddo
 if(ok==0)then
  write(6,*) 'wrong solutions'
  stop
 end if
 return
 end subroutine bouncefree
 
 subroutine assign(wma,ii,vlist)
 use configuration
 implicit none 
 integer(4)ii,i1,i2,i3,i4,i5,i6,en,ex,vt(0:5),li,vlist(6)
 real(8)tot,wma(6,6)  

  do  li=1,6
   vt(li-1)=vlist(li) 
  end do
 
  do en=0,5
   
   tot=sum(wma(en+1,:)) 
   ex=0
   prob(ii,vt(en),en,ex)=wma(en+1,ex+1)/tot 
   !write(6,*)ii,vt(en),en,ex,wma(en+1,ex+1)/tot 
   do ex=1,5     
     prob(ii,vt(en),en,ex)=prob(ii,vt(en),en,ex-1)+wma(en+1,ex+1)/tot
    ! write(6,*)ii,vt(en),en,ex,wma(en+1,ex+1)/tot
   enddo 
  end do

 end subroutine assign
subroutine probupdates !_loopupdate_XXZ
 
 use configuration
 implicit none
 integer :: i,j,k,ok,update,el,vlist(6)
 real(kind=8)ma(6,6),mach_prec,tot,smallest,tt,tpt,Jzt
 

 allocate(vtex(nb,nver),prob(nb,nver,0:5,0:5))

 smallest=10000000000.0
! prob(bond, vertex type( bond p), entrance leg, exit leg)

 do i=1,nb
   
  if(mod(i,2)==1)then
   tt=t2
   tpt=tp2 
   Jzt=Jz2
  elseif(mod(i,2)==0)then
   tt=t
   tpt=tp
   Jzt=Jz
  end if 
  ! calculation of the vertices 

  ! diagonal vertices

  do j=1,8
   vtex(i,j)=cadd+ 0.5d0/z*(epsdis(bsites(1,i))*legspin(0,j)+ &
             epsdis(bsites(2,i))*legspin(1,j)+ epsdis(bsites(3,i))*legspin(2,j)     )-&
             0.25*Jzt*( legspin(0,j)*legspin(1,j)+legspin(0,j)*legspin(2,j) + legspin(1,j)*legspin(2,j) )
  end do

  ! hopping vertices
  do j=9,20
   vtex(i,j)=tt/2.0
  end do

  ! S+S+ S-S- vertices
  do j=21,32
   vtex(i,j)=tpt/2.0
  end do

  do j=1,nver
 
   if(j<9)then 
    if(vtex(i,j)<smallest)then
     smallest=vtex(i,j) 
    end if        
   end if
   if(vtex(i,j)<0)then
    write(6,*) 'negative vertex i, j , V(i,j)', i,j, vtex(i,j)       
    stop 
   end if        

  end do

  
  ! Solution of directed loop equations
  do j=1,nver
   
   ! finds out the closed sets of vertices
   do el=0,5
    vlist(el+1)=newvtx(0,el,j)     
   end do  
   
   ! determine the solutions 
   call bouncefree(ma,i,vlist,ok)
   
   ! assign the solutions to the probabilities
   call assign(ma,i,vlist)
 
  end do

 
 end do

 write(6,*) '+++++++++++++++ Smallest diagonal VTEX=====',smallest
 
 end subroutine probupdates !_loopupdate_XXZ


! ----------- probabilities according to the loop update (not directed loop
! update)
 subroutine probupdates_heat_bath
 
 use configuration
 implicit none
 integer :: i,j,k,m
 real(8) deno,tt,tpt,Jzt
 allocate(vtex(nb,nver),prob(nb,nver,0:5,0:5))

 prob=0.0d0
! prob(bond, vertex type( bond p), entrance leg, exit leg)

 do i=1,nb
  
 if(mod(i,2)==1)then
   tt=t2
   tpt=tp2
   Jzt=Jz2
  elseif(mod(i,2)==0)then
   tt=t
   tpt=tp
   Jzt=Jz
  end if
  
 ! diagonal vertices
  do j=1,8
   vtex(i,j)=cadd+ 0.5d0/z*(epsdis(bsites(1,i))*legspin(0,j)+ &
             epsdis(bsites(2,i))*legspin(1,j)+ epsdis(bsites(3,i))*legspin(2,j)     )-&
             0.25*Jzt*( legspin(0,j)*legspin(1,j)+legspin(0,j)*legspin(2,j) + legspin(1,j)*legspin(2,j) )
  
  end do

  ! hopping vertices
  do j=9,20
   vtex(i,j)=tt/2.0  
  end do

  ! S+S+ S-S- vertices
  do j=21,32
   vtex(i,j)=tpt/2.0
  end do 

 end do
 

 do i=1,nb
     ! generation of the probabilities
 
  do j=1,nver
   do k=0,5 ! loop entrance leg 

    !write(6,*)'evertex,entrance k, exit m, overtex' 
    deno=0.0d0
    do m=0,5 ! loop exit leg
     !newvtx(en,ex,i)
     deno=deno+vtex(i,newvtx(k,m,j))
     !write(6,*)'evertex,entrance k, exit m, overtex'
    ! write(6,*)j,k,m,newvtx(k,m,j)
    end do

    prob(i,j,k,0)=vtex(i,newvtx(k,0,j))/deno
    do m=1,5
     prob(i,j,k,m)=prob(i,j,k,m-1)+vtex(i,newvtx(k,m,j))/deno
    end do
    
   end do

  end do

 end do

 end subroutine probupdates_heat_bath



 subroutine newvertex
 
 use configuration
 implicit none
 integer i,j,k,li,le,vertex,en,ex
 integer,allocatable :: legt(:)

 nver=32
 allocate(newvtx(0:5,0:5,nver),optype(nver),legspin(0:5,nver),legt(0:5))


 ! diagonal vertices

 !vertex type 1    0 0 0             3 4 5   (names of the legs)
 !                 =====             =====
 !                 0 0 0             0 1 2 

 optype(1)=0 

 legspin(0,1)=-1
 legspin(1,1)=-1
 legspin(2,1)=-1
 legspin(3,1)=-1
 legspin(4,1)=-1
 legspin(5,1)=-1
 


!vertex type 2    X 0 0
!                 =====
!                 X 0 0

 optype(2)=0

 legspin(0,2)=1
 legspin(1,2)=-1
 legspin(2,2)=-1
 legspin(3,2)=1
 legspin(4,2)=-1
 legspin(5,2)=-1


!vertex type 3   0 X 0 
!                =====
!                0 X 0 

 optype(3)=0

 legspin(0,3)=-1
 legspin(1,3)=1
 legspin(2,3)=-1
 legspin(3,3)=-1
 legspin(4,3)=1
 legspin(5,3)=-1


!vertex type 4   0 0 X 
!                =====
!                0 0 X 

 optype(4)=0

 legspin(0,4)=-1
 legspin(1,4)=-1
 legspin(2,4)=1
 legspin(3,4)=-1
 legspin(4,4)=-1
 legspin(5,4)=1
 

!vertex type 5   X X 0 
!                =====
!                X X 0 

 optype(5)=0

 legspin(0,5)=1
 legspin(1,5)=1
 legspin(2,5)=-1
 legspin(3,5)=1
 legspin(4,5)=1
 legspin(5,5)=-1


!vertex type 6   X 0 X 
!                =====
!                X 0 X 

 optype(6)=0

 legspin(0,6)=1
 legspin(1,6)=-1
 legspin(2,6)=1
 legspin(3,6)=1
 legspin(4,6)=-1
 legspin(5,6)=1


!vertex type 7   0 X X 
!                =====
!                0 X X 

 optype(7)=0

 legspin(0,7)=-1
 legspin(1,7)=1
 legspin(2,7)=1
 legspin(3,7)=-1
 legspin(4,7)=1
 legspin(5,7)=1


!vertex type 8   X X X 
!                =====
!                X X X 

 optype(8)=0

 legspin(0,8)=1
 legspin(1,8)=1
 legspin(2,8)=1
 legspin(3,8)=1
 legspin(4,8)=1
 legspin(5,8)=1



! S+S- S-S+ vertices


!vertex type 9   0 X 0  
!                =====
!                X 0 0 

 optype(9)=1

 legspin(0,9)=1
 legspin(1,9)=-1 
 legspin(2,9)=-1
 legspin(3,9)=-1
 legspin(4,9)=1
 legspin(5,9)=-1 


!vertex type 10  X 0 0  
!                =====
!                0 X 0 
      
 optype(10)=1
 
 legspin(0,10)=-1
 legspin(1,10)=1
 legspin(2,10)=-1
 legspin(3,10)=1
 legspin(4,10)=-1
 legspin(5,10)=-1



!vertex type 11  0 0 X  
!                =====
!                0 X 0 
      
 optype(11)=1
 
 legspin(0,11)=-1
 legspin(1,11)=1
 legspin(2,11)=-1
 legspin(3,11)=-1
 legspin(4,11)=-1
 legspin(5,11)=1


!vertex type 12  0 X 0  
!                =====
!                0 0 X 

 optype(12)=1
 
 legspin(0,12)=-1
 legspin(1,12)=-1
 legspin(2,12)=1
 legspin(3,12)=-1
 legspin(4,12)=1
 legspin(5,12)=-1


!vertex type 13  0 0 X  
!                =====
!                X 0 0 

 optype(13)=1

 legspin(0,13)=1
 legspin(1,13)=-1
 legspin(2,13)=-1
 legspin(3,13)=-1
 legspin(4,13)=-1
 legspin(5,13)=1


!vertex type 14  X 0 0  
!                =====
!                0 0 X 

 optype(14)=1

 legspin(0,14)=-1
 legspin(1,14)=-1
 legspin(2,14)=1
 legspin(3,14)=1
 legspin(4,14)=-1
 legspin(5,14)=-1


!vertex type 15  0 X X  
!                =====
!                X 0 X 

 optype(15)=1

 legspin(0,15)=1
 legspin(1,15)=-1
 legspin(2,15)=1
 legspin(3,15)=-1
 legspin(4,15)=1
 legspin(5,15)=1


!vertex type 16  X 0 X  
!                =====
!                0 X X 

 optype(16)=1

 legspin(0,16)=-1
 legspin(1,16)=1
 legspin(2,16)=1
 legspin(3,16)=1
 legspin(4,16)=-1
 legspin(5,16)=1


!vertex type 17  X 0 X  
!                =====
!                X X 0 

 optype(17)=1

 legspin(0,17)=1
 legspin(1,17)=1
 legspin(2,17)=-1
 legspin(3,17)=1
 legspin(4,17)=-1
 legspin(5,17)=1


!vertex type 18  X X 0  
!                =====
!                X 0 X 

 optype(18)=1

 legspin(0,18)=1
 legspin(1,18)=-1
 legspin(2,18)=1
 legspin(3,18)=1
 legspin(4,18)=1
 legspin(5,18)=-1


!vertex type 19  X X 0  
!                =====
!                0 X X 

 optype(19)=1

 legspin(0,19)=-1
 legspin(1,19)=1
 legspin(2,19)=1
 legspin(3,19)=1
 legspin(4,19)=1
 legspin(5,19)=-1


!vertex type 20  0 X X  
!                =====
!                X X 0 

 optype(20)=1

 legspin(0,20)=1
 legspin(1,20)=1
 legspin(2,20)=-1
 legspin(3,20)=-1
 legspin(4,20)=1
 legspin(5,20)=1


! S+S+ S-S- vertices

!vertex type 21  0 0 0  
!                =====
!                X X 0 

 optype(21)=1

 legspin(0,21)=1
 legspin(1,21)=1
 legspin(2,21)=-1
 legspin(3,21)=-1
 legspin(4,21)=-1
 legspin(5,21)=-1

 
!vertex type 22  0 0 0  
!                =====
!                X 0 X 

 optype(22)=1

 legspin(0,22)=1
 legspin(1,22)=-1
 legspin(2,22)=1
 legspin(3,22)=-1
 legspin(4,22)=-1
 legspin(5,22)=-1


!vertex type 23  0 0 0  
!                =====
!                0 X X 

 optype(23)=1

 legspin(0,23)=-1
 legspin(1,23)=1
 legspin(2,23)=1
 legspin(3,23)=-1
 legspin(4,23)=-1
 legspin(5,23)=-1



!vertex type 24  0 0 X  
!                =====
!                X X X 

 optype(24)=1

 legspin(0,24)=1
 legspin(1,24)=1
 legspin(2,24)=1
 legspin(3,24)=-1
 legspin(4,24)=-1
 legspin(5,24)=1


!vertex type 25  0 X 0  
!                =====
!                X X X 

 optype(25)=1

 legspin(0,25)=1
 legspin(1,25)=1
 legspin(2,25)=1
 legspin(3,25)=-1
 legspin(4,25)=1
 legspin(5,25)=-1


!vertex type 26  X 0 0  
!                =====
!                X X X 

 optype(26)=1

 legspin(0,26)=1
 legspin(1,26)=1
 legspin(2,26)=1
 legspin(3,26)=1
 legspin(4,26)=-1
 legspin(5,26)=-1



!vertex type 27  X X X  
!                =====
!                0 0 X 

 optype(27)=1

 legspin(0,27)=-1
 legspin(1,27)=-1
 legspin(2,27)=1
 legspin(3,27)=1
 legspin(4,27)=1
 legspin(5,27)=1


!vertex type 28  X X X  
!                =====
!                0 X 0 

 optype(28)=1

 legspin(0,28)=-1
 legspin(1,28)=1
 legspin(2,28)=-1
 legspin(3,28)=1
 legspin(4,28)=1
 legspin(5,28)=1


!vertex type 29  X X X  
!                =====
!                X 0 0 

 optype(29)=1

 legspin(0,29)=1
 legspin(1,29)=-1
 legspin(2,29)=-1
 legspin(3,29)=1
 legspin(4,29)=1
 legspin(5,29)=1


!vertex type 30  X X 0  
!                =====
!                0 0 0 

 optype(30)=1

 legspin(0,30)=-1
 legspin(1,30)=-1
 legspin(2,30)=-1
 legspin(3,30)=1
 legspin(4,30)=1
 legspin(5,30)=-1



!vertex type 31  X 0 X  
!                =====
!                0 0 0 

 optype(31)=1

 legspin(0,31)=-1
 legspin(1,31)=-1
 legspin(2,31)=-1
 legspin(3,31)=1
 legspin(4,31)=-1
 legspin(5,31)=1


!vertex type 32  0 X X  
!                =====
!                0 0 0 

 optype(32)=1

 legspin(0,32)=-1
 legspin(1,32)=-1
 legspin(2,32)=-1
 legspin(3,32)=-1
 legspin(4,32)=1
 legspin(5,32)=1


!write(6,*) 'newvertex table'
!write(6,*) 'entrance leg, exit leg, input vertex, output vertex '
do i=1,nver
 do en=0,5
   do ex=0,5
    legt=legspin(:,i)
    legt(en)=-legt(en)
    legt(ex)=-legt(ex)
    call searchvtx(vertex,legt)
    newvtx(en,ex,i)=vertex
    !write(6,*) en,ex,i,vertex
   end do 
 end do 
end do

end subroutine newvertex

subroutine searchvtx(vertex,legt)
use configuration
implicit none
integer(4) vertex, i ,j,k,legt(0:5),sumi

search: do i=1,nver

 sumi=0
 do j=0,5
  sumi=sumi+abs(legt(j)-legspin(j,i))
 end do
 if(sumi==0)then
  vertex=i 
  exit search
 end if
end do search

end subroutine searchvtx


!-----------------------------!
 subroutine adjustcutoff(step)
!-----------------------------!
 use configuration
 implicit none

 integer, allocatable :: stringcopy(:),vtypecopy(:)
 integer :: mmnew,step

! mmnew=nh+nh/3 !serial version

  mmnew=nh0+nh0/3

 if (mmnew<=mm) return

 allocate(stringcopy(mm))
 allocate(vtypecopy(mm))

 stringcopy(:)=sm(:)
 vtypecopy(:)=vtype(:)

 deallocate(sm,vtype)

 allocate(sm(mmnew),vtype(mmnew))
 sm(1:mm)=stringcopy(:)
 sm(mm+1:mmnew)=0
 vtype(1:mm)=vtypecopy(:)
 vtype(mm+1:mmnew)=0

 deallocate(stringcopy,vtypecopy)

 mm=mmnew

 deallocate (vtx,vtypel)
 
 !deallocate(A) 
 allocate(vtx(0:6*mm-1),vtypel(mm))
 !allocate(A(mm,L,2))
! open(unit=10,file='results.txt',status='replace')
 
! close(10)

 end subroutine adjustcutoff

!------------------------------------!
 subroutine writeresults(msteps,bins,size,rank,nbins,iseedd)
!------------------------------------!
 use configuration;use mpi; use measurementdata; implicit none

 integer :: i,msteps,bins,shift,size,rank,ierr,j,nbins,qv,k,iseedd
 integer(4) stats(MPI_STATUS_SIZE)

 real(8) :: wdata1(8),wdata2(8), nk0,nk01,errnk0, dii,value ! ,wdatgf(L,L),wdat2gf(L,L)
 complex(8)valuec

 enrg1=enrg1/msteps
 enrg2=enrg2/msteps
 stiff=stiff/msteps
 den=den/msteps
 comp=comp/msteps
 comp2=comp2/msteps
 kinetic=-kinetic/(msteps*beta)
 gfden=gfden/dble(msteps)
 enrg2=(enrg2-enrg1*(enrg1+1.d0))/dble(L)
 enrg1=-enrg1/(beta*dble(L))+cadd*dble(nb)/dble(L)
 stiff=stiff/(beta)
 loc=loc/dble(msteps) ! local density
 sq=sq/dble(sqmsteps)
 rdsq=rdsq/dble(sqmsteps)
 idsq=idsq/dble(sqmsteps)
 !fourq=fourq/dble(fqmsteps)
 real_bond=real_bond/dble(fqmsteps)
 !fourb=fourb/dble(fqmsteps)
 !write(6,*) 'sqmsteps',sqmsteps 
 tri=tri/msteps 
 
 token(1)=enrg1
 token(2)=enrg2
 token(3)=den
 token(4)=kinetic
 token(5)=gfden
 token(6)=stiff
 token(7)=comp
 token(8)=comp2

 call mpi_reduce(token,token0,8,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(gf,gf0,indist,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(tri,tri0,2,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)


! gfm is memory expensive (not used anymore)
! cN is cpu intensive   (not used anymore)

! call mpi_reduce(gfm,gfm0,L*L,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)   
! iN=aimag(cN)
! rN=real(cN)
! call  mpi_reduce(iN,iN0,nqv*nbase*nbase,mpi_real8,mpi_sum,0,mpi_comm_world,ierr) 
! call  mpi_reduce(rN,rN0,nqv*nbase*nbase,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

 call mpi_reduce(loc,loc0,L,mpi_real8,mpi_sum,0,mpi_comm_world,ierr) 
 call mpi_reduce(sq,sq0,nqv*nbase*nbase,mpi_real8,mpi_sum,0,mpi_comm_world,ierr) 
 call mpi_reduce(rdsq,rdsq0,nqv*nbase,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(idsq,idsq0,nqv*nbase,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
! call mpi_reduce(fourq,fourq0,6*6*nqv,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
! call mpi_reduce(fourb,fourb0,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(real_bond,real_bond0,6*6*indist,mpi_real8,mpi_sum,0,mpi_comm_world,ierr) 


 if(rank==0)then ! rank 0

  token0=token0/dble(size)
  gf0=gf0/dble(size)
  !gfm0=dble(size)
  loc0=loc0/dble(size)
  !iN0=iN0/dble(size)
  !rN0=rN0/dble(size) 
  sq0=sq0/dble(size)
  rdsq0=rdsq0/dble(size) 
  idsq0=idsq0/dble(size)
  !fourq0=fourq0/dble(size)
  !fourb0=fourb0/dble(size)
  real_bond0=real_bond0/dble(size)
  tri0=tri0/dble(size)

  ! accumulate the measurements
  data1(1)=data1(1)+token0(1)
  data1(2)=data1(2)+token0(2)
  data1(3)=data1(3)+token0(3)
  data1(4)=data1(4)+token0(4)
  data1(6)=data1(6)+token0(6)
  data1(5)=data1(5)+token0(5)
  data1(7)=data1(7)+token0(7) 
  data1(8)=data1(8)+token0(8)
  
  datagf=datagf+gf0/dble(nl*msteps)
 ! dataiN=dataiN+iN0/dble(nl*msteps)
 ! datarN=datarN+rN0/dble(nl*msteps)
 ! datagfm=datagfm+gfm0/dble(nl*msteps)
  dataloc=dataloc+loc0  
  datasq=datasq+sq0 
  datardsq=datardsq+rdsq0 
  dataidsq=dataidsq+idsq0
!  datafourq=datafourq+fourq0
  !datafourb=datafourb+fourb0 
  datareal_bond=datareal_bond+real_bond0

  datatri=datatri+tri0

   
  ! and the square of them
  data2(1)=data2(1)+token0(1)**2
  data2(2)=data2(2)+token0(2)**2
  data2(3)=data2(3)+token0(3)**2
  data2(4)=data2(4)+token0(4)**2
  data2(6)=data2(6)+token0(6)**2
  data2(5)=data2(5)+token0(5)**2
  data2(7)=data2(7)+token0(7)**2
  data2(8)=data2(8)+token0(8)**2
  data2gf=data2gf+(gf0/dble(nl*msteps))**2
 ! data2iN=data2iN+(iN0/dble(nl*msteps))**2
 ! data2rN=data2rN+(rN0/dble(nl*msteps))**2 
  !data2gfm=data2gfm+(gfm0/dble(nl*msteps))**2 
  dataloc2=dataloc2+loc0**2
  data2sq=data2sq+sq0**2
  data2rdsq=data2rdsq+rdsq0**2 
  data2idsq=data2idsq+idsq0**2 
!  data2fourq=data2fourq+fourq0**2
!  data2fourb=data2fourb+fourb0**2
  data2real_bond=data2real_bond+real_bond0**2

  data2tri=data2tri+tri0**2
 

  do i=1,8
    wdata1(i)=data1(i)/bins
    wdata2(i)=data2(i)/bins
    wdata2(i)=sqrt(abs(wdata2(i)-wdata1(i)**2)/bins)
  end do


  nk0=0.0d0
  errnk0=0.0d0
  gf0=data2gf/dble(bins)
  gf=datagf/dble(bins)
  gf0=sqrt(abs(gf0 -( gf )**2 )/dble(bins))

  tri=datatri/dble(bins)
  tri0=data2tri/dble(bins)
  tri0=sqrt(abs(tri0 -( tri )**2 )/dble(bins))
  


! structure factor, FT is measured already
sq=datasq/dble(bins)
sq0=data2sq/dble(bins)
sq0=sqrt(abs(sq0-sq**2)/dble(bins))

loc=dataloc/bins
loc0=dataloc2/bins 
loc0=sqrt(abs(loc0 -(loc)**2 )/dble(bins))

rdsq=datardsq/dble(bins)
rdsq0=data2rdsq/dble(bins)
rdsq0=sqrt(abs(rdsq0-rdsq**2)/dble(bins))

idsq=dataidsq/dble(bins)
idsq0=data2idsq/dble(bins)
idsq0=sqrt(abs(idsq0-idsq**2)/dble(bins))

open(unit=18,file='isq.dat',form='formatted',status='unknown')
write(18,*) 'sum rule structure factor (should be the density (if nqv=ndim!))', &
 (sum(sq(:,1,1))+sum(sq(:,3,3))+sum(sq(:,3,3)))/dble(nbase)
do i=1,nqv
 write(18,256) q(i,1),q(i,2), idsq(i,1),idsq(i,2),idsq(i,3)
end do

close(18)

open(unit=17,file='sq.dat',form='formatted',status='unknown')

do i=1,nqv
 write(17,256) q(i,1),q(i,2),sq(i,1,1)+sq(i,2,2)+sq(i,3,3), sq0(i,1,1)+sq0(i,2,2)+sq0(i,3,3), &
  rdsq(i,1)*rdsq(i,1)+rdsq(i,2)*rdsq(i,2)+rdsq(i,3)*rdsq(i,3)  
end do

close(17)

256 format(5ES16.8)
!************


!*** compressibility bootcomp()
call bootcomp(wdata1,wdata2)

!************************************************************
  open(10,file='results.txt',status='replace')
  write(10,*)' Cut-off L : ',mm
  write(10,*)' Number of bins completed : ',bins
  write(10,*)' ========================================='
  write(10,10)' -E/N       : ',wdata1(1),wdata2(1)
  write(10,10)'  C/N       : ',wdata1(2),wdata2(2)
  write(10,10)'  <Sz>      : ',wdata1(3),wdata2(3)
  write(10,10)'  <K>       : ',wdata1(4),wdata2(4)
  write(10,10)'  <N_gf>    : ',wdata1(5),wdata2(5)
  write(10,10)'  rho_s     : ',wdata1(6),wdata2(6)
  write(10,10)'  k         : ',compress,ecomp
  write(10,10)'  n1^2      : ',tri(1),tri0(1)
  write(10,10)'  n2^2      : ',tri(2),tri0(2)
  !write(10,10)'  fourb     : ',fourb*(2.0/tp)**2,fourb0*(2.0/tp)**2
  write(10,*)' ========================================='
  10 format(1x,a,2f14.8)
  close(10)

open(unit=55,file='data.dat',form='unformatted',status='unknown')
rewind(55)
write(55)bins,data1,data2,datagf,data2gf,dataloc,dataloc2,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,datareal_bond,&
         data2real_bond,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,lx,ly,beta,mu,t,tp,Jz,iseedd,maxloop,nl
close(55)

write(12,*) token0(1),token0(3),token0(6)
flush(12)
rewind(13)
write(13,*)bins

 
token=0.0d0
token0=0.0d0


end if ! rank 0 

 enrg1=0.d0
 enrg2=0.d0
 stiff=0.d0
 gf=0.0d0 
 gf0=0.0d0
 den=0.0d0
 comp=0.0d0
 comp2=0.0d0 
 kinetic=0.0d0
 gfden=0.0d0
 loc=0.0d0
 loc0=0.0d0
 tri=0.0d0
 tri0=0.0d0
 !iN=0.0d0
 !rN=0.0d0
 !iN0=0.0d0
 !rN0=0.0d0 
 !cN=cmplx(0.0d0,0.0d0)
 sq=0.0d0
 sq0=0.0d0
 rdsq=0.0d0
 idsq=0.0d0
! fourq=0.0d0
! fourq0=0.0d0
 !fourb=0.0
 !fourb0=0.0d0
 real_bond=0.0d0
 real_bond0=0.0d0


 rewind(37)
 write(37)mm,nh
 write(37)sm,vtype,spin
 flush(37)
  

! gfm=0.0d0
! gfm0=0.0d0 
   
 end subroutine writeresults
!---------------------------!

subroutine bootcomp(wdata1,wdata2)
use configuration
use measurementdata
implicit none
integer(kind=4)qq,ll,i,j,k
real(kind=8) u,v,x,y,nboot,n2boot
real(kind=8) drand1,wdata1(8),wdata2(8)

compress=0.0d0
ecomp=0.0d0

do qq=1,boots

  u=drand1()
  v=drand1()
  x=dsqrt(-2.0d0*log(u))*dcos(2*pi*v)
  y=dsqrt(-2.0d0*log(u))*dsin(2*pi*v)
  nboot=wdata1(7)+x*wdata2(7)
  
  u=drand1()
  v=drand1()
  x=dsqrt(-2.0d0*log(u))*dcos(2*pi*v)
  y=dsqrt(-2.0d0*log(u))*dsin(2*pi*v)
  n2boot=wdata1(8)+y*wdata2(8)
 
  compress=compress+ beta*(n2boot-nboot**2)
  ecomp=ecomp+(beta*(n2boot-nboot**2))**2

   
end do

compress=compress/boots
ecomp=ecomp/boots
ecomp=sqrt(abs(compress**2-ecomp))

end subroutine bootcomp

subroutine check(mc)
use configuration
use measurementdata
implicit none


integer :: i,j,b,op,s1,s2,am,jj(0:2),opt,kkk,kk,nt,kin,btype,i1,i2,a1,a2,ax1,ax2,ay1,ay2,counter,mc,modi,mea
real(8)current,compi,comp2i

 enrg1=enrg1+dfloat(nh)
 enrg2=enrg2+dfloat(nh)**2
 
 gfden=gfden+densgf/dble(nl)
 loc=loc+(1.0d0 + dble(spin))/2.0d0
 rhoi=(1.0d0 + dble(spin))/2.0d0

! structure factor

if(mod(mc,sqmeas)==0)then
 counter=1
 !call DGEMV('N',nqv,L,1.0d0,cost,nqv,rhoi,1,0.0d0,rnq,1)
 !call DGEMV('N',nqv,L,1.0d0,sint,nqv,rhoi,1,0.0d0,inq,1)
 cnq=0.0d0
 snq=0.0d0
 do j=1,nqv
  do i=1,L
   modi=mod(i,nbase)
   if(modi==0)modi=nbase
   cnq(j,modi)=cnq(j,modi)+cost(j,i)*rhoi(i)
   snq(j,modi)=snq(j,modi)+sint(j,i)*rhoi(i)
  end do 
 end do

 sqt=0.0d0
 do kk=1,nqv 
  do i=1,nbase
   do j=1,nbase 
    sqt(kk,i,j)=cnq(kk,i)*cnq(kk,j)+snq(kk,i)*snq(kk,j)
   enddo
  enddo
 end do

  
 rdsq=rdsq+cnq
 idsq=idsq+snq 
 sq=sq+sqt
 
 sqmsteps=sqmsteps+1


 
end if

!spintt=spin for debugging purposes 

NNN=0.0d0
do kkk=1,L
  NNN=NNN+(dble(spin(kkk)))/2.d0
end do



current=NNN
NNN=0.0d0
compi=0.0d0
comp2i=0.0d0


jj(:)=0
kin=0
mea=0


 do i=1,mm
    opt=sm(i)
    op=vtype(i)
    b=sm(i)/2
    btype=mod(b,2) 
    

    ! the usual measurement of superfluid stiffness from S+S- terms
    !**************************************************************
    
    if(op==9.or.op==15) then
      if(btype==1)then
       jj(0)=jj(0)+1
      else
       jj(0)=jj(0)+1
      end if    
    elseif(op==10.or.op==16)then
      if(btype==1)then
       jj(0)=jj(0)-1
      else
       jj(0)=jj(0)-1
      endif
    elseif(op==11.or.op==17)then
       if(btype==1)then
       jj(0)=jj(0)-1
       jj(1)=jj(1)+1
      else
       jj(1)=jj(1)-1
      endif
    elseif(op==12.or.op==18)then
      if(btype==1)then
       jj(0)=jj(0)+1
       jj(1)=jj(1)-1
      else
       jj(1)=jj(1)+1
      endif  
    elseif(op==13.or.op==20)then
      if(btype==1)then
       jj(1)=jj(1)+1
      else
       jj(1)=jj(1)-1
       jj(0)=jj(0)+1
      endif 
    elseif(op==14.or.op==19)then
      if(btype==1)then
       jj(1)=jj(1)-1
      else
       jj(1)=jj(1)+1
       jj(0)=jj(0)-1
      endif      
    end if

    !******************************************************** 

    
    if(op>8)then
      b=opt/2

      ! ** structure factor (perfect statistics (expensive)) **
      !rnq(:)=rnq(:)+(cost( :,bsites(1,b))*( -dble(spin(bsites(1,b))) ) + &
      !               cost( :,bsites(2,b))*(- dble(spin(bsites(2,b))))         )
      
      ! inq(:)=inq(:)+(sint( :,bsites(1,b))*( -dble(spin(bsites(1,b))) ) + &
      !                sint( :,bsites(2,b))*(- dble(spin(bsites(2,b) )))         )      
      !sqt=sqt+rnq**2+inq**2
      !counter=counter+1 
      ! ********************** 

      ! ** This bit does not work anymore with the triangular plaquettes **
      !spin(bsites(1,b))=-spin(bsites(1,b))
      !spin(bsites(2,b))=-spin(bsites(2,b))
      !********************************************************************

      kin=kin+1.0d0
    end if

!    if(sum((1.0d0-dble(spintt))/2.0d0)- sum((1.0d0-dble(spin))/2.0d0)/=0)write(6,*)'wrong' DEBUGGING

    if(mod(mc,fqmeas)==0)then
     !if(op>20)then
     ! call bq(i,op)
     !end if 
     if(op>8.and.op<21)then
       call bqplusminus(i,op) 
     end if        
     mea=1  
    end if 
     
   if(op>20.and.op<27)then
     current=current-2.0d0
   end if
   if(op>26)then
     current=current+2.0d0
   end if 

   compi=compi+(current/dble(L)+0.5d0)
   comp2i=comp2i+(current/dble(L)+0.5d0)**2
   NNN=NNN+current

  ! write(6,*)'NNN,current', NNN,current/dble(L) 

  
 end do

if(mea==1)then
 fqmsteps=fqmsteps+1
 mea=0
end if

NNN=NNN/(dble(L)*dble(mm))
kinetic=kinetic+dble(kin)/dble(L)
den=den+NNN


comp=comp+compi/dble(mm)
comp2=comp2+comp2i/dble(mm)

!sqt=sqt/dble(counter)

!sq=sq+sqt

!stiff=stiff+ dble(jj(0))**2/(2*lx)**2*0 +dble(jj(1))**2/(2*ly)**2 !+dble(jj(2))**2/lz**2

stiff=stiff+ dble(jj(0))**2/(3*lx*ly) +dble(jj(1))**2/(3*ly*lx)


call comptri()
tri=tri+trir
!do i=1,mm
!write(6,*)vtype(i)
!end do
!write(6,*)fourb
!stop
end subroutine check


subroutine comptri()
use configuration
use measurementdata
implicit none
integer(4)i,b1,ii,k
real(8)temp

trir=0.0d0

do i=1,nb

 ii=mod(i,2)

 if(ii==0)then
  temp=0.0d0
  do k=1,3
   temp=temp+rhoi(bsites(k,i))
  end do 
  temp=temp/3
  trir(1)=trir(1)+temp**2
 elseif(ii==1)then
  temp=0.0d0
  do k=1,3
   temp=temp+rhoi(bsites(k,i))
  end do
  temp=temp/3
  trir(2)=trir(2)+temp**2  
 end if
end do

trir=2.0d0*trir/(dble(nb))



end subroutine comptri

subroutine bqplusminus(i,op)
use configuration
use measurementdata
implicit none
integer :: i,b1,b2,btype1,btype2,s(2,2),qv,bn,alph,bet,subalp,subbet
integer :: op,opnext,inext,break,optn,opn,o1,p1,q1,m1,n1

b1=sm(i)/2

!******* find if the next non unit operator next to op *******
if(i==mm)then
 inext=1
 break=1
      !optn=sm(1)
      !bn=sm(1)/2
      !opn=vtype(1)
else
inext=i+1
break=1
      !optn=sm(i+1)
      !bn=sm(i+1)/2
      !opn=vtype(i+1)
end if

loopb: do while(break==1)
 optn=sm(inext)
 bn=sm(inext)/2
 opn=vtype(inext)
 if(opn==0)then
   if(inext==mm)then
   inext=1
  else
   inext=inext+1
  end if
 else
  break=0
 end if
end do loopb

opnext=opn
b2=bn
!********************************************************

  !! hopping vertices
  !do j=9,20
  ! vtex(i,j)=tt/2.0
  !end do
if(opnext<9.or.opnext>20)then
 return
else

 call gsitesplusminus(op,opnext,b1,b2,s)
! o1=1
! p1=2
! m1=2
! n1=7
! if(s(1,1)==o1.and.s(1,2)==p1.and.s(2,1)==m1.and.s(2,2)==n1)then
!   fourb=fourb+dble(nh-1)/(beta)**2
   !write(6,*) b1,b2,op,opnext
! end if
! if(s(1,1)==p1.and.s(1,2)==o1.and.s(2,1)==m1.and.s(2,2)==n1)then
!  fourb=fourb+1!dble(nh-1)/(beta)**2
   !write(6,*) b1,b2,op,opnext
! end if
! if(s(1,1)==o1.and.s(1,2)==p1.and.s(2,1)==n1.and.s(2,2)==m1)then
!  fourb=fourb+1!dble(nh-1)/(beta)**2
   !write(6,*) b1,b2,op,opnext
! end if
! if(s(1,1)==p1.and.s(1,2)==o1.and.s(2,1)==n1.and.s(2,2)==m1)then
!  fourb=fourb+1!dble(nh-1)/(beta)**2
  !write(6,*) b1,b2,op,opnext
! end if

 alph=site_to_bondt(s(1,1),s(1,2),1)
 bet=site_to_bondt(s(2,1),s(2,2),1)
 subalp=site_to_bondt(s(1,1),s(1,2),2)
 subbet=site_to_bondt(s(2,1),s(2,2),2)

 real_bond(alph,bet,site_to_d(subalp,subbet))=real_bond(alph,bet,site_to_d(subalp,subbet))+dble((nh-1))/(beta)**2

! do qv=1,nqv

!  fourq(alph,bet,qv)=fourq(alph,bet,qv)+cos(2.0*pi*( q(qv,1)*( ord(subalp,1)-ord(subbet,1) ) +&
!                                                     q(qv,2)*( ord(subalp,2)-ord(subbet,2) ) ) )*&
!                      dble((nh-1))/(beta)**2/dble(L)

! end do

end if
end subroutine bqplusminus

subroutine gsitesplusminus(opf,opn,btf,btn,s)
use configuration
use measurementdata
integer :: opf,opn,s(2,2),bt,btn,btf,optt,bobo

optt=opf
bt=btf
if(optt==9)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt)
elseif(optt==10)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt) 
elseif(optt==11)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt) 
elseif(optt==12)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==13)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==14)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt) 
elseif(optt==15)then 
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt) 
elseif(optt==16)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt)
elseif(optt==17)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==18)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==19)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==20)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt)
endif
optt=opn
bt=btn
if(optt==9)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==10)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==11)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==12)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==13)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==14)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==15)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==16)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==17)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==18)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==19)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==20)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
endif
end subroutine 



subroutine bq(i,op)
use configuration
use measurementdata
implicit none
integer(4)i,b1,b2,btype1,btype2,s(2,2),qv,bn,alph,bet,subalp,subbet
integer(4)op,opnext,inext,break,optn,opn,o1,p1,q1,m1,n1

b1=sm(i)/2

!******* find if the next non unit operator next to op *******
if(i==mm)then
 inext=1
 break=1
      !optn=sm(1)
      !bn=sm(1)/2
      !opn=vtype(1)
else
inext=i+1
break=1
      !optn=sm(i+1)
      !bn=sm(i+1)/2
      !opn=vtype(i+1)
end if

loopb: do while(break==1)
 optn=sm(inext)
 bn=sm(inext)/2
 opn=vtype(inext)
 if(opn==0)then
  if(inext==mm)then
   inext=1
  else
   inext=inext+1
  end if
 else
  break=0
 end if
end do loopb

opnext=opn
b2=bn
!********************************************************
if(opnext<21)then
 return
else

 call gsites(op,opnext,b1,b2,s)
! o1=1
! p1=2
! m1=2
! n1=7
! if(s(1,1)==o1.and.s(1,2)==p1.and.s(2,1)==m1.and.s(2,2)==n1)then
!   fourb=fourb+dble(nh-1)/(beta)**2
   !write(6,*) b1,b2,op,opnext
! end if
! if(s(1,1)==p1.and.s(1,2)==o1.and.s(2,1)==m1.and.s(2,2)==n1)then
!  fourb=fourb+1!dble(nh-1)/(beta)**2
   !write(6,*) b1,b2,op,opnext
! end if
! if(s(1,1)==o1.and.s(1,2)==p1.and.s(2,1)==n1.and.s(2,2)==m1)then
!  fourb=fourb+1!dble(nh-1)/(beta)**2
   !write(6,*) b1,b2,op,opnext
! end if
! if(s(1,1)==p1.and.s(1,2)==o1.and.s(2,1)==n1.and.s(2,2)==m1)then
!  fourb=fourb+1!dble(nh-1)/(beta)**2
  !write(6,*) b1,b2,op,opnext
! end if
 
 alph=site_to_bondt(s(1,1),s(1,2),1) 
 bet=site_to_bondt(s(2,1),s(2,2),1)
 subalp=site_to_bondt(s(1,1),s(1,2),2)
 subbet=site_to_bondt(s(2,1),s(2,2),2)

 real_bond(alph,bet,site_to_d(subalp,subbet))=real_bond(alph,bet,site_to_d(subalp,subbet))+dble((nh-1))/(beta)**2

! do qv=1,nqv
 
!  fourq(alph,bet,qv)=fourq(alph,bet,qv)+cos(2.0*pi*( q(qv,1)*( ord(subalp,1)-ord(subbet,1) ) +&
!                                                     q(qv,2)*( ord(subalp,2)-ord(subbet,2) ) ) )*&
!                      dble((nh-1))/(beta)**2/dble(L)           

! end do

end if
end subroutine bq

subroutine gsites(opf,opn,btf,btn,s)
use configuration
use measurementdata
integer(4)opf,opn,s(2,2),bt,btn,btf,optt

optt=opf
bt=btf
if(optt==21)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt)
elseif(optt==22)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt) 
elseif(optt==23)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt) 
elseif(optt==24)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt)
elseif(optt==25)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==26)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt) 
elseif(optt==27)then 
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt) 
elseif(optt==28)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==29)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==30)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(2,bt)
elseif(optt==31)then
 s(1,1)=bsites(1,bt)
 s(1,2)=bsites(3,bt)
elseif(optt==32)then
 s(1,1)=bsites(2,bt)
 s(1,2)=bsites(3,bt)
endif
optt=opn
bt=btn
if(optt==21)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==22)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==23)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==24)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==25)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==26)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==27)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==28)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==29)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==30)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(2,bt)
elseif(optt==31)then
 s(2,1)=bsites(1,bt)
 s(2,2)=bsites(3,bt)
elseif(optt==32)then
 s(2,1)=bsites(2,bt)
 s(2,2)=bsites(3,bt)
endif

end subroutine


subroutine coord(i,a1,a2)
use configuration
use measurementdata
implicit none
integer(4)s1,a1,a2,k,kk,i


if(mod(i,3*ly)==0)then
  k=3*ly
 else
  k=mod(i,3*ly)
 end if

 !write(*,*)i, ceiling(dble(k)/dble(3))-1,mod(i,3) !! with ceiling... and mod i,3 one can get the a2 coordinate

 if(mod(i,3)==0)then
  kk=3
 else
  kk=mod(i,3)
 end if


! write(*,*)i, ceiling(dble(k)/dble(3))-1, ceiling(dble(i)/dble(3*ly))-1,kk

 if(kk==1)then
  a2=2*(ceiling(dble(k)/dble(3))-1)
  a1=2*(ceiling(dble(i)/dble(3*ly))-1)
 elseif(kk==2)then
  a2=2*(ceiling(dble(k)/dble(3))-1)
  a1=2*(ceiling(dble(i)/dble(3*ly))-1) +1
 elseif(kk==3)then
  a2=2*(ceiling(dble(k)/dble(3))-1)+1
  a1=2*(ceiling(dble(i)/dble(3*ly))-1)
 end if

end subroutine coord

subroutine cdist()
use configuration
use measurementdata
implicit none
integer(4) i,j,k,d,ci,cit(1),cnt,modi,found
real(8) dist9(9),aa1,aa2,candidate(2),diff,eps


eps=1.0d-12

allocate(distances(nbase*L,dd),cdistances(nbase*L))

! kagome 
ba1(1)=2.0d0*dble(lx)
ba1(2)=0.0d0
ba2(1)=2.0d0*dble(ly)*dcos(pi/3.0d0)
ba2(2)=2.0d0*dble(ly)*dsin(pi/3.0d0)


d=1

do i=1,nbase !,L
 do j=1,L ! ,L
   
   aa1=0.0d0
   aa2=0.0d0 
   dist9(1)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ) )**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=1.0d0
   aa2=0.0d0 
   dist9(2)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  ) 
  
   aa1=1.0d0
   aa2=1.0d0 
   dist9(3)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=0.0d0
   aa2=1.0d0  
   dist9(4)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=-1.0d0
   aa2=1.0d0
   dist9(5)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=-1.0d0
   aa2=0.0d0
   dist9(6)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=-1.0d0
   aa2=-1.0d0
   dist9(7)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
   
   aa1=0.0d0
   aa2=-1.0d0
   dist9(8)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )

   aa1=1.0d0
   aa2=-1.0d0
   dist9(9)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
   
   
   cit=minloc(dist9)
   ci=cit(1)
   

   if(ci==1)then
    aa1=0.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)1,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==2)then
    aa1=1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)2,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==3)then
    aa1=1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)3,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==4)then
    aa1=0.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)4,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==5)then
    aa1=-1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)5,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==6)then
    aa1=-1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)6,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==7)then
    aa1=-1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)7,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==8)then
    aa1=0.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)8,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==9)then
    aa1=1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)9,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   endif 
   !write(6,*) '=================================' 
   !write(6,*)i,j,ci,dist9  
   !write(6,*) 'distance',candidate 
   !write(6,*) '================================'
    
   
    distances(d,:)=candidate
    d=d+1 
   
 end do
end do



indist=nbase*L
!write(6,*)'indist',indist,d-1 
!do i=1,nbase*L
! write(6,*)i,distances(i,1), distances(i,2),sqrt(distances(i,1)**2+ distances(i,2)**2)
!end do

allocate(site_to_d(L,L))

cdistances=0

do i=1,L
 do j=1,L

   modi=mod(i,nbase)
   if(modi==0)modi=nbase
   found=0 
   loopkk: do k=1,L

    aa1=0.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !write(6,*)'00', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk 
    end if  

    aa1=1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'10', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if
  
    aa1=1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'11', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if

    aa1=0.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'01', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if

    aa1=-1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'-11', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:) 
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if

    aa1=-1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'-10', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k) +1
     found=1
     exit loopkk
    end if     
 
    aa1=-1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))

     !write(6,*)'-1-1', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k) +1
     found=1
     exit loopkk
    end if

    aa1=0.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'0-1', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if
    aa1=1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'11', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)=cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if
 
   end do loopkk
  ! write(6,*) 'found', i,j,found
 
 end do
end do  


!write(6,*) 'distances'
!do i=1,L
! do j=1,L
!  write(6,*)i,j,site_to_d(i,j),distances(site_to_d(i,j),:)
! end do
!enddo


!write(6,*)'cdistances'
!do i=1,nbase*L
!write(6,*) i, cdistances(i),distances(i,:)
!end do

end subroutine cdist
subroutine restarts(iseedd,rank)
use configuration
use measurementdata
implicit none
integer lxt,lyt,iseedd,iseeddr,i,nf,rank
real(8)betat,mut,tt,tpt,Jzt,eps

eps=0.0000000001

 if(rank==0)then
  open(unit=12,file='fort.12',form='formatted',status='unknown')
  open(unit=13,file='fort.13',form='formatted',status='unknown')
  !open(unit=12,file='fort.12',form='unformatted',status='old',position='append')
  !open(unit=13,file='fort.13',form='unformatted',status='old')
  rewind(12)
  rewind(13)
  read(13,*)nf
  do i=1,nf
   read(12,*)
  end do
  write(6,*)'bins read in unit 12', nf
 end if 
    
  open(unit=55,file='data.dat',form='unformatted',status='old')
  rewind(55) 
  read(55)kstart,data1,data2,datagf,data2gf,dataloc,dataloc2,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,datareal_bond,&
          data2real_bond,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,lxt,lyt,betat,mut,tt,tpt,Jzt,iseeddr,maxloop,nl
  close(55)       
  kstart=kstart+1
  allocate(tryloop(nl,4))

   ! reading configurations for each specific rank
  WRITE(filename, fmt = '(A3,I0,A4)')'conf',rank,'.dat'
  !print *, trim(filename)
  open(unit=37,file=filename,form='unformatted',status='unknown')
  rewind(37)
  read(37)mm,nh
  allocate(vtypel(mm),vtx(0:6*(mm)-1))
  allocate(sm(mm),vtype(mm))
  read(37)sm,vtype,spin 
 
 
 ! write(6,*) 'when restarting you need a different seed for theQMC',iseedd,iseeddr 
 !stop
  if(lxt.ne.lx)then
   write(6,*) 'check restart parameters lx',lxt,lx
   stop
  end if
  if(lyt.ne.ly)then
   write(6,*) 'check restart parameters ly',lyt,ly
   stop
  end if
  if(abs(betat-beta).ge.eps)then
   write(6,*) 'check restart parameters beta',betat,beta
   stop
  end if
  if(abs(mut-mu).ge.eps)then
   write(6,*) 'check restart parameters mu',mut,mu
   stop
  end if
  if(abs(tt-t).ge.eps)then
   write(6,*) 'check restart parameters t',tt,t
   stop
  end if
  if(abs(tpt-tp).ge.eps)then
   write(6,*) 'check restart parameters tp',tpt,tp
   stop
  end if
  if(abs(Jzt-Jz).ge.eps)then
   write(6,*) 'check restart parameters Jz',Jzt,Jz
   stop
  end if
!  if(iseeddr==iseedd)then
!   write(6,*) 'when restarting you need a different seed for the QMC',iseedd,iseeddr
!   write(6,*) 'code selects one for you'
!   CALL SYSTEM('echo $RANDOM>random.dat')
!   open(unit=100,file='random.dat',form='formatted',status='unknown')
!   read(100,*) iseedd
!   close(100)
!   CALL SYSTEM('rm random.dat')
!  endif
 

end subroutine restarts
