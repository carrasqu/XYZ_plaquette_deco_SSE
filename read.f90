
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
integer :: maxloop ! maximum length of the loop operators
integer :: mloopc ! coefficient for maximum loop length maxloop=mloopc*<n>
integer :: termal ! termalization call sign
integer :: totvisit
integer :: looplength
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
real(8) :: t,t2     ! "hopping" matrix element for S^+S^- +  S^-S^+
real(8) :: tp,tp2   ! "hopping" matrix element for S^+S^+ +  S^-S^-
real(8) :: Jz,Jz2    ! SzSz term (or NN  interaction V for hardcore bosons )
real(8) :: sup   ! amplitude of the superlattice potential is present
real(8) :: vo    ! strength of the optical trap potential vo*(r)**2, r distance from the center of the trap
real(8) :: timeli,timelo,timeme,timew

real(8), allocatable :: epsdis(:) ! espdis(i)=mu+ epsilon_i : disorder configuration plus overall chemical potential
real(8), allocatable :: prob(:,:,:,:) ! probabilities for the loop updates
                       !prob(bond,vertex type(st bond p), entrance leg, exit leg) en principio.  

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

integer, allocatable :: spin(:) ! spin (-1 or 1 ) or (boson 0 or 1 respectively) configuration
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
 real(kind=8),allocatable :: fourq(:,:,:)    ! real part of Bq
 real(kind=8),allocatable :: fourq0(:,:,:)
 real(kind=8),allocatable :: datafourq(:,:,:)      
 real(kind=8),allocatable :: data2fourq(:,:,:) 
 real(kind=8),allocatable :: afourq(:,:,:)    ! real part of Bq
 real(kind=8),allocatable :: afourq2(:,:,:)  

! real(kind=8),allocatable :: fourb    ! real part of Bq
! real(kind=8),allocatable :: fourb0
! real(kind=8),allocatable :: datafourb
! real(kind=8),allocatable :: data2fourb

real(kind=8),allocatable :: real_bond(:,:,:)
real(kind=8),allocatable :: real_bond0(:,:,:) 
real(kind=8),allocatable :: datareal_bond(:,:,:)
real(kind=8),allocatable :: data2real_bond(:,:,:)
real(kind=8),allocatable :: databootreal_bond(:,:,:)

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
real(kind=8),allocatable :: ardsq2(:,:)
real(kind=8),allocatable :: ardsq(:,:)
real(kind=8),allocatable :: databootrdsq(:,:)
complex(8),allocatable :: rhosz(:)
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

 integer :: nbins,msteps,isteps,iseedd,k,m,irstart,mstepsq

 !parallelization variables
 integer :: size,rank,ierr

 real(8) :: drand1,avt,avn, time1,time,avn0,avt0

 !parallel initialization


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

 fqmeas=sqmeas
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
 
 restart=1
 call restarts()
 
 call writeresults(msteps,kstart-1,size,rank,nbins)


 !---------------------------
 call mpi_finalize(ierr)
 stop

 end program hcb_sse
!================================!


subroutine makelattice() 
!------------------------------------------------------------------------------
! Constructs the list of sites bsites(1,b) and bsites(2,b)  bsites(3,b)  in plaquette b
!------------------------------------------------------------------------------

use configuration
use measurementdata
implicit none

integer is,x1,x2,y1,y2,z1,z2,i,bond,j,qv,typeb,wh,sub,k,kk
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
snq(nqv,nbase),ardsq2(nqv,nbase),ardsq(nqv,nbase),databootrdsq(nqv,nbase),rhosz(ndim) )


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

allocate(fourq(6,6,nqv),fourq0(6,6,nqv),datafourq(6,6,nqv),data2fourq(6,6,nqv),afourq(6,6,nqv),afourq2(6,6,nqv))

fourq=0.0d0
fourq0=0.0d0
datafourq=0.0d0
data2fourq=0.0d0

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

do j=1,nqv
 read(27,*)q(j,1),q(j,2)
end do

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

allocate(cf(indist),cf0(indist),datacf(indist),data2cf(indist),databootcf(indist))

allocate(real_bond(6,6,indist),real_bond0(6,6,indist),datareal_bond(6,6,indist),data2real_bond(6,6,indist),&
         databootreal_bond(6,6,indist) )

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




!gfm=0.0d0
!datagfm=0.0d0
!data2gfm=0.0d0



! disorder configuration in the local chemical potential
 call rand_init(disseed)

      
 epsdis=mu


!stop



! organizes the bonds

allocate(bbond(2,2*L))
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
   
end do
write(6,*)ndim,sub-1

!do i=1,L
! do j=1,L
! if(site_to_bondt(i,j)>0)then
! write(6,*) i,j,site_to_bondt(i,j)  
! endif
! end do
!end do

end subroutine makelattice 


!------------------------------------!
 subroutine writeresults(msteps,bins,size,rank,nbins)
!------------------------------------!
 use configuration;use mpi; use measurementdata; implicit none

 integer :: i,msteps,bins,shift,size,rank,ierr,j,nbins,qv,k,ii
 integer(4) stats(MPI_STATUS_SIZE)

 real(8) :: wdata1(8),wdata2(8), nk0,nk01,errnk0, dii,value ! ,wdatgf(L,L),wdat2gf(L,L)
 complex(8)valuec
 character(len=1024) :: filename,str2
 write(6,*) bins

 if(rank==0)then ! rank 0

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


! *** sq **** 

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

! estimates the error in <n_{q}><n_{-q}>

call bootrdsq()

open(unit=18,file='bisq.dat',form='formatted',status='unknown')

write(18,*) 'sum rule structure factor (should be the density (if nqv=ndim!))', &
 (sum(sq(:,1,1))+sum(sq(:,3,3))+sum(sq(:,3,3)))/dble(nbase)
do i=1,nqv
 write(18,256) q(i,1),q(i,2), idsq(i,1),idsq(i,2),idsq(i,3)
end do

close(18)

open(unit=17,file='bootsq.dat',form='formatted',status='unknown')

do i=1,nqv
 write(17,256) q(i,1),q(i,2),sq(i,1,1)+sq(i,2,2)+sq(i,3,3), sq0(i,1,1)+sq0(i,2,2)+sq0(i,3,3), &
  ardsq(i,1)+ ardsq(i,2)+ardsq(i,3), ardsq2(i,1)+ ardsq2(i,2)+ardsq2(i,3)
end do
close(17)


open(unit=17,file='bootsq11.dat',form='formatted',status='unknown')

do i=1,nqv
 write(17,256) q(i,1),q(i,2),sq(i,1,1), sq0(i,1,1), ardsq(i,1), ardsq2(i,1)
end do
close(17)


open(unit=17,file='bootsq22.dat',form='formatted',status='unknown')

do i=1,nqv
 write(17,256) q(i,1),q(i,2),sq(i,2,2), sq0(i,2,2), ardsq(i,2), ardsq2(i,2)
end do
close(17)

open(unit=17,file='bootsq33.dat',form='formatted',status='unknown')

do i=1,nqv
 write(17,256) q(i,1),q(i,2),sq(i,3,3), sq0(i,3,3), ardsq(i,3), ardsq2(i,3)
end do
close(17)

256 format(6ES16.8)

!************


! real space diagonal correlation function ****************************

call realspaceDiag()

!**********************************************************************


!************** momentum distribution from one-body density matrix  *******************
 call nok( )

! write(6,*)'sum rule with a basis',  ((sum(rnk(:,1,1))) + (sum(rnk(:,2,2))) + (sum(rnk(:,3,3)))) ! should be the number of "bosons" 

 open(unit=20,file='rn_k.dat',form='formatted',status='unknown')
 open(unit=21,file='in_k.dat',form='formatted',status='unknown')

  write(21,*)'sum rule is the number of bosons if nqv=ndim',  ((sum(rnk(:,1,1))) + (sum(rnk(:,2,2))) + (sum(rnk(:,3,3))))
 do i=1,nqv

  write(20,255)q(i,1),q(i,2),  rnk(i,1,1)+rnk(i,2,2)+rnk(i,3,3) , rnk2(i,1,1)+rnk2(i,2,2)+rnk2(i,3,3) ! rnk(i,1,2),rnk(i,1,3),rnk(i,2,1),rnk(i,2,2),rnk(i,2,3),rnk(i,3,1),rnk(i,3,2),rnk(i,3,3)
  write(21,255)q(i,1),q(i,2),  ink(i,1,1)+ink(i,2,2)+ink(i,3,3) , ink2(i,1,1)+ink2(i,2,2)+ink2(i,3,3)
  
 end do 
do i=1,3

 write(96,*)(rnk(1,i,j),j=1,3)

enddo

write(96,*)
do i=1,3

 write(96,*)(rnk2(1,i,j),j=1,3)

enddo

 close(20)
 close(21) 

255 format(6ES16.8)


! real space single particle density matrix
call real_space_nij()

!************************************************************

!*** compressibility bootcomp()
call bootcomp(wdata1,wdata2)

!************************************************************

! ********* four point correlation function *****************

real_bond=datareal_bond/dble(bins)
real_bond0=data2real_bond/dble(bins)
real_bond0=sqrt(abs(real_bond0-real_bond**2)/dble(bins))
!nqv=0
!do i=1,indist
! do j=1,6
!  do k=1,6
!   write(6,*)j,k,i ,real_bond(j,k,i)
!   if(real_bond(j,k,i)>0)nqv=nqv+1
!  enddo
! enddo
!end do
!write(6,*)"cuantos", nqv
!stop
call fftfour()

!***real space bond bond***

call real_space_bond()

!**************************

open(unit=20,file='fourq.dat',form='formatted',status='unknown')
do i=1,nqv
 write(20,255)q(i,1),q(i,2),(fourq(1,1,i)+fourq(2,2,i)+fourq(3,3,i)+fourq(4,4,i)+fourq(5,5,i)+fourq(6,6,i))*(2.0/t)**2,&
                            (fourq0(1,1,i)+fourq0(2,2,i)+fourq0(3,3,i)+fourq0(4,4,i)+fourq0(5,5,i)+fourq0(6,6,i))*(2.0/t)**2
end do
close(20)


i=2
!write(filen,'("idensit_t_",E10.4)')time
 !write (filename, "(A5,I2)") "fourq", 


do i=1,6
 do j=1,6
  write(filename,'("fourq",I1)')i  
  filename=trim(filename)
  write(str2,"(I1)")j
  str2=trim(str2)
  filename=trim(filename)//str2
  filename=trim(filename)
  filename=trim(filename)//".dat" 
  filename=trim(filename)
!  write(6,*) i,j, 'file=',trim(filename)
  open(unit=20,file=trim(filename),form='formatted',status='unknown')
  do ii=1,nqv
   write(20,255)q(ii,1),q(ii,2),(fourq(i,j,ii))*(2.0/t)**2,&
                            (fourq0(i,j,ii))*(2.0/t)**2
  end do
  close(20)
 enddo
enddo




!************************************************************
  
  open(10,file='xBootresults.txt',status='replace')
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
 ! write(10,10)'  fourb     : ',fourb*(2.0/tp)**2,fourb0*(2.0/tp)**2
  write(10,*)' ========================================='
  10 format(1x,a,2f14.8)
  close(10)

 
 if(densiono==1)then
  open(unit=18,file='localdensity.dat',form='formatted',status='unknown')
  loc=sqrt(abs(dataloc2/bins -(dataloc/bins)**2 )/bins)
  do i=1,L 
    write(18,*) ord(i,1),ord(i,2), dataloc(i)/dble(bins),loc(i) 
  !write(18)loc
  enddo
  close(18)
 end if


 end if ! rank 0 

  
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
subroutine nok()
use configuration
use measurementdata
implicit none
integer(kind=4)qq,ll,i,j,k,modj,modk
real(kind=8) u,v,x,y
real(kind=8) drand1,fact,dv(dd),mag,eps
complex(8)ii
ii=cmplx(0,1.0d0)
 
 eps=0.00000000001
 rnk2=0.0d0
 arnk=0.0d0
 ink2=0.0d0
 aink=0.0d0


do qq=1,boots

 databootgf=0.0d0
 do i=1,indist
  u=drand1()
  v=drand1()
  x=dsqrt(-2.0d0*log(u))*dcos(2*pi*v)
  y=dsqrt(-2.0d0*log(u))*dsin(2*pi*v)
  databootgf(i)=gf(i)+x*gf0(i)
 enddo

 rnk=0.0d0
 ink=0.0d0

   
  do i=1,nqv
   modj=1 
   do j=1,indist
    modk=mod(j,nbase)
    if(modk==0)modk=nbase

    mag=sqrt(distances(j,1)**2+distances(j,2)**2)
    
    if(mag<eps)then
     fact=1.0d0
    else
     fact=0.5
    end if
      
    dv(:)=distances(j,:) 
   !write(6,*)'q',i,j,modj,modk 
    rnk(i,modj,modk)=rnk(i,modj,modk)+ cos(2.0d0*pi*( dv(1)*q(i,1)+&
         dv(2)*q(i,2) ))*(databootgf(j)/dble(cdistances(j)))*dble(L)*fact
        
        !write(6,*) exp(2.0d0*pi*ii*(dv(1))*q(i,1)+&
        ! (dv(2))*q(i,2)) 
   ink(i,modj,modk)=ink(i,modj,modk)+ sin(2.0d0*pi*( (dv(1))*q(i,1)+&
         (dv(2))*q(i,2) ) )*(databootgf(j)/dble(cdistances(j)))*dble(L)*fact

         ! /dble(cdistances(site_to_d(j,k)) Because the rho_ij has been summed over all cdistances equivalent distances
         ! *L because one measures rho_ij/L during the loopupdate
         ! *fact: Because the actual content of rho_ij= <b_i^\dag b_j> + <b_j^\dag b_i> i\=j
         ! /ndim comes form the definition of the fourier transform
    if(mod(j,L)==0)then
    modj=modj+1
   endif         
  
   end do 
  end do
 rnk2=rnk2+rnk**2
 arnk=arnk+rnk 

 ink2=ink2+ink**2
 aink=aink+ink  

end do

rnk2=rnk2/dble(boots)
ink2=ink2/dble(boots)

arnk=arnk/dble(boots)
aink=aink/dble(boots)

rnk2=sqrt(abs(arnk**2-rnk2))
ink2=sqrt(abs(aink**2-ink2))

rnk=arnk
ink=aink


end subroutine nok


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


write(6,*)'cdistances'
do i=1,nbase*L
write(6,*) i, cdistances(i),distances(i,:)
end do

end subroutine cdist
subroutine restarts()
use configuration
use measurementdata
implicit none
integer lxt,lyt,iseedd
real(8)betat,mut,tt,tpt,Jzt,eps

eps=0.0000000001
if(restart==0)then
  kstart=1
 elseif(restart>0)then

  open(unit=55,file='data.dat',form='unformatted',status='old')
  rewind(55) 
  read(55)kstart,data1,data2,datagf,data2gf,dataloc,dataloc2,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,datareal_bond,&
         data2real_bond,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,lxt,lyt,betat,mut,tt,tpt,Jzt
         kstart=kstart+1
 
!  rewind(55)
!write(55)bins,data1,data2,datagf,data2gf,dataloc,dataloc2,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,datareal_bond,&
!         data2real_bond,datasq,data2sq,datardsq,data2rdsq,dataidsq,data2idsq,lx,ly,beta,mu,t,tp,Jz,iseedd,maxloop,nl
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
 close(55)
 end if
end subroutine restarts

subroutine fftfour()
use configuration
use measurementdata
implicit none
integer(kind=4)qq,ll,i,j,k,modj,modk,alp,bet,cc
real(kind=8) u,v,x,y,dv(dd)
real(kind=8) drand1,fact

afourq2=0.0d0
afourq=0.0d0

do qq=1,boots

 databootreal_bond=0.0d0
 do i=1,indist
  do alp=1,6
   do bet=1,6
  
    u=drand1()
    v=drand1()
    x=dsqrt(-2.0d0*log(u))*dcos(2*pi*v)
    y=dsqrt(-2.0d0*log(u))*dsin(2*pi*v)
    databootreal_bond(alp,bet,i)=real_bond(alp,bet,i)+x*real_bond0(alp,bet,i)
   enddo
  enddo
 enddo

 fourq=0.0d0
  do i=1,nqv
   !cc=0
   modj=1  
   do j=1,indist !L
    !do k=1,L
     
     dv(:)=distances(j,:)
     modk=mod(j,nbase)
     if(modk==0)modk=nbase  
     !write(6,*)i,j,modj,modk
 
     if(modj==1.and.modk==1)then
      do alp=1,6
       do bet=1,6
         
        fourq(alp,bet,i)=fourq(alp,bet,i)+databootreal_bond(alp,bet,j)*&
                                        cos(2.0d0*pi*( q(i,1)*( dv(1) ) +&
                                                       q(i,2)*( dv(2) ) ) )/dble(cdistances(j)*nbase)
       cc=cc+1
       end do 
      end do
     end if
     if(mod(j,L)==0)then
      modj=modj+1
     endif
    !end do
   end do
   !write(6,*) 'cc',cc,6**2*ndim 
  end do

 afourq2=afourq2+fourq**2
 afourq=afourq+fourq
end do

afourq2=afourq2/dble(boots)
afourq=afourq/dble(boots)


afourq2=sqrt(abs(afourq**2-afourq2))
fourq=afourq
fourq0=afourq2

end subroutine fftfour

subroutine bootrdsq()
use configuration
use measurementdata
implicit none
integer(kind=4)qq,ll,i,j,k,modj,modk,alp,bet
real(kind=8) u,v,x,y
real(kind=8) drand1,fact

ardsq2=0.0d0
ardsq=0.0d0

do qq=1,boots

 databootrdsq=0.0d0
 do i=1,nqv
  do alp=1,nbase
    u=drand1()
    v=drand1()
    x=dsqrt(-2.0d0*log(u))*dcos(2*pi*v)
    y=dsqrt(-2.0d0*log(u))*dsin(2*pi*v)
    databootrdsq(i,alp)=rdsq(i,alp) +x*rdsq0(i,alp)
  enddo
 enddo

 ardsq=ardsq+databootrdsq*databootrdsq
 ardsq2=ardsq2+(databootrdsq*databootrdsq)**2

end do

ardsq=ardsq/dble(boots)
ardsq2=ardsq2/dble(boots)
ardsq2=sqrt(abs(ardsq**2-ardsq2 ))

end subroutine bootrdsq

subroutine real_space_nij()
use configuration
use measurementdata
implicit none
integer(kind=4)qq,ll,i,j,k,modj,modk
real(kind=8) u,v,x,y
real(kind=8) drand1,fact,dv(dd),mag,eps

eps=0.00000000001

open(unit=108,file='realnij.dat',form='formatted',status='unknown')
open(unit=109,file='T1realnij.dat',form='formatted',status='unknown')
open(unit=110,file='T2realnij.dat',form='formatted',status='unknown')

do j=1,indist
 modk=mod(j,nbase)
 if(modk==0)modk=nbase

 mag=sqrt(distances(j,1)**2+distances(j,2)**2)
 write(6,*) j, distances(j,:)
 if(mag<eps)then
  fact=1.0d0
 else
  fact=0.5
 end if

 write(108,*) mag,databootgf(j)/dble(cdistances(j))*dble(L)*fact 

end do

do j=1,indist
 modk=mod(j,nbase)
 mag=sqrt(distances(j,1)**2+distances(j,2)**2)
 if(mag<eps)then
  fact=1.0d0
 else
  fact=0.5
 end if
 if(modk==0)modk=nbase
 if(modk==1)then
  if(abs(distances(j,1))<eps)then
   write(109,*) distances(j,2),gf(j)/dble(cdistances(j))*dble(L)*fact,gf0(j)/dble(cdistances(j))*dble(L)*fact  
  end if
  if(abs(distances(j,2))<eps)then
   write(110,*) distances(j,1),gf(j)/dble(cdistances(j))*dble(L)*fact, gf0(j)/dble(cdistances(j))*dble(L)*fact
  end if
 end if

end do



close(108)
close(109)
close(110)

end subroutine real_space_nij

subroutine realspaceDiag()
use configuration
use measurementdata
implicit none
integer(4)i,j,k,qv,dist,counter,modj
real(8)eps
complex(8)iii
iii=cmplx(0.0d0,1.0d0)
eps=0.00000000001

open(unit=108,file='Diag_real_sz_r_11.dat',form='formatted',status='unknown')
write(6,*) 'aqui'


counter=1
do i=1,1
 do j=1,L
  modj=mod(j,nbase)
  dist=i*j
  if(modj==1)then
   rhosz(counter)=0.0d0
   write(6,*)'counter ndim',counter,ndim
   do qv=1,ndim
    rhosz(counter)=rhosz(counter)+cos( 2.0d0*pi*( q(qv,1)*distances(dist,1)+q(qv,2)*distances(dist,2) ))*&
    (sq(qv,1,1)-ardsq(qv,1))+&
       iii*sin( 2.0d0*pi*( q(qv,1)*distances(dist,1)+q(qv,2)*distances(dist,2) ))*&
      (sq(qv,1,1)-ardsq(qv,1))
   end do
   if(abs(distances(dist,2) )<eps)then
    write(108,*) distances(dist,1),distances(dist,2),real(rhosz(counter)),aimag(rhosz(counter))
   end if
   counter=counter+1
  end if
 end do
end do
write(6,*)'despues'
close(108)

end subroutine realspaceDiag

subroutine real_space_bond()
use configuration
use measurementdata
implicit none
integer(kind=4)qq,ll,i,j,k,modj,modk
real(kind=8) u,v,x,y
real(kind=8) drand1,fact,dv(dd),mag,eps

eps=0.00000000001

open(unit=110,file='T1realbond11.dat',form='formatted',status='unknown')
open(unit=111,file='T1realbond22.dat',form='formatted',status='unknown')
open(unit=112,file='T1realbond33.dat',form='formatted',status='unknown')
open(unit=113,file='T1realbond44.dat',form='formatted',status='unknown')
open(unit=114,file='T1realbond55.dat',form='formatted',status='unknown')
open(unit=115,file='T1realbond66.dat',form='formatted',status='unknown')
open(unit=116,file='T1realbond12.dat',form='formatted',status='unknown')
open(unit=117,file='T1realbond13.dat',form='formatted',status='unknown')
open(unit=118,file='T1realbond14.dat',form='formatted',status='unknown')
open(unit=119,file='T1realbond15.dat',form='formatted',status='unknown')
open(unit=120,file='T1realbond16.dat',form='formatted',status='unknown')



do j=1,indist
 modk=mod(j,nbase)
 mag=sqrt(distances(j,1)**2+distances(j,2)**2)
 if(mag<eps)then
  fact=1.0d0
 else
  fact=0.5
 end if
 if(modk==0)modk=nbase
 if(modk==1)then
  if(abs(distances(j,1))<eps)then
 !  write(109,*) distances(j,2),gf(j)/dble(cdistances(j))*dble(L)*fact,gf0(j)/dble(cdistances(j))*dble(L)*fact  
  end if
  if(abs(distances(j,2))<eps)then
   
   if(real_bond0(1,1,j)>eps)write(110,*) distances(j,1),real_bond(1,1,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(1,1,j)/dble(cdistances(j)*nbase)
   if(real_bond0(2,2,j)>eps)write(111,*) distances(j,1),real_bond(2,2,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(2,2,j)/dble(cdistances(j)*nbase)
   if(real_bond0(3,3,j)>eps)write(112,*) distances(j,1),real_bond(3,3,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(3,3,j)/dble(cdistances(j)*nbase)
   if(real_bond0(4,4,j)>eps)write(113,*) distances(j,1),real_bond(4,4,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(4,4,j)/dble(cdistances(j)*nbase)
   if(real_bond0(5,5,j)>eps)write(114,*) distances(j,1),real_bond(5,5,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(5,5,j)/dble(cdistances(j)*nbase)
   if(real_bond0(6,6,j)>eps)write(115,*) distances(j,1),real_bond(6,6,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(6,6,j)/dble(cdistances(j)*nbase)
   if(real_bond0(1,2,j)>eps)write(116,*) distances(j,1),real_bond(1,2,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(1,2,j)/dble(cdistances(j)*nbase)
   if(real_bond0(1,3,j)>eps)write(117,*) distances(j,1),real_bond(1,3,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(1,3,j)/dble(cdistances(j)*nbase)
   if(real_bond0(1,4,j)>eps)write(118,*) distances(j,1),real_bond(1,4,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(1,4,j)/dble(cdistances(j)*nbase)
   if(real_bond0(1,5,j)>eps)write(119,*) distances(j,1),real_bond(1,5,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(1,5,j)/dble(cdistances(j)*nbase)
   if(real_bond0(1,6,j)>eps)write(120,*) distances(j,1),real_bond(1,6,j)/dble(cdistances(j)*nbase)&
     ,real_bond0(1,6,j)/dble(cdistances(j)*nbase)
  end if
 end if

end do



close(110)
close(111)
close(112)
close(113)
close(114)
close(115)
close(116)
close(117)
close(118)
close(119)
close(120)

end subroutine real_space_bond

 
 
