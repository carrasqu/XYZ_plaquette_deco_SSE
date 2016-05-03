ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c written by Federico Becca
c  Given two primitive vectors and the two translations defining the box
c  the program finds:
c  1) all the (rigid) spatial symmetries of the infinite lattice: 
c  2) the coordinates of all the lattice points and the table of distances
c  3) all the point symmetries preserved by the chosen cluster
c
c     Each 2x2 symmetry matrix is defined by an angle theta:
c
c     Rotations:                                
c                   cos(theta)      sin(theta)
c                  -sin(theta)      cos(theta)
c
c     Reflections:
c                   cos(theta)      sin(theta)
c                   sin(theta)     -cos(theta) 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program lattice
      implicit none
      INTEGER(4) nbase,nsite,ntotalsite,nsimmax,nmom
      INTEGER(4) i,j,k,irot,iref,isg,itest,i60,numax,nu12,nrot,nref,ndim
     >,n1,n2,id2,iflag,itest1,itest2,idt2,ii,it,jj,jt,l,m,nind
     >,nindm,nindms,irotc,irefc,iflagi,maxmulti,maxmultis
     >,ib1,ib2,iflagibc,ik1,ik2,is,ikm1,ikm2,ikm1b,ikm2b,ikm,ik
     >,mm,l1,l2,k1,k2,n1min,n1max,n2min,n2max,iflagref,jshift
     >,idists,ibroken,icount,input1,jsav,mt
     >,iclass,index,iii
      CHARACTER(5) oper1,oper2,oper3,oper4
      REAL(8) a1,a2,phi,rho,pi,small,axx,ayy
     >,ctgp,st,sinp,t,ct,test11,test22,test12,test21,beta,sb,test
     >,test11b,test21b,xi,yi,xij,yij,d2,xt,yt,d2t,cth,xj,yj,det,sth,yb1
     >,dx,dy,testb1,testb2,xb1,dij,d,test1,test2
     >,t1t2,t1,t2,rr,rt1,rt2,rt1x,rt2x,rrx,cost1,cost2
     >,phiKINpp_r,phiKINpp_i,phiKINmm_r,phiKINmm_i
     >,phiKINss_r,phiKINss_i,phiKINaa_r,phiKINaa_i
     >,uno,eig1,eig2,qx,qy,disc,disc1,disc2,check,checkx

      INTEGER(4), dimension(:), allocatable :: i1,i2
      INTEGER(4), dimension(:), allocatable :: isymc
      INTEGER(4), dimension(:), allocatable :: imulti,imultis
      INTEGER(4), dimension(:), allocatable :: imultis_scra
      INTEGER(4), dimension(:), allocatable :: imark
      INTEGER(4), dimension(2,2) :: enne
      INTEGER(4), dimension(:,:), allocatable :: inewb0,inewb1,inewb2
      INTEGER(4), dimension(:,:), allocatable :: boundary
      INTEGER(4), dimension(:,:), allocatable :: isymt
      INTEGER(4), dimension(:,:), allocatable :: imark2
      INTEGER(4), dimension(:,:,:), allocatable :: sym
      INTEGER(4), dimension(:,:,:), allocatable :: ivic,ivict,ivicns
      INTEGER(4), dimension(:,:,:), allocatable :: ivicns_scra
      REAL(8), dimension(:), allocatable :: xb,yb
      REAL(8), dimension(:), allocatable :: point
      REAL(8), dimension(:), allocatable :: dindip
      REAL(8), dimension(2,2) :: a,b,enne1,amat
      REAL(8), dimension(:,:), allocatable :: x,y,ord
      REAL(8), dimension(:,:), allocatable :: dist

      open(unit=90,file='geometry.d',status='old')               !input file

      pi=dacos(-1.d0)
      small=1.d-5

c     ENTRIES      
c     a1=modulus of the first primitive vector
c     a2=modulus of the second primitive vector
c     phi=angle between the two primitive vectors
c
c     enne(i,j) s a 2x2 integer matrix defining the two translations which 
c     identify the cluster:
c     T_i = sum_j enne(i,j)*a_j   
c     We assume that the leading components are enne(1,1)>0 and enne(2,2)>0 

ccccccccccccccccccccccccc  READING PART ccccccccccccccccccccccccc
      read(90,*) oper1,a1,oper2,a2,phi
      if(oper1.eq.'sqrt') then
       a1=dsqrt(a1)
      elseif(oper1.ne.'') then
       write(6,*)'wrong operator'
       stop
      endif
      if(oper2.eq.'sqrt') then
       a2=dsqrt(a2)
      elseif(oper2.ne.'') then
       write(6,*)'wrong operator'
       stop
      endif
      phi=phi*pi/180.d0

      read(90,*) nbase                      !number of elements of the basis

      ALLOCATE(xb(nbase))
      ALLOCATE(yb(nbase))

      if(nbase.ne.1)then
       do i=1,nbase
        read(90,*) oper3,axx,oper4,ayy
        if(oper3.eq.'sqrt') axx=dsqrt(axx)
        if(oper4.eq.'sqrt') ayy=dsqrt(ayy)
        xb(i)=axx
        yb(i)=ayy
       enddo
      else
       xb(1)=0.d0
       yb(1)=0.d0
      endif

      read(90,*) enne(1,1),enne(1,2),enne(2,1),enne(2,2)
      read(90,*) ib1,ib2

      if(enne(1,1).le.0.or.enne(2,2).le.0)then
       write(6,*)'box positioning'
       stop
      endif

      nsite=abs(enne(1,1)*enne(2,2)-enne(1,2)*enne(2,1))

c     Primitive vectors
      a(1,1)=a1                !a_1 x component
      a(2,1)=0.d0              !a_1 y component
      a(1,2)=a2*dcos(phi)      !a_2 x component
      a(2,2)=a2*dsin(phi)      !a_2 y component
c     Reciprocal vectors
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      b(1,1)=a(2,2)/det
      b(2,1)=-a(1,2)/det
      b(1,2)=-a(2,1)/det
      b(2,2)=a(1,1)/det
ccccccccccccccccccccccccc  END READING PART ccccccccccccccccccccccccc

ccccccccccccccccccccccccc  INFINITE LATTICE PART ccccccccccccccccccccccccc
      rho=a1/a2

      nsimmax=20
      ALLOCATE(sym(2,2,nsimmax))
      ALLOCATE(point(nsimmax))

      do i=1,nsimmax
       point(i)=0.d0
       do j=1,2
        do k=1,2
         sym(j,k,i)=0
        enddo
       enddo
      enddo

      ctgp=1.d0/dtan(phi)
      sinp=dsin(phi)

c  ROTATIONS
c  point(irot)= rotation angle in units of 2*pi
c  sym(i,j,k)=matrix nu_{ij} defining the symmetry k

c  Identity
 
      irot=1
      point(1)=0.d0
      sym(1,1,1)=1
      sym(2,2,1)=1

c  Inversion
 
      irot=2
      point(2)=0.5d0
      sym(1,1,2)=-1
      sym(2,2,2)=-1

c  90 degrees

      do isg=-1,1,2
       t=isg*pi/2.d0
       ct=dcos(t)
       st=dsin(t)
       test11=ct+st*ctgp 
       if(itest(test11).eq.0) then
        test22=ct-st*ctgp 
        if(itest(test22).eq.0) then
         test12=-rho*st/sinp
         if(itest(test12).eq.0) then
          test21=st/sinp/rho
          if(itest(test21).eq.0) then
           irot=irot+1
           point(irot)=isg*0.25d0 
           sym(1,1,irot)=nint(test11)
           sym(1,2,irot)=nint(test12)
           sym(2,1,irot)=nint(test21)
           sym(2,2,irot)=nint(test22)
          endif
         endif
        endif
       endif
      enddo

c  60 degrees

      do isg=-1,1,2
       do i60=1,2
        t=isg*i60*pi/3.d0
        ct=dcos(t)
        st=dsin(t)
        test11=ct+st*ctgp 
        if(itest(test11).eq.0) then
         test22=ct-st*ctgp 
         if(itest(test22).eq.0) then
          test12=-rho*st/sinp
          if(itest(test12).eq.0) then
           test21=st/sinp/rho
           if(itest(test21).eq.0) then
            irot=irot+1
            point(irot)=dfloat(isg*i60)/6.d0
            sym(1,1,irot)=nint(test11)
            sym(1,2,irot)=nint(test12)
            sym(2,1,irot)=nint(test21)
            sym(2,2,irot)=nint(test22)
           endif
          endif
         endif
        endif
       enddo
      enddo

      nrot=irot              !total number of rotations

c  REFLECTIONS

      numax=nint(abs(rho/sinp))
      iref=0

      do nu12=-numax,numax
       sb=dfloat(nu12)*sinp/rho
       test=abs(sb)-1.d0
       if(test.le.0.d0) then
        beta=dasin(sb)
        test11=-dsin(beta-phi)/sinp
        if(itest(test11).eq.0) then
         test21=-dsin(beta-2.d0*phi)/sinp/rho
         if(itest(test21).eq.0) then
          iref=iref+1
          point(nrot+iref)=beta/(2.d0*pi)
          sym(1,1,nrot+iref)=nint(test11)
          sym(1,2,nrot+iref)=nu12
          sym(2,1,nrot+iref)=nint(test21)
          sym(2,2,nrot+iref)=-nint(test11)
         endif
        endif
        beta=pi-beta
        test11b=-dsin(beta-phi)/sinp
        if(itest(test11b).eq.0) then
         test21b=-dsin(beta-2.d0*phi)/sinp/rho
         if(itest(test21b).eq.0) then
          if(abs(test11-test11b).gt.small.or.
     &       abs(test21-test21b).gt.small) then
           iref=iref+1
           point(nrot+iref)=beta/(2.d0*pi)
           sym(1,1,nrot+iref)=nint(test11b)
           sym(1,2,nrot+iref)=nu12
           sym(2,1,nrot+iref)=nint(test21b)
           sym(2,2,nrot+iref)=-nint(test11b)
          endif
         endif
        endif
       endif
      enddo

      nref=iref
      write(6,*)
      write(6,333) nrot,nref
      write(6,*)
333   format(1x,'NROT=',i4,5x,'NREF=',i4)

      if(nrot+nref.gt.nsimmax)then
       write(6,*)'nsimmax too small',nrot+nref
       stop
      endif
ccccccccccccccccccccccccc  END INFINITE LATTICE PART ccccccccccccccccccccccccc

ccccccccccccccccccccccccc  CLUSTER PART ccccccccccccccccccccccccc
      ndim=0

      ALLOCATE(x(nsite,nbase))
      ALLOCATE(y(nsite,nbase))
      ALLOCATE(i1(nsite))
      ALLOCATE(i2(nsite))

c     do n1=0,enne(1,1)+abs(enne(2,1))
c      do n2=0,enne(2,2)+abs(enne(1,2))
c       id2=n1**2+n2**2
c       iflag=0
c       do isg=-1,1,2
c        itest1=n1-enne(1,1)*isg 
c        itest2=n2-enne(1,2)*isg 
c        idt2=itest1**2+itest2**2
c        if(itest1.ge.0.and.itest2.ge.0.and.idt2.lt.id2) iflag=1
c        itest1=n1-enne(2,1)*isg 
c        itest2=n2-enne(2,2)*isg 
c        idt2=itest1**2+itest2**2
c        if(itest1.ge.0.and.itest2.ge.0.and.idt2.lt.id2) iflag=1
c        itest1=n1-(enne(1,1)+enne(2,1))*isg 
c        itest2=n2-(enne(1,2)+enne(2,2))*isg 
c        idt2=itest1**2+itest2**2
c        if(itest1.ge.0.and.itest2.ge.0.and.idt2.lt.id2) iflag=1
c        itest1=n1-(enne(1,1)-enne(2,1))*isg 
c        itest2=n2-(enne(1,2)-enne(2,2))*isg 
c        idt2=itest1**2+itest2**2
c        if(itest1.ge.0.and.itest2.ge.0.and.idt2.lt.id2) iflag=1
c       enddo
c       if(iflag.eq.0)then
c        ndim=ndim+1
c        do j=1,nbase
c  real space coordinates
c         x(ndim,j)=xb(j)+dfloat(n1)*a(1,1)+dfloat(n2)*a(1,2)
c         y(ndim,j)=yb(j)+dfloat(n1)*a(2,1)+dfloat(n2)*a(2,2)
c        enddo
c  components along the primitive vectors
c        i1(ndim)=n1
c        i2(ndim)=n2
c       endif
c      enddo
c     enddo

      n1min=min(0,enne(1,1),enne(2,1),enne(1,1)+enne(2,1))
      n1max=max(0,enne(1,1),enne(2,1),enne(1,1)+enne(2,1))
      n2min=min(0,enne(1,2),enne(2,2),enne(1,2)+enne(2,2))
      n2max=max(0,enne(1,2),enne(2,2),enne(1,2)+enne(2,2))

      t1t2=0.d0
      t1=0.d0
      t2=0.d0
      do i=1,2
       do j=1,2
        amat(i,j)=0.d0
        do k=1,2
         amat(i,j)=amat(i,j)+a(k,i)*a(k,j)       !a_i scalar a_j
        enddo
        t1t2=t1t2+enne(1,i)*enne(2,j)*amat(i,j)
        t1=t1+enne(1,i)*enne(1,j)*amat(i,j)
        t2=t2+enne(2,i)*enne(2,j)*amat(i,j)
       enddo
      enddo

      do n1=n1min,n1max
       do n2=n2min,n2max
        iflag=0
        if(n1.eq.enne(1,1).and.n2.eq.enne(1,2)) iflag=1     !boundary
        if(n1.eq.enne(2,1).and.n2.eq.enne(2,2)) iflag=1     !boundary
        if(n1.eq.(enne(1,1)+enne(2,1)).and.n2.eq.(enne(1,2)+enne(2,2)))
     >                                          iflag=1     !boundary
        rr=n1**2*amat(1,1)+2.d0*n1*n2*amat(1,2)+n2**2*amat(2,2) 
        rt1=n1*enne(1,1)*amat(1,1)+n1*enne(1,2)*amat(1,2)+
     >      n2*enne(1,1)*amat(2,1)+n2*enne(1,2)*amat(2,2)
        rt2=n1*enne(2,1)*amat(1,1)+n1*enne(2,2)*amat(1,2)+
     >      n2*enne(2,1)*amat(2,1)+n2*enne(2,2)*amat(2,2)
        disc1=(t1*rr-rt1**2)
        disc2=(t2*rr-rt2**2)
        disc=(t1*rr-rt1**2)*(t2*rr-rt2**2)
        if(abs(disc1).lt.small.or.abs(disc2).lt.small)then
         check=t1t2*rr-rt1*rt2
        else
         check=t1t2*rr-rt1*rt2+dsqrt(disc)
        endif
        if(abs(check).gt.small)iflag=1
        rt1x=t1+t1t2-rt1
        rt2x=t2+t1t2-rt2
        rrx=rr+t1+t2-2.d0*(rt1+rt2-t1t2)
        disc1=rrx
        disc2=t1
        disc=rrx*t1
        if(abs(disc1).lt.small.or.abs(disc2).lt.small)then
         cost1=rt1x
        else
         cost1=rt1x-dsqrt(disc)
        endif
        disc1=rrx
        disc2=t2
        disc=rrx*t2
        if(abs(disc1).lt.small.or.abs(disc2).lt.small)then
         cost2=rt2x
        else
         cost2=rt2x-dsqrt(disc)
        endif
        if(abs(cost1).lt.small)iflag=1
        if(abs(cost2).lt.small)iflag=1
        disc1=(t1*rrx-rt1x**2)
        disc2=(t2*rrx-rt2x**2)
        disc=(t1*rrx-rt1x**2)*(t2*rrx-rt2x**2)
        if(abs(disc1).lt.small.or.abs(disc2).lt.small)then
         checkx=t1t2*rrx-rt1x*rt2x
        else
         checkx=t1t2*rrx-rt1x*rt2x+dsqrt(disc)
        endif
        if(abs(checkx).gt.small)iflag=1
        if(iflag.eq.0)then
         ndim=ndim+1
         do j=1,nbase
          x(ndim,j)=xb(j)+dfloat(n1)*a(1,1)+dfloat(n2)*a(1,2)
          y(ndim,j)=yb(j)+dfloat(n1)*a(2,1)+dfloat(n2)*a(2,2)
         enddo
         i1(ndim)=n1
         i2(ndim)=n2
        endif
       enddo
      enddo

      if(nsite.ne.ndim)then
       write(6,*)'wrong site number',nsite,ndim
       stop
      endif

      write(6,*)'Nsite=',ndim
      write(6,*)'Nbase=',nbase
      allocate(ord(ndim*nbase,2))
      write(6,*)
      write(6,*)'CLUSTER'
      write(6,*)
      do i=1,ndim
       do j=1,nbase
        write(19,*)(i-1)*nbase+j,x(i,j),y(i,j)
        ord((i-1)*nbase+j,1)=x(i,j)
        ord((i-1)*nbase+j,2)=y(i,j)
       enddo
      enddo

c  Distances

      ntotalsite=nsite*nbase

      ALLOCATE(dist(ntotalsite,ntotalsite))
      ALLOCATE(boundary(ntotalsite,ntotalsite))

c     write(6,*)
c     write(6,*)'DISTANCES'
c     write(6,*)

      do i=1,ndim
       do ii=1,nbase
        it=(i-1)*nbase+ii
        xi=x(i,ii)
        yi=y(i,ii)
        do j=1,ndim
         do jj=1,nbase
          jt=(j-1)*nbase+jj
          xj=x(j,jj)
          yj=y(j,jj)
          xij=xi-xj
          yij=yi-yj
          d2=xij**2+yij**2
          dist(it,jt)=d2
          iflagibc=0
          do isg=-1,1,2
           xt=xij+isg*(enne(1,1)*a(1,1)+enne(1,2)*a(1,2)) 
           yt=yij+isg*(enne(1,1)*a(2,1)+enne(1,2)*a(2,2)) 
           d2t=xt**2+yt**2
           if(d2t.lt.dist(it,jt)) then
            dist(it,jt)=d2t
            if(ib1.eq.1) iflagibc=1
           endif
           xt=xij+isg*(enne(2,1)*a(1,1)+enne(2,2)*a(1,2)) 
           yt=yij+isg*(enne(2,1)*a(2,1)+enne(2,2)*a(2,2)) 
           d2t=xt**2+yt**2
           if(d2t.lt.dist(it,jt)) then
            dist(it,jt)=d2t
            iflagibc=0
            if(ib2.eq.1) iflagibc=1
           endif
           xt=xij+isg*((enne(1,1)+enne(2,1))*a(1,1)+
     &                 (enne(1,2)+enne(2,2))*a(1,2)) 
           yt=yij+isg*((enne(1,1)+enne(2,1))*a(2,1)+
     &                 (enne(1,2)+enne(2,2))*a(2,2)) 
           d2t=xt**2+yt**2
           if(d2t.lt.dist(it,jt)) then
            dist(it,jt)=d2t
            iflagibc=0
            if(ib1.eq.1.and.ib2.eq.0) iflagibc=1
            if(ib1.eq.0.and.ib2.eq.1) iflagibc=1
           endif
           xt=xij+isg*((enne(1,1)-enne(2,1))*a(1,1)+
     &                 (enne(1,2)-enne(2,2))*a(1,2)) 
           yt=yij+isg*((enne(1,1)-enne(2,1))*a(2,1)+
     &                 (enne(1,2)-enne(2,2))*a(2,2)) 
           d2t=xt**2+yt**2
           if(d2t.lt.dist(it,jt)) then
            dist(it,jt)=d2t
            iflagibc=0
            if(ib1.eq.1.and.ib2.eq.0) iflagibc=1
            if(ib1.eq.0.and.ib2.eq.1) iflagibc=1
           endif
          enddo
          dist(it,jt)=dsqrt(dist(it,jt))
          boundary(it,jt)=1
          if(iflagibc.eq.1) boundary(it,jt)=-1
c         write(6,*) it,jt,dist(it,jt)
         enddo
        enddo
       enddo
      enddo

c  independent distances

      ALLOCATE(dindip(ntotalsite))

      do i=1,nsite*nbase
       dindip(i)=0.d0
      enddo

      dindip(1)=dist(1,1)
      k=1
      do j=2,ndim*nbase
       do l=1,k
        test=dist(1,j)-dindip(l)
        if(abs(test).lt.small) goto 21
        if(test.lt.0.d0) goto 11
       enddo
       l=k+1
11     continue            ! j is at l-th position
       if(l.le.k)then
        do m=k,l,-1
         dindip(m+1)=dindip(m) 
        enddo
       endif
       dindip(l)=dist(1,j)
       k=k+1
21     continue            ! j is already in the list
      enddo

      do i=2,k
       dindip(i-1)=dindip(i)
      enddo
      dindip(k)=0.d0

      nind=k-1               ! number of independent distances

      write(6,*)
      write(6,*) 'TOTAL INDEPENDENT DISTANCES   =',nind
      write(6,*)
      do i=1,nind
       write(6,*) i,dindip(i)
      enddo

c  SYMMETRIES PRESERVED IN THE CLUSTER

c  inverse N matrix
      det=enne(1,1)*enne(2,2)-enne(1,2)*enne(2,1)
      enne1(1,1)=enne(2,2)/det
      enne1(1,2)=-enne(1,2)/det
      enne1(2,1)=-enne(2,1)/det
      enne1(2,2)=enne(1,1)/det

      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      b(1,1)=a(2,2)/det
      b(2,1)=-a(1,2)/det
      b(1,2)=-a(2,1)/det
      b(2,2)=a(1,1)/det

c  Rotations

      ALLOCATE(isymc(nsimmax))

      ALLOCATE(inewb0(nsimmax,nbase))
      ALLOCATE(inewb1(nsimmax,nbase))
      ALLOCATE(inewb2(nsimmax,nbase))

      irotc=0

      do m=1,nrot
       test11=0.d0
       do i=1,2
        do j=1,2
         test11=test11+enne(1,i)*sym(i,j,m)*enne1(j,1)
        enddo
       enddo
       if(itest(test11).eq.0)then
        test12=0.d0
        do i=1,2
         do j=1,2
          test12=test12+enne(1,i)*sym(i,j,m)*enne1(j,2)
         enddo
        enddo
        if(itest(test12).eq.0)then
         test21=0.d0
         do i=1,2
          do j=1,2
           test21=test21+enne(2,i)*sym(i,j,m)*enne1(j,1)
          enddo
         enddo 
         if(itest(test21).eq.0)then
          test22=0.d0
          do i=1,2
           do j=1,2
            test22=test22+enne(2,i)*sym(i,j,m)*enne1(j,2)
           enddo
          enddo
          if(itest(test22).eq.0)then
           if(nbase.ne.1)then              !check symmetries inside the basis 
            cth=dcos(point(m)*2.d0*pi)
            sth=dsin(point(m)*2.d0*pi)
            iflag=0
            do i=1,nbase
             xb1=cth*xb(i)+sth*yb(i) 
             yb1=-sth*xb(i)+cth*yb(i) 
             iflagi=0
             do j=1,nbase
              dx=xb1-xb(j)
              dy=yb1-yb(j)
              testb1=dx*b(1,1)+dy*b(2,1)
              testb2=dx*b(1,2)+dy*b(2,2)
              if(itest(testb1).eq.0.and.itest(testb2).eq.0)then
               iflagi=1
               inewb0(irotc+1,i)=j
               inewb1(irotc+1,i)=nint(testb1)
               inewb2(irotc+1,i)=nint(testb2)
              endif
             enddo
             iflag=iflag+iflagi
            enddo
           else
            iflag=1
            inewb0(irotc+1,1)=1
            inewb1(irotc+1,1)=0
            inewb2(irotc+1,1)=0
           endif
           if(iflag.eq.nbase)then 
            irotc=irotc+1
            isymc(irotc)=m      
c  the symmetries of the infinite lattice preserved in the cluster
           endif
          endif
         endif
        endif
       endif
      enddo

c  Reflections

      irefc=0

      do m=nrot+1,nrot+nref
       test11=0.d0
       do i=1,2
        do j=1,2
         test11=test11+enne(1,i)*sym(i,j,m)*enne1(j,1)
        enddo
       enddo
       if(itest(test11).eq.0)then
        test12=0.d0
        do i=1,2
         do j=1,2
          test12=test12+enne(1,i)*sym(i,j,m)*enne1(j,2)
         enddo
        enddo
        if(itest(test12).eq.0)then
         test21=0.d0
         do i=1,2
          do j=1,2
           test21=test21+enne(2,i)*sym(i,j,m)*enne1(j,1)
          enddo
         enddo
         if(itest(test21).eq.0)then
          test22=0.d0
          do i=1,2
           do j=1,2
            test22=test22+enne(2,i)*sym(i,j,m)*enne1(j,2)
           enddo
          enddo
          if(itest(test22).eq.0)then
           if(nbase.ne.1)then              !check symmetries inside the basis
            cth=dcos(point(m)*2.d0*pi)
            sth=dsin(point(m)*2.d0*pi)
            iflag=0
            do i=1,nbase
             xb1=cth*xb(i)+sth*yb(i)
             yb1=sth*xb(i)-cth*yb(i)
             iflagi=0
             do j=1,nbase
              dx=xb1-xb(j)
              dy=yb1-yb(j)
              testb1=dx*b(1,1)+dy*b(2,1)
              testb2=dx*b(1,2)+dy*b(2,2)
              if(itest(testb1).eq.0.and.itest(testb2).eq.0)then
               iflagi=1
               inewb0(irotc+irefc+1,i)=j
               inewb1(irotc+irefc+1,i)=nint(testb1)
               inewb2(irotc+irefc+1,i)=nint(testb2)
              endif
             enddo
             iflag=iflag+iflagi
            enddo
           else
            iflag=1
            inewb0(irotc+irefc+1,1)=1
            inewb1(irotc+irefc+1,1)=0
            inewb2(irotc+irefc+1,1)=0
           endif
           if(iflag.eq.nbase)then
            irefc=irefc+1
            isymc(irotc+irefc)=m    
c  the symmetries of the infinite lattice preserved in the cluster
           endif
          endif
         endif
        endif
       endif
      enddo

      write(6,*)
      write(6,*)'CLUSTER ROTATIONS   =',irotc
      write(6,*)'CLUSTER REFLECTIONS =',irefc
      write(6,*)

c  Allowed momenta (in units of 2*pi)

      nmom=0
      n1max=enne(1,1)+abs(enne(2,1))
      n2max=enne(2,2)+abs(enne(1,2))

      small=1.d-9
      uno=1.d0-small
      do n1=-n1max,n1max
       do n2=-n2max,n2max
        eig1=enne1(1,1)*n1+enne1(1,2)*n2
        eig2=enne1(2,1)*n1+enne1(2,2)*n2
        if(eig1.ge.0.d0.and.eig1.lt.uno.and.
     >     eig2.ge.0.d0.and.eig2.lt.uno)then
         qx=eig1*b(1,1)+eig2*b(1,2)
         qy=eig1*b(2,1)+eig2*b(2,2)
         nmom=nmom+1
         write(27,*) qx,qy
        endif
       enddo
      enddo

      write(6,*)'Independent momenta=',nmom

      if(nmom.ne.ndim)then
       write(6,*)'problems in determining the momenta'
       stop
      endif

       
      do i=1,ndim
       write(28,*) x(i,1),y(i,1)
      enddo
      
ccccccccccccccccccccccccc  END CLUSTER PART ccccccccccccccccccccccccc

ccccccccccccccccccccccccc  TABLE OF NEIGHBORS ccccccccccccccccccccccccc

      nindm=1
      ALLOCATE(imulti(nindm))

      do k=1,nindm
       d=dindip(k)
       imulti(k)=0
       do j=1,ndim*nbase
        dij=dist(j,1)
        if(abs(dij-d).lt.small) imulti(k)=imulti(k)+1
       enddo 
      enddo

      maxmulti=0
      do k=1,nindm
       maxmulti=max(maxmulti,imulti(k))
      enddo

      ALLOCATE(ivic(ntotalsite,maxmulti,nindm))

      do k=1,nindm
       do i=1,ndim*nbase
        do j=1,maxmulti
         ivic(i,j,k)=0
        enddo
       enddo
      enddo

      do k=1,nindm
       d=dindip(k)
       do i=1,ndim*nbase
        ii=0
        do j=1,ndim*nbase
         dij=dist(i,j)
         if(abs(dij-d).lt.small) then
          ii=ii+1
          ivic(i,ii,k)=j
          if(boundary(i,j).lt.0) ivic(i,ii,k)=-ivic(i,ii,k)
         endif
        enddo
       enddo
      enddo

c  In this table the order of neighbors for each site is given by the 
c  lexicographic order of the sites
c     write(6,*) 
c     write(6,*) 'Table of neighbors'
c     do k=1,nindm
c      write(6,*) 'At distance k=',k
c      do i=1,ndim*nbase
c       write(6,334) i,(ivic(i,j,k),j=1,imulti(k))
c      enddo
c     enddo

      ALLOCATE(isymt(ntotalsite,ndim))  ! table of translational symmetries

      do i=1,ndim            ! sites
       it=0
       do k=1,ndim           ! translations
        ik1=i1(i)+i1(k)
        ik2=i2(i)+i2(k)
        it=it+1
        iflag=0
        do l=1,ndim
         l1=ik1-i1(l)
         l2=ik2-i2(l)
         test1=enne1(1,1)*l1+enne1(2,1)*l2
         test2=enne1(1,2)*l1+enne1(2,2)*l2
         if(itest(test1).eq.0.and.itest(test2).eq.0) then
          iflag=iflag+1
          ik=l
         endif
        enddo
        if(iflag.ne.1) then
         write(6,*)'wrong reduction'
         stop
        endif
        do ii=1,nbase
         j=(i-1)*nbase+ii
         jj=(ik-1)*nbase+ii
         isymt(j,it)=jj
        enddo
       enddo
      enddo

      ALLOCATE(ivict(ntotalsite,maxmulti,nindm))
      ALLOCATE(imark(ndim))

      do k=1,nindm
       do i=1,ndim*nbase
        do j=1,maxmulti
         ivict(i,j,k)=0
        enddo
       enddo
      enddo

      do i=1,ndim
       imark(i)=0
      enddo

      IF(nbase.ne.1) THEN ! we do not use inversion symmetry
       do k=1,nindm
        do i=1,nbase
         do j=1,imulti(k) 
          ii=abs(ivic(i,j,k))
          do it=1,nsite   ! traslations
           m=isymt(i,it)
           l=isymt(ii,it)
           ivict(m,j,k)=l
           if(boundary(m,l).lt.0) ivict(m,j,k)=-ivict(m,j,k)
          enddo
         enddo
        enddo
       enddo
      ELSE                ! we order neighbors using the inversion symmetry 
       do k=1,nindm
        do ii=1,nsite
         imark(ii)=0
        enddo
        jshift=imulti(k)/2
        jj=0
        do j=1,imulti(k) 
         i=abs(ivic(1,j,k))
         if(imark(i).eq.0) then
          imark(i)=1  ! mark the site...
          iflag=0
          ik1=i1(i)-i1(1)
          ik2=i2(i)-i2(1)
          ikm1=sym(1,1,2)*ik1+sym(2,1,2)*ik2+i1(1)
          ikm2=sym(1,2,2)*ik1+sym(2,2,2)*ik2+i2(1)
          do l=1,ndim              !reduction
           l1=ikm1-i1(l)
           l2=ikm2-i2(l)
           test1=enne1(1,1)*l1+enne1(2,1)*l2
           test2=enne1(1,2)*l1+enne1(2,2)*l2
           if(itest(test1).eq.0.and.itest(test2).eq.0) then
            iflag=iflag+1
            ikm=l         ! inverted site
           endif
          enddo
          if(iflag.ne.1) then
           write(6,*) 'wrong reduction'
           stop
          endif
          if(i.ne.ikm) then ! there is an inverted site
           jj=jj+1
           iflagref=0
           imark(ikm)=1  ! ...and its inverted
          elseif(i.eq.ikm) then ! the inverted site is the site itself
           jj=jj+1
           iflagref=1
          endif
          do it=1,nsite   ! traslations
           m=isymt(1,it)
           l=isymt(i,it)
           if(iflagref.eq.0) then
            ivict(m,jj,k)=l
            if(boundary(m,l).lt.0) ivict(m,jj,k)=-ivict(m,jj,k)
            l=isymt(ikm,it)
            ivict(m,jj+jshift,k)=l
            if(boundary(m,l).lt.0) 
     >      ivict(m,jj+jshift,k)=-ivict(m,jj+jshift,k)
           elseif(iflagref.eq.1) then
            ivict(m,jj,k)=l
            if(boundary(m,l).lt.0)
     >      ivict(m,jj,k)=-ivict(m,jj,k)
           endif
          enddo
         endif
        enddo
       enddo
      ENDIF

c  In this table the order of neighbors for each site is given by the 
c  order of the first site
      write(6,*) 
      write(6,*) 'Ordered Table of neighbors'
      do k=1,nindm
       write(6,*) 'At distance k=',k
       do i=1,ndim*nbase
        write(6,334) i,(ivict(i,j,k),j=1,imulti(k))
       enddo
      enddo

334   format(i4,4x,100(i5,2x))

ccccccccccccccccccccccccc  WRITE INFORMATION ccccccccccccccccccccccccc

      rewind(9)
      write(9) ndim,nbase,nindm,maxmulti
      write(9) ivict,ord

ccccccccccccccccccccccccc  END WRITE INFORMATION ccccccccccccccccccccccccc

      stop
      end

      integer*4 function itest(test)
      implicit none
      integer*4 ii
      real*8 test,check,small
      parameter(small=1.d-10)
      
      ii=nint(test)
      check=test-dfloat(ii)
      if(abs(check).lt.small)then
       itest=0      
      else
       itest=1
      endif

      return
      end
