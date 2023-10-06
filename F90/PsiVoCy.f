c
c   	aps 2021 03 26
c	double precision by compiler option
c
	program TestPsiSIVC
c
c--------------- input part from here
c
	use mem
c
        real, allocatable :: vzK(:),vrK(:),vzBo(:)
	real, allocatable :: vrBo(:),vzBr(:),vrBr(:)
c
	real psiout
	real nslope
c
	integer maxiter, niter, model,status
c
	character*10 s(3),inpstr
c
	logical update
c
c       open files
c
	OPEN(UNIT= 5,FORM='FORMATTED',STATUS='UNKNOWN',FILE='inpa.dat')
c
c       some numbers
c	
	pih = 2.*atan(1.)
	pi  = 4.*atan(1.)
c
c       read configutration inpa.dat
c
c       read number of singularities N
c
	read(5,*) inpstr,NP
c	write(*,*)inpstr,NP
c
c       set up case - propeller ct > 0
c
	read(5,*) inpstr,ct
c	write(*,*)inpstr,ct
c	
c       computational model 1: Bontempo, 2: van Kuik, 3: Branlard
c
	s(1) = 'Bontempo '
	s(2) = 'van Kuik '
	s(3) = 'Brandlard'
c
	read(5,*) inpstr,zsivc
c	write(*,*)inpstr,zsivc
c
c       some numbers concerning case
c
	vzinf = sqrt(1. + ct)
	vz0   = 0.5*(1. + vzinf)
	vzsinf= vz0
	rinf  = sqrt(vz0/vzinf)
c
	a     = 0.5*(sqrt(1.+ct)-1.)
	cpm   = 4.*a*(1+a)**2
c
	rsivc = rinf
c
c      allocate memory
c
	allocate(rsl  (0:np),stat=status)
	allocate(zsl  (0:np),stat=status)
	allocate(ga   (0:np),stat=status)
	allocate(rsln (0:np),stat=status)
	allocate(zsln (0:np),stat=status)
	allocate(vz   (0:np),stat=status)
	allocate(vr   (0:np),stat=status)
c
	allocate(vrK  (0:np),stat=status)
	allocate(vrBo (0:np),stat=status)
	allocate(vrBr (0:np),stat=status)
	allocate(vzK  (0:np),stat=status)
	allocate(vzBo (0:np),stat=status)
	allocate(vzBr (0:np),stat=status)
c
c	write(*,*)
c	write(*,*)'From axial momentum theory:'
c	write(*,*)'==========================='
c	write(*,'(a20,3f12.4)')'cT vzdisk cp ',ct,1.+a,cpm
c	write(*,*)
c
c       ini wake shape
c
	rsl(0) = 1.
	zsl(0) = 0.
c
	do i =0,np
	 rsln(i) = 0.
	 zsln(i) = 0.
	end do
c
c       ini stength 
c
	gasivc = 1.-sqrt(1.+ct)
c
	do i=1,np
	   ga(i) = gasivc
	end do
c
c	define and stretch zsl
	
	dphi  = pih/real(np)
c
	do i=1,np
	   phi = i*dphi
	   zsl(i) = zsivc*(1.-cos(phi))
c	   write(*,*)cos(phi),zsl(i)
	end do

c
	if (nwt.eq.0)then
	   do i=1,np
	      rsl(i) = 1.
	    end do
        else
	   do i=1,np
	      rsl(i) = rinf+(1.-rinf)*exp(-al*zsl(i))
	   end do
	endif
c
	do i=1,np
	   vz(i) = 0.
	   vr(i) = 0.
	end do
c
c-------------------------------- 
c
	zanf =  zsivc-5.
	zend =  zsivc+5.
	nz   = 100
c
	ranf = 0.01
	rend = 1.
	nr   =20
c
	dz = (zend-zanf)/float(nz)
	dr = (rend-ranf)/float(nr)
c	write(*,*)'zanf zend dz',zanf, zend, dz
c	write(*,*)'ranf rend dr',ranf, rend, dr
c	write(*,*)
c  
c	write(*,*)'Cylinder at z r = ',zsivc,rcyl
c	write(*,*)'strength',gasivc
c
c	write(*,100)'z','r','psicyl'
c
	do i =0,nz
	   z = zanf+i*dz
	   do j = 0, nr
	      r = ranf+j*dr
	         zz = z
                 rr = r
   	         write(*,101)z,r,psiBrCyl(zz,rr)
          end do
	  write(*,*)
	end do
c
c---------------------------------------------------------------------
c
100	format(3a18)
101	format(3e18.6)
c
c-----------------------------------------------------------------------
	end
c-----------------------------------------------------------------------
