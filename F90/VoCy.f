c
c   	aps 2020 11 13
c	double precision by compiler option
c
	program TestInduction
c
c--------------- input part from here
c
	use mem
c
        real, allocatable :: vzK(:),vrK(:),vzBo(:)
	real, allocatable :: vrBo(:),vzBr(:),vrBr(:)
c
	character*10 s(3),inpstr
c
	integer status
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
	read(5,*) inpstr,restart
c
c       read number of singularities N
c
	read(5,*) inpstr,NP
	write(*,*)inpstr,NP
c
c       set up case - propeller ct > 0
c
	read(5,*) inpstr,ct
	write(*,*)inpstr,ct
c	
c       computational model 1: Bontempo, 2: van Kuik, 3: Branlard
c
	s(1) = 'Bontempo '
	s(2) = 'van Kuik '
	s(3) = 'Brandlard'
c
	read(5,*) inpstr,zsivc
	write(*,*)inpstr,zsivc
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
	write(*,*)
	write(*,*)'From axial momentum theory:'
	write(*,*)'==========================='
	write(*,'(a20,3f12.4)')'cT vzdisk cp ',ct,1.+a,cpm
	write(*,*)
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
	zanf = 0.
	zend = zsivc+1.
c
	rcyl = 0.5
c
	nn  = 10000
	dz = (zend-zanf)/float(nn)
c  
	write(*,*)'Cylinder at z r = ',zsivc,rcyl
	write(*,*)'strength',gasivc
	write(101,*)'Cylinder at z r = ',zsivc,rcyl
c
	write(101,100)'z','vzK','vrK','vzBo','vrBo','vzBr','vrBr',
     +                'psicyl'
c

	do i =0,nn
	   z = zanf+i*dz
c	   write(*,*)'z ',z
	   write(101,101)z,vzsivcK (z,rcyl),vrsivcK (z,rcyl),
     +                     vzsivcBo(z,rcyl),vrsivcBo(z,rcyl),
     +                     vzsivcBr(z,rcyl),vrsivcBr(z,rcyl),
     +                     psiBrCyl(z,rcyl)
	end do
c
crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
c
	write(*,*)
        write(*,*)'VORTEX RINGS'
	write(*,*)'velocities at slipstream position'
        write(100,*)'VORTEX RINGS'
	write(100,100)'z', 'r','strength','vzK ','vrK ', 
     +                                    'vzBo','vrBo',
     +                                    'vzBr','vrBr'
c
	do i=1,np
	   vzK(i)  = 0.
	   vrK(i)  = 0.
	   vzBo(i) = 0.
	   vrBo(i) = 0.
	   vzBr(i) = 0.
	   vrBr(i) = 0.
	end do
c	 
	do n = 1,np
c
c          induction from all m's
c          1. Branlard
c
	   do m = 1,np
	     if (m.eq.n)then
 	       vzBr(n) = vzBr(n) + vzself(n)
               vrBr(n) = vrBr(n) + vrself(n)
             else
	        vzBr(n) = vzBr(n) + vznmBr(n,m)
	        vrBr(n) = vrBr(n) + vrnmBr(n,m)
             end if
	   end do
c
c          induction from all m's
c          2. von Kuik
c
	   do m = 1,np
	     if (m.eq.n)then
 	       vzK(n) = vzK(n) + vzself(n)
               vrK(n) = vrK(n) + vrself(n)
             else
	        vzK(n) = vzK(n) + vznmK(n,m)
	        vrK(n) = vrK(n) + vrnmK(n,m)
             end if
	   end do
c
c          induction from all m's
c          3. Bontempo
c
	   do m = 1,np
	     if (m.eq.n)then
 	       vzBo(n) = vzBo(n) + vzself(n)
               vrBo(n) = vrBo(n) + vrself(n)
             else
	        vzBo(n) = vzBo(n) + vznmBo(n,m)
	        vrBo(n) = vrBo(n) + vrnmBo(n,m)
             end if
	   end do
c
	  write(100,101)zsl(n),rsl(n),ga(n),vzK(n),vrK(n),
     +                  vzBo(n),vrBo(n),vzBr(n),vrBr(n)
	end do
c
	write(*,*)'End'
c
c---------------------------------------------------------------------
c
100	format(10a18)
101	format(10e18.6)
102	format(a15,2f12.4)
c
c-----------------------------------------------------------------------
	end
c-----------------------------------------------------------------------
