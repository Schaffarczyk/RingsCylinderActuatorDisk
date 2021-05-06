c
c
c----------------------------------------------------
	subroutine indvel
c----------------------------------------------------
c
	use mem
c
	do i=1,np
	      vz(i) = 0.
	      vr(i) = 0.
	   end do
c
c          main loop over all vortex strips for INDUCED velocities
c
	   do n = 1,np
c
	      do m = 1,np
	         if (m.eq.n)then
	            if (selfind)then
c
c                      self induction Bontempo Eq (6) and (7)
c
                       vz(n)   = vz(n) + vzself(n)
                       vr(n)   = vr(n) + vrself(n)
	             end if
                 else  
c
c                induction from all (other) m's
c
             	    select case(model)
             	    case(1)
	               vz(n) = vz(n) + vznmBo(n,m)
	               vr(n) = vr(n) + vrnmBo(n,m)
	            case(2)
	               vz(n) = vz(n) + vznmK(n,m)
	               vr(n) = vr(n) + vrnmK(n,m)
	            case(3)
	               vz(n) = vz(n) + vznmBr(n,m)
	               vr(n) = vr(n) + vrnmBr(n,m)
	            case Default
	            end select 	
                 end if
	      end do
c	      end m loop
c
c             induction from cylinder to sniplet n
c
	      z = 0.5*(zsl(n)+zsl(n-1))
	      r = 0.5*(rsl(n)+rsl(n-1))  
c       
              select case(model)
              case(1)
	         vz(n) = vz(n) + vzsivcBo(z,r)
	         vr(n) = vr(n) + vrsivcBo(z,r)
	      case(2)
	         vz(n) = vz(n) + vzsivcK(z,r) 
	         vr(n) = vr(n) + vrsivcK(z,r) 
	      case(3)
	         vz(n) = vz(n) + vzsivcBr(z,r) 
	         vr(n) = vr(n) + vrsivcBr(z,r) 
	      case Default
	      end select      
c
	end do
c       end n loop   
c	write(*,*)'end n loop'
c
	return
	end
c
c----------------------------------------------------------------------
	subroutine output
c----------------------------------------------------------------------
c
	use mem
c
	if(restart)then
	   OPEN(UNIT= 3,FORM='FORMATTED',STATUS='REPLACE',FILE='sls.DAT')
	else
	   OPEN(UNIT= 3,FORM='FORMATTED',STATUS='unknown',FILE='sls.DAT')
	end if
c
        write(3,99)'i','z','r','vz','vr','ga','psi-Ri','psi-Cyl',
     +             'vz-s','vr-s'
c
	do i=1,np
	   zpsi   = 0.5*(zsl(i)+zsl(i-1))
	   rpsi   = 0.5*(rsl(i)+rsl(i-1))
	   write(3,100)i,zsl(i),rsl(i),1.+vz(i),vr(i),ga(i),
     +     psiBr(zpsi,rpsi),psiBrCyl(zpsi,rpsi),vzself(i),vrself(i)
	end do
c
	write(*,*)
c
c**********************************************************************
c
c       calculate cP from Bontempo's Eq. (19) and print out vz at disk
c
	cpmom  = cpm
	cpcode = cpp(cp)
c
	write(*,*)'--------------  Finally  -------------------------'
		dev = (cpcode-cpmom)/cpmom
	write(*,105)'cpmom     cp      dev = ',cpmom,cpcode,dev
	write(*,*)'-----------------------------------------------'
c	                  
c**********************************************************************
c
c       print out vz at z = 0 (disk)
c
	write(1,*)
	write(1,*)
	nd = 100
	drd = 1./float(nd)
	do i=0,nd-1
	   rd = 0.5*drd+drd*real(i)
	   vzout = 1.+vzd(rd)
	   vzoutc= vzdcyl(rd)
	   write(1,102)rd,vzout,vzoutc
	end do
c	                  
c**********************************************************************
c
c       print out vr at z = 0 (disk)
c
	write(1,*)
	write(1,*)
	nd = 1000
	drd = 1./float(nd)
	do i=0,nd-1
	   rd = 0.5*drd+drd*real(i)
	   vrout = vrd(rd)
	   vroutc= vrdcyl(rd)
	   write(11,102)rd,vrout,vroutc
	end do
c
c**********************************************************************
c
c       print out stream function at z = 0 (disk)
c
	nd = 100
	drd = 1./float(nd)
	do i=0,nd-1
	   rd = 0.5*drd+drd*real(i)
	   write(2,102)rd,psiBr(0.,rd),psiBrCyl(0.,rd)
	end do
c
c**********************************************************************
c
c       print out psi in defined range 
c
	nz = 100
	nr = 100
c
	za = 0.
	ze = 5.
	ra = 1.e-3
	re = rsivc
c
	dr = (re-ra)/float(nr)
        dz = (ze-za)/float(nz)
	do i=0, nr
           rp = ra + i*dr
	   do j = 0, nz
             zp = za+j*dz
	     psipr = psiBr(zp,rp)+psiBrCyl(zp,rp)
	     write(12,107)rp,zp,psipr 
	   end do
           write(12,*)
	end do
c
c**********************************************************************
c
c       print out vz at r = 0 (axis)
c
	write(4,*)
	write(4,*)
	na = 500
	dz = zsivc/float(na)
	r  = rsivc
	do i=0,na
	   za = dz*real(i)
	   write(4,102)za,vza(za),vzsivcBr(za,1.e-6),vzacylan(za)
	end do
c
99	format(10a12) 
100	format(i12,9f12.6)
102	format(4f18.8)
105	format(a25,2f12.6,e10.2)
107	format(3f12.6)
c
	return
	end
c
c----------------------------------------------------------------
	subroutine initialize
c----------------------------------------------------------------------
c
	use mem
c
	integer values(8), status
c
	real stime,etime
	real d(4)
	real al
c
	character*14 s(3),inpstr,is(4),slt(0:3)
	character*8 date
	character*10 time
	character*5 zone
c
	OPEN(UNIT= 5,FORM='FORMATTED',STATUS='UNKNOWN',FILE='inpa.dat')

c       some numbers
c	
	pih = 2.*atan(1.)
	pi  = 4.*atan(1.)
c
c       some strings
c
	s(1) = 'Bontempo '
	s(2) = 'van Kuik '
	s(3) = 'Brandlard'
c
	is(1) = 'Bontempo'
	is(2) = 'vK:norma'
	is(3) = 'vK: psi '
c
	slt(0) = ' r = 1    '
	slt(1) = ' exp simp '
	slt(2) = ' exp adv  '
	slt(3) = ' Eq. (D.9)'
c
	call date_and_time(date,time,zone,values)
c
	write(*,*)
	write(*,*)'================================================ '
	write(*,'(a6,a4,a1,a2,a1,a2,a6,a4)')
     +  'date ',date(1:4),' ',date(5:6),' ',date(7:8),' time ',time
        Write(*,*)'  Vortex Code for constantly loaded Actuator Disk'
	write(*,*)'  INITIALIZE                                     '
	write(*,*)'================================================ '
c
c       read configuration file inpa.dat
c       --------------------------------
c
c       restart from previosu case ?
c
	read(5,*) inpstr,restart
	write(*,*)inpstr,restart
c
c       read number of singularities NP
c
	read(5,*) inpstr,NP
	write(*,700)inpstr,NP
700	format(a10,i6)
c
c       set up case - cT thrust coefficient, propeller ct > 0
c
	read(5,*) inpstr,ct
	write(*,701)inpstr,ct
701	format(a10,f12.6)
c
	if(ct.gt.0.)then
	   write(*,*)'PROPELLER case'
	else
	   write(*,*)'TURBINE case'
	end if	
c
c	cylinder
c
	read(5,*) inpstr,zsivc
	write(*,701)inpstr,zsivc
c
c       include Bontempo's self induction ?
c
c       logical selfind
c
	read(5,*) inpstr,selfind
	write(*,*)inpstr,selfind
c
c       slip stream iteration parameters
c
c       logical update
c
	read(5,*) inpstr,update
c
c       maximum number of iterations
c
	read(5,*) inpstr,maxiter
c
c       some under-relaxation for sls 
c
	read(5,*) inpstr,under
c
c       epsilon to finish iteration
c
	read(5,*) inpstr,epsiter
c
c       computational model for induced velocities
c	1: Bontempo, 2: van Kuik, 3: Branlard
c
	read(5,*) inpstr,model
	write(*,*)inpstr,model, s(model)
704	format(a10,f12.6,a10)
c
c       ini wake shape by an assumed exp 
c
c       nwt 0 -> rsl = 1 
c	nwt 1 -> rsl = exp
c
	read(5,*) inpstr,nwt
c
c       slipstream saturation length (1 to 3)
c
	read(5,*) inpstr,al
	write(*,701)inpstr,al
c
c       model for sls iteration see line 240 ff
c       ---------------------------------------
c	1: Bontempo, 
c       2: van Kuik: normal to stream lines 
c       3: van Kuik Psi
c
	read(5,*) inpstr,modelit
	write(*,*)
c
c*****************************************************************
c       some numbers concerning case from 1D momentum theory
c	
c	total unpertubed velocity for all case = 1
c
c	total velocity far downstream
	vzinf = sqrt(1. + ct)
c	velocity at disk
	vz0   = 0.5*(1. + vzinf)
c	radial extension of slipstream far downstream
	rinf  = sqrt(vz0/vzinf)
c	axial induction at disk
	a     = 0.5*(sqrt(1.+ct)-1.)
	cpm   = 4.*a*(1+a)**2
c
c	stream function far downstream at slip stream
c       i) inflow inlcuded
c
	psiwake = 0.25*(1.+sqrt(1.+ct))
c
c       ii) only induced velocity
c
c	psiwake = ct/(4.*sqrt(1.+ct))
c
	rsivc = rinf
c
c       open more files
c
	OPEN(UNIT= 1,FORM='FORMATTED',STATUS='UNKNOWN',FILE='vzd.DAT')
	OPEN(UNIT= 2,FORM='FORMATTED',STATUS='UNKNOWN',FILE='psid.DAT')
	OPEN(UNIT= 3,FORM='FORMATTED',STATUS='UNKNOWN',FILE='sls.DAT')
	OPEN(UNIT= 4,FORM='FORMATTED',STATUS='UNKNOWN',FILE='vza.DAT')

	OPEN(UNIT=10,FORM='FORMATTED',STATUS='UNKNOWN',FILE='conv.DAT')
	OPEN(UNIT=11,FORM='FORMATTED',STATUS='UNKNOWN',FILE='vrd.DAT')
	OPEN(UNIT=12,FORM='FORMATTED',STATUS='UNKNOWN',FILE='streaml.DAT')
	OPEN(UNIT=13,FORM='FORMATTED',STATUS='UNKNOWN',FILE='waked.DAT')
	OPEN(UNIT=14,FORM='FORMATTED',STATUS='UNKNOWN',FILE='psidev.DAT')
c
c       allocate memory
c
	allocate(rsl  (0:np),stat=status)
	allocate(zsl  (0:np),stat=status)
	allocate(ga   (0:np),stat=status)
	allocate(rsln (0:np),stat=status)
	allocate(zsln (0:np),stat=status)
	allocate(vz   (0:np),stat=status)
	allocate(vr   (0:np),stat=status)
c
c       ini wake shape
c
	rsl(0) = 1.
	zsl(0) = 0.
c
c       initialize vortex strength of cylinder
c
        gasivc = 1.-sqrt(1.+ct)
c
	write(*,108)'gasivc=  ',gasivc
	write(*,108)'psiwake= ',psiwake
	write(*,*)
c
c       give all rings the same strength as the cylinder
c
	if(restart)then
	   read(3,*)
	   do l=1,np
	      read(3,*)id,d1,rsl(l),da,d3,ga(l),(d(ll),ll=1,4)
           end do
	   write(*,*)'done read sls.DAT for restart'
	   close(3)
	else   
	   do i=1,np
	      ga(i) = gasivc
	   end do
	end if
c
c	define  zsl distribution
c	
	dphi  = pih/real(np)
c
c       simple cosine stretch
c
	do i=1,np
	   phi = i*dphi
	   zsl(i) = zsivc*(1.-cos(phi))
c	   write(*,*)cos(phi),zsl(i)
	end do
c
c       initialize slipstream shape
c
c       nwt = 0 -> rsl = 1 (cylinder)
c      	nwr = 1 -> simple exponential blending
c      	nwr = 2 -> more general exponential blending
c      	nwr = 3 -> according to vK Eq. (D.9) 
c
           select case(nwt)
           case(1)
	      do i=1,np
	         rsl(i) = rinf+(1.-rinf)*exp(-al*zsl(i))
	      end do
	   case(2)
c
c          from gnuplot fit of vK data cT = 1 and cT = -8./9. only
c
	      if (cT.gt.0.)then
                 exsl = 0.757037
                 a2   = 1.7743
	      else
	         exsl = 0.707145
                 a2   = 0.734338
              endif
c
	      do i=1,np
	         rsl(i) = rinf+(1.-rinf)*exp(-a2*zsl(i)**exsl)
	      end do
c
	   case(3)
c
c          accoring vz_ind(r=0) accoriding van Kuik's Eq. (D.9)
c
	   do i=1,np
	      u0    = 1.-0.5*gasivc
	      r0 = i*(rsivc-1.)/np + 1.
	      uwake = 1.-0.5*gasivc*(1.+zsl(i)/(sqrt(r0*r0+zsl(i)**2)))
	      rsl(i)= sqrt(u0/uwake)
	   end do
c
	      case Default
	         do i=1,np
	            rsl(i) = 1.
	         end do	
              end select 
c
c
c          initialize new slip stream coordiantes from iteration
c
	do i =0,np
	   rsln(i) = 0.
	   zsln(i) = 0.
	end do
c
c*****************************************************************
c       output some values from inpa.dat
c
	write(*,'(a20,4f12.4)')'a vz0 vzinf rinf',a, vz0, vzinf, rinf
	write(*,*)
	write(*,'(a22,2f12.6)')'Momentum Theory: cP= ',cpm
	write(*,*)
c
	write(*,  *)'update  ',update
	write(*,702)'maxiter ',maxiter
	write(*,703)'under   ',under
	write(*,703)'epsiter ',epsiter
	write(*,706)'sl model',nwt,slt(nwt)
	write(*,705)'itersc  ',modelit,is(modelit)
	write(*,*)
702	format(a10,i6)
703	format(a10,e12.4)
705     format(a10,i10,x,a10)
706     format(a10,i6,a10)
c
	close(5)
c
108	format(a20,f12.8)
c
	return
	end
c
c**********************************************************************
c
c       calculate cP from Bontempo Eq. (19) 
c
	real function cpp(cp)
c
	use mem
c
	integer nd
c
	cp = 0.
	a  = 0.
	nd = 2000
	drd = 1./float(nd)
c
c       Simpson's rule
c
c       Factor 2: Bontempo Eq (18)
c
	do i  = 2,nd-2,2
	   r  = drd*real(i)
	   f1 = 2.*vzd(r)
	   f2 = 2.*vzd(r+drd)
           a  = a + 2.*r*f1 + 4.*(r+drd)*f2
        end do
c   
	f1   = 2.*vzd(1. )
	fdrd = 2.*vzd(drd)
	a = (a + f1 + 4.*drd*fdrd)*drd/3.

	cp = ct*(1. + a)
c
	return
	end
c
c--------------------------------------------------------------
c
	real function vzsivcBo(z,r)
c
c       induction of v_z at P = z from cylinder at zsivc
c	Bontempo Eq (1)
c
	USE MEM
c
	real k,ks,n,delta,eps, rs,zs, f1, f2
	integer isl
c
	eps = 1.e-7
c
	vzsivcBo = 0.
c
	zs = (z-zsivc)/rsivc
	rs = r/rsivc
c
	ks = 4.*rs/(zs*zs+(1.+rs)*(1.+rs))
	k  = sqrt(ks)
	n  = 4.*rs/((1.+rs)**2)
c
	delta = pi
	if(r.gt.rsivc)delta = 0.
c
	if(rs.gt.eps) then
	   f1       = zs/(sqrt(zs**2+(1.+rs)**2))
	   f2       = (rs-1.)/(rs+1.)
	   f2       = ellk(ks)-f2*ellp(n,ks)
	   vzsivcBo = -(delta+f1*f2)*gasivc/(2.*pi)
	else
           f1 = 0.25 + ellk(k)*zs/(2.*pi*sqrt(4.+zs*zs))
	   vzsivcBo = -gasivc*f1
        endif
c
c	write(*,*) 'r vzsivc',r,rsivc,vzsivcBo
c
	return
	end
c
c----------------------------------------------------------------------------------
c
	real function vzsivcK(z,r)
c
c	van Kuik Eqs (D.6 and 7)
c
	use mem
c
	real k, n, ns,ks, theta, delta,eps
c
	dz = zsivc-z
	dr = rsivc-r
c
	eps = 1.e-8
c
	if(r.gt.rsivc)         delta=0.
	if(r.lt.rsivc)         delta=2.*pi
	if(abs(r-rsivc).lt.eps)delta=pi
c
        A = dz**2 + r**2 + rsivc**2
	B = -2.*r*rsivc
	rho1 = sqrt(dz*dz+dr*dr)
	rho2 = sqrt(dz**2+(r+rsivc)**2)	
	k    = sqrt(1.-(rho1/rho2)**2) 
	ks = k*k
	n    = 2.*(sqrt(r*rsivc))/(r+rsivc)
	ns = n*n
c
	f2 = (r-rsivc)/(r+rsivc)
c	
c       2021 04 121:  n -> n**2, k -> k**2
c
	f2 = ellk(ks)-f2*ellp(ns,ks)
	f2 = 2.*f2*(z-zsivc)/rho2
	theta = delta + f2  
c
	vzsivcK = -gasivc*theta/(4.*pi)
c
	return
	end
c
c----------------------------------------------------------------------------------
c
	real function vzsivcBr(z,r)
c
c	Brandlard Eq. (36.70)
c
	use mem
c
	real k,ks,n,f1,f2,f3,zc,rc,eps
c
c       important !!
c
	eps = 1.e-12
c
	rc = r
c
	zc = z - zsivc
	f3 = rsivc-rc
c
	ks = 4.*rc*rsivc/((rsivc+rc)**2+zc**2)
	k  = sqrt(ks)
	n  = 4.*rc*rsivc/(rsivc+r)**2
c
	f1 = ellk(ks) + ((rsivc-rc)/(rsivc+rc))*ellp(n,ks)
	f2 = zc*k/(2.*pi*sqrt(rc*rsivc))
c
c       fixed 2021 04 12
c
	if (abs(f3).lt.eps) then
	   f3 = 0.5
	else
           f3 = (f3 + abs(f3))/(2.*abs(f3))
        endif
c
	vzsivcBr = -0.5*gasivc*(f3+f2*f1)
c
	write(222,'(7f12.6)')z,zc,r,f1,f2,f3,vzsivcBr
c	vzsivcBr = 0.
c
	return
	end
c
c----------------------------------------------------------------------------------
c
	real function vrsivcBr(z,r)
c
	use mem
c
	real zc,rc,k,ks,f1,eps
c
	eps = 1.e-12
	rc  = r
	zc  = z-zsivc
c	
	if (rc.gt.eps) then
c
c  	   Brandlard Eq. (36.69)
c
	   ks  = 4.*rc*rsivc/((r+rsivc)**2+zc**2)
           k   = sqrt(ks)  
c
	   f1 = (2.-ks)*ellk(ks)-2.*elle(ks)
	   f1 = f1/k                  
c               
	   vrsivcBr = -gasivc*f1*sqrt(rsivc/rc)/(2.*pi)
c
	else
	   vrsivcBr = 0.
	end if
c
	return
	end
c----------------------------------------------------------------------------------
c
	real function vrsivcK(z,r)
c
c	van Kuik Eq (D.8)
c
	use mem
c
	real ks,k,n
	REAL A, B, rho1, rho2, f1
c
	dz = z-zsivc
	dr = r-rsivc
c
        A    = dz**2 + r**2 + rsivc**2
	B    = -2.*r*rsivc
	rho1 = sqrt(dz*dz+dr*dr)
	rho2 = sqrt(dz*dz+(rsivc+r)**2)	
	ks    = 1.-(rho1/rho2)**2 
	k    = sqrt(ks) 
c
	f1 = (2./k-k)*ellk(ks)-2./k*elle(ks)
c
	vrsivcK = -gasivc*sqrt(rsivc/r)*f1/(2.*pi)
c
	return
	end
c
c----------------------------------------------------------------------------------
c
	real function vrsivcBo(z,r)
c
c	Bontempo Eq (2) axial induction form vortex cylinder
c
	use mem
	REAL k, ks, zs, rs, f0, f1, f2
c	
	zs  = (z-zsivc)/rsivc
	rs  = r/rsivc
c
	ks  = 4.*rs/(zs**2+(1.+rs)**2)
	k  = sqrt(ks)
c
	f0 = sqrt(zs**2+(1.+rs)**2)
c
c       Bontempo:  argument of ellliptic functions: k 
c	Br and vK: ks
c
	f1 = elle(ks)-(1.-0.5*ks)*ellk(ks)
	f2 = sqrt(zs**2+(1.+rs)**2)
c
c       Bontempo' prefactor 
c	in conflict with van Kuik Eq (D.8) from Branlard (2016) Eq (36.69)
c
c       sign changed to "+" (Bontempo gives "-")
c
	vrsivcBo = 2.*gasivc*f1/(pi*ks*f0)
c
c       2021 03 22: this is off around zsivc by some faktor
c
c
c       again Branlard
c
	f0        = (2.-ks)*ellk(ks)/k - 2.*elle(ks)/k
	vrsivcBo  = -gasivc*f0/(2.*pi*sqrt(rs))
c
	return
	end
c
c------------------------------------------------------------------------------------
c
	real function vznmBo(npt,n)
c
c	Bontempo Eq (4) induced axial velocity from vortex rings n 
c       to point P (npt)
c
	use mem
c
	integer npt,n
	real k, ks
c
	zp = 0.5*(zsl(npt)+zsl(npt-1))
        rp = 0.5*(rsl(npt)+rsl(npt-1))
c	
	zn = 0.5*(zsl(n)+zsl(n-1))
        rn = 0.5*(rsl(n)+rsl(n-1))
c
        ztil = (zp-zn)/rn
	rtil = rp/rn
c
	dz = zsl(n)-zsl(n-1)
	dr = rsl(n)-rsl(n-1)
	ds = sqrt(dz*dz+dr*dr)
c
	k = sqrt(4.*rtil/(ztil**2+(1.+rtil)**2))
	ks = k*k
c
c       2 pi (Bontempo)  
c
	f1 = -ga(n)*ds/(2.*pi*rn*sqrt(ztil**2+(1.+rtil)**2))
        f2 =  1.+2.*(rtil-1.)/(ztil**2+(rtil-1.)**2)
c
	vznmBo = f1*(ellk(ks)-f2*elle(ks))
c
	return
	end
c
c------------------------------------------------------------------------------------
c
	real function vznmBr(pt,n)
c
c	Branlard induced axial velocity from vortex ring n 
c       to point P (pt)
c	rn  -> r0 
c	rp  -> r
c
	use mem
c
	integer pt,n
	real ks, k, f1, f2
c
	zp = 0.5*(zsl(pt)+zsl(pt-1))
        rp = 0.5*(rsl(pt)+rsl(pt-1))
	zn = 0.5*(zsl(n)  +zsl(n-1))
        rn = 0.5*(rsl(n)  +rsl(n-1))	
c
	dz = zsl(n)-zsl(n-1)
	dr = rsl(n)-rsl(n-1)
	ds = sqrt(dz*dz+dr*dr)
c
	ks = 4.*rp*rn/((rn+rp)**2 + (zp-zn)**2)
c
c       Branlard Eq (35.8)
c
	f1 = rn/(2.*rp)*ks/(1.-ks)-(2.-ks)/(2.*(1.-ks))
	f1 = f1*elle(ks) + ellk(ks)
	f2 = -ga(n)*ds*sqrt(ks)*sqrt(rn/rp)/(4.*pi*rn)
c
	vznmBr = f2*f1
c
	return
	end
c
c------------------------------------------------------------------------------------
c
	real function vznmK(npt,n)
c
c	van Kuik Eq (D.1) induced axial velocity from vortex rings n 
c       to point P (npt)
c
c      replace k to ks
c
	use mem
c
	integer npt,n
	real k, ks, I1, I2
c
        zp = 0.5*(zsl(npt)+zsl(npt-1))
        rp = 0.5*(rsl(npt)+rsl(npt-1))
        zn = 0.5*(zsl(n)+zsl(n-1))
        rn = 0.5*(rsl(n)+rsl(n-1))	
c
	dz = zp-zn
	dr = rp-rn
c
        dzp= zsl(n)-zsl(n-1)
        drp= rsl(n)-rsl(n-1)
	dsp = sqrt(dzp**2+drp**2)
c
        A = dz**2 + rp**2 + rn**2
	B = -2.*rp*rn
	rho1 = sqrt(dz*dz+dr*dr)
	rho2 = sqrt(dz**2+(rp+rn)**2)	
	ks   = 1.-(rho1/rho2)**2
	I1   = (4./rho2)*ellk(ks)
	I2   = (4.*elle(ks))/((1.-ks)*rho2**3)
c
	f1 = (rn+rp*A/B)*I2
	f2 = rp*I1/B
c
c       changed sign
c       dsp
c
c	vK -> 4 pi

	vznmK = -ga(n)*dsp*rn*(f1-f2)/(4.*pi)
c
	return
	end
c
c------------------------------------------------------------------------------------
c
	real function vrnmBo(npt,n)
c
c	Bontempo Eq (5) calculated induced radial velocity 
c	at a Point P (npt) form all other sheets (n ne npt)
c
	use mem
c
	integer n,npt
	real k,ks
c	
  	zp = 0.5*(zsl(npt)+zsl(npt-1))
        rp = 0.5*(rsl(npt)+rsl(npt-1))
c
 	zn = 0.5*(zsl(n)+zsl(n-1))
        rn = 0.5*(rsl(n)+rsl(n-1))	
c
        ztil = (zp-zn)/rn
	rtil = rp/rn
c
	dz = zsl(n)-zsl(n-1)
	dr = rsl(n)-rsl(n-1)
	ds = sqrt(dz*dz+dr*dr)
c
	ks = 4.*rtil/(ztil**2+(1.+rtil)**2)
	k  = sqrt(ks)
c
	f1 = ga(n)*ds*ztil/rtil
c
        f1 = f1/(2.*pi*rsl(n)*sqrt(ztil**2+(1.+til)**2))
	f2 = ellk(ks)-elle(ks)*(1.+2.*rtil/(ztil**2+(rtil-1.)**2))
c
	vrnmBo =f1*f2
c
	return
	end
c
c------------------------------------------------------------------------------------
c
	real function vrnmBr(npt,n)
c
c       induction of radial velocity from n to npt
c       n   ->      r0
c	npt -> P -> r
c
	use mem
c
	integer n,npt
	real ks,k,f1, f2
c	
	zp = 0.5*(zsl(npt)+zsl(npt-1))
        rp = 0.5*(rsl(npt)+rsl(npt-1))
	zn = 0.5*(zsl(n  )+zsl(n  -1))
        rn = 0.5*(rsl(n  )+rsl(n  -1))	
c
	dz = zsl(n)-zsl(n-1)
	dr = rsl(n)-rsl(n-1)
	ds = sqrt(dz*dz+dr*dr)
c
	ks = 4.*rp*rn/((rn+rp)**2+(zp-zn)**2)
	k  = sqrt(ks)
c
c       Branlard Eq (35.7) 
c
	f1 = (2.-ks)/(2.*(1.-ks))
        f1 = f1*elle(ks) - ellk(ks)
	f2 = -ga(n)*ds*k*(zp-zn)*((rn/rp))**1.5/(4.*pi*rn*rn)
c
	vrnmBr = f1*f2
c
	return
	end
c
c---------------------------------------------------------------------------------------
c
	real function vrnmK(npt,n)
c
c	van Kuik Eq (D.1) calculated induced radial velocity 
c	at a Point P (npt) form all other sheets (n ne npt)
c
	use mem
c
	integer n,npt
	real k,den,I1,I2
c	
  	zp = 0.5*(zsl(npt)+zsl(npt-1))
        rp = 0.5*(rsl(npt)+rsl(npt-1))
 	zn = 0.5*(zsl(n)  +zsl(n-1))
        rn = 0.5*(rsl(n)  +rsl(n-1)) 
c
	dzp = zsl(n)-zsl(n-1)
	drp = rsl(n)-rsl(n-1)
	dsp = sqrt(dzp**2+drp**2)
c
	dz = zp-zn
	dr = rp-rn
c
        A = dz**2 + rp**2 + rn**2
	B = -2.*rp*rn
	rho1 = sqrt(dz**2+dr**2)
	rho2 = sqrt(dz**2+(rp+rn)**2)	
	k    = sqrt(1.-(rho1/rho2)**2) 
	I1 = 4.*ellk(k)/rho2
	I2 = 4.*elle(k)/((1.-k)*rho2**3)
c
	f1 = dz/B
	f2 = I1-A*I2
	vrnmK = -ga(n)*dsp*rn*(f1*f2)/(2.*pi)
c
	return
	end
c
c---------------------------------------------------------------------------
c
	real function vzself(m)
c
c	Bontempo Eq (6) axial self induction
c
	use mem
c
	real f1,f2,dz,dr,ds,cosb,rslm
	real betam,betamm1,betamp1
	real k
c
	dz = zsl(m)-zsl(m-1)
	dr = rsl(m)-rsl(m-1)
	ds = sqrt(dz*dz+dr*dr)
c
	betam = dr/dz
c
	if (m.eq.np) then
  	   dz = zsl(m)-zsl(m-1)
	   dr = rsl(m)-rsl(m-1)
	else
  	   dz = zsl(m+1)-zsl(m)
	   dr = rsl(m+1)-rsl(m)
	end if
	betamp1 = dr/dz
c	
	if (m.eq.1) then
	   betamm1 = betam
	else
	   dz = zsl(m-1)-zsl(m-2)
	   dr = rsl(m-1)-rsl(m-2)
	   betamm1 = dr/dz
	endif
c
	cosb = 1./(sqrt(1.+betam*betam))
	rslm = 0.5*(rsl(m)+rsl(m-1))
c
	f1 = -cosb*(betamp1-betamm1)/(8.*pi)
	f2 = -ds/(4.*pi*rslm)*(-0.25+log(8.*pi*rsl(m)/ds))
c
        vzself = (f1 + f2)*ga(m)
c
	return
	end
c---------------------------------------------------------------------------
c
	real function vrself(m)
c
c	Bontempo Eq (7) radial self induction
c
	use mem
c
	real betam,betamm1,betamp1,dz,dr
	integer m
c
	dz = zsl(m)-zsl(m-1)
	dr = rsl(m)-rsl(m-1)
c
	betam = dr/dz
c
	if (m.eq.np) then
  	   dz = zsl(m)-zsl(m-1)
	   dr = rsl(m)-rsl(m-1)
	else
  	   dz = zsl(m+1)-zsl(m)
	   dr = rsl(m+1)-rsl(m)
	end if
	betamp1 = dr/dz
c	
	if (m.eq.1) then
	   betamm1 = betam
	else
	   dz = zsl(m-1)-zsl(m-2)
	   dr = rsl(m-1)-rsl(m-2)
	   betamm1 = dr/dz
	endif
c
	sinb =betam/(sqrt(1.+betam*betam))
c
	vrself = -(betamp1-betamm1)*ga(m)*sinb/(8.*pi)	
c
	return
	end
c
c-----------------------------------------------------------------------
c
	real function vzd(r)
c
c	induction of all vortex rings  to points 
c	at the disk z = 0  to evaluate power
c
	use mem
c
	real k,vzdd,r,ks,f,f1,f2,I1,I2
c
	vzd  = 0.
	vzdd = 0.	
    	zp   = 0. 
        rp   = r
c
	do n = 1,np
           select case(model)
              case(1)
		zn = 0.5*(zsl(n)+zsl(n-1))
        	rn = 0.5*(rsl(n)+rsl(n-1))
c
        	ztil = (zp-zn)/rn
		rtil = rp/rn
c
		dz = zsl(n)-zsl(n-1)
		dr = rsl(n)-rsl(n-1)
		ds = sqrt(dz*dz+dr*dr)
c
		k = sqrt(4.*rtil/(ztil**2+(1.+rtil)**2))
c
		f1 = -ga(n)*ds/(2.*pi*rn*sqrt(ztil**2+(1.+rtil)**2))
        	f2 =  1.+2.*(rtil-1.)/(ztil**2+(rtil-1.)**2)
c
                vzdd = vzdd +  f1*(ellk(k)-f2*elle(k))
c	
	      case(2)
        	zn = 0.5*(zsl(n)+zsl(n-1))
        	rn = 0.5*(rsl(n)+rsl(n-1))	
c
		dz = zp-zn
		dr = rp-rn
c
        	dzp= zsl(n)-zsl(n-1)
        	drp= rsl(n)-rsl(n-1)
		dsp = sqrt(dzp**2+drp**2)
c
        	A = dz**2 + rp**2 + rn**2
		B = -2.*rp*rn
		rho1 = sqrt(dz*dz+dr*dr)
		rho2 = sqrt(dz**2+(rp+rn)**2)	
		k    = sqrt(1.-(rho1/rho2)**2) 
		I1 = 4./rho2*ellk(k)
		I2 = 4.*elle(k)/((1.-k)*rho2**3)
c
		f1 = (rn+rp*A/B)*I2
		f2 = rp*I1/B
c
		vzdd = vzdd-ga(n)*dsp*rn*(f1-f2)/(4.*pi)    
c           
	      case(3)
		zp = 0.
 	        rp = r
		zn = 0.5*(zsl(n)  +zsl(n-1))
  	        rn = 0.5*(rsl(n)  +rsl(n-1))	
c
		dz = zsl(n)-zsl(n-1)
		dr = rsl(n)-rsl(n-1)
		ds = sqrt(dz*dz+dr*dr)
c
		ks = 4.*rp*rn/((rn+rp)**2 + (zp-zn)**2)
c
c       	 Branlard Eq (35.8)  
c
		f1 = rn/(2.*rp)*ks/(1.-ks)-(2.-ks)/(2.*(1.-ks))
		f1 = f1*elle(ks) + ellk(ks)
		f2 = -ga(n)*ds*sqrt(ks)*sqrt(rn/rp)/(4.*pi*rn)
c
		vzdd = vzdd + f2*f1
c	               
	    end select 	
	end do	
c        
	vzd    = vzdd + vzdcyl(r)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	real function vrd(r)
c
c	induction of all vortex rings  to points 
c	at the disk z = 0  to evaluate power
c
	use mem
c
	real k,ks,vzdd,f1,f2
c
c	use mem
c
	vrd  = 0.
	vrdd = 0.	
c
	do n = 1,np
c
	   zn = 0.5*(zsl(n)+zsl(n-1))
	   rn = 0.5*(rsl(n)+rsl(n-1))
	   zp = 0.
           rp = r
c
           ztil = -zn/rn
	   rtil = r/rn
c
	   dz = zsl(n)-zsl(n-1)
	   dr = rsl(n)-rsl(n-1)
	   ds = sqrt(dz*dz+dr*dr)
c
c	   write(*,*)'model ',model
c
           select case(model)
              case(1)
	         k  = sqrt(4.*rtil/(ztil**2+(1.+rtil)**2))
  	         f1 = ga(n)*ds*ztil/rtil
c
c                Bontempo 2 pi probally wrong (4 pi, 2021 jan 25)
c
                 f1 = f1/(4.*pi*rsl(n)*sqrt(ztil**2+(1.+til)**2))
	         f2 = ellk(k)-elle(k)*(1.+2.*rtil/(ztil**2+(rtil-1.)**2))
c
	         vrdd = vrdd + f1*f2 
	      case(2)
c                
	      case(3)
		ks = 4.*rp*rn/((rn+rp)**2+(zp-zn)**2)
		k  = sqrt(ks)
c
c       	Branlard Eq (35.7) 
c
		f1 = (2.-ks)/(2.*(1.-ks))
	        f1 = f1*elle(ks) - ellk(ks)
		f2 = -ga(n)*ds*k*(zp-zn)*(rn/rp)**1.5/(4.*pi*rn*rn)
c
		vrdd = vrdd + f1*f2
c	               
	    end select 	
	end do	
c
	vrd    = vrdd + vrdcyl(r)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	real function vzdcyl(r)
c
c	Bontempo Eq (2) induction of vortex cylinder 
c	at the disk z = 0  to evaluate power
c
	use mem
c
	real k,vzd,rc,ks,f,f1,f2
c
	vzdcyl = 0.	
	rc = r
c
c       Contribution from cylinder
c
        select case(model)
           case(1)
              vzdcyl   = vzsivcBo(0.,rc)
           case(2)
              vzdcyl   = vzsivcK (0.,rc)
           case(3)
              vzdcyl   = vzsivcBr(0.,rc)
        end select 	
c
c	write(*,*)'r vz',r, vzdcyl
c
	return
	end
c-----------------------------------------------------------------------
c
	real function vrdcyl(r)
c
c	Bontempo Eq  induction of vortex cylinder 
c	at the disk z = 0  
c
	use mem
c
	real k,vrd,r,ks,f,f1,f2
c	real rsl(0:np), zsl(0:np), ga(np)
c	real vz(np), vr(np)
c
	vrdcyl = 0.	
c
c       Contribution from cylinder
c
        select case(model)
           case(1)
              vrd   = vrsivcBo(0.,r)
           case(2)
              vrd   = vrsivcK(0.,r)
           case(3)
              vrd   = vrsivcBr(0.,r)
        end select 	
c
	vrdcyl = vrd
c
	return
	end
c
c-----------------------------------------------------------------------
c
	real function vza(z)
c
c       calculates induced axial velocity at r = 0 (axis)
c
	use mem
c
	real k, vzad
c
	vzad = 0.
	zcyl = z
c
	do n = 1,np
c
	   zref = 0.5*(zsl(n)+zsl(n-1))
	   rref = 0.5*(rsl(n)+rsl(n-1))

           ztil = (zcyl-zref)/rref
	   rtil = 0.
c
	   dz = zsl(n)-zsl(n-1)
	   dr = rsl(n)-rsl(n-1)
	   ds = sqrt(dz*dz+dr*dr)
c
  	   k = 0.
c
	   f1 = -ga(n)*ds/(2.*pi*rsl(n)*sqrt(ztil**2+(1.+rtil)**2))
	   f2 = 1.+2.*(rtil-1.)/(ztil**2+(rtil-1.)**2)
	   f2 = ellk(k)-f2*elle(k)
c
	   vzad = vzad + f1*f2
c
	end do
c
	vza  = vzad
c	vzy = vzad
c
	return
	end
c
c-----------------------------------------------------------------------
c
	real function psiBr(z,r)
c
c	 calculates streamfunction
c
	use mem
c
	real k, ks, f1, f2, rp, zp, ds,zr,zrel
c
	psiBr = 0.	
c
        rp = r
	zp = z
c
	do n = 1,np
c
	   zr = zsl(n)
	   rr = rsl(n)
c
	   zrel = zr-zp
c
	   dz = zsl(n)-zsl(n-1)
	   dr = rsl(n)-rsl(n-1)
	   ds = sqrt(dz*dz+dr*dr)
c
  	   ks = 4.*rp*rr/((rp+rr)**2+zrel**2)
	   k  = sqrt(ks)
c
c	   Eq(35.5) branlard: PSI = r * psi_theta ! 
c
	   f1 = -ga(n)*ds*sqrt(rp*rr)/(2.*pi)
	   f2 = (2./k-k)*ellk(ks)-(2./k)*elle(ks)
c
	   psiBr =  psiBr + f1*f2
c
	end do
c
c       contribution from inflow
c
	psiBr =  psiBr + 0.5*rp*rp
c
	return
	end
c
c-----------------------------------------------------------------------
c
	real function psiRingBr(n)
c
c	 calculates streamfunction of ring n at
c
	use mem
c
	real k, ks, f1, f2, rp, zp, ds,zr,zrel
	integer n
c
	zr = 0.5*(zsl(n)+zsl(n-1))
	rr = 0.5*(rsl(n)+rsl(n-1))
c
        dz = zsl(n)-zsl(n-1)
	dr = rsl(n)-rsl(n-1)
	ds = sqrt(dz*dz+dr*dr)
	rp = 0.5*ds
c
  	ks = 4.*rp*rr/((rr+rp)**2+zr**2)
	k  = sqrt(ks)
c
c	Eq(35.5) branlard: PSI = r * psi ! 
c
	f1 = -ga(n)*ds*sqrt(rp*rr)/(2.*pi)
	f2 = (2./k-k)*ellk(ks)-(2./k)*elle(ks)
c
	psiRingBr = f1*f2
c	write(*,*)n,psiringbr
c
	return
	end
c
c----------------------------------------------------------------
c
	real function psiBrCyl(z,r)
c
c	calculates streamfunction for a semi infinite cylinder
c	we use Eq. (36.39) from Emmanuel, E-mail of 2021 02 08
c
	use mem
c
	real m,m0,a,f1,f2,f3,f4,f
	real rc, zc, gam
c
	gam = -gasivc
c
	f = 0.	
	psiBrcyl = 0.
c
        zc  = z-zsivc
	rc  = r
c
        m   = 4.*rc*rsivc/((rc+rsivc)**2+(zc**2))
        m0  = 4.*rc*rsivc/((rc+rsivc)*(rc+rsivc))
c
c       psi = r * PSI_theta
c
	a  = gam*zc*sqrt(m*rsivc*r)/(2.*pi*m0)

	f1 = (1.-m0+(m0/m))*ellk(m)
c
	f2 = -(m0/m)*elle(m)	
c
        f3 = (m0-1.)*ellp(m0,m)
c
	f  = a*(f1 +f2 +f3)  
c
	if(rc.gt.rsivc)then
	   f4 = 0.25*gam*rsivc**2
        else
           f4 = 0.25*gam*r**2
        end if
c
	psiBrcyl = f + f4
c
c
c------------------------------------------------------------------
c       (debug)
c	write(*,*)'1 und 2',psiBrcyl2, psiBrcyl1
c
c	write(777,778)'z r m m0 f1 f2 f3 f4',
c     +                z,r,'**',m,m0,'**',f1,f2,'**',
c     +                f3,f4,'**',psiBrCyl
c778	format(a22,4(2e12.4,a2),e12.4)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	real function psivK(z,r)
c
c	 calculates streamfunction from all rings
c	 we use vK Eq. (D.1.), 3rd
c
	use mem
c
	real k, ks, zr, rr, rho1, rho2, psir, rc, zc
c
	psivK = 0.	
	psir  = 0.
c
	rc = r
	zc = z
c
	do n = 1,np
	   zr = 0.5*(zsl(n)+zsl(n-1))
	   rr = 0.5*(rsl(n)+rsl(n-1))
c
	   dz = zsl(n)-zsl(n-1)
	   dr = rsl(n)-rsl(n-1)
	   ds = sqrt(dz*dz+dr*dr)
c
	   rho1 = sqrt((zc-zr)**2+(rc-rr)**2)
	   rho2 = sqrt((zc-zr)**2+(rc+rr)**2) 	
           ks   = 1.-(rho1/rho2)**2
           k    = sqrt(ks)
           psir  = ga(n)*ds*sqrt(rc*rr)/(2.*pi)
	   psir  = psir*((2./k -k)*ellk(ks)-2./k*elle(ks))

	   psivK = psivK + psir
c
	end do
c
	return
	end
c
c------------------------------------------------------------------------------------
c
	real function vzacylan(z)
c
c       induced axial velocity at axis, r = 0 
c       from a cylinder with gasivc at  z = 0
c	radius of cylinder = rsivc
c	Eq. (D.9) van Kuik
c
	use mem
	real f1
c
	f1 = sqrt(rsivc**2 + z**2)
	f1 = z/f1
	vzacylan = 0.5*gasivc*(1. + f1)
c
	return
	end
c
cEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
c--------------------------------------------------------------------------------------
cEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
c
	real function ellk(k)
	integer ie1
	real k, ks
c                       2
c       k(k) = rfl(0,1-k ,1)
c
	ks = k
c
	ellk = rfl(0.,1.-k,1.,ie1)
	if(ie1.ne.0)write(99,*)'error ellk'
c
	return
	end
c--------------------------------------------------------------------------------------
c
	real function elle(k)
	integer ie1,ie2
	real k, ks
c                       2           2         2
c       e(k) = rfl(0,1-k ,1) - 1/3 k rdl(0,1-k ,1)
c
	ks = k
c
        elle = rfl(0.,1.-ks,1.,ie1) - ks/3.*rdl(0.,1.-ks ,1.,ie2)	 
	if(ie1.ne.0.or.ie2.ne.0)write(99,*)'error elle'
c
	return
	end
c
c------------------------------------------------------------------------------------------
c
	real function ellp(n,k)
	real k, n
	integer ie1,ie2
c
c                  legendre form of elliptic integral of 3rd kind
c                  ----------------------------------------------
c
c
c                               phi         2         -1
c                  p(phi,k,n) = int (1+n sin (theta) )   *
c                                0
c
c                                      2    2         -1/2
c                                 *(1-k  sin (theta) )     d theta
c
c
c                                               2          2   2
c                            = sin (phi) rfl(cos (phi), 1-k sin (phi),1)
c
c                                     3             2         2   2
c                             -n/3 sin (phi) rjl(cos (phi),1-k sin (phi),
c
c                                      2
c                             1,1+n sin (phi))
c
c	2021 02 11: Branlard Eq (C.66) n -> (-n)
c

        ellp = rfl(0.,1.-k,1.,ie1) + n/3.*rjl(0.,1.-k ,1.,1.-n,ie2)
	if(ie1.ne.0.or.ie2.ne.0)write(99,*)'error ellp'
c
	return
	end
c
cLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
c
c
CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+              Lawrence Livermore National Laboratory              +CC
CC+                                                                  +CC
CC+     Livermore Computing Center  Mathematical Software Library    +CC
CC+    Mathematical Software Service Library (MSSL) -- Class Three   +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+  Each Class Three routine is guaranteed only to meet minimal     +CC
CC+  documentation and programming standards.  Although the quality  +CC
CC+  of a particular routine may be high, this is not assured.  A    +CC
CC+  Class Three routine is recommended only in cases where no       +CC
CC+  suitable routine exists in other libraries (e.g., MSSL or       +CC
CC+  SLATEC).                                                        +CC
CC+                                                                  +CC
CC+  These routines are distributed exclusively for use in support   +CC
CC+  of LLNL programs.  Check with the LLNL Code Release Center or   +CC
CC+  the LC Client Services HotLine, (510)422-4531, before moving    +CC
CC+  this source code to a non-LLNL system.                          +CC
CC+                                                                  +CC
CC+  Support for Class Three routines is not guaranteed, but is      +CC
CC+  available in many cases.                                        +CC
CC+                                                                  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +                        N O T I C E                         +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +  This report was prepared as an account of work sponsored  +  +CC
CC+  +  by the United States government.  Neither the United      +  +CC
CC+  +  States government nor any of their employees, nor any of  +  +CC
CC+  +  their contractors, subcontractors, or their employees,    +  +CC
CC+  +  makes any warranty, expressed or implied, or assumes any  +  +CC
CC+  +  legal liability or responsibility for the accuracy,       +  +CC
CC+  +  completeness or usefulness of any information, apparatus, +  +CC
CC+  +  product or process disclosed, or represents that its use  +  +CC
CC+  +  would not infringe privately-owned rights.                +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+                                                                  +CC
CC+  Please refer questions regarding this software to the LC        +CC
CC+  Client Services Hotline, (510)422-4531, or use the NMG "mail"   +CC
CC+  command.                                                        +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC


      real function rfl(x,y,z,ier) 
c
c***begin prologue  rfl
c***date written   790801   (yymmdd)
c***revision date  831228   (yymmdd)
c***category no.  c14
c***keywords  duplication theorem,elliptic integral,incomplete,complete,
c             integral of the first kind,taylor series
c***author  carlson, b.c., ames laboratory-doe
c             iowa state university, ames, iowa  50011
c           notis, e.m., ames laboratory-doe
c             iowa state university, ames, iowa  50011
c           pexton, r.l., lawrence livermore national laboratory
c             livermore, california  94550
c***purpose  incomplete or complete elliptic integral of the 1st kind.
c            for x, y, and z nonnegative and at most one of them zero,
c            rfl(x,y,z) = integral from zero to infinity of
c                                -1/2     -1/2     -1/2
c                      (1/2)(t+x)    (t+y)    (t+z)    dt,
c            if one of them is zero, the integral is complete.
c***description
c
c   1.     rfl
c          evaluate an incomplete (or complete) elliptic integral
c          of the first kind
c          standard fortran function routine
c          single precision version
c          the routine calculates an approximation result to
c          rfl(x,y,z) = integral from zero to infinity of
c
c                               -1/2     -1/2     -1/2
c                     (1/2)(t+x)    (t+y)    (t+z)    dt,
c
c          where x, y, and z are nonnegative and at most one of them
c          is zero.  if one of them is zero, the integral is complete.
c          the duplication theorem is iterated until the variables are
c          nearly equal, and the function is then expanded in taylor
c          series to fifth order.
c
c   2.     calling sequence
c          rfl( x, y, z, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - single precision, nonnegative variable
c
c          y      - single precision, nonnegative variable
c
c          z      - single precision, nonnegative variable
c
c
c
c          on return     (values assigned by the rfl routine)
c
c          rfl     - single precision approximation to the integral
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine.  it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c          x, y, z are unaltered.
c
c
c   3.    error messages
c
c         value of ier assigned by the rfl routine
c
c                  value assigned         error message printed
c                  ier = 1                amin1(x,y,z) .lt. 0.0e0
c                      = 2                amin1(x+y,x+z,y+z) .lt. lolim
c                      = 3                amax1(x,y,z) .gt. uplim
c
c
c
c   4.     control parameters
c
c                  values of lolim,uplim,and errtol are set by the
c                  routine.
c
c          lolim and uplim determine the valid range of x, y and z
c
c          lolim  - lower limit of valid arguments
c
c                   not less than 5 * (machine minimum).
c
c          uplim  - upper limit of valid arguments
c
c                   not greater than (machine maximum) / 5.
c
c
c                     acceptable values for:   lolim      uplim
c                     ibm 360/370 series   :   3.0e-78     1.0e+75
c                     cdc 6000/7000 series :   1.0e-292    1.0e+321
c                     univac 1100 series   :   1.0e-37     1.0e+37
c                     cray                 :   2.3e-2466   1.09e+2465
c                     vax 11 series        :   1.5e-38     3.0e+37
c
c
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the routine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c
c
c          errtol - relative error due to truncation is less than
c                   errtol ** 6 / (4 * (1-errtol)  .
c
c
c
c              the accuracy of the computed approximation to the inte-
c              gral can be controlled by choosing the value of errtol.
c              truncation of a taylor series after terms of fifth order
c              introduces an error less than the amount shown in the
c              second column of the following table for each value of
c              errtol in the first column.  in addition to the trunca-
c              tion error there will be round-off error, but in prac-
c              tice the total error from both sources is usually less
c              than the amount given in the table.
c
c
c
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    3.0e-19
c                           3.0e-3    2.0e-16
c                           1.0e-2    3.0e-13
c                           3.0e-2    2.0e-10
c                           1.0e-1    3.0e-7
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c***long description
c
c   rfl special comments
c
c
c
c          check by addition theorem: rfl(x,x+z,x+w) + rfl(y,y+z,y+w)
c          = rfl(0,z,w), where x,y,z,w are positive and x * y = z * w.
c
c
c          on input:
c
c          x, y, and z are the variables in the integral rfl(x,y,z).
c
c
c          on output:
c
c
c          x, y, and z are unaltered.
c
c
c
c          ********************************************************
c
c          warning: changes in the program may improve speed at the
c          expense of robustness.
c
c
c
c   special functions via rfl
c
c
c                  legendre form of elliptic integral of 1st kind
c                  ----------------------------------------------
c
c
c                                             2         2   2
c                  f(phi,k) = sin(phi) rfl(cos (phi),1-k sin (phi),1)
c
c
c                                  2
c                  k(k) = rfl(0,1-k ,1)
c
c
c
c                  bulirsch form of elliptic integral of 1st kind
c                  ----------------------------------------------
c
c
c                                          2 2    2
c                  el1(x,kc) = x rfl(1,1+kc x ,1+x )
c
c
c
c
c                  lemniscate constant a
c                  ---------------------
c
c
c                       1      4 -1/2
c                  a = int (1-s )    ds = rfl(0,1,2) = rfl(0,2,1)
c                       0
c
c
c    -------------------------------------------------------------------
c          subroutines or functions needed
c              - errpex
c              - r1pext
c              - fortran abs,amax1,amin1,sqrt
c***references  carlson, b.c. and notis, e.m.
c                 algorithms for incomplete elliptic integrals
c                 acm transactions on mathematical software,vol.7,no.3,
c                 sept, 1981, pages 398-403
c               carlson, b.c.
c                 computing elliptic integrals by duplication
c                 numer. math. 33, (1979), 1-16
c               carlson, b.c.
c                 elliptic integrals of the first kind
c                 siam j. math. anal. 8 (1977), 231-242
c***routines called  r1pext,errpex
c***end prologue  rfl
      integer ier,itodo
      real lolim, uplim, epslon, errtol
      real c1, c2, c3, e2, e3, lamda
      real mu, s, x, xn, xndev
      real xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
c
      character msg*80
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     the original routine is designed for cft
c
c
c     for  civic  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(10)
c
c
c     for  chat  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(8)
c
c
c      other changes are listed subsequently.
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c
c
c
      data itodo/1/
c
c
c
c***first executable statement  rfl
      if(itodo.eq.1)then
c
c
      errtol=(4.0*r1pext(3))**(1.0/6.0)
c
c
      lolim = 5.0e0 * (r1pext(1) )
c
      uplim = (r1pext(2) ) / 5.0e0
	uplim = 1.e38
c	
c      write(*,*)'*****************************'
c      write(*,*)'lolim uplim ',lolim, uplim
c      write(*,*)'*****************************'

c
      c1=1.0e0/24.0e0
c
      c2=3.0e0/44.0e0
c
      c3=1.0e0/14.0e0
c
c
      itodo=0
c
      end if
c
c
c         call error handler if necessary.
c
c
    5 rfl=0.0e0
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     for civic users insert the following three statements:
c
c
c          do 9191 igogo=1,10
c
c          msg(igogo)='        '
c
c9191 continue
c
c
c     for chat users insert the following three statements:
c
c
c          do 9191 igogo=1,8
c
c          msg(igogo)='          '
c
c9191 continue
c
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c  anf ab
c
      if( amin1(x,y,z).lt.0.0e0) then
      ier=1
      msg= 'rfl - error: amin1(x,y,z).lt.0.0e0 (where x=r1, y=r2, z=r3)'
      call errpex (msg,ier,x,y,z,0.0,0.0)
      return
      end if
      if (amin1(x+y,x+z,y+z).lt.lolim) then
	write(99,*)'lolim ',lolim
      ier=2
      msg='rfl-error:amin1(x+y,x+z,y+z).lt.lolim(where x=r1,y=r2,z=r3)'
      call errpex (msg,ier,x,y,z,0.0,0.0)
      return
      end if
      if (amax1(x,y,z).gt.uplim) then
      ier=3
	write(*,*)'uplim ',uplim
      msg='rfl - error: amax1(x,y,z).gt.uplim (where x=r1, y=r2, z=r3)'
      call errpex (msg,ier,x,y,z,0.0,0.0)
      return
      end if
c
c   end ab
c
   20 ier = 0
      xn = x
      yn = y
      zn = z
c
c
c
   30 mu = (xn+yn+zn)/3.0e0
      xndev = 2.0e0 - (mu+xn)/mu
      yndev = 2.0e0 - (mu+yn)/mu
      zndev = 2.0e0 - (mu+zn)/mu
      epslon = amax1( abs(xndev), abs(yndev), abs(zndev))
      if (epslon.lt.errtol) go to 40
      xnroot =  sqrt(xn)
      ynroot =  sqrt(yn)
      znroot =  sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      zn = (zn+lamda)*0.250e0
      go to 30
c
c
c
   40 e2 = xndev*yndev - zndev*zndev
      e3 = xndev*yndev*zndev
      s = 1.0e0 + (c1*e2-0.10e0-c2*e3)*e2 + c3*e3
      rfl = s/ sqrt(mu)
c
   50 return
c
c
c
      end
      subroutine errpex(messg,ierr,r1,r2,r3,r4,r5)    
c
c
      character messg*80
c
c
c
c***first executable statement  errpex
c
c
      write (77,101) messg
  101 format (a)
c
      write(77,102) ierr,r1,r2,r3,r4,r5
  102 format (2x,i3,3x,2e16.7,/,3x,3e16.7)
      return
      end
c
c     machine constants
c
      real function r1pext(i)      
c
c
      character msg*80
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
c
      real rmach(5)
c
      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))
c
c     machine constants as given
c
       rmach(1)  =2.0**(-149)
       rmach(2)  =2.0**(127)
c
c
c     machine constants for the ketos, real*8, 2020 nov 20
c
       rmach(1)  =2.0**(-1023)
       rmach(2)  =2.0**( 1020)
c
       rmach(3)  =2.0**(-33)
       rmach(4)  =2.0**(-32)
       rmach(5)  = 0.30102999566
c
c
c     machine constants for the cdc 6000/7000 series.
c
c     data rmach(1) / 00014000000000000000b /
c     data rmach(2) / 37767777777777777777b /
c     data rmach(3) / 16404000000000000000b /
c     data rmach(4) / 16414000000000000000b /
c     data rmach(5) / 17164642023241175720b /
c
c     machine constants for the cray 1
c
c      data rmach(1) / 200004000000000000000b /
c      data rmach(2) / 577777777777777777777b /
c      data rmach(3) / 377214000000000000000b /
c      data rmach(4) / 377224000000000000000b /
c      data rmach(5) / 377774642023241175720b /
c
      if (i .lt. 1  .or.  i .gt. 5)  then
         msg=  'r mach -- i out of bounds'
         call errpex (msg , i, 0.0, 0.0, 0.0, 0.0, 0.0)
      end if
c
      r1pext = rmach(i)
      return
c
      end
c
CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+              Lawrence Livermore National Laboratory              +CC
CC+                                                                  +CC
CC+     Livermore Computing Center  Mathematical Software Library    +CC
CC+    Mathematical Software Service Library (MSSL) -- Class Three   +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+  Each Class Three routine is guaranteed only to meet minimal     +CC
CC+  documentation and programming standards.  Although the quality  +CC
CC+  of a particular routine may be high, this is not assured.  A    +CC
CC+  Class Three routine is recommended only in cases where no       +CC
CC+  suitable routine exists in other libraries (e.g., MSSL or       +CC
CC+  SLATEC).                                                        +CC
CC+                                                                  +CC
CC+  These routines are distributed exclusively for use in support   +CC
CC+  of LLNL programs.  Check with the LLNL Code Release Center or   +CC
CC+  the LC Client Services HotLine, (510)422-4531, before moving    +CC
CC+  this source code to a non-LLNL system.                          +CC
CC+                                                                  +CC
CC+  Support for Class Three routines is not guaranteed, but is      +CC
CC+  available in many cases.                                        +CC
CC+                                                                  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +                        N O T I C E                         +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +  This report was prepared as an account of work sponsored  +  +CC
CC+  +  by the United States government.  Neither the United      +  +CC
CC+  +  States government nor any of their employees, nor any of  +  +CC
CC+  +  their contractors, subcontractors, or their employees,    +  +CC
CC+  +  makes any warranty, expressed or implied, or assumes any  +  +CC
CC+  +  legal liability or responsibility for the accuracy,       +  +CC
CC+  +  completeness or usefulness of any information, apparatus, +  +CC
CC+  +  product or process disclosed, or represents that its use  +  +CC
CC+  +  would not infringe privately-owned rights.                +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+                                                                  +CC
CC+  Please refer questions regarding this software to the LC        +CC
CC+  Client Services Hotline, (510)422-4531, or use the NMG "mail"   +CC
CC+  command.                                                        +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC
      real function rdl(x,y,z,ier)  
c***begin prologue  rdl
c***date written   790801   (yymmdd)
c***revision date  831228   (yymmdd)
c***category no.  c14
c***keywordls  duplication theorem,elliptic integral,incomplete,complete,
c             integral of the second kind,taylor series
c***author  carlson, b.c., ames laboratory-doe
c             iowa state university, ames, iowa  50011
c           notis, e.m., ames laboratory-doe
c             iowa state university, ames, iowa  50011
c           pexton, r.l., lawrence livermore national laboratory
c             livermore, california  94550
c***purpose  incomplete or complete elliptic integral of the 2nd kind.
c            for x and y nonnegative, x+y and z positive,
c            rdl(x,y,z) = integral from zero to infinity of
c                                -1/2     -1/2     -3/2
c                      (3/2)(t+x)    (t+y)    (t+z)    dt.
c            if x or y is zero, the integral is complete.
c***description
c
c   1.     rdl
c          evaluate an incomplete (or complete) elliptic integral
c          of the second kind
c          standardl fortran function routine
c          single precision version
c          the routine calculates an approximation result to
c          rdl(x,y,z) = integral from zero to infinity of
c                              -1/2     -1/2     -3/2
c                    (3/2)(t+x)    (t+y)    (t+z)    dt,
c          where x and y are nonnegative, x + y is positive, and z is
c          positive.  if x or y is zero, the integral is complete.
c          the duplication theorem is iterated until the variables are
c          nearly equal, and the function is then expanded in taylor
c          series to fifth order.
c
c   2.     calling sequence
c
c          rdl( x, y, z, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - single precision, nonnegative variable
c
c          y      - single precision, nonnegative variable
c
c                   x + y is positive
c
c          z      - real, positive variable
c
c
c
c          on return     (values assigned by the rdl routine)
c
c          rdl     - real approximation to the integral
c
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine.  it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c
c          x, y, z are unaltered.
c
c   3.    error messages
c
c         value of ier assigned by the rdl routine
c
c                  value assigned         error message printed
c                  ier = 1                amin1(x,y) .lt. 0.0e0
c                      = 2                amin1(x + y, z ) .lt. lolim
c                      = 3                amax1(x,y,z) .gt. uplim
c
c
c   4.     control parameters
c
c                  values of lolim,uplim,and errtol are set by the
c                  routine.
c
c          lolim and uplim determine the valid range of x, y, and z
c
c          lolim  - lower limit of valid arguments
c
c                    not less  than 2 / (machine maximum) ** (2/3).
c
c          uplim  - upper limit of valid arguments
c
c                    not greater than (0.1e0 * errtol / machine
c                    minimum) ** (2/3), where errtol is described below.
c                    in the following table it is assumed that errtol
c                    will never be chosen smaller than 1.0e-5.
c
c
c                    acceptable values for:   lolim      uplim
c                    ibm 360/370 series   :   6.0e-51     1.0e+48
c                    cdc 6000/7000 series :   5.0e-215    2.0e+191
c                    univac 1100 series   :   1.0e-25     2.0e+21
c                    cray                 :   3.0e-1644   1.69e+1640
c                    vax 11 series        :   1.0e-25     4.5e+21
c
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the routine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c          errtol    relative error due to truncation is less than
c                    3 * errtol ** 6 / (1-errtol) ** 3/2.
c
c
c
c              the accuracy of the computed approximation to the inte-
c              gral can be controlled by choosing the value of errtol.
c              truncation of a taylor series after terms of fifth ordler
c              introduces an error less than the amount shown in the
c              second column of the following table for each value of
c              errtol in the first column.  in addition to the trunca-
c              tion error there will be round-off error, but in prac-
c              tice the total error from both sources is usually less
c              than the amount given in the table.
c
c
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    4.0e-18
c                           3.0e-3    3.0e-15
c                           1.0e-2    4.0e-12
c                           3.0e-2    3.0e-9
c                           1.0e-1    4.0e-6
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c***long description
c
c   rdl special comments
c
c
c
c          check: rdl(x,y,z) + rdl(y,z,x) + rdl(z,x,y)
c          = 3 /  sqrt(x * y * z), where x, y, and z are positive.
c
c
c          on input:
c
c          x, y, and z are the variables in the integral rdl(x,y,z).
c
c
c          on output:
c
c
c          x, y, and z are unaltered.
c
c
c
c          ********************************************************
c
c           warning: changes in the program may improve speed at the
c          expense of robustness.
c
c
c
c    -------------------------------------------------------------------
c
c
c   special functions via rdl and rfl
c
c
c                  legendre form of elliptic integral of 2nd kind
c                  ----------------------------------------------
c
c
c                                             2         2   2
c                  e(phi,k) = sin(phi) rfl(cos (phi),1-k sin (phi),1) -
c
c                       2   3             2         2   2
c                  -1/3k sin (phi) rdl(cos (phi),1-k sin (phi),1)
c
c
c                                  2          2         2
c                  e(k) = rfl(0,1-k ,1) - 1/3k rdl(0,1-k ,^1)
c
c
c
c                  bulirsch form of elliptic integral of 2nd kind
c                  ----------------------------------------------
c
c                                               2 2    2
c                  el2(x,kc,a,b) = ax rfl(1,1+kc x ,1+x ) +
c
c                                                       2 2    2
c                                 +1/3(b-a) x rdl(1,1+kc x ,1+x )
c
c
c
c                  legendre form of alternative elliptic integral of 2nd
c                  -----------------------------------------------------
c                        kind
c                        ----
c
c                            q     2       2   2  -1/2
c                  d(q,k) = int sin p  (1-k sin p)     dp
c                            0
c
c
c
c                                   3           2     2   2
c                  d(q,k) = 1/3 (sin q)  rdl(cos q,1-k sin q,1)
c
c
c
c
c
c                  lemniscate constant b
c                  ---------------------
c
c
c
c                       1    2    4 -1/2
c                  b = int  s (1-s )    ds
c                       0
c
c
c                  b = 1/3 rdl (0,2,1)
c
c
c
c
c                  heuman's lambda function
c                  ------------------------
c
c
c
c                  pi/2 lambda0(a,b) =
c
c                                        2              2
c                     = sin(b) (rfl(0,cos (a),1)-1/3 sin (a) *
c
c                                2               2         2       2
c                      *rdl(0,cos (a),1)) rfl(cos (b),1-cos (a) sin (b),1)
c
c                             2       3             2
c                     -1/3 cos (a) sin (b) rfl(0,cos (a),1) *
c
c                              2         2       2
c                      *rdl(cos (b),1-cos (a) sin (b),1)
c
c
c
c                  jacobi zeta function
c                  --------------------
c
c
c                             2                 2       2   2
c                  z(b,k) = (k/3) sin(b) rfl(cos (b),1-k sin (b),1)
c
c
c                                       2             2
c                             *rdl(0,1-k ,1)/rfl(0,1-k ,1)
c
c                               2       3           2       2   2
c                            -(k /3) sin (b) rdl(cos (b),1-k sin (b),1)
c
c
c    -------------------------------------------------------------------
c          subroutines or functions needed
c              - errpex
c              - r1pext
c              - fortran abs,amax1,amin1,sqrt
c***references  carlson, b.c. and notis,e .m.
c                 algorithms for incomplete elliptic integrals
c                 acm transactions on mathematical software,vol.7,no.3,
c                 sept, 1981, pages 398-403
c               carlson, b.c.
c                 computing elliptic integrals by duplication
c                 numer. math. 33, (1979), 1-16
c               carlson, b.c.
c                 elliptic integrals of the first kind
c                 siam j. math. anal. 8 (1977), 231-242
c***routines called  r1pext,errpex
c***end prologue  rdl
      integer ier,itodo
      real lolim, uplim, epslon, errtol
      real c1, c2, c3, c4, ea, eb, ec, ed, ef, lamda
      real mu, power4, sigma, s1, s2, x, xn, xndev
      real xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
c
c
      character msg*80
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     the original routine is designed for cft
c
c
c     for  civic  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(10)
c
c     for  chat  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(8)
c
c
c      other changes are listed subsequently.
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
      data itodo/1/
c
c
c
c***first executable statement  rdl
      if(itodo.eq.1)then
c
c
      errtol=(r1pext(3)/3.0e0)**(1.0e0/6.0e0)
c
c
      lolim = 2.0e0/(r1pext(2))**(2.0e0/3.0e0)
c
      uplim=r1pext(1)**(1.0e0/3.0e0)
      uplim=(0.10e0*errtol)**(1.0e0/3.0e0)/uplim
      uplim=uplim**2.0e0
c
c
      c1 = 3.0e0/14.0e0
      c2 = 1.0e0/6.0e0
      c3 = 9.0e0/22.0e0
      c4 = 3.0e0/26.0e0
c
c
      itodo=0
c
      end if
c
c
c
c         call error handler if necessary.
c
c
c
    5 rdl=0.0e0
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     for civic users insert the following three statements:
c
c
c          do 9191 igogo=1,10
c
c          msg(igogo)='        '
c
c9191 continue
c
c
c     for chat users insert the following three statements:
c
c
c          do 9191 igogo=1,8
c
c          msg(igogo)='          '
c
c9191 continue
c
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
      if( amin1(x,y).lt.0.0e0) then
      ier=1
      msg=  'rdl - error: amin1(x,y).lt.0.0 (where x=r1, y=r2)'
      call errpex (msg,ier,x,y,0.0,0.0,0.0)
      return
      end if
      if (amin1(x+y,z).lt.lolim) then
      ier=2
      msg=  'rdl - error: amin1(x+y,z).lt.lolim (where x=r1, y=r2)'
      call errpex (msg,ier,x,y,0.0,0.0,0.0)
      msg=  'rdl - error: amin1(x+y,z).lt.lolim (where z=r1, lolim=r2)'
      call errpex (msg,ier,z,lolim,0.0,0.0,0.0)
      return
      end if
      if (amax1(x,y,z).gt.uplim) then
      ier=3
      msg=  'rdl - error: amax1(x,y,z).gt.uplim (where x=r1, y=r2)'
      call errpex (msg,ier,x,y,0.0,0.0,0.0)
      msg=  'rdl - error: amax1(x,y,z).gt.uplim (where z=r1, uplim=r2)'
      call errpex (msg,ier,z,lolim,0.0,0.0,0.0)
      return
      end if
c
   20 ier = 0
      xn = x
      yn = y
      zn = z
      sigma = 0.0e0
      power4 = 1.0e0
c
c
c
   30 mu = (xn+yn+3.0e0*zn)*0.20e0
      xndev = (mu-xn)/mu
      yndev = (mu-yn)/mu
      zndev = (mu-zn)/mu
      epslon = amax1( abs(xndev), abs(yndev), abs(zndev))
      if (epslon.lt.errtol) go to 40
      xnroot =  sqrt(xn)
      ynroot =  sqrt(yn)
      znroot =  sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      sigma = sigma + power4/(znroot*(zn+lamda))
      power4 = power4*0.250e0
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      zn = (zn+lamda)*0.250e0
      go to 30
c
c
c
   40 ea = xndev*yndev
      eb = zndev*zndev
      ec = ea - eb
      ed = ea - 6.0e0*eb
      ef = ed + ec + ec
      s1 = ed*(-c1+0.250e0*c3*ed-1.50e0*c4*zndev*ef)
      s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea))
      rdl = 3.0e0*sigma + power4*(1.0e0+s1+s2)/(mu* sqrt(mu))
c
   50 return
c
c
c
      end

      real function rjl(x,y,z,p,ier)       
c***begin prologue  rjl
c***date written   790801   (yymmdd)
c***revision date  831229   (yymmdd)
c***category no.  c14
c***keywords  duplication theorem,elliptic integral,incomplete,complete,
c             integral of the third kind,taylor series
c***author  carlson, b.c., ames laboratory-doe
c             iowa state university, ames, iowa  50011
c           notis, e.m., ames laboratory-doe
c             iowa state university, ames, iowa  50011
c           pexton, r.l., lawrence livermore national laboratory
c             livermore, california  94550
c***purpose  incomplete or complete elliptic integral of the 3rd kind.
c            for x, y, and z nonnegative, at most one of them zero, and
c            p positive, rjl(x,y,z,p) = integral from zero to infinity of
c                                  -1/2     -1/2     -1/2     -1
c                        (3/2)(t+x)    (t+y)    (t+z)    (t+p)  dt,
c            if x or y or z is 0, the integral is complete.
c***description
c
c   1.     rjl
c          standard fortran function routine
c          single precision version
c          the routine calculates an approximation result to
c          rjl(x,y,z,p) = integral from zero to infinity of
c
c                                -1/2     -1/2     -1/2     -1
c                      (3/2)(t+x)    (t+y)    (t+z)    (t+p)  dt,
c
c          where x, y, and z are nonnegative, at most one of them is
c          zero, and p is positive.  if x or y or z is zero, the
c          integral is complete.  the duplication theorem is iterated
c          until the variables are nearly equal, and the function is
c          then expanded in taylor series to fifth order.
c
c
c   2.     calling sequence
c          rjl( x, y, z, p, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - single precision, nonnegative variable
c
c          y      - single precision, nonnegative variable
c
c          z      - single precision, nonnegative variable
c
c          p      - single precision, positive variable
c
c
c          on  return     (values assigned by the rjl routine)
c
c          rjl     - single precision approximation to the integral
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine.  it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c
c          x, y, z, p are unaltered.
c
c
c   3.    error messages
c
c         value of ier assigned by the rjl routine
c
c                  value assigned        error message printed
c                  ier = 1               amin1(x,y,z) .lt. 0.0e0
c                      = 2               amin1(x+y,x+z,y+z,p) .lt. lolim
c                      = 3               amax1(x,y,z,p) .gt. uplim
c
c
c
c   4.     control parameters
c
c                  values of lolim,uplim,and errtol are set by the
c                  routine.
c
c
c          lolim and uplim determine the valid range of x y, z, and p
c
c          lolim is not less than the cube root of the value
c          of lolim used in the routine for rcl.
c
c          uplim is not greater than 0.3 times the cube root of
c          the value of uplim used in the routine for rcl.
c
c
c                     acceptable values for:   lolim      uplim
c                     ibm 360/370 series   :   2.0e-26     3.0e+24
c                     cdc 6000/7000 series :   5.0e-98     3.0e+106
c                     univac 1100 series   :   5.0e-13     6.0e+11
c                     cray                 :   1.32e-822   1.4e+821
c                     vax 11 series        :   2.5e-13     9.0e+11
c
c
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the routine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c
c
c
c          relative error due to truncation of the series for rjl
c          is less than 3 * errtol ** 6 / (1 - errtol) ** 3/2.
c
c
c
c              the accuracy of the computed approximation to the inte-
c              gral can be controlled by choosing the value of errtol.
c              truncation of a taylor series after terms of fifth order
c              introduces an error less than the amount shown in the
c              second column of the following table for each value of
c              errtol in the first column.  in addition to the trunca-
c              tion error there will be round-off error, but in prac-
c              tice the total error from both sources is usually less
c              than the amount given in the table.
c
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    4.0e-18
c                           3.0e-3    3.0e-15
c                           1.0e-2    4.0e-12
c                           3.0e-2    3.0e-9
c                           1.0e-1    4.0e-6
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c***long description
c
c   rjl special comments
c
c
c          check by addition theorem: rjl(x,x+z,x+w,x+p)
c          + rjl(y,y+z,y+w,y+p) + (a-b) * rjl(a,b,b,a) + 3 / sqrt(a)
c          = rjl(0,z,w,p), where x,y,z,w,p are positive and x * y
c          = z * w,  a = p * p * (x+y+z+w),  b = p * (p+x) * (p+y),
c          and b - a = p * (p-z) * (p-w).  the sum of the third and
c          fourth terms on the left side is 3 * rcl(a,b).
c
c
c          on input:
c
c          x, y, z, and p are the variables in the integral rjl(x,y,z,p).
c
c
c          on output:
c
c
c          x, y, z, and p are unaltered.
c
c          ********************************************************
c
c          warning: changes in the program may improve speed at the
c          expense of robustness.
c
c ------------------------------------------------------------
c
c
c   special functions via rjl,rdl, and rfl
c
c
c                  legendre form of elliptic integral of 3rd kind
c                  ----------------------------------------------
c
c
c                               phi         2         -1
c                  p(phi,k,n) = int (1+n sin (theta) )   *
c                                0
c
c                                      2    2         -1/2
c                                 *(1-k  sin (theta) )     d theta
c
c
c                                               2          2   2
c                            = sin (phi) rfl(cos (phi), 1-k sin (phi),1)
c
c                                     3             2         2   2
c                             -n/3 sin (phi) rjl(cos (phi),1-k sin (phi),
c
c                                      2
c                             1,1+n sin (phi))
c
c
c
c                  bulirsch form of elliptic integral of 3rd kind
c                  ----------------------------------------------
c
c
c                                            2 2    2
c                  el3(x,kc,p) = x rfl(1,1+kc x ,1+x ) +
c
c                                          3           2 2    2     2
c                               +1/3(1-p) x  rjl(1,1+kc x ,1+x ,1+px )
c
c
c                                            2
c                  cel(kc,p,a,b) = a rfl(0,kc ,1) +
c
c
c                                 +1/3(b-pa) rjl(0,kc ,1,p)
c
c
c
c
c                  heuman's lambda function
c                  ------------------------
c
c
c                                 2                     2      2    1/2
c                  l(a,b,p) = (cos(a)sin(b)cos(b)/(1-cos (a)sin (b))   )
c
c                                            2         2       2
c                            *(sin(p) rfl(cos (p),1-sin (a) sin (p),1)
c
c                                 2       3            2       2
c                            +(sin (a) sin (p)/(3(1-cos (a) sin (b))
c
c                                    2         2       2
c                            *rjl(cos (p),1-sin (a) sin (p),1,1-
c
c                                2       2          2       2
c                            -sin (a) sin (p)/(1-cos (a) sin (b))))
c
c
c
c
c                  pi/2 lambda0(a,b) =l(a,b,pi/2) =
c
c                                        2              2
c                     = sin(b) (rfl(0,cos (a),1)-1/3 sin (a) *
c
c                               2               2         2       2
c                     *rdl(0,cos (a),1)) rfl(cos (b),1-cos (a) sin (b),1)
c
c                             2       3             2
c                     -1/3 cos (a) sin (b) rfl(0,cos (a),1) *
c
c                             2         2       2
c                     *rdl(cos (b),1-cos (a) sin (b),1)
c
c
c
c                  jacobi zeta function
c                  --------------------
c
c
c                             2                     2   2    1/2
c                  z(b,k) = (k/3) sin(b) cos(b) (1-k sin (b))
c
c
c                                       2      2   2                 2
c                             *rjl(0,1-k ,1,1-k sin (b)) / rfl (0,1-k ,1)
c
c
c    -------------------------------------------------------------------
c
c          subroutine or functions needed
c              - rcl
c              - errpex
c              - r1pext
c              - fortran abs,amax1,amin1,sqrt
c***references  carlson, b.c. and notis, e.m.
c                 algorithms for incomplete elliptic integrals
c                 acm transactions on mathematical software,vol.7,no.3,
c                 sept, 1981, pages 398-403
c               carlson, b.c.
c                 computing elliptic integrals by duplication
c                 numer. math. 33, (1979), 1-16
c               carlson, b.c.
c                 elliptic integrals of the first kind
c                 siam j. math. anal. 8 (1977), 231-242
c***routines called  r1pext,rcl,errpex
c***end prologue  rjl
      integer ier,itodo
      real alfa, beta, c1, c2, c3, c4, ea, eb, ec, e2, e3
      real lolim, uplim, epslon, errtol, etolrcl
      real lamda, mu, p, pn, pndev
      real power4, rcl, sigma, s1, s2, s3, x, xn, xndev
      real xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
c
c
      character msg*80
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     the original routine is designed for cft
c
c
c     for  civic  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(10)
c
c
c     for  chat  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(8)
c
c
c      other changes are listed subsequently.
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c
c
      data itodo/1/
c
c
c
c***first executable statement  rjl
      if(itodo.eq.1)then
c
c
      errtol=(r1pext(3)/3.0e0)**(1.0e0/6.0e0)
c
c
      lolim =( 5.0e0 * r1pext(1))**(1.0e0/3.0e0)
c
      uplim = 0.30e0*( r1pext(2) / 5.0e0)**(1.0e0/3.0e0)
c
c
c
      c1 = 3.0e0/14.0e0
      c2 = 1.0e0/3.0e0
      c3 = 3.0e0/22.0e0
      c4 = 3.0e0/26.0e0
c
c
c
      itodo=0
c
      end if
c
c
c
c         call error handler if necessary.
c
c
    5 rjl=0.0e0
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     for civic users insert the following three statements:
c
c
c          do 9191 igogo=1,10
c
c          msg(igogo)='        '
c
c9191 continue
c
c
c     for chat users insert the following three statements:
c
c
c          do 9191 igogo=1,8
c
c          msg(igogo)='          '
c
c9191 continue
c
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
      if( amin1(x,y,z).lt.0.0e0) then
      ier=1
      msg='rjl-error: amin1(x,y,z).lt.0.0e0 (where x=r1, y=r2, z=r3)'
      call errpex (msg,ier,x,y,z,p,uplim)
      return
      end if
      if (amin1(x+y,x+z,y+z,p).lt.lolim) then
      ier=2
      msg= 'rjl-error: amin1(x+y,x+z,y+z,p).lt.lolim (where x=r1, y=r2,
     *z=r3, p=r4, lolim=r5)'
      call errpex (msg,ier,x,y,z,p,lolim)
      return
      end if
      if (amax1(x,y,z,p).gt.uplim) then
      ier=3
      msg= 'rjl-error: amax1(x,y,z,p).gt.uplim (where x=r1, y=r2,
     *z=r3, p=r4, uplim=r5)'
      call errpex (msg,ier,x,y,z,p,uplim)
      return
      end if
c
c
c
   20 ier = 0
      xn = x
      yn = y
      zn = z
      pn = p
      sigma = 0.0e0
      power4 = 1.0e0
c
c
c
   30 mu = (xn+yn+zn+pn+pn)*0.20e0
      xndev = (mu-xn)/mu
      yndev = (mu-yn)/mu
      zndev = (mu-zn)/mu
      pndev = (mu-pn)/mu
      epslon = amax1( abs(xndev), abs(yndev), abs(zndev), abs(pndev))
      if (epslon.lt.errtol) go to 40
      xnroot =  sqrt(xn)
      ynroot =  sqrt(yn)
      znroot =  sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      alfa = pn*(xnroot+ynroot+znroot) + xnroot*ynroot*znroot
      alfa = alfa*alfa
      beta = pn*(pn+lamda)*(pn+lamda)
      sigma = sigma + power4*rcl(alfa,beta,ier)
      power4 = power4*0.250e0
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      zn = (zn+lamda)*0.250e0
      pn = (pn+lamda)*0.250e0
      go to 30
c
c
c
   40 ea = xndev*(yndev+zndev) + yndev*zndev
      eb = xndev*yndev*zndev
      ec = pndev*pndev
      e2 = ea - 3.0e0*ec
      e3 = eb + 2.0e0*pndev*(ea-ec)
      s1 = 1.0e0 + e2*(-c1+0.750e0*c3*e2-1.50e0*c4*e3)
      s2 = eb*(0.50e0*c2+pndev*(-c3-c3+pndev*c4))
      s3 = pndev*ea*(c2-pndev*c3) - c2*pndev*ec
      rjl = 3.0e0*sigma + power4*(s1+s2+s3)/(mu* sqrt(mu))
c
c
   50 return
c
c
c
      end

      real function rcl(x,y,ier)    
      integer ier,itodo
      real c1, c2, errtol, lamda, lolim
      real mu, s, sn, uplim, x, xn, y, yn
c
c***begin prologue  rcl
c***date written   790801   (yymmdd)
c***revision date  830622   (yymmdd)
c***category no.  c14
c***keywords  duplication theorem,elementary functions,
c             elliptic integrals,taylor series
c
c
c***author  carlson, b.c., ames laboratory-doe
c             iowa state university, ames, iowa  50011
c           notis, e.m., ames laboratory-doe,
c             iowa state university, ames, iowa  50011
c           pexton, r.l., lawrence livermore national laboratory
c             livermore, california  94550
c***purpose  the routine calculates an approximation result to
c            rcl(x,y) = integral from zero to infinity of
c                              -1/2     -1
c                    (1/2)(t+x)    (t+y)  dt,
c            where x is nonnegative and y is positive.
c***description
c
c   1.     rcl
c          standard fortran function routine
c          single precision version
c          the routine calculates an approximation result to
c          rcl(x,y) = integral from zero to infinity of
c
c                              -1/2     -1
c                    (1/2)(t+x)    (t+y)  dt,
c
c          where x is nonnegative and y is positive.  the duplication
c          theorem is iterated until the variables are nearly equal,
c          and the function is then expanded in taylor series to fifth
c          order.  logarithmic, inverse circular, and inverse hyper-
c          bolic functions can be expressed in terms of rcl.
c
c
c   2.     calling sequence
c          rcl( x, y, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - single precision, nonnegative variable
c
c          y      - single precision, positive variable
c
c
c
c          on return  (values assigned by the rcl routine)
c
c          rcl     - single precision approximation to the integral
c
c          ier    - integer to indicate normal or abnormal termination.
c
c                     ier = 0 normal and reliable termination of the
c                             routine.  it is assumed that the requested
c                             accuracy has been achieved.
c
c                     ier > 0 abnormal termination of the routine
c
c          x and y are unaltered.
c
c
c   3.    error messages
c
c         value of ier assigned by the rcl routine
c
c                  value assigned         error message printed
c                  ier = 1                x.lt.0.0e0.or.y.le.0.0e0
c                      = 2                x+y.lt.lolim
c                      = 3                amax1(x,y) .gt. uplim
c
c
c   4.     control parameters
c
c                  values of lolim,uplim,and errtol are set by the
c                  routine.
c
c          lolim and uplim determine the valid range of x and y
c
c          lolim  - lower limit of valid arguments
c
c                   not less  than 5 * (machine minimum)  .
c
c          uplim  - upper limit of valid arguments
c
c                   not greater than (machine maximum) / 5 .
c
c
c                     acceptable values for:   lolim       uplim
c                     ibm 360/370 series   :   3.0e-78     1.0e+75
c                     cdc 6000/7000 series :   1.0e-292    1.0e+321
c                     univac 1100 series   :   1.0e-37     1.0e+37
c                     cray                 :   2.3e-2466   1.09e+2465
c                     vax 11 series        :   1.5e-38     3.0e+37
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the routine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c
c          errtol  - relative error due to truncation is less than
c                    16 * errtol ** 6 / (1 - 2 * errtol).
c
c
c              the accuracy of the computed approximation to the inte-
c              gral can be controlled by choosing the value of errtol.
c              truncation of a taylor series after terms of fifth order
c              introduces an error less than the amount shown in the
c              second column of the following table for each value of
c              errtol in the first column.  in addition to the trunca-
c              tion error there will be round-off error, but in prac-
c              tice the total error from both sources is usually less
c              than the amount given in the table.
c
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    2.0e-17
c                           3.0e-3    2.0e-14
c                           1.0e-2    2.0e-11
c                           3.0e-2    2.0e-8
c                           1.0e-1    2.0e-5
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c***long description
c
c   rcl special comments
c
c
c
c
c                  check: rcl(x,x+z) + rcl(y,y+z) = rcl(0,z)
c
c                  where x, y, and z are positive and x * y = z * z
c
c
c          on input:
c
c          x and y are the variables in the integral rcl(x,y).
c
c          on output:
c
c          x and y are unaltered.
c
c
c
c                    rcl(0,1/4)=rcl(1/16,1/8)=pi=3.14159...
c
c                    rcl(9/4,2)=ln(2)
c
c
c
c          ********************************************************
c
c          warning: changes in the program may improve speed at the
c          expense of robustness.
c
c
c   --------------------------------------------------------------------
c
c   special functions via rcl
c
c
c
c                  ln x                x .gt. 0
c
c                                             2
c                  ln(x) = (x-1) rcl(((1+x)/2)  , x )
c
c
c   --------------------------------------------------------------------
c
c                   arcsin x            -1 .le. x .le. 1
c
c                                        2
c                   arcsin x = x rcl (1-x  ,1 )
c
c   --------------------------------------------------------------------
c
c                   arccos x            0 .le. x .le. 1
c
c
c                                      2        2
c                   arccos x = sqrt(1-x ) rcl(x  ,1 )
c
c   --------------------------------------------------------------------
c
c                   arctan x            -inf .lt. x .lt. +inf
c
c                                         2
c                   arctan x = x rcl(1,1+x  )
c
c   --------------------------------------------------------------------
c
c                   arccot x            0 .le. x .lt. inf
c
c                                   2   2
c                   arccot x = rcl(x  ,x +1 )
c
c   --------------------------------------------------------------------
c
c                   arcsinh x           -inf .lt. x .lt. +inf
c
c                                        2
c                   arcsinh x = x rcl(1+x  ,1 )
c
c   --------------------------------------------------------------------
c
c                   arccosh x           x .ge. 1
c
c                                     2         2
c                   arccosh x = sqrt(x -1) rcl(x  ,1 )
c
c   --------------------------------------------------------------------
c
c                   arctanh x           -1 .lt. x .lt. 1
c
c                                          2
c                   arctanh x = x rcl(1,1-x  )
c
c   --------------------------------------------------------------------
c
c                   arccoth x           x .gt. 1
c
c                                    2   2
c                   arccoth x = rcl(x  ,x -1 )
c
c   --------------------------------------------------------------------
c
c          subroutines or functions needed
c              - errpex
c              - r1pext
c              - fortran abs,amax1,sqrt
c***references  carlson, b.c. and notis, e.m.
c                 algorithms for incomplete elliptic integrals
c                 acm transactions on mathematical software,vol.7,no.3,
c                 sept, 1981, pages 398-403
c               carlson, b.c.
c                 computing elliptic integrals by duplication
c                 numer. math. 33, (1979), 1-16
c               carlson, b.c.
c                 elliptic integrals of the first kind
c                 siam j. math. anal. 8 (1977), 231-242
c***routines called  r1pext,errpex
c
      character msg*80
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     the original routine is designed for cft
c
c
c     for  civic  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(10)
c
c
c     for  chat  users replace the above character statement
c
c
c        '  character msg*80  '   with the following :
c
c
c           dimension msg(8)
c
c
c      other changes are listed subsequently.
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
      data itodo/1/
c***first executable statement  rcl
      if(itodo.eq.1)then
c
c
      errtol=(r1pext(3)/16.0)**(1.0/6.0)
c
c
      lolim = 5.0e0 * r1pext(1)
      uplim = r1pext(2) / 5.0e0
c
c
      c1 = 1.0e0/7.0e0
      c2 = 9.0e0/22.0e0
c
c
c
      itodo=0
c
      end if
c
c
c         call error handler if necessary.
c
c
    5 rcl=0.0e0
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
c     for civic users insert the following three statements:
c
c
c          do 9191 igogo=1,10
c
c          msg(igogo)='        '
c
c9191 continue
c
c
c     for chat users insert the following three statements:
c
c
c          do 9191 igogo=1,8
c
c          msg(igogo)='          '
c
c9191 continue
c
c
c
c
c###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
c
c
      if (x.lt.0.0e0.or.y.le.0.0e0) then
      ier=1
      msg='rcl - error: x.lt.0.0e0.or.y.le.0.0e0 (where x=r1, y=r2)'
      call errpex (msg,ier,x,y,0.0,0.0,0.0)
      return
      end if
      if (x+y.lt.lolim)   then
      ier=2
      msg='rcl - error: x+y.lt.lolim (where x=r1, y=r2)'
      call errpex (msg,ier,x,y,0.0,0.0,0.0)
      msg='rcl - error: x+y.lt.lolim (where lolim=r1)'
      call errpex (msg,ier,lolim,0.0,0.0,0.0,0.0)
      return
      end if
      if (amax1(x,y).gt.uplim) then
      ier=3
      msg='rcl - error: amax1(x,y).gt.uplim (where x=r1, y=r2)'
      call errpex (msg,ier,x,y,0.0,0.0,0.0)
      msg='rcl - error: amax1(x,y).gt.uplim (where uplim=r1)'
      call errpex (msg,ier,uplim,0.0,0.0,0.0,0.0)
      return
      end if
c
c
c
   20 ier = 0
      xn = x
      yn = y
c
c
c
   30 mu = (xn+yn+yn)/3.0e0
      sn = (yn+mu)/mu - 2.0e0
      if ( abs(sn).lt.errtol) go to 40
      lamda = 2.0e0* sqrt(xn)* sqrt(yn) + yn
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      go to 30
c
   40 s = sn*sn*(0.30e0+sn*(c1+sn*(0.3750e0+sn*c2)))
      rcl = (1.0e0+s)/ sqrt(mu)
c
c
c
   50 return
c
      end

