c
c	---------------------------------
c
c   	aps start:    2020 11 13
c	
c       try more F90:  2021 02 16
c
	program VortexCode
c
	use mem
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
	OPEN(UNIT= 1,FORM='FORMATTED',STATUS='UNKNOWN',FILE='vzd.DAT')
	OPEN(UNIT= 2,FORM='FORMATTED',STATUS='UNKNOWN',FILE='psid.DAT')
	OPEN(UNIT= 3,FORM='FORMATTED',STATUS='UNKNOWN',FILE='sls.DAT')
	OPEN(UNIT= 4,FORM='FORMATTED',STATUS='UNKNOWN',FILE='vza.DAT')
	OPEN(UNIT= 5,FORM='FORMATTED',STATUS='UNKNOWN',FILE='inpa.dat')
c	OPEN(UNIT=10,FORM='FORMATTED',STATUS='UNKNOWN',FILE='slsit.DAT')
	OPEN(UNIT=11,FORM='FORMATTED',STATUS='UNKNOWN',FILE='vrd.DAT')
c
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
c       read configuration file inpa.dat
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
c	start of cylinder
c
	read(5,*) inpstr,zsivc
	write(*,701)inpstr,zsivc
c
c       slip stream iteration parameters
c
c
c       logical update:
c
	read(5,*) inpstr,update
	write(*,*)inpstr,update
c
c       maximum number of iterations
c
	read(5,*) inpstr,maxiter
	write(*,702)inpstr,maxiter
702	format(a10,i6)
c
c       some under-relaxation for sls 
c
	read(5,*) inpstr,under
	write(*,703)inpstr,under
703	format(a10,e12.4)
c
c       epsilon to finish iteration
c
	read(5,*) inpstr,epsiter
	write(*,703)inpstr,epsiter
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
	write(*,702)inpstr,nwt
c
c      slipstream saturation length (around 1)
c
	read(5,*) inpstr,al
	write(*,701)inpstr,al
c
c       model for sls iteration see line 240 ff
c
	read(5,*) inpstr,modelit
	write(*,702)inpstr,modelit
c
c       some numbers concerning case
c
	vzinf = sqrt(1. + ct)
	vz0   = 0.5*(1. + vzinf)
	rinf  = sqrt(vz0/vzinf)
c
	a     = 0.5*(sqrt(1.+ct)-1.)
	cpm   = 4.*a*(1+a)**2
c
	psiwake = 0.25*(1.+sqrt(1.+ct))
c
	rsivc = rinf
c
c      output some values from inpa.dat
c
	write(*,*)'======================================='
	write(*,'(a10,i6,a10,f10.2)')'NP = ',NP,' zsivc ',zsivc
	write(*,'(a20,4f12.4)')'ct vz0 vzinf rinf',ct, vz0, vzinf, rinf
	write(*,*)
	write(*,'(a22,2f12.4)')'Momentum Theory: cP= ',cpm
	write(*,*)
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
c       initialize vortex strengths 
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
c       simple cosine stretch
c
	do i=1,np
	   phi = i*dphi
	   zsl(i) = zsivc*(1.-cos(phi))
c	   write(*,*)cos(phi),zsl(i)
	end do
c
c       initialize slipstream shape
c       nwt = 0 -> rsl = 1 (like a cylinder)
c      	nwr = 1 -> simple exponetial blending
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
c***********************************************************
c       slipstream iteration
c***********************************************************
c
	niter   = 0
        erri    = 1.
	errcp   = 1.
c
	do while 
     +     (niter.lt.maxiter)
c
c	   write(*,'(a20,i3,e10.3)')'niter errcp',niter, errcp
c	   write(*,*)
c
	   do i=1,np
	      vz(i) = 0.
	      vr(i) = 0.
	   end do
c
c          main loop for INDUCED velocities
c
c	   loop over all vortex strips
c
	   do n = 1,np
c
c             induction from all (other) m's
c
	      do m = 1,n
	         if (m.eq.n)then
c
c                   self induction Bontempo Eq (6) and (7)
c
                       vz(n) = vz(n) + vzself(n)
                       vr(n) = vr(n) + vrself(n)
                 else  
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
c
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
c          end n loop
c	      
	end do
c
c	update vortex strengths according to Bontempo's Eq. (11)
c       
	do i = 1,np
	   vm    = sqrt(vr(i)**2+vz(i)**2)
	   vm = 1. + vm
	   ga(i) = -ct/(2.*vm)
	end do
c
c	update wake shape 
c
	zsln(0) = 0.
	rsln(0) = 1.
c
	if(update)then
	   do mm = 1,np
	     dz   = zsl(mm)-zsl(mm-1)
	     dr   = rsl(mm)-rsl(mm-1)
c
	     ds   = sqrt(dz*dz+dr*dr)
c
             vz2  = (vz(mm))**2 
             vzz  = sqrt(vz2)

	     vr2  = vr(mm)**2
             vrr  = sqrt(vr2)
	     vav  = sqrt(vz2+vr2)
c
	     select case(modelit)
c
c----------------------------------------------------------------------------
c
c               (1) From Bontempo Eqs 13 and 14
c                z should also change, but this does not work
c	    	     
	     case (1)
 	        dzn  = vzz*ds/vav
                drn  = vrr*ds/vav 	
	     case (2)
c
c	        own approaches: 
c	        (2) incline sniplet towards vel-vector of total flow
c
                alfav = atan2(vr(mm),vz(mm))
                alfas = atan2(dr,dz)
                dalfa = alfav-alfas
                drn = ds*sin(dalfa) 
                dzn = ds*cos(dalfa)
c
	     case (3)
c
c	        (3) shift normal to existing shape
c
                sslope = vr(mm)/(vz(mm))
                nslope = 1./sslope
                cossl  = 1.    /sqrt(1. + nslope**2)
                sinsl  = nslope/sqrt(1. + nslope**2)
                drn    = -ds*sinsl
                dzn    = ds*cossl
c
	     case (4)
c
c	        (4) from van Kuik Eq (D.11) and u_z = 1/r d psi/dr -> dr = d psi/(r u_z)
c               we set d psi = psi-wake - psi-sls = psi-inf - psid
c               gives problems close  to z= 0
c
	        zw    = zsl(mm)
	        rw    = rsl(mm)
	        vzw   = 1.+vz(mm)
c
	        dpsi  = psiBr(zw,rw) - psiwake
c
	        drn   = dpsi/(vzw*rw)
	        dzn   = 0.
c
c	        write(*,*)'z r vzw dpsi drn',zw,rw,vzw,dpsi,drn
c
	    case default
c
	    end select
c
c----------------------------------------------------------------------------------------------------
c            now set new slipstream coordinates         
c   
c            approach Bontempo: 
c	     new slipstream directely from z=0 without ref to older (exp) shape
c    
c 	      zsln(mm) = zsln(mm-1) + under*dzn
c             rsln(mm) = rsln(mm-1) + under*drn
c
c            approach APS 2020 12 18: change with ref to old shape
c
              zsln(mm) = zsl(mm) + under*dzn
              rsln(mm) = rsl(mm) + under*drn
c
c	      if(rsln(mm).lt.rsivc)rsln(mm)=rsivc
c
c         ende do mm
          end do
c
c         check convergence
c
          erri= 0.
	  do i=1,np
	     dz  = zsl(i)-zsln(i)
	     dr  = rsl(i)-rsln(i)
             erri = erri + dz*dz + dr*dr
	  enddo
c
	  cpw   = cp(ct,model)
	  errcp = abs((cpw-cpm)/cpm)
	  if(errcp.lt.epsiter)niter=maxiter
c
c         output
c
	  write(*,'(a20,i3,2f10.5)')'niter cp rsl ',
     +             niter,cpw,rsl(410)
c
c	  write out new and old slipstream
c
c	  do i = 1,np
c	    write(10,'(4f12.6)')zsl(i),rsl(i),zsln(i),rsln(i)
c	  end do
c	  write(10,*)
c
c        now update slip stream
c 
	  do i = 1,np
            zsl(i) =  zsln(i)
            rsl(i) =  rsln(i)
	  end do
c
	  rsivc = rsl(np)
	  zsivc = zsl(np)
c
	  niter = niter + 1
c
c       end if update
c
        end if
c
        end do
c
c       end do while niter.lt.maxiter
c
	write(*,*)'end iteration'
c
c**********************************************************************
c	output
c**********************************************************************

	write(*,'(a10,E13.3)')'errsl ',erri
	write(*,'(a10,E13.3)')'errcp ',errcp
	write(*,'(a10,i5)')'   niter ',niter
c
        write(3,99)'i','z','r','vz','vr','ga','psi-vK','psi-Br'
c
c       rings
c
	do i=1,np
	   zpsi   = 0.5*(zsl(i)+zsl(i-1))
	   rpsi   = 0.5*(rsl(i)+rsl(i-1))
	   write(3,100)i,zsl(i),rsl(i),vz(i),vr(i),ga(i),
     +     psivK(zpsi,rpsi),psiBr(zpsi,rpsi)
	end do
c
c       cylinder
c
	dz   = zsivc-zsl(np)
	dr   = rsivc-rsl(np)
	ds   = sqrt(dz*dz+dr*dr)
c
        select case(model)
           case(1)
              vzs  = vzsivcBo(zsivc,rsivc)
	      vrs  = vrsivcBo(zsivc,rsivc)
	   case(2)
              vzs  = vzsivcK(zsivc,rsivc) 
	      vrs  = vrsivcK(zsivc,rsivc)	     
	   case(3)
	      vzs  = vzsivcBr(zsivc,rsivc) 
	      vrs  = vrsivcBr(zsivc,rsivc)
	   case Default
	end select      
c
c       stream function
c
c	psicyl = psiBrcyl(zsivc,rsivc)
c	write(3,100)np+1,zsivc,rsivc,vzs,vrs,gasivc,psicyl,psicyl
c
	write(*,*)
c
c**********************************************************************
c
c       calculate cP from Bontempo's Eq. (19) and print out vz at disk
c
	cpr = cp(ct,model)
c
	write(*,*)'--------------  Finally  -------------------------'
		abw = (cpm-cpr)/cpm
	write(*,105)'cpmom     cp      dev = ',cpm,cpr,abw
	abw = 1000.*(rsl(410)-rinf)/rinf
	write(*,105)'rinv  rsl(410) dev(pm)= ',rinf,rsl(410),abw
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
	   vzout = 1.+vzd(rd,model)
	   vzoutc= vzdcyl(rd,model)
	   write(1,102)rd,vzout,vzoutc
	end do
c
c	                  
c**********************************************************************
c
c       print out vr at z = 0 (disk)
c
	write(1,*)
	write(1,*)
	nd = 100
	drd = 1./float(nd)
	do i=0,nd-1
	   rd = 0.5*drd+drd*real(i)
	   vrout = vrd(rd,model)
	   vroutc= vrdcyl(rd,model)
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
	   write(2,102)rd,psiBr(0.,rd),psiBrCyl(0.,rd),psivK(0.,rd)
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
	   write(4,102)za,vza(za),vzsivcBr(za,0.),vzacylan(za,r)
	end do
c
c---------------------------------------------------------------------
c
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(11)
c
99	format(8a12) 
100	format(i12,7f12.6)
101	format(a15,4f10.4)
102	format(4f18.8)
103     format(a20,5F9.3)
104	format(a20,2f12.6)
105	format(a25,2f12.6,e10.2)
106	format(a)
c
c-----------------------------------------------------------------------
	end

