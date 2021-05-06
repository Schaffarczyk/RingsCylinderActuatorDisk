c       ---------------------------------
c
c       Alois Peter Schaffarczyk
c	Kiel University of Applied Sciences
c	Kiel, Germany
c
c   	start coding:  2020 11 13
c	
c       try more F90:  2021 02 16
c
c	polish:        2021 05 03
c
	program VortexCode
c
	use mem
c
	call cpu_time(stime)
	call initialize
c
c***********************************************************
c       slipstream iteration
c***********************************************************
c
	niter   = 0
c
	write(*,*)
	write(*,*)'================================================ '
	write(*,*)'start slip stream iteration'
	write(*,*)'================================================ '
	write(*,'(4a10)')'niter', 'psidev', 'cp', 'cpdev'
c
	do while (niter.le.maxiter)
c
c       calculate induced velocities
c
	call indvel
c
c	update vortex strengths according to Bontempo's Eq. (11)
c       
	do i = 1,np
c
	   vm = sqrt(vr(i)**2+((1.0+vz(i))**2))
	   ga(i) = -ct/(2.*vm)
c
	end do
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
c            start slip stream iteation
c
	     select case(modelit)
c
c----------------------------------------------------------------------------
c    	     
	     case (1)
c
c               (1) From Bontempo Eqs 13 and 14
c                z should also change, but this does not work
c
 	        dzn  = vz(mm)*ds/vav
                drn  = vr(mm)*ds/vav 	

	     case (2)
c
c	        (2) shift normal to existing shape
c	        
                sslope = vr(mm)/(vz(mm))
                nslope = -1./sslope
                cossl  = 1.    /sqrt(1. + nslope**2)
                sinsl  = nslope/sqrt(1. + nslope**2)
                drn    = under*ds*sinsl
                dzn    = under*ds*cossl
c
	     case (3)
c
c	        (3) from van Kuik Eq (D.11) using
c               u_z = 1/r dpsi/dr -> dr = dpsi/(r u_z)
c               we set dpsi = psi-wake - psi-sls 
c               gives problems close  to z= 0
c
	        zw    = 0.5*(zsl(mm)+zsl(mm-1))
	        rw    = 0.5*(rsl(mm)+rsl(mm-1))
	        vzw   = 1. + 0.5*(vz(mm)+vz(mm-1)) 
c
	        psir  = psiBr(zw,rw) + psiBrCyl(zw,rw)
	        dpsi  = psiwake - psir
c
                drn   = under*dpsi/(vzw*rw)
	        dzn   = 0.
c
	    case default
c
	    end select
c	    select case(modelit)
c
c           DEBUG
c
c	    if(zw.gt.2.4.and.zw.lt.2.45)then
c	       write(*,888)'zw rw vz dpsi dr',zw,rw,vzw,dpsi,drn
c888	       format(a24,5f9.3)
c	    end if
c
c----------------------------------------------------------------------------------------------------
c            change to new slipstream coordinates         
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
          end do
c         ende do mm
c
c         check convergence
c
          errsl= 0.
	  do i=1,np
	     zw = .5*(zsl(i)+zsl(i-1))
	     rw = .5*(rsl(i)+rsl(i-1))
	     dr =  rsl(i)-rsl(i-1)
	     zr =  zsl(i)-zsl(i-1)
	     ds = sqrt(dr*dr+dz*dz)
             dpsi = abs(psiBr(zw,rw)+ psiBrCyl(zw,rw)-psiwake)
	     dpsi = dpsi*ds
             errsl = errsl  + dpsi
	  enddo
	  errsl = errsl/zsivc
c
	  cpw   = cpp(cp)
	  errcp = abs((cpw-cpm)/cpm)
	  if(errcp.lt.epsiter)niter=maxiter
c
c         output
c
	  write(*,'(i10,e10.2,f10.6,e10.2)')niter,errsl,cpw,errcp
c
c	  write out cp for convergence check

	 write(10,'(i5,f12.8,e10.3)')niter,cpw,errcp
c
c        now update slip stream
c 
	  do i = 1,np
            zsl(i) =  zsln(i)
            rsl(i) =  rsln(i)
	  end do
c
c         DEBUG: wake development during iteration
c	
	  do i=1,np
  	    zp   = 0.5*(zsl(i)+zsl(i-1))
	    rp   = 0.5*(rsl(i)+rsl(i-1))
	    write(13,110),zp,rp
	    write(14,110),zp,psiBr(zp,rp)+psiBrCyl(zp,rp)
          end do
	  write(13,*)
	  write(14,*)
c
        end if
c       end if update
c
	niter = niter + 1
c
c	write(*,*) 'niter 2: ',niter
c
        end do
c
c       end do while niter.lt.maxiter
c
	write(*,*)'end slip stream iteration'
	write(*,'(a10,E13.3)')'err PSI ',errsl
c
	call output
c
c---------------------------------------------------------------------
c
	call cpu_time(etime)
	write(*,*)
	write(*,109)'elapsed time (sec): ', etime-stime,' secs'
	write(*,*)
c
c	close used files
c
	close(1)
	close(2)
	close(3)
	close(4)
	close(10)
	close(11)
	close(12)
	close(13)
c
c       format specifiers
c
101	format(a15,4f10.4)
103     format(a20,5F9.3)
104	format(a20,2f12.6)
106	format(a)
109	format(a20,f12.2,a6)
110	format(2f12.6)
c
	end

