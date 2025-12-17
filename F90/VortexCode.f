c       --------------------------------------------------------------
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
c	V2:	       2021 
c
c       V3:            Start: 2020 August 22
c
	program VortexCode
c
	use mem
c
	real nslope
c
        logical inccp
c
	call cpu_time(stime)
	call initialize
c
c***********************************************************
c       slipstream iteration
c***********************************************************
c
	niter    = 0
        errcpnew = 1.
c
	write(*,*)
	write(*,*)'================================================ '
	write(*,*)'start slip stream iteration'
	write(*,*)'================================================ '
	write(*,'(4a10,a15)')'niter','cperr','psierr','fferr','cP'
c
	do while (niter.le.maxiter)
c
	errcpold = errcpnew
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
c-------------------------------------------------------------------------------------------
c    	     
	     case (1)
c
c               (1) From Bontempo Eqs 13 and 14
c                z should also change, but this does not work
c
                dzz   = zsl(mm)+zsl(mm-1)

 	        dzn  = amin1(dzz,1.)*under*(1.+vz(mm))*ds/vav
                drn  = amin1(dzz,1.e-3)*under*vr(mm)*ds/vav 	
c
	     case (2)
c
c	        (2) rotate slip stream pieces - modification of van Kuik's (D.13)
c	        
                dzz   = zsl(mm)+zsl(mm-1)
c
	        beta   = atan2(dr,dz)
	        vrr    = vr(mm)
	        vzz    = 1. + vz(mm)
	        alfa   = atan2(vrr,vzz)
c
c               rotate in direction of v
c
                dphi   = alfa-beta
c
                dzn = dz*cos(dphi)-dr*sin(dphi)
                drn = dz*sin(dphi)+dr*cos(dphi)
c
	     case (3)
c
c	        (3) from van Kuik Eq (D.12) using
c               u_z = 1/r dpsi/dr -> dr = dpsi/(r u_z)
c               we use dpsi = psi-wake - psi-sls 
c
                dzz   = .5*(zsl(mm)+zsl(mm-1))
	        zw    = dzz
	        rw    = rsl(mm)
	        vzw   = 1. + 0.5*(vz(mm)+vz(mm-1)) 
	        vrw   =      0.5*(vr(mm)+vr(mm-1)) 
c
	        psir  = psiBr(zw,rw) + psiBrCyl(zw,rw)
	        dpsi  = psiwake - psir
c
c                drn   = amin1(dzz,.1)*under*dpsi/(vzw*rw)
                drn   = amin1(dzz,.1)*under*dpsi
c
c               2022 08 22: underrelaxation depends on z-ring
c
c               debug
c
c               write(*,'(5a12)')'zw','psiwake','psir','rw','drn'
	        if(zw.gt.1..and.zw.lt.5.)then
                   write(222,'(5e12.4)')zw,psiwake,psir,rw,drn
                endif
c
c---------------------------------------------------------------------------------------
c              add change of dz as well
c
                dzn   = dpsi/(vrw*rw)
	        dzn   = amin1(dzz,1.)*under*dzn
                if (dzn.lt.0.) then
                   dzn = amax1(dzn,-1.e-2)
                else
                   dzn = amin1(dzn, 1.e-2)
                endif
	        dzn = 0.  
c
	    case default
c
	    end select
c	    end select case(modelit)
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
	  errff= 0.
	  do i=1,np
	     zw = .5*(zsl(i)+zsl(i-1))
	     rw = .5*(rsl(i)+rsl(i-1))
	     dr =  rsl(i)-rsl(i-1)
	     zr =  zsl(i)-zsl(i-1)
	     ds = sqrt(dr*dr+dz*dz)
             dpsi = abs(psiBr(zw,rw)+ psiBrCyl(zw,rw)-psiwake)
	     dpsi = dpsi*ds
             errsl = errsl  + dpsi
c
	     vm  = sqrt(vr(i)**2+((1.0+vz(i))**2))
	     dff = (vm*ga(i) + 0.5*ct)*ds
	     errff = errff + dff
	  enddo
c
	  errsl = errsl/zsivc
	  errff = errff/zsivc
c
	  cpw      = cpp(cp)
	  errcp    = abs((cpw-cpm)/cpm)
c
c         new: 2022 08 22: prevent iteration if errcp rizes
c
          errcpold = errcpnew
          errcpnew = errcp
c***********************************************************************
c
c         quit iteration if accuracy reached or cp rises again
c
c
c          write(*,'(2e10.3)')errcp,epsiter
          inccp = .true.
 	  if(errcp.lt.epsiter.or.
     +       errsl.lt.epsiter.or.
     +       errcpnew.gt.errcpold)
     +    GOTO 321
C***************************************************************************
c         output
c
	  write(*,'(i10,3e10.2,f15.7)')niter,errcp,errsl,errff,cpw
c
c	  write out cp for convergence check

	 write(10,'(i5,f12.8,e10.3)')niter,cpw,errcp
c
c       save data from last iteration
c       2022 08 22
c
	   do i = 1,np
	      gao(i)  = ga(i)
	      rslo(i) = rsl(i)
	      zslo(i) = zsl(i)
              vzo(i)  = vz(i)
              vro(i)  = vr(i)
	   end do
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
c-----------------------------------------------------------
        end if
c       end if update
c
	niter = niter + 1
c
        end do
c
c       end do while niter.lt.maxiter
c
321	CONTINUE
c
c       restore data from last iteration
c       2022 08 22
c
        if(inccp)then
           write(*,*)'load last result'
	   do i = 1,np
	      ga (i) = gao( i)
	      rsl(i) = rslo(i)
	      zsl(i) = zslo(i)
              vz (i) = vzo (i)
              vr (i) = vro (i)
	   end do
        endif
c
	write(*,*)'========================='
	write(*,*)'end slip stream iteration'
	write(*,*)'========================='
c
	call output
c
c---------------------------------------------------------------------
c
	call cpu_time(etime)
	write(*,*)
        wsec = etime-stime
        imin = int(wsec/60.)
        ihrs = int(imin/60 )
        imin = imin - 60 *ihrs
        wsec = wsec - 3600.*ihrs - 60.*imin
        write(*,*)'wsec',wsec
	write(*,109)'elapsed time:'
        write(*,111)ihrs ,'hrs',imin,'min',int(wsec),'sec'
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
	close(14)
	close(15)
c
c       format specifiers
c
101	format(a15,4f10.4)
103     format(a20,5F9.3)
104	format(a20,2f12.6)
106	format(a)
109	format(a15)
110	format(2f12.6)
111     format(3(i3,a5))
c
	end

