c       --------------------------------------------------------------
c
c       VortexCode-V4:  structural change of the discretization
c                       (Claude, 2026 07 27)
c
c       Following van Kuik & Lignarolo, Wind Energy 19 (2016),
c       Appendix A+B, but with a spacing-proportional vortex kernel:
c
c       - N vortex rings with circulations Gam_i sit AT NODES
c         (x_i, r_i); ring 1 is FIXED at the disc edge (0,1).
c       - vortex kernel: for evaluation distances rho < 0.5*beta*ds_i
c         the smoothed (Marshall-type) expressions are used.
c         beta = exp(-1/2) = 0.607 matches the analytic strip
c         integral (a FIXED delta as in the paper biases cP when the
c         cell size varies -- tested with the python prototype).
c       - semi-infinite vortex tube at x = LtubeFac*Rwake, radius
c         Rwake, strength gamma_wake (from momentum theory).
c       - convergence scheme (single phase; the paper's v_n
c         "fine-tuning" phase was tested and is NOT stable/needed):
c            Gam_i <- Gam_i + dG*(Gam*_i - Gam_i),
c                     Gam*_i = -cT/(2 v_i) * ds_i     (force free)
c            r_i   <- r_i + dA*(psi_wake - psi_i)/(r_i*(1+vz_i))
c         with UNIFORM damping dA = dG = 0.05 -- stable because of
c         the kernel regularization (no saw-tooth mode).
c
c       Validation (python prototype adcode_v4.py, Betz cT = -8/9):
c         N =  400:  cP = -0.5920118  (err +5.8e-4), maxdpsi 3e-13
c         N =  800:  cP = -0.5923623  (err +2.3e-4), maxdpsi 2e-12
c         N = 1600:  cP = -0.5925434  (err +4.9e-5), maxdpsi 8e-12
c       i.e. a UNIQUE fixed point (no cP band, no drift) and cP
c       errors well below the V3 midpoint scheme at comparable N.
c
c       Files: link with mem.f and SubA-V3.f (elliptic integrals and
c       the cylinder functions with the CORRECTED v_r sign).
c
c       --------------------------------------------------------------
c
	module memv4
c
	real, allocatable :: dsc(:), psv(:), gnew(:)
	real beta, dampa, dampg, epspsi, c1, c2, ltfac
	integer ncpout
c
	end module memv4
c
c       --------------------------------------------------------------
c
	program VortexCodeV4
c
	use mem
	use memv4
c
	real vzp, vrp, psp
c
	call cpu_time(stime)
	call initv4
c
	write(*,*)
	write(*,*)'================================================'
	write(*,*)'start V4 iteration (node collocation + kernel)'
	write(*,*)'================================================'
	write(*,'(a10,2a13,a15,a13)')
     +      'niter','maxdpsi','rmsdpsi','cP','cperr'
c
	niter = 0
	dpmax = 1.
c
	do while (niter.le.maxiter)
c
c          fields at all ring positions
c
	   do i = 1,np
	      call ringsv4(zsl(i),rsl(i),vzp,vrp,psp)
	      vz(i)  = vzp + vzsivcBr(zsl(i),rsl(i))
	      vr(i)  = vrp + vrsivcBr(zsl(i),rsl(i))
	      psv(i) = psp + psiBrCyl(zsl(i),rsl(i))
     +                     + 0.5*rsl(i)*rsl(i)
	   end do
c
c          circulation update (force free), uniform damping
c
	   do i = 1,np
	      vm    = sqrt((1.+vz(i))**2 + vr(i)**2)
	      gstar = -ct/(2.*vm)*dsc(i)
	      ga(i) = ga(i) + dampg*(gstar-ga(i))
	   end do
c
c          position update from psi condition; ring 1 pinned
c
	   dpmax = 0.
	   dprms = 0.
	   do i = 2,np
	      dpsi  = psiwake - psv(i)
	      dpmax = amax1(dpmax,abs(dpsi))
	      dprms = dprms + dpsi*dpsi
	      dr    = dampa*dpsi/(rsl(i)*(1.+vz(i)))
	      if(dr.gt. 0.01) dr =  0.01
	      if(dr.lt.-0.01) dr = -0.01
	      rsl(i) = rsl(i) + dr
	   end do
	   dprms = sqrt(dprms/float(np-1))
c
c          refresh centered arclength weights
c
	   call dscupd
c
c          output & stop
c
	   if(mod(niter,ncpout).eq.0)then
	      cpw   = cppv4()
	      errcp = abs((cpw-cpm)/cpm)
	      write(* ,'(i10,2e13.3,f15.8,e13.3)')
     +           niter,dpmax,dprms,cpw,errcp
	      write(10,'(i10,2e13.3,f15.8,e13.3)')
     +           niter,dpmax,dprms,cpw,errcp
	   endif
c
	   if(dpmax.lt.epspsi) GOTO 321
c
	   niter = niter + 1
	end do
c
321	CONTINUE
c
	write(*,*)'========================='
	write(*,*)'end V4 iteration'
	write(*,*)'========================='
c
	cpw   = cppv4()
	errcp = abs((cpw-cpm)/cpm)
	write(*,113)'cP  (final)       ',cpw,errcp
	write(*,113)'cP  (momentum)    ',cpm,0.
113	format(a20,f14.9,e12.3)
c
c       final sheet output: x, r, gamma (=Gam/ds), vn, dpsi
c
	do i = 1,np
	   if(i.eq.1)then
	      tx = zsl(2)-zsl(1)
	      tr = rsl(2)-rsl(1)
	   else if(i.eq.np)then
	      tx = zsl(np)-zsl(np-1)
	      tr = rsl(np)-rsl(np-1)
	   else
	      tx = zsl(i+1)-zsl(i-1)
	      tr = rsl(i+1)-rsl(i-1)
	   endif
	   tn = sqrt(tx*tx+tr*tr)
	   vn = (-tr*(1.+vz(i))+tx*vr(i))/tn
	   write(3,110) zsl(i),rsl(i),ga(i)/dsc(i),vn,
     +                  psiwake-psv(i)
	end do
110	format(5e16.7)
c
	call cpu_time(etime)
	write(*,*)'elapsed seconds: ',etime-stime
c
	close(3)
	close(10)
c
	end
c
c       --------------------------------------------------------------
	subroutine initv4
c       --------------------------------------------------------------
c
	use mem
	use memv4
c
	character*20 inpstr
	real gg, phi, r0l, u0, uw
c
	pih = 2.*atan(1.)
	pi  = 4.*atan(1.)
c
	OPEN(UNIT= 5,FORM='FORMATTED',STATUS='UNKNOWN',
     +       FILE='inpa-V4.dat')
c
	read(5,*) inpstr,ct
	read(5,*) inpstr,np
	read(5,*) inpstr,ltfac
	read(5,*) inpstr,beta
	read(5,*) inpstr,dampa
	read(5,*) inpstr,dampg
	read(5,*) inpstr,maxiter
	read(5,*) inpstr,epspsi
	read(5,*) inpstr,c1
	read(5,*) inpstr,c2
	read(5,*) inpstr,ncpout
	close(5)
c
	write(*,*)'***** VortexCode V4 *****'
	write(*,*)'cT      ',ct
	write(*,*)'N       ',np
	write(*,*)'LtubeFac',ltfac
	write(*,*)'beta    ',beta
	write(*,*)'dampA/G ',dampa,dampg
	write(*,*)'maxiter ',maxiter
	write(*,*)'epsPSI  ',epspsi
c
c       momentum theory values
c
	vzinf   = sqrt(1.+ct)
	vz0     = 0.5*(1.+vzinf)
	rsivc   = sqrt(vz0/vzinf)
	gasivc  = 1.-sqrt(1.+ct)
	psiwake = 0.25*(1.+sqrt(1.+ct))
	a       = 0.5*(sqrt(1.+ct)-1.)
	cpm     = 4.*a*(1.+a)**2
c
c       the tube starts at zsivc = ltfac*Rwake
c
	zsivc   = ltfac*rsivc
c
	allocate(zsl(0:np))
	allocate(rsl(0:np))
	allocate(ga (0:np))
	allocate(vz (0:np))
	allocate(vr (0:np))
	allocate(dsc (np))
	allocate(psv (np))
	allocate(gnew(np))
c
c       node grid: stretched cosine, ring 1 exactly at the edge
c
	do i = 1,np
	   phi = float(i-1)*pi/(c1*float(np-1))
	   gg  = (1.-cos(phi))/(1.-cos(pi/c1))
	   zsl(i) = zsivc*gg**c2
	end do
	zsl(1) = 0.
c
c       initial shape: van Kuik Eq. (D.10)-like
c
	rsl(1) = 1.
	u0 = 1.-0.5*gasivc
	do i = 2,np
	   r0l = 1. + float(i-1)/float(np-1)*(rsivc-1.)
	   uw  = 1.-0.5*gasivc*(1.+zsl(i)/
     +           sqrt(r0l*r0l+zsl(i)*zsl(i)))
	   rsl(i) = sqrt(u0/uw)
	end do
c
	call dscupd
c
c       initial circulations from the far-wake strength
c
	do i = 1,np
	   ga(i) = gasivc*dsc(i)
	end do
c
	OPEN(UNIT= 3,FORM='FORMATTED',STATUS='UNKNOWN',
     +       FILE='sls-V4.DAT')
	OPEN(UNIT=10,FORM='FORMATTED',STATUS='UNKNOWN',
     +       FILE='conv-V4.DAT')
c
	return
	end
c
c       --------------------------------------------------------------
	subroutine dscupd
c       --------------------------------------------------------------
c
c       centered arclength weights ds_i; the last ring gets half the
c       spacing to its neighbour plus half the gap to the tube start
c
	use mem
	use memv4
c
	real seg, segm, gap
c
	do i = 1,np
	   if(i.eq.1)then
	      seg = sqrt((zsl(2)-zsl(1))**2+(rsl(2)-rsl(1))**2)
	      dsc(1) = 0.5*seg
	   else if(i.eq.np)then
	      segm = sqrt((zsl(np)-zsl(np-1))**2
     +                   +(rsl(np)-rsl(np-1))**2)
	      gap  = sqrt((zsivc-zsl(np))**2+(rsivc-rsl(np))**2)
	      dsc(np) = 0.5*segm + 0.5*gap
	   else
	      seg  = sqrt((zsl(i+1)-zsl(i))**2
     +                   +(rsl(i+1)-rsl(i))**2)
	      segm = sqrt((zsl(i)-zsl(i-1))**2
     +                   +(rsl(i)-rsl(i-1))**2)
	      dsc(i) = 0.5*(seg+segm)
	   endif
	end do
c
	return
	end
c
c       --------------------------------------------------------------
	subroutine ringsv4(zp,rp,vzs,vrs,pss)
c       --------------------------------------------------------------
c
c       vz, vr, psi at the point (zp,rp) induced by ALL np rings
c       (circulations ga(i) at nodes zsl(i),rsl(i)), with the
c       spacing-proportional vortex kernel (Marshall-type smoothed
c       expressions inside rho < 0.5*beta*ds_i).
c       Branlard ring formulas (35.5/35.7/35.8) with Gam = gamma*ds.
c
	use mem
	use memv4
c
	real zp,rp,vzs,vrs,pss
	real ks,ak,f1,g1,rho,dloc,ee,ek
c
	vzs = 0.
	vrs = 0.
	pss = 0.
c
	do j = 1,np
c
	   rho  = sqrt((zp-zsl(j))**2+(rp-rsl(j))**2)
	   dloc = beta*dsc(j)
c
	   if(rho.lt.0.5*dloc)then
c
c             smoothed (kernel) expressions
c
	      vzs = vzs - ga(j)/(4.*pi*rsl(j))
     +                    *(log(16.*rsl(j)/dloc)-0.25)
	      pss = pss - ga(j)*rsl(j)/(2.*pi)
     +                    *(log(16.*rsl(j)/dloc)-1.5)
	   else
c
	      ks = 4.*rp*rsl(j)/((rsl(j)+rp)**2+(zp-zsl(j))**2)
	      if(ks.gt.1.-1.e-13) ks = 1.-1.e-13
	      if(ks.lt.1.e-14)    ks = 1.e-14
	      ek = ellk(ks)
	      ee = elle(ks)
	      ak = sqrt(ks)
c
	      f1 = rsl(j)/(2.*rp)*ks/(1.-ks)-(2.-ks)/(2.*(1.-ks))
	      vzs = vzs - ga(j)*ak*sqrt(rsl(j)/rp)
     +                    /(4.*pi*rsl(j))*(f1*ee+ek)
c
	      g1 = (2.-ks)/(2.*(1.-ks))*ee - ek
	      vrs = vrs - ga(j)*ak*(zp-zsl(j))
     +                    *(rsl(j)/rp)**1.5
     +                    /(4.*pi*rsl(j)*rsl(j))*g1
c
	      pss = pss - ga(j)*sqrt(rp*rsl(j))/(2.*pi)
     +                    *((2./ak-ak)*ek-(2./ak)*ee)
	   endif
c
	end do
c
	return
	end
c
c       --------------------------------------------------------------
	real function cppv4()
c       --------------------------------------------------------------
c
c       cP = cT*(1 + mean(vz at disc)), Simpson over 3000 points
c
	use mem
	use memv4
c
	integer nd
	real vzp,vrp,psp,rr,f1,f2,acc
c
	nd  = 3000
	drd = 1./float(nd)
	acc = 0.
c
	do i = 2,nd-2,2
	   rr = drd*float(i)
	   call ringsv4(0.,rr,vzp,vrp,psp)
	   f1 = 2.*(vzp + vzsivcBr(0.,rr))
	   call ringsv4(0.,rr+drd,vzp,vrp,psp)
	   f2 = 2.*(vzp + vzsivcBr(0.,rr+drd))
	   acc = acc + 2.*rr*f1 + 4.*(rr+drd)*f2
	end do
c
	call ringsv4(0.,1.,vzp,vrp,psp)
	f1 = 2.*(vzp + vzsivcBr(0.,1.))
	call ringsv4(0.,drd,vzp,vrp,psp)
	f2 = 2.*(vzp + vzsivcBr(0.,drd))
	acc = (acc + f1 + 4.*drd*f2)*drd/3.
c
	cppv4 = ct*(1. + acc)
c
	return
	end
