c
	program testElliptic
c
	use mem
c
	real k, n, a1, a2
	real m, m0, a, f1, f2, f3, f4
	real rc, zc
c
	pi     = 4.*atan(1.)
	ct     = 1.
	gasivc = 1.-sqrt(1.+ct)
c
c----------------------------------------------------------------------------------------
c
	write(*,*)
	write(*,*)'************* Test Elliptic Integrals ************'
c
	write(*,*)'single values ? ja = 1'
	read(*,*)ising
c
	if (ising.ne.1)goto 777	
c
999	write(*,*)'input k n'
	read(*,*) k, n
	a1 = k
	a2 = n
c
	write(*,100)'k   elle',k,elle(a1)
	write(*,100)'k   ellk',k,ellk(a1)
	write(*,101)'k n ellp',k,n,ellp(a1,a2)
	write(*,*) 'weiter ? (=1)'
	read(*,*)iw
	if(iw.eq.1)goto 999
c
777	write(*,*)
	write(*,*)'************* Test stream function SIVC ************'
c
	val = 0.	
c
        rsivc = 1.
	write(*,*)'z r ?'
	read(*,*)z,r
c
        m  = 4.*r*rsivc/((r+rsivc)*(r+rsivc)+z*z)
        m0 = 4.*r*rsivc/((r+rsivc)*(r+rsivc))
c
c       psi = r * PSI_theta
c
	a = gasivc*sqrt(rsivc*r)*z/(2.*pi*m0)

	f1 = (1.-m0+(m0/m))*ellk(m)
c
	f2 = -(m0/m)*elle(m)	
c
        f3 = (m0-1.)*ellp(m0,m)
c
	val = a*(f1 +f2 +f3)  
c
	if(r.gt.rsivc)then
	   f4 = 0.25*gasivc*rsivc**2
        else
           f4 = 0.25*gasivc*r**2
        end if
c
	val =val + f4
c
	write(*,778)'z r m m0 f1 f2 f3 f4',
     +                z,r,'**',m,m0,'**',f1,f2,'**',
     +                f3,f4,'**',psiBrCyl
778	format(a22,4(2e12.4,a2),e12.4)
c
	write(*,*)'r z val ',r,z, val
	write(*,*)'weiter ? ja = 1'
	read(*,*) iw
	if (iw.eq.1) goto 777
c

	write(*,*)'Ende'
c
100	format(a10, f20.3,f16.10)
101	FORMAT(a10,2f10.3, f16.10)
c
	end
