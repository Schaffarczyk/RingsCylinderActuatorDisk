c
c	module mem instead of "evil" common
c
	module mem
c
	real, allocatable :: rsl(:),  zsl(:), ga(:)
c
c       rsl,zsl: radial and axial position of vortex rings
c	ga:      strength of vortex rings
c
	real, allocatable :: rsln(:), zsln(:)
c
c       new resl and zsl
c
	real, allocatable :: vz(:), vr(:)
c
c       vz, vr: axial and radial velocities
c
	integer np,maxiter,niter,modelit,model
c
c       number of vortex rings
c
	real pi, pih, rsivc, zsivc, gasivc, ct, cpm, psiwake
	real under
c
c       pi, 0.5*pi
c	rSIVC, zSIVC, gaSIVC
c	radius, start and strength of Semi-Infite Vortex Cylinder
c
	logical update, selfind, restart	
c
	end module mem
