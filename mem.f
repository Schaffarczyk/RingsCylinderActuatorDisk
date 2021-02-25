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
	integer np
c
c       number of vortex rings
c
	real pi, pih, rsivc, zsivc, gasivc
c
c       pi, 0.5*pi
c	zSIVC, gaSIVC: start and strength of Semi-Infite Vortex Cylinder
c	
	end module mem
