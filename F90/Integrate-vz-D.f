c
c	aps 2021 April 8
c
c       Integrate vzd to check psi
c
c
	real vz(1000), r(1000)
c
	OPEN(UNIT= 1,FORM='FORMATTED',STATUS='UNKNOWN',FILE='vzd.DAT')
c	
	i = 1
c
100	read(1,*,err=101,end=101) r(i),vz(i),du
c	write(*,*)i,r(i),vz(i)
	i = i + 1
	goto 100
c
101	i = i - 1
	write(*,*)'n= ',i
	write(*,*)
c
	psi = 0.
	do j = 2,i
	   rint = 0.5*(r(j)+r(j-1))
           vint = 0.5*(vz(j)+vz(j-1)-2.)
	   dr   = r(j)-r(j-1)
	   dpsi = rint*vint*dr
c	   write(*,*)rint, vint, dr
	   psi = psi + dpsi
	end do
c
	write(*,*)'psi disk = ',psi	
c
	end

