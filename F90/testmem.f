c
c
	program testmem
c
	real, allocatable :: test(:)
	integer status
c
111	write(*,*)'input np'
	if(np.lt.0)goto 112

	read (*,*) np
c
	call system('free -m')
c
	allocate(test(0:np),stat=status)
c
	call system('free -m')
c
	goto 111
112	end
c
