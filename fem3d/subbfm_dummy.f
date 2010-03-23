c
c $Id: subbfm_dummy.f,v 1.1 2008-04-29 15:31:32 georg Exp $
c
c bfm dummy module
c
c revision log :
c
c 29.04.2008    ggu     bfm model integrated in main branch
c
c**************************************************************

	subroutine bfm_module(it,dt)

c administers bfm ecological model

	implicit none

	integer it
	real dt

	real getpar

	integer ibfm
	save ibfm
	data ibfm / 0 /

        if( ibfm .lt. 0 ) return

        if( ibfm .eq. 0 ) then
          ibfm = nint(getpar('ibfm'))
          if( ibfm .le. 0 ) ibfm = -1
          if( ibfm .lt. 0 ) return
	end if

	write(6,*) 'BFM module has not been linked - ibfm = ',ibfm

	stop 'error stop bfm_module: ibfm'

	end

c**************************************************************
