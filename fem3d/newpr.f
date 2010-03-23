c
c $Id: newpr.f,v 1.24 2010-02-22 15:38:36 georg Exp $
c
c administration of external output files
c
c contents :
c
c subroutine prious(mode)				writes OUT file
c function wrrc77(iunit,nvers,it,knausm,knaus,u,v,z)	writes EXT file
c subroutine locous
c
c revision log :
c
c 03.11.1993	ggu	$$cmplerr	compiler warnings hydro
c 03.11.1993	ggu	$$ITRES	never write at time itmres but idtres later
c 11.11.1993	ggu	$$N66		write to unit 66
c 05.04.1994	ggu	$$WRIT8		old writes commented out
c 24.06.1994	ggu	$$RST		open rst file only if desired
c 31.10.1995	ggu	$$OLDOUS	old ous file format
c 24.08.1998	ggu	new write to terminal
c 04.02.2000	ggu	no prious, dobefor/after3d,
c 04.02.2000	ggu	wrrc77 to subext
c 11.10.2002	ggu	write OUV file to variable unit (not 88 anymore)
c 21.08.2003	ggu	deleted locspc
c 25.11.2004	ggu	locous not called anymore
c 22.02.2010	ccf	new routine for tidal pot. (tideforc), locaus deleted
c
c************************************************************

	subroutine dobefor3d

c to do befor time step (3D specific)

        implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	call dobefor

	!call tidenew(it)	!tidal potential
        call tideforc(it)       !tidal potential !ccf

	end

c************************************************************

	subroutine doafter3d

c to do after time step (3D specific)

        implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	call doafter

	end

c************************************************************
