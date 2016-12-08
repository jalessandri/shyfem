c
c $Id: subwat.f,v 1.15 2009-05-21 09:24:00 georg Exp $
c
c water volume, surface and flux routines
c
c contents :
c
c subroutine connod(k,dz,con,coe)	changes concentration in node k
c
c subroutine volz3(k,dvol)                      inputs water volume
c
c subroutine volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)
c       computes volume and total mass of concentration in column of node k
c subroutine volno0(k,lmax,nlvddi,s,dvol,dcon)
c       inputs concentration in finite volume (node) k
c
c revision log :
c
c 30.04.1998    ggu     routines from subflx.f merged into this file
c 20.06.1998    ggu     two dimensional common blocks regolarized (nen3v,...)
c 21.08.1998    ggu     xv eliminated
c 22.10.1999    ggu     file cleaned, added 3D routines
c 20.11.2001    ggu     input concentration with ambient value ($AMB0)
c 25.03.2003    ggu     new routine zrise
c 28.04.2009    ggu     links re-structured 
c 
c*****************************************************************

	subroutine connod(k,dz,con,coe)

c changes concentration according to a water level rise dz with conz. con
c in node k on variable coe
c
c k	node where to input
c dz	water level change 
c con	concentraion of water injected
c coe	variable to change
c
c written 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
c revised 19.01.94 by ggu   $$conz - implementation of concentration
c revised 20.01.94 by ggu   $$lumpc - evaluate conz at node
c revised 04.12.97 by ggu   concentration adjusted in own routine
c 16.01.2001 by ggu   new concentration unique for node

	use mod_geom
	use mod_hydro
	use evgeom
	use basin

	implicit none

c arguments
	integer k
	real con,dz
	real coe(3,1)
c local
	integer ie,i,ii,nl
	integer ibase
	integer elems(maxlnk)
	real depth
	real area,vol,dvol,voltot,cnew
	real massold,massnew
	real dvoltot

	integer ithis

	call get_elems_around(k,maxlnk,nl,elems)

	voltot = 0.
	dvoltot = 0.
	cnew = 0.
	massold = 0.
	massnew = 0.

	do i=1,nl
	  ie=elems(i)
	  if(ie.le.0) stop 'error stop connod: internal error'
	  ii = iikthis(k,ie)
	  area = get_finvol_area_of_element(ie)
	  depth=hm3v(ii,ie)+zenv(ii,ie)-dz		!$$lumpc
	  vol = area * depth
	  dvol = area * dz
	  dvoltot = dvoltot + dvol
	  voltot = voltot + vol + dvol
	  massold = massold + coe(ii,ie)*vol
	  cnew = cnew + coe(ii,ie)*vol + con*dvol
	end do

	do i=1,nl
	  ie=elems(i)
	  ii = iikthis(k,ie)
	  coe(ii,ie) = cnew / voltot
	end do

	massnew = cnew
c	write(6,*) 'ggu0: ',dvoltot
c	write(6,*) 'ggu1: ',massold,massnew,massnew-massold

	massnew = 0.
	do i=1,nl
	  ie=elems(i)
	  ii = iikthis(k,ie)
	  area = get_finvol_area_of_element(ie)
	  depth=hm3v(ii,ie)+zenv(ii,ie)		!$$lumpc
	  vol = area * depth
	  massnew = massnew + coe(ii,ie)*vol
	end do
c	write(6,*) 'ggu2: ',massold,massnew,massnew-massold

	end

c*****************************************************************

	subroutine volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)

c computes volume and total mass of concentration in column of node k
c + volume of upper layer
c new version that computes only up to layer lmax (lmax > 0)

	use levels

	implicit none

c arguments
	integer k		!node defining column			(in)
	integer lmax		!maximum level to which compute		(in)
	integer nlvddi		!vertical dimension			(in)
	real s(nlvddi,1)	!variable (temperature, salinity,...)	(in)
	real area		!area of column 			(out)
	real vol		!total volume of column			(out)
	real vol0		!volume of first layer			(out)
	real svol		!total mass of s in column		(out)
c local
	integer l,ilevel,nlev
	integer mode
	real volume
c functions
	real volnode,areanode

	nlev = lmax
	if( lmax .le. 0 ) nlev = nlvddi		!all levels
	ilevel = min(nlev,ilhkv(k))

	mode = +1	!new time level
	vol=0.
	svol=0.

	do l=1,ilevel
	  volume = volnode(l,k,mode)
	  vol = vol + volume
	  svol = svol + volume * s(l,k)
	end do

	vol0 = volnode(1,k,mode)
	area = areanode(1,k)

	end

c******************************************************

	subroutine volno0(k,lmax,nlvddi,s,dvol,dcon)

c inputs concentration in finite volume (node) k
c ( 3d version ) -> only up to layer lmax

	use levels

	implicit none

c arguments
	integer k		!node defining volume
	integer lmax		!maximum level to which introduce
	integer nlvddi		!vertical dimension
	real s(nlvddi,1)	!variable (temperature, salinity,...)
	real dvol		!change in water volume (whole time step)
	real dcon		!concentration of new volume dvol
c local
	logical debug
	integer l,nlev
	real area,vol,vol0,svol
	real alpha,beta,gamma

	debug=.false.

        if( dcon .le. -990. ) return !input with ambient concentration $AMB0

	call volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)

	nlev = ilhkv(k)
	if( lmax .gt. 0 .and. nlev .gt. lmax ) nlev = lmax

	alpha=dvol/(vol+dvol)
	beta=vol0/(vol0+dvol)
	gamma=dvol*svol/(vol*(vol0+dvol))
	gamma=dvol*dcon+alpha*(svol-dcon*vol)
	gamma=gamma/(vol0+dvol)

	do l=1,nlev
	  if(l.eq.1) then
	    s(l,k)=beta*s(l,k)+alpha*beta*(dcon-s(l,k))+gamma
	  else
	    s(l,k)=s(l,k)+alpha*(dcon-s(l,k))
	  end if
	  if(debug) write(6,*) k,l,s(l,k)
	end do

	end

c******************************************************

