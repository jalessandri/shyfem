c
c $Id: subuti.f,v 1.11 2009-05-21 09:24:00 georg Exp $
c
c utility routines for 2d/3d model
c
c contents :
c
c getxy(k,x,y)					coordinates x/y for node k
c getexy(ie,x,y)				coordinates x/y for element ie
c
c real function areael(ie)			area for element ie
c real function areavl(k)			area for finite volume k
c
c function flxnod(k)            		discharge around node k
c
c subroutine energ(ielem,kenerg,penerg)		kinetic & potential energy
c
c subroutine stmima(a,nkn,nlvdi,ilhkv,amin,amax)
c                                       computes min/max of 3d field
c subroutine n2ebar(cn,ce)
c		copies concentrations from node value to element value
c
c revision log :
c
c 19.08.1998	ggu	new routines volco0, volno0
c 26.08.1998	ggu	subroutine stmima transferred from newbcl0
c 18.12.1999	ggu	/iweich/ -> /iwegv/ (bug)
c 29.03.2000	ggu	new routine getxy and getexy
c 16.05.2000	ggu	routine volel removed
c 28.04.2009    ggu     links re-structured
c 08.06.2010    ggu     new routine for computing 3D kin/pot energy
c 07.07.2011    ggu     bug fix in areael (unstable computation)
c 29.04.2015    ggu     energy now in Joule
c
c******************************************

	subroutine getxy0(k,x,y)

c gets coordinates x/y for node k

	use basin

	implicit none

	integer k
	real x,y

	x = xgv(k)
	y = ygv(k)

	end

c******************************************

	subroutine getexy0(ie,x,y)

c gets coordinates x/y for element ie

	use basin

	implicit none

	integer ie
	real x(3), y(3)

	integer k,ii,n
	integer kn(3)

	x = 0.
	y = 0.
	call basin_get_vertex_nodes(ie,n,kn)

	do ii=1,n
	  k = kn(ii)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	end

c******************************************

        subroutine baric0(ie,x,y)

c finds baricentre of element
c
c ie            number of element
c x,y           coordinates of baricentre (return value)

        use basin

        implicit none

        integer ie
        real x,y

        call basin_element_average2(ie,xgv,ygv,x,y)

        end

c******************************************

	function areael(ie)

c area for element ie
c
c double precision version - bug fix 07.07.2011

	use basin

	implicit none

c arguments
	real areael
	integer ie
c local
	integer n,kn(3)
	double precision x1,x2,x3,y1,y2,y3
	double precision, parameter :: half = 0.5

	call basin_get_vertex_nodes(ie,n,kn)

	if( n /= 3 ) then
	  stop 'error stop areael: not ready'	!FIXME_GGU
	end if

	x1=xgv(kn(1))
	y1=ygv(kn(1))
	x2=xgv(kn(2))
	y2=ygv(kn(2))
	x3=xgv(kn(3))
	y3=ygv(kn(3))

	areael = half * ( (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )

	end

c******************************************

	function areavl(k)

c area for finite volume k

	use mod_geom
	use evgeom
	use basin

	implicit none

c arguments
	real areavl
	integer k
c local
	logical blink
	integer ie,ii,i,nl
	integer elems(maxlnk)
	real area

	blink = .true.
	area=0.

	if( blink ) then

	  call get_elems_around(k,maxlnk,nl,elems)

          do i=1,nl
            ie=elems(i)
	    area = area + get_finvol_area_of_element(ie)
          end do

	else

	  do ie=1,nel
	    do ii=1,3
	      if( link_is_k(ii,ie,k) ) then
	        area = area + get_finvol_area_of_element(ie)
	      end if
	    end do
	  end do

	end if

	areavl = area

	end

c*****************************************************************

        function flxnod(k)

c computes discharge (m**3/sec) into finite volume around node k
c ... value flxnod has to be multiplied by dt to obtain total volume
c
c depending on value of blink uses link structure or not
c
c discharge into node n:     Q = 12 * aj * ( b(n)*U + c(n)*V )
c formula for 1D:            Q =  6 * aj *   b(n)*U
c volume difference:         dV = dt * Q

	use mod_geom
	use mod_hydro_baro
	use evgeom
	use basin

        implicit none

c arguments
        real flxnod
        integer k
c local
        integer i,nl,ie,ii
	integer elems(maxlnk)
        logical blink
        real flux,area,b,c

        blink = .true.
        flux=0.

        if( blink ) then

	  call get_elems_around(k,maxlnk,nl,elems)

          do i=1,nl
            ie=elems(i)
            ii=iikthis(k,ie)
	    area = get_total_area_of_element(ie)
	    call get_derivative_of_vertex(ii,ie,b,c)
            flux = flux + area * (unv(ie)*b+vnv(ie)*c)
          end do

        else

          do ie=1,nel
            do ii=1,3
	      if( link_is_k(ii,ie,k) ) then
	        area = get_total_area_of_element(ie)
                flux = flux + area * (unv(ie)*b+vnv(ie)*c)
              end if
            end do
          end do

        end if

        flxnod = flux

        end

c*****************************************************************

	subroutine energ(ielem,kenergy,penergy)

c computation of kinetic & potential energy [Joule]
c
c	pot = (g/2) * rho * area * z*z
c	kin = (1/2) * rho * area * (U*U+V*V)/H

	use mod_depth
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use basin

	implicit none

	integer ielem		!element (0 for total energy - all elements)
	real kenergy		!kinetic energy (return)
	real penergy		!potential energy relative to z=0 (return)

	include 'pkonst.h'

	integer ie,ii,ie1,ie2,n
	double precision area,pot,kin,z,zz

	if(ielem.gt.0.and.ielem.le.nel) then
	  ie1 = ielem
	  ie2 = ielem
	else
	  ie1 = 1
	  ie2 = nel
	end if

	kin=0.
	pot=0.
	do ie=ie1,ie2
	  call get_vertex_area_of_element(ie,n,area)
	  area = n * area	!total area of element
	  zz=0.
	  do ii=1,n
	    z = zenv(ii,ie)
	    zz = zz + z*z
	  end do
	  zz = zz / n
	  pot=pot+area*zz
	  kin=kin+area*(unv(ie)**2+vnv(ie)**2)/hev(ie)
	end do

        penergy = 0.5*grav*pot
        kenergy = 0.5*kin

	end

c***************************************************************

	subroutine energ3d(kenergy,penergy,ia_ignore)

c computation of kinetic & potential energy [Joule]
c
c	pot = (g/2) * rho * area * z*z
c	kin = (1/2) * rho * area * (U*U+V*V)/H

	use mod_layer_thickness
	use mod_ts
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	real kenergy		!kinetic energy (return)
	real penergy		!potential energy relative to z=0 (return)
	integer ia_ignore	!area code to be ignored

	include 'pkonst.h'

	integer ie,ii,l,lmax,ia,k,n
	double precision area,pot,kin,z,zz
	double precision h,uu,vv,rho

	kin=0.
	pot=0.

	do ie=1,nel

	  ia = iarv(ie)
	  if( ia .eq. ia_ignore ) cycle

	  call get_vertex_area_of_element(ie,n,area)
	  area = n * area	!total area of element

	  zz=0.
	  do ii=1,n
	    z = zenv(ii,ie)
	    zz = zz + z*z
	  end do
	  rho = rowass + basin_element_average(nlvdi,1,ie,rhov)
	  zz = zz / n
          pot = pot + area * rho * zz

	  lmax = ilhv(ie)
	  do l=1,lmax
	    rho = rowass + basin_element_average(nlvdi,l,ie,rhov)
	    h = hdenv(l,ie)
	    uu = utlnv(l,ie)
	    vv = vtlnv(l,ie)
	    kin = kin + area * rho * (uu*uu + vv*vv) / h
	  end do

	end do

        penergy = 0.5*grav*pot
        kenergy = 0.5*kin

	end

c***************************************************************

        subroutine stmima(a,nkn,nlvddi,ilhkv,amin,amax)

c computes min/max of 3d field

        implicit none

c arguments
        integer nkn,nlvddi
        real a(nlvddi,nkn)
        integer ilhkv(nkn)
        real amin,amax
c local
        integer lmax,k,l

        amin=a(1,1)
        amax=amin

        do k=1,nkn
          lmax=ilhkv(k)
          do l=1,lmax
            if(a(l,k).gt.amax) amax=a(l,k)
            if(a(l,k).lt.amin) amin=a(l,k)
          end do
        end do

        end

c**********************************************************************

        subroutine n2ebar(cn,ce)

c copies concentrations from node value to element value (only wet areas)

	use mod_geom_dynamic
	use evgeom
	use basin

        implicit none

        real cn(nkn)
        real ce(3,nel)

        integer ie,ii,k,n
	integer kn(3)

        do ie=1,nel
          if( iwegv(ie) > 0 ) cycle
	  call basin_get_vertex_nodes(ie,n,kn)
          do ii=1,n
            k = kn(ii)
            ce(ii,ie) = cn(k)
          end do
	  if( n == 2 ) then
	    ce(3,ie) = 0.5 * ( ce(1,ie) + ce(2,ie) )
	  end if
        end do

        end

c******************************************************

