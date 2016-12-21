c
c $Id: suplin1.f,v 1.1 2009-11-18 17:16:00 georg Exp $
c
c utilitiy routines for section plot (velocities)
c
c revision log :
c
c 13.10.2009    ggu     routines written from scratch
c
c*******************************************************************

	subroutine prepare_vel(pp3)

	use mod_hydro_print
	use mod_hydro_plot
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real pp3(nlvdi,nkn)

	integer k,l
	real u,v,w
	real href

	call make_vertical_velocity	!compute wlnv from utlnv,vtlnv

	href = 0.
	call mkht3(nlvdi,het3v,href)

	call make_vel_from_tra(het3v)
	call vel_to_node
	
	do k=1,nkn
	  do l=1,nlv
	    u = uprv(l,k)
	    v = vprv(l,k)
	    w = wprv(l,k)
	    pp3(l,k) = sqrt( u*u + v*v + w*w )
	  end do
	end do

	end

c*******************************************************************

	subroutine make_vel_from_tra(het3v)

c from transports to velocities (on elements)

	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        real het3v(nlvdi,nel)		!layer depth at elements

	integer ie,l,lmax
	real h,rh

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    h = het3v(l,ie)
	    if( h .le. 0. ) goto 99
	    rh = 1. / h
	    ulnv(l,ie) = utlnv(l,ie) * rh
	    vlnv(l,ie) = vtlnv(l,ie) * rh
	  end do
	end do

	return
   99	continue
	write(6,*) 'ie,l,h: ',ie,l,h
	stop 'error stop make_vel_from_tra: zero depth'
	end

c*******************************************************************

	subroutine vel_to_node

c transfers velocities at elements to nodes 
c and vertical velocities to center of layer

	use mod_hydro_print
	use mod_hydro_vel
	use evgeom
	use levels
	use basin

	implicit none

	integer ie,ii,k,l,lmax,n
	integer kn(3)
	real area

        uprv=0.
        vprv=0.
        wprv=0.

        do ie=1,nel
	  call get_vertex_area_of_element(ie,n,kn,area)
	  lmax = ilhv(ie)
          do ii=1,n
            k=kn(ii)
	    do l=1,lmax
              wprv(l,k)=wprv(l,k)+area
              uprv(l,k)=uprv(l,k)+area*ulnv(l,ie)
              vprv(l,k)=vprv(l,k)+area*vlnv(l,ie)
            end do
	  end do
	end do

	where( wprv > 0 )
          uprv=uprv/wprv
          vprv=vprv/wprv
	end where

        do k=1,nkn
          do l=1,nlv
            wprv(l,k)=0.5*(wlnv(l,k)+wlnv(l-1,k))
          end do
        end do

	end

c*******************************************************************

