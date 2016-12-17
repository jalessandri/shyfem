
!==================================================================
        module average
!==================================================================

	implicit none

        INTERFACE		 element_to_node
        MODULE PROCEDURE 
     +				 element_to_node_2d
     +				,element_to_node_3d
     +				,element_to_node_2d_minmax
     +				,element_to_node_3d_minmax
        END INTERFACE

        INTERFACE		 vertex_to_element
        MODULE PROCEDURE 
     +				 vertex_to_element_2d
        END INTERFACE

        INTERFACE		 vertex_to_node
        MODULE PROCEDURE 
     +				 vertex_to_node_2d
        END INTERFACE

        INTERFACE		 node_to_vertex
        MODULE PROCEDURE 
     +				 node_to_vertex_2d
     +				,node_to_vertex_2d_mask
        END INTERFACE

        INTERFACE		 node_to_element
        MODULE PROCEDURE 
     +				 node_to_element_2d
     +				,node_to_element_3d
        END INTERFACE

        INTERFACE		 create_node_indicator
        MODULE PROCEDURE 
     +				 create_node_indicator_r
        END INTERFACE

        INTERFACE		 extract_on_vertices
        MODULE PROCEDURE 
     +				 extract_on_vertices_3d
        END INTERFACE

	PRIVATE ::
     +				 element_to_node_2d
     +				,element_to_node_3d
     +				,element_to_node_2d_minmax
     +				,element_to_node_3d_minmax
	PRIVATE ::
     +				 vertex_to_element_2d
	PRIVATE ::
     +				 vertex_to_node_2d
	PRIVATE ::
     +				 node_to_vertex_2d
     +				,node_to_vertex_2d_mask
	PRIVATE ::
     +				 node_to_element_2d
     +				,node_to_element_3d
	PRIVATE ::
     +				 create_node_indicator_r
	PRIVATE ::
     +				 extract_on_vertices_3d

!==================================================================
        contains
!==================================================================

	pure subroutine element_to_node_2d(elv,nov)

c transforms element values to nodal values (weights are area)
c
c (2D version)

	use evgeom
	use basin

	implicit none

c arguments
        real, intent(in)  :: elv(nel)     !array with element values
        real, intent(out) :: nov(nkn)     !array with nodal values

c local
        integer k,ie,ii,n
	integer kn(3)
        real area,value
        real aux(nkn)     !aux array (nkn)

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

        nov = 0.
        aux = 0.

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
	  call get_vertex_area_of_element(ie,n,kn,area)
          value = elv(ie)
          do ii=1,n
            k = kn(ii)
            nov(k) = nov(k) + area*value
            aux(k) = aux(k) + area
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

	where ( aux > 0. ) nov = nov / aux

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine element_to_node_2d

c******************************************************************

	pure subroutine element_to_node_3d(elv,nov)

c transforms element values to nodal values (weights are area)
c
c (3D version)

	use evgeom
	use levels
	use basin

	implicit none

c arguments
        real, intent(in)  :: elv(nlvdi,nel)     !array with element values
        real, intent(out) :: nov(nlvdi,nkn)     !array with nodal values

c local
        integer k,ie,ii,l,lmax,n
	integer kn(3)
        real area,value
        real aux(nlvdi,nkn)     !aux array (nkn)

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

        nov = 0.
        aux = 0.

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
	  call get_vertex_area_of_element(ie,n,kn,area)
	  lmax = ilhv(ie)
	  do l=1,lmax
            value = elv(l,ie)
            do ii=1,n
              k = kn(ii)
              nov(l,k) = nov(l,k) + area*value
              aux(l,k) = aux(l,k) + area
	    end do
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

	where ( aux > 0. ) nov = nov / aux

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine element_to_node_3d

c******************************************************************

	pure subroutine element_to_node_2d_minmax(mode,elv,nov)

c transforms element values to nodal values (no weights - use min/max)
c
c (2D version)

	use basin

	implicit none

c arguments
	integer, intent(in) :: mode	     !<0:min, >0:max, =0:aver
        real, intent(in)    :: elv(nel)      !array with element values
        real, intent(out)   :: nov(nkn)      !array with nodal values

c local
        integer k,ie,ii,n
	integer kn(3)
        real rinit,value

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

	if( mode < 0 ) then
	  nov = 1.e+30
	else if( mode > 0 ) then
	  nov = -1.e+30
	else
	  call element_to_node_2d(elv,nov)
	  return
	end if

c-----------------------------------------------------------
c get min/max values
c-----------------------------------------------------------

        do ie=1,nel
	  call basin_get_vertex_nodes(ie,n,kn)
          value = elv(ie)
          do ii=1,n
            k = kn(ii)
	    if( mode > 0 ) then
              nov(k) = max(nov(k),value)
	    else
              nov(k) = min(nov(k),value)
	    end if
	  end do
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine element_to_node_2d_minmax

c******************************************************************

	pure subroutine element_to_node_3d_minmax(mode,elv,nov)

c transforms element values to nodal values (no weights - use min/max)
c
c (3D version)

	use levels
	use basin

	implicit none

c arguments
	integer, intent(in) :: mode	 	   !<0:min, >0:max, =0:aver
        real, intent(in)    :: elv(nlvdi,nel)      !array with element values
        real, intent(out)   :: nov(nlvdi,nkn)      !array with nodal values

c local
        integer k,ie,ii,l,lmax,n
	integer kn(3)
        real rinit,value

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

	if( mode < 0 ) then
	  nov = 1.e+30
	else if( mode > 0 ) then
	  nov = -1.e+30
	else
	  call element_to_node_3d(elv,nov)
	  return
	end if

c-----------------------------------------------------------
c get min/max values
c-----------------------------------------------------------

        do ie=1,nel
	  call basin_get_vertex_nodes(ie,n,kn)
	  lmax = ilhv(ie)
	  do l=1,lmax
            value = elv(l,ie)
            do ii=1,n
              k = kn(ii)
	      if( mode > 0 ) then
                nov(l,k) = max(nov(l,k),value)
	      else
                nov(l,k) = min(nov(l,k),value)
	      end if
	    end do
          end do
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine element_to_node_3d_minmax

c******************************************************************
c******************************************************************
c******************************************************************

        pure subroutine vertex_to_element_2d(el3v,elv)

c transforms vertex values to element values
c
c (2D version)

	use basin

        implicit none

        real, intent(in)  :: el3v(3,nel)  !array with vertex values
        real, intent(out) :: elv(nel)     !array with element values

        integer ie

c-----------------------------------------------------------
c convert values
c-----------------------------------------------------------

        do ie=1,nel
          elv(ie) = basin_vertex_average(ie,el3v)
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine vertex_to_element_2d

c******************************************************************

        pure subroutine vertex_to_node_2d(el3v,nov)

c transforms vertex values to nodal values
c
c (2D version)

	use basin
	use evgeom

        implicit none

        real, intent(in)  :: el3v(3,nel)  !array with vertex values
        real, intent(out) :: nov(nkn)     !array with node values

        integer ie,ii,k,n
	integer kn(3)
	real area,value
	real aux(nkn)

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

        nov = 0.
        aux = 0.

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
	  call get_vertex_area_of_element(ie,n,kn,area)
          do ii=1,n
            k = kn(ii)
            value = el3v(ii,ie)
            nov(k) = nov(k) + area*value
            aux(k) = aux(k) + area
	  end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

	where ( aux > 0. ) nov = nov / aux

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine vertex_to_node_2d

c******************************************************************

	pure subroutine node_to_vertex_2d(nov,el3v)

c transforms nodal values to element values
c
c (2D version)

	use basin

        implicit none

        real, intent(in)  :: nov(nkn)     !array with nodal values
        real, intent(out) :: el3v(3,nel)  !array with vertex values

        integer ie,ii,k,n
	integer kn(3)

c-----------------------------------------------------------
c convert values
c-----------------------------------------------------------

        do ie=1,nel
	  call basin_get_vertex_nodes(ie,n,kn)
	  do ii=1,n
	    k = kn(ii)
	    el3v(ii,ie) = nov(k)
	  end do
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine node_to_vertex_2d

c******************************************************************

	pure subroutine node_to_vertex_2d_mask(nov,mask,el3v)

c transforms nodal values to element values
c
c (2D version)

	use basin

        implicit none

        real, intent(in)     :: nov(nkn)      !array with nodal values
        logical, intent(in)  :: mask(nel)     !array with mask
        real, intent(out)    :: el3v(3,nel)   !array with vertex values

        integer ie,ii,k,n
	integer kn(3)

c-----------------------------------------------------------
c convert values
c-----------------------------------------------------------

        do ie=1,nel
	  if( .not. mask(ie) ) cycle
	  call basin_get_vertex_nodes(ie,n,kn)
	  do ii=1,n
	    k = kn(ii)
	    el3v(ii,ie) = nov(k)
	  end do
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine node_to_vertex_2d_mask

c******************************************************************

        pure subroutine node_to_element_2d(nov,elv)

c transforms nodal values to element values
c
c (2D version)

	use basin

        implicit none

        real, intent(in)  :: nov(nkn)     !array with nodal values
        real, intent(out) :: elv(nel)     !array with element values

        integer ie

c-----------------------------------------------------------
c convert values
c-----------------------------------------------------------

        do ie=1,nel
          elv(ie) = basin_element_average(ie,nov)
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine node_to_element_2d

c******************************************************************

        pure subroutine node_to_element_3d(nov,elv)

c transforms nodal values to element values
c
c (3D version)

	use levels
	use basin

        implicit none

c arguments
        real, intent(in)  :: nov(nlvdi,nkn)	!array with nodal values
        real, intent(out) :: elv(nlvdi,nel)	!array with element values

c local
        integer ie,l,lmax

c-----------------------------------------------------------
c convert values
c-----------------------------------------------------------

        do ie=1,nel
          lmax = ilhv(ie)
          do l = 1,lmax
	    elv(l,ie) = basin_element_average(nlvdi,l,ie,nov)
          end do
	end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine node_to_element_3d

c******************************************************************

	pure subroutine average_over_nodes(val,aver)

c computes average of val (defined on nodes) over total basin
c
c 2D version

	use evgeom
	use basin

	implicit none

        real, intent(in)  :: val(nkn)     !array with element values
	real, intent(out) :: aver	  !average

        integer k,ie,ii,n
	integer kn(3)
        double precision area,value
        double precision accum,area_tot

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

	accum = 0.
	area_tot = 0.

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
	  call get_vertex_area_of_element(ie,n,kn,area)
	  do ii=1,n
	    k = kn(ii)
	    value = val(k)
	    accum = accum + value * area
	    area_tot = area_tot + area
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

	if ( area_tot > 0. ) accum = accum / area_tot

	aver = accum

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine average_over_nodes

c******************************************************************

	pure subroutine create_node_indicator_r(mask,ind)

	use basin

	implicit none

        logical, intent(in)  :: mask(nel)     !mask on elements
	real, intent(out)    :: ind(nkn)      !indicator

        integer k,ie,ii,n
	integer kn(3)

        ind = 0

        do ie=1,nel
	  call basin_get_vertex_nodes(ie,n,kn)
          if( mask(ie) ) then
            do ii=1,n
              k = kn(ii)
              ind(k) = 1
            end do
          end if
        end do

	end subroutine create_node_indicator_r

c******************************************************************

	pure subroutine extract_on_vertices_3d(l,ie,nov,vert)

	use basin
	use levels

	implicit none

	integer, intent(in)  :: l,ie
	real, intent(in)  :: nov(nlvdi,nkn)      !nodal value
        real, intent(out) :: vert(3)             !vertex values

        integer k,ii,n
	integer kn(3)

	call basin_get_vertex_nodes(ie,n,kn)

        do ii=1,n
          k = kn(ii)
	  vert(ii) = nov(l,k)
	end do
	
	end subroutine extract_on_vertices_3d

!==================================================================
        end module average
!==================================================================

