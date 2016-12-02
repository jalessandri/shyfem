
!==================================================================
        module average
!==================================================================

	implicit none

        INTERFACE element_to_node
        MODULE PROCEDURE 
     +				 element_to_node_2d
     +				,element_to_node_3d
     +				,element_to_node_3d_minmax
        END INTERFACE

        INTERFACE node_to_element
        MODULE PROCEDURE 
     +				 node_to_element_2d
     +				,node_to_element_3d
        END INTERFACE

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
	  call get_vertex_area_of_element(ie,n,area)
          value = elv(ie)
          do ii=1,n
            k = nen3v(ii,ie)
            nov(k) = nov(k) + area*value
            aux(k) = aux(k) + area
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

        nov = nov / aux

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
	  call get_vertex_area_of_element(ie,n,area)
	  lmax = ilhv(ie)
	  do l=1,lmax
            value = elv(l,ie)
            do ii=1,n
              k = nen3v(ii,ie)
              nov(l,k) = nov(l,k) + area*value
              aux(l,k) = aux(l,k) + area
	    end do
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

	where ( aux > 0. ) 
	  nov = nov / aux
	end where

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end subroutine element_to_node_3d

c******************************************************************

	pure subroutine element_to_node_3d_minmax(mode,elv,nov)

c transforms element values to nodal values (no weights - use min/max)
c
c (3D version)

	use levels
	use basin

	implicit none

c arguments
	integer, intent(in) :: mode	 	   !min (-1),  max (+1) or aver
        real, intent(in)    :: elv(nlvdi,nel)      !array with element values
        real, intent(out)   :: nov(nlvdi,nkn)      !array with nodal values

c local
        integer k,ie,ii,l,lmax,n
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
	  n = basin_get_vertex_of_element(ie)
	  lmax = ilhv(ie)
	  do l=1,lmax
            value = elv(l,ie)
            do ii=1,n
              k = nen3v(ii,ie)
	      if( mode .eq. 1 ) then
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

        pure subroutine node_to_element_2d(nov,elv)

c transforms nodal values to element values
c
c (2D version)

	use basin

        implicit none

        real, intent(in)  :: nov(nkn)     !array with nodal values (in)
        real, intent(out) :: elv(nel)     !array with element values (out)

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
	  call get_vertex_area_of_element(ie,n,area)
	  do ii=1,n
	    k = nen3v(ii,ie)
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

!==================================================================
        end module average
!==================================================================

