c
c $Id: bdist.f,v 1.5 2009-05-21 09:24:00 georg Exp $
c
c distance routines
c
c contents :
c
c subroutine shdist(rdist)
c                       makes distance array from open boundaries
c subroutine mkdist(nkn,idist,rdist)
c                       makes distance array from given nodes
c
c revision log :
c
c 20.08.2003    ggu     new routine wr2dn()
c 05.01.2005    ggu     routines for writing nos file into subnsoa.f
c 07.01.2005    ggu     documentation for shdist
c 28.04.2009    ggu     links re-structured
c
c****************************************************************

        subroutine shdist(rdist)

c makes distance array from open boundaries
c
c rdist is contained between 0 and 1
c if nadist (=d) is not given rdist = 1 (default)
c otherwise the first d2=d/2 rows of nodes have rdist = 0
c and then the next ones have rdist = i/d with i = d2+1, d2+d
c the rest has again rdist = 1
c
c example: nadist = d = 4,   d2 = d/2 = 2
c
c   row i:   1   2   3   4   5   6   7   8   ...
c   rdist:   0   0  1/4 2/4 3/4  1   1   1   ...

	implicit none

	include 'param.h'

        real rdist(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,neldim)			!element index
        common /nen3v/nen3v

c local variables

        integer idist(nkndim)

        integer i,k,kk
        integer nadist
        integer ibc,n,itype,nk
        real r,d,d2

	integer iapini,ipint
        integer nbnds,itybnd,nkbnds,kbnds
        real getpar

c-----------------------------------------------------------------
c get parameters
c-----------------------------------------------------------------

        do k=1,nkn
          rdist(k) = 1.
          idist(k) = 0
        end do

        nadist = nint(getpar('nadist'))
        if( nadist .le. 0 ) return

c-----------------------------------------------------------------
c gather open boundary nodes
c-----------------------------------------------------------------

        n = 0

        do ibc=1,nbc
          itype = itybnd(ibc)
          if( itype .eq. 1 .or. itype .eq. 2 ) then
            nk = nkbnds(ibc)
            do i=1,nk
              k = kbnds(ibc,i)
              idist(k) = 1
            end do
          end if
        end do

c-----------------------------------------------------------------
c make distance
c-----------------------------------------------------------------

        write(6,*) 'Making distance rdist'
        call mkdist(nkn,idist,rdist)

c-----------------------------------------------------------------
c adjust distance
c-----------------------------------------------------------------

        d = nadist
        d2 = d / 2.

        do k=1,nkn
          r = rdist(k)
          r = r - d2            !first rows no adv terms
          r = r / d             !slowly introduce
          r = max(0.,r)
          r = min(1.,r)
          !if( rdist(k) .eq. 1 ) write(6,*) 'bdist ',rdist(k),d,d2,r
          rdist(k) = r
        end do

c-----------------------------------------------------------------
c write dist (nos) file
c-----------------------------------------------------------------
 
        call wrnos2d('dist','distance from boundary nodes',rdist)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c*************************************************************************** 
                         
        subroutine mkdist(nkn,idist,rdist)

c makes distance array from given nodes
c
c rdist of open boundary nodes is 1
c other nodes are > 1 (integer)
c example: neibors of rdist=1 nodes have rdist=2 etc.

        implicit none

        integer nkn
        integer idist(1)
        real rdist(1)

	include 'links.h'

        integer k,kk,i
        integer n
        integer idact,idnew,nfound

c----------------------------------------------------------
c initialize with first level
c----------------------------------------------------------

c       idist is already initialized with first level

c----------------------------------------------------------
c loop on levels
c----------------------------------------------------------

        idact = 1
        nfound = nkn

        do while( nfound .gt. 0 )

          idnew = idact + 1
          nfound = 0

          do k=1,nkn
            if( idist(k) .eq. idact ) then

	      call set_node_links(k,n)

              do i=1,n
                kk = lnk_nodes(i)
                if( idist(kk) .eq. 0 ) then
                  idist(kk) = idnew
                  nfound = nfound + 1
                end if
              end do
            end if
          end do

          idact = idnew

        end do

c----------------------------------------------------------
c set real value
c----------------------------------------------------------

        do k=1,nkn
          rdist(k) = idist(k)
        end do

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

        end

c*************************************************************************** 
