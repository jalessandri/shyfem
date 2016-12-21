c
c $Id: links.f,v 1.8 2009-05-21 09:24:00 georg Exp $
c
c geometric routines (static)
c
c contents :
c
c revision log :
c
c 10.08.2003    ggu	new routines for ieltv, kantv
c 18.10.2005    ggu	more debug output, new routines get_link, get_lenk
c 06.04.2009    ggu	changed nlidim to nlkddi
c 30.03.2011    ggu	bverbose introduced
c 30.05.2015    ggu	new subroutine mklenkii to compute lenkiiv
c 02.12.2015    ggu	routines get_link() and get_lenk() deleted
c
c****************************************************************
c****************************************************************
c****************************************************************

        subroutine mklenk(nlkddi,ilinkv,lenkv)

c sets up vector with element links and a pointer to it
c
c only array nen3v is needed, no aux array is needed
c
c ilinkv    pointer to links
c lenkv     link to element numbers
c
c number of links of node n : nl = ilinkv(n+1)-ilinkv(n)
c links of node n           : ( lenkv ( ilinkv(n)+i ), i=1,nl )

	use basin

        implicit none

c arguments
        integer nlkddi
        integer ilinkv(nkn+1)
        integer lenkv(nlkddi)
c local
        logical binside,bdebug
        integer nloop
        integer ie,ii,i,n,k,nv,kk
	integer kn(3)
        integer ip,ip0,ip1,ipe
        integer ie0,ie1,k0,k1
        integer nfill
	integer n1d,n2d
	integer ie1d(ngr),ie2d(ngr)

	ip = 0
	ip1 = 0

c-------------------------------------------------------------------
c first the total number of elements (links) for each node is established
c-------------------------------------------------------------------

	bdebug = .false.

        ilinkv=1           !one more link -> will be corrected later

        do ie=1,nel
	  call basin_get_vertex_nodes(ie,nv,kn)
          do ii=1,nv
            k=kn(ii)
            ilinkv(k)=ilinkv(k)+1
          end do
        end do

c-------------------------------------------------------------------
c now create pointer into array
c-------------------------------------------------------------------

        n=0
        do k=1,nkn
          i=ilinkv(k)
          ilinkv(k)=n
          n=n+i
        end do
        ilinkv(nkn+1)=n

        if( n .gt. nlkddi ) goto 98

c-------------------------------------------------------------------
c now enter the element numbers of the links
c-------------------------------------------------------------------

        lenkv=0

        do ie=1,nel
	  call basin_get_vertex_nodes(ie,nv,kn)
          do ii=1,nv
            k=kn(ii)
	    if( k .gt. nkn ) goto 97
            ip=ilinkv(k)+1      !first entry
            ip1=ilinkv(k+1)     !last possible entry
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              ip=ip+1
            end do
            if(ip.gt.ip1) goto 99       !error 1
            lenkv(ip)=ie
          end do
        end do

c-------------------------------------------------------------------
c sort element entries
c-------------------------------------------------------------------

        do k=1,nkn

	  !bdebug = ( k == 13 )

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)

	  do while( lenkv(ip1) == 0 )
	    ip1 = ip1 - 1
	  end do

	  if( ip1 < ip0 ) goto 96

	  if( bdebug ) then
	    write(6,*) 'start node = ',k
	    do ip=ip0,ip1
		ie = lenkv(ip)
		write(6,*) ie,(nen3v(ii,ie),ii=1,3)
	    end do
	    write(6,*) 'end node = ',k
	  end if

c         ----------------------------------------------------------
c         copy 1d elements to end of list
c         ----------------------------------------------------------

	  if( basin_has_1d() ) then
	    n1d = 0
	    n2d = 0
	    ie1d = 0
	    ie2d = 0
	    do ip=ip0,ip1
	      ie = lenkv(ip)
	      if( basin_element_is_1d(ie) ) then
	        n1d = n1d + 1
		ie1d(n1d) = ie
	      else
	        n2d = n2d + 1
		ie2d(n2d) = ie
	      end if
	    end do
	    lenkv(ip0:ip0+n2d-1) = ie2d(1:n2d)
	    lenkv(ip0+n2d:ip0+n2d+n1d-1) = ie1d(1:n1d)
	    !from here on check
	    ip = ip0
	    ipe = ip1			!this is end of total list
	    do while( ip <= ipe )
	      if( basin_element_is_1d(lenkv(ip)) ) exit
	      ip = ip + 1
	    end do
	    ip1 = ip - 1		!this is end of 2d list
	    do while( ip <= ipe )
	      if( .not. basin_element_is_1d(lenkv(ip)) ) goto 95
	      ip = ip + 1
	    end do
	  end if

c         ----------------------------------------------------------
c         if k is boundary node, find first element (in anti-clockwise sense)
c         ----------------------------------------------------------

          binside = .true.
          nloop = 0
          ip=ip0

	  if( bdebug ) write(6,*) 'starting: ',ip0,ip1,ipe

          do while( binside .and. nloop .le. ip1-ip0 )
            ipe=ip
	    ie = lenkv(ipe)
            k1=kknext(k,ie)
            ip=ip0
	    if( bdebug ) write(6,*) 'start: ',ipe,ie,k1,ip,ip1
            do while(ip.le.ip1)
	      ie = lenkv(ip)
	      kk = kknext(k1,ie)
	    if( bdebug ) write(6,*) 'comparing: ',ie,ip,k,k1,kk
              if(kk.eq.k) goto 1
              ip=ip+1
            end do
            binside = .false.
    1       continue
	    if( bdebug ) write(6,*) 'end: ',ipe,ie,k1,ip,nloop,binside
            nloop = nloop + 1
          end do

c         --------------------------------------------------------
c         at ipe is first element --> now swap
c         --------------------------------------------------------

	  if( bdebug ) then
	    write(6,*) 'first element: ',ip0,ipe
	    write(6,*) 'first element: ',lenkv(ip0),lenkv(ipe)
	  end if

          if( .not. binside ) then      !boundary element
	    call swap_element_i(ip0,ipe,lenkv)
          end if

c         ----------------------------------------------------------
c         sort next elements
c         ----------------------------------------------------------

          do while(ip0.le.ip1-1.and.lenkv(ip0).ne.0)
            k1=kkbhnd(k,lenkv(ip0))
            ip=ip0+1
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              if(kknext(k,lenkv(ip)).eq.k1) exit
              ip=ip+1
            end do

            ip0=ip0+1

            if(ip.le.ip1.and.lenkv(ip).ne.0) then  !swap here
	      call swap_element_i(ip0,ip,lenkv)
            end if
          end do

        end do  !loop over all nodes

c-------------------------------------------------------------------
c compact structure of non boundary nodes
c
c for elements we could theoretically compact all entries,
c but we want to reuse ilinkv also for the node pointer,
c and therefore we need one more entry for boundary nodes
c-------------------------------------------------------------------

        nfill = 0

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( lenkv(ip1) .eq. 0 ) ip1 = ip1 - 1

          ie0 = lenkv(ip0)
          ie1 = lenkv(ip1)
	  !write(6,*) 'a: ',ie0,ie1,ip0,ip1
          k0 = kknext(k,ie0)
          k1 = kkbhnd(k,ie1)
	  !write(6,*) 'b: ',k0,k1
          if( k0 .ne. k1 ) ip1 = ip1 + 1        !boundary node

          ilinkv(k) = nfill
          do ip=ip0,ip1
            nfill = nfill + 1
            lenkv(nfill) = lenkv(ip)
          end do

        end do

        ilinkv(nkn+1) = nfill

c-------------------------------------------------------------------
c check structure
c-------------------------------------------------------------------

	call checklenk(ilinkv,lenkv)

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

        return
   95   continue
	write(6,*) '2d elements after 1d elements'
	call link_element_info(k,ilinkv,lenkv)
        stop 'error stop mklenk : internal error (4)'
   96   continue
        write(6,*) 'ip0,ip1: ',ip0,ip1
        stop 'error stop mklenk : internal error (3)'
   97   continue
        write(6,*) 'node: ',k,'  nkn: ',nkn
        stop 'error stop mklenk : internal error (2)'
   98   continue
        write(6,*) n,nlkddi
        write(6,*) nkn,nel,2*nkn+nel
        stop 'errro stop mklenk: nlkddi'
   99   continue
        !write(6,*) k,ilinkv(k),ip,ip1
        write(6,*) 'node: ',k
        write(6,*) 'element: ',ie
        write(6,*) 'first entry: ',ilinkv(k)+1
        write(6,*) 'last entry: ',ip1
        write(6,*) 'actual pointer: ',ip
        stop 'error stop mklenk : internal error (1)'
        end

c****************************************************************

        subroutine checklenk(ilinkv,lenkv)

c checks structure of lenkv

	use basin

        implicit none

c arguments
        integer ilinkv(nkn+1)
        integer lenkv(*)
c local
	integer nbnd,nnull
	integer k,k0,k1,ie,ie0,ie1
	integer ip,ip0,ip1,ip2d
	integer i
	integer ipp,ii
	logical bverbose,bbound,bhas1d,bhas2d

	bverbose = .false.

        nbnd = 0        !total number of boundary nodes
        nnull = 0       !total number of 0 entries

c-------------------------------------------------------------------
c loop over nodes
c-------------------------------------------------------------------

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)

	  bbound = .false.
          if( lenkv(ip1) .eq. 0 ) then
	    ip1 = ip1 - 1
	    bbound = .true.
	  end if

	  do ip=ip0,ip1
	    if( basin_element_is_1d(lenkv(ip)) ) exit
	  end do
	  ip2d = ip - 1
	  bhas1d = ( ip2d /= ip1 )
	  bhas2d = ( ip2d >= ip0 )

	  do ip=ip2d+1,ip1
	    if( .not. basin_element_is_1d(lenkv(ip)) ) then
	      write(6,*) '2d after 1d elements found: '
	      call link_element_info(k,ilinkv,lenkv)
              stop 'error stop checklenk: structure of lenkv (4)'
	    end if
	  end do

	  if( bhas2d .and. bhas1d .and. .not. bbound ) then
	    write(6,*) 'has 1d elements but is not boundary'
	    call link_element_info(k,ilinkv,lenkv)
            stop 'error stop checklenk: structure of lenkv (5)'
	  end if

          do ip=ip0,ip1
            ie = lenkv(ip)
            if( ie .eq. 0 ) then	!0 in index found
              write(6,*) 'Node (internal) k = ',k
              write(6,*) k,ip0,ip1
              write(6,*) (lenkv(i),i=ip0,ip1)
              stop 'error stop checklenk: structure of lenkv (2)'
            end if
          end do

          do ip=ip0+1,ip2d
            ie0 = lenkv(ip-1)
            ie1 = lenkv(ip)
            k0 = kkbhnd(k,ie0)
            k1 = kknext(k,ie1)
            if( k0 .ne. k1 ) then	!something is wrong
              write(6,*) 'Node (internal) k = ',k
              write(6,*) ip0,ip1,ip
              write(6,*) ie0,ie1,k0,k1
              write(6,*) ie0,(kiithis(i,ie0),i=1,3)
              write(6,*) ie1,(kiithis(i,ie1),i=1,3)
              write(6,*) (lenkv(i),i=ip0,ip1)
	      write(6,*) 'element list: '
	      do ipp=ip0,ip1
		ie = lenkv(ipp)
		write(6,*) ie,(nen3v(ii,ie),ii=1,3)
	      end do
              stop 'error stop checklenk: structure of lenkv (3)'
            end if
          end do

          ie0 = lenkv(ip0)
          ie1 = lenkv(ip1)
          k0 = kknext(k,ie0)
          k1 = kkbhnd(k,ie1)
          if( k0 .ne. k1 ) nbnd = nbnd + 1

          if( lenkv(ilinkv(k+1)) .eq. 0 ) nnull = nnull + 1

        end do

	if( bverbose ) then
          write(6,*) 'checklenk: ',nnull,nbnd,nkn
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

        end

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine mklenkii(nlkddi,ilinkv,lenkv,lenkiiv)

c sets up vector lenkiiv
c
c only array nen3v is needed, no aux array is needed
c
c ilinkv    pointer to links
c lenkv     link to element numbers
c
c number of links of node n : nl = ilinkv(n+1)-ilinkv(n)
c links of node n           : ( lenkv ( ilinkv(n)+i ), i=1,nl )

	use basin

        implicit none

c arguments
        integer nlkddi
        integer ilinkv(nkn+1)
        integer lenkv(nlkddi)
        integer lenkiiv(nlkddi)
c local
        integer ie,i,n,k,ii,ibase

c-------------------------------------------------------------------
c compute the vertex number for elements connected to node k
c-------------------------------------------------------------------

	lenkiiv = 0

        do k=1,nkn
	  n = ilinkv(k+1)-ilinkv(k)
	  ibase = ilinkv(k)
	  if( lenkv(ibase+n) .le. 0 ) n = n - 1
	  do i=1,n
	    ie = lenkv(ibase+i)
	    ii = iikthis(k,ie)
	    if( ii == 0 ) goto 99
	    lenkiiv(ibase+i) = ii
	  end do
        end do

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	return
   99	continue
	write(6,*) 'k,ie,n,ibase ',k,ie,n,ibase
	stop 'error stop mklenkii: internal error (1)'
	end

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine mklink(ilinkv,lenkv,linkv)

c sets up vector with node links
c
c ilinkv and lenkv must have already been set up
c
c only array nen3v is needed, no aux array is needed
c
c ilinkv    pointer to links
c lenkv     link to element numbers
c linkv     link to node numbers

	use basin

        implicit none

c arguments
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer linkv(*)
c local
	logical bdebug
        integer ie,k
        integer ip,ip0,ip1,ip2d

c-------------------------------------------------------------------
c loop over nodes
c-------------------------------------------------------------------

	bdebug = .false.

        do k=1,nkn

	  bdebug = k == 4007

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( lenkv(ip1) .eq. 0 ) ip1 = ip1 - 1

	  do ip=ip0,ip1
	    if( basin_element_is_1d(lenkv(ip)) ) exit
	  end do
	  ip2d = ip-1

	  if( bdebug ) write(6,*) ip0,ip1,ip2d

          do ip=ip0,ip2d
            ie = lenkv(ip)
            linkv(ip) = kknext(k,ie)
          end do

c         ----------------------------------------------------------
c         if boundary node insert last node
c         ----------------------------------------------------------

          ip=ilinkv(k+1)
          if( lenkv(ip) .eq. 0 ) then
            linkv(ip2d+1) = kkbhnd(k,ie)
          end if

c         ----------------------------------------------------------
c         if 1d nodes insert now
c         ----------------------------------------------------------

          do ip=ip2d+1,ip1
            ie = lenkv(ip)
            linkv(ip+1) = kknext(k,ie)
	  end do

	  if( bdebug ) then
	    write(6,*) 'start node list'
            ip1=ilinkv(k+1)
	    do ip=ip0,ip1
	      write(6,*) ip,linkv(ip)
	    end do
	    write(6,*) 'end node list'
	    call link_element_info(k,ilinkv,lenkv)
	  end if

        end do

c-------------------------------------------------------------------
c check structure
c-------------------------------------------------------------------

	call checklink(ilinkv,lenkv,linkv)

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

        end

c****************************************************************

        subroutine checklink(ilinkv,lenkv,linkv)

c checks structure of linkv

	use basin

        implicit none

c arguments
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer linkv(*)
c local
	integer k,k1,i
	integer ip,ip0,ip1
	integer ipk,ipk0,ipk1
	logical bverbose

	bverbose = .false.

	k1 = 0				!compiler warnings
	ipk0 = 0
	ipk1 = 0

c-------------------------------------------------------------------
c loop over nodes
c-------------------------------------------------------------------

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( linkv(ip1) .eq. 0 ) then	!0 in index found
            write(6,*) 'Node (internal) k = ',k
            write(6,*) ip0,ip1
            write(6,*) (linkv(ip),ip=ip0,ip1)
	    call link_element_info(k,ilinkv,lenkv)
            stop 'error stop checklink: internal error (1)'
          end if

	  do ip=ip0,ip1
	    k1 = linkv(ip)
            ipk0=ilinkv(k1)+1
            ipk1=ilinkv(k1+1)
	    ipk = ipk0
	    do while( ipk .le. ipk1 .and. linkv(ipk) .ne. k )
	      ipk = ipk + 1
	    end do
	    if( ipk .gt. ipk1 ) then	!node not found
              write(6,*) 'Node (internal) k = ',k
              write(6,*) ip0,ip1
              write(6,*) (linkv(i),i=ip0,ip1)
              write(6,*) k1,ipk0,ipk1
              write(6,*) (linkv(i),i=ipk0,ipk1)
              stop 'error stop checklink: internal error (2)'
	    end if
	  end do
        end do

	if( bverbose ) then
          write(6,*) 'checklink: ',nkn
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

        end

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine mkkant(ilinkv,lenkv,linkv,kantv)

c makes vector kantv

	use basin

        implicit none

c arguments
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer linkv(*)
        integer kantv(2,nkn)
c local
	integer k
	integer ip,ip0,ip1

        do k=1,nkn
          kantv(1,k)=0
          kantv(2,k)=0
        end do

        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( lenkv(ip1) .eq. 0 ) then		!boundary node
            kantv(1,k) = linkv(ip0)
            kantv(2,k) = linkv(ip1)
          end if
        end do

c-------------------------------------------------------------------
c check structure
c-------------------------------------------------------------------

	call checkkant(kantv)

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

        end

c****************************************************************

        subroutine checkkant(kantv)

c checks structure of kantv

	use basin

        implicit none

c arguments
        integer kantv(2,nkn)
c local
	integer k,k1,k2
	integer nbnd,nint
	logical bverbose

	bverbose = .false.

	nbnd = 0
	nint = 0

c-------------------------------------------------------------------
c loop over nodes
c-------------------------------------------------------------------

        do k=1,nkn
	  k1 = kantv(1,k)
	  k2 = kantv(2,k)
	  if( k1 .gt. 0 .and. k2 .gt. 0 ) then
	    nbnd = nbnd + 1
	    if( k .ne. kantv(2,k1) .or. k .ne. kantv(1,k2) ) then
              write(6,*) 'Node (internal) k = ',k
	      write(6,*) 'k1,k2: ',k1,k2
	      write(6,*) 'backlink: ',kantv(2,k1),kantv(1,k2)
	      stop 'error stop checkkant: structure of kantv (2)'
	    end if
	  else if( k1 .eq. 0 .and. k2 .eq. 0 ) then
	    nint = nint + 1
	  else
            write(6,*) 'Node (internal) k = ',k
	    write(6,*) 'k1,k2: ',k1,k2
	    stop 'error stop checkkant: structure of kantv (1)'
	  end if
        end do

	if( bverbose ) then
	  write(6,*) 'checkkant: ',nkn,nint,nbnd
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

        end

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine mkielt(ilinkv,lenkv,linkv,ieltv)

c makes vector ieltv (without open boundary nodes)

	use basin

        implicit none

c arguments
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer linkv(*)
        integer ieltv(3,nel)
c local
	integer k,ie,ii
	integer ip,ip0,ip1
	integer inext

c-------------------------------------------------------------------
c initialize
c-------------------------------------------------------------------

        do ie=1,nel
          do ii=1,3
            ieltv(ii,ie) = 0
          end do
        end do

c-------------------------------------------------------------------
c loop over links
c-------------------------------------------------------------------

        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1) - 1

          do ip=ip0,ip1
            ie=lenkv(ip)
            ieltv(inext(k,ie),ie)=lenkv(ip+1)
          end do

c	  -------------------------------------------------------------------
c	  last element
c	  -------------------------------------------------------------------

          ie=lenkv(ip1+1)
          if(ie.gt.0) then                              !is internal node
            ieltv(inext(k,ie),ie)=lenkv(ip0)
          end if

        end do

c-------------------------------------------------------------------
c check structure
c-------------------------------------------------------------------

	call checkielt(ieltv)

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c****************************************************************

        subroutine checkielt(ieltv)

c checks structure of ieltv

	use basin

        implicit none

c arguments
        integer ieltv(3,nel)
c local
	integer k,ie,ii
	integer kn,inn,ien,ienn
        integer nbnd,nobnd,nintern
	integer inext,knext,kthis
	logical bverbose

c-------------------------------------------------------------------
c initialize
c-------------------------------------------------------------------

	bverbose = .false.

        nbnd = 0
        nobnd = 0
        nintern = 0

c-------------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------------

        do ie=1,nel
          do ii=1,3
            ien = ieltv(ii,ie)
            if( ien .gt. 0 ) then
            	nintern = nintern + 1
                if( ien .gt. nel ) goto 99
                k = kthis(ii,ie)
                kn = knext(k,ie)
                inn = inext(kn,ien)
                ienn = ieltv(inn,ien)
                if( ie .ne. ienn ) goto 98
            else if( ien .eq. 0 ) then
                nbnd = nbnd + 1
            else if( ien .eq. -1 ) then
                nobnd = nobnd + 1
            else
                goto 99
            end if
          end do
        end do

	if( bverbose ) then
          write(6,*) 'checkielt is ok'
          write(6,*) '  internal sides =      ',nintern
          write(6,*) '  boundary sides =      ',nbnd
          write(6,*) '  open boundary sides = ',nobnd
          write(6,*) '  total sides =         ',nintern+nbnd+nobnd
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

        return
   98   continue
        write(6,*) 'Element (internal) ie = ',ie
        write(6,*) 'ii,ien,nel: ',ii,ien,nel
        write(6,*) 'k,kn,inn,ienn: ',k,kn,inn,ienn
        write(6,*) 'nen3v: ',ie,(kthis(ii,ie),ii=1,3)
        write(6,*) 'nen3v: ',ien,(kthis(ii,ien),ii=1,3)
        write(6,*) 'ieltv: ',ie,(ieltv(ii,ie),ii=1,3)
        write(6,*) 'ieltv: ',ien,(ieltv(ii,ien),ii=1,3)
        stop 'error stop checkielt: corrupt data structure of ieltv (1)'
   99   continue
        write(6,*) 'Element (internal) ie = ',ie
        write(6,*) 'ii,ien,nel: ',ii,ien,nel
        stop 'error stop checkielt: corrupt data structure of ieltv (2)'
	end

c****************************************************************

        subroutine update_ielt(ibound,ieltv)

c updates vector ieltv with open boundary nodes

	use basin

        implicit none

c arguments
        integer ibound(nkn)	! >0 => open boundary node
        integer ieltv(3,nel)
c local
	integer k,ie,ii,i
	integer ip,ip0,ip1
	integer knext,ibhnd,kthis

c-------------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------------

        do ie=1,nel
          do ii=1,3
	    k = kthis(ii,ie)
	    if( ibound(k) .gt. 0 .and. ibound(knext(k,ie)) .gt. 0 ) then
	      i = ibhnd(k,ie)
              ieltv(i,ie) = -1
	    end if
          end do
        end do

c-------------------------------------------------------------------
c check structure
c-------------------------------------------------------------------

	call checkielt(ieltv)

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine mkdxy(dxv,dyv)

c initializes dxv,dyv

	use basin

        implicit none

c arguments
	real dxv(nkn),dyv(nkn)
c local
	integer k,ie,ii,i
	integer ip,ip0,ip1
	integer knext,ibhnd,kthis

c-------------------------------------------------------------------
c loop over nodes
c------------------------------------------------------------------

        do k=1,nkn
	  dxv(k) = 0.
	  dyv(k) = 0.
        end do

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine swap_element_i(i1,i2,array)

	implicit none

	integer i1,i2
	integer array(*)

	integer iaux

        iaux=array(i1)
        array(i1)=array(i2)
        array(i2)=iaux

	end

c****************************************************************

	subroutine link_element_info(k,ilinkv,lenkv)

	use basin

	implicit none

	integer k
        integer ilinkv(nkn+1)
        integer lenkv(*)

	integer ie,ip0,ip1,ip,n,ii

        ip0=ilinkv(k)+1
        ip1=ilinkv(k+1)
        if( lenkv(ip1) .eq. 0 ) ip1 = ip1 - 1

	write(6,*) 'node = ',k
	write(6,*) 'ip0,ip1: ',ip0,ip1
	do ip=ip0,ip1
	  ie = lenkv(ip)
	  if( ie == 0 ) then
	    write(6,*) ip,ie,0
	  else
	    n = basin_get_vertex_of_element(ie)
	    write(6,*) ip,ie,n,(nen3v(ii,ie),ii=1,n)
	  end if
	end do

	end

c****************************************************************

