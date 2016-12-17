
!***************************************************************

	subroutine grd_to_basin

	use grd
	use basin

	implicit none

	integer nk,ne,nl,nne,nnl,nt
	integer k,ie,ii,ibase,n,i,nv,il
	integer kn(3)
	integer nhn,nhe,nhl,nht
	real flag
	logical bdebug

!--------------------------------------------------------
! initialize
!--------------------------------------------------------

	flag = -999.
	bdebug = .false.

	call grd_get_params(nk,ne,nl,nne,nnl)
	nt = ne + nl
	call basin_init(nk,nt)

!--------------------------------------------------------
! copy params
!--------------------------------------------------------

	dcorbas = dcor_grd
	dirnbas = dirn_grd
	descrr = title_grd

!--------------------------------------------------------
! copy arrays
!--------------------------------------------------------

	xgv = xv
	ygv = yv
	ipv = ippnv
	iarnv = ianv
	ipev(1:ne) = ippev(1:ne)
	iarv(1:ne) = iaev(1:ne)
	ipev(ne+1:ne+nl) = ne + ipplv(1:nl)
	iarv(ne+1:ne+nl) = ialv(1:nl)

!--------------------------------------------------------
! copy element index
!--------------------------------------------------------

	nen3v = 0

	do ie=1,ne
	  ibase = ipntev(ie-1)
	  n = ipntev(ie) - ipntev(ie-1)
	  if( n .ne. 3 ) then
	    write(6,*) 'element is not triangle: ',ie,ippev(ie)
	    stop 'error stop grd_to_basin: not a triangle'
	  end if
	  do ii=1,3
	    nen3v(ii,ie) = inodev(ibase+ii)
	  end do
	end do

!--------------------------------------------------------
! copy channel index
!--------------------------------------------------------

	do il=1,nl
	  ie = ne + il
	  ibase = ipntlv(il-1)
	  n = ipntlv(il) - ipntlv(il-1)
	  if( n .ne. 2 ) then
	    write(6,*) 'line is not channel: ',il,ipplv(il)
	    stop 'error stop grd_to_basin: not a channel'
	  end if
	  do ii=1,2
	    nen3v(ii,ie) = inodlv(ibase+ii)
	  end do
	end do

!--------------------------------------------------------
! check depth and copy
!--------------------------------------------------------

	nhe = 0
	do ie=1,ne
	  if( hhev(ie) .ne. flag ) nhe = nhe + 1
	end do

	nhl = 0
	do il=1,nl
	  if( hhlv(il) .ne. flag ) nhl = nhl + 1
	end do

	nhn = 0
	do k=1,nk
	  if( hhnv(k) .ne. flag ) nhn = nhn + 1
	end do

	nht = nhe + nhl

	if( nht > 0 .and. nhn > 0 ) then
	  write(6,*) 'nht,nhn: ',nht,nhn
	  if( nht == nel .and. nhn == nkn ) then
	    write(6,*) 'can use both depths...'
	    write(6,*) '... using element values'
	    nhn = 0
	  else if( nht == nel ) then
	    write(6,*) 'depth on nodes incomplete...'
	    write(6,*) '... using element values'
	    nhn = 0
	  else if( nhn == nkn ) then
	    write(6,*) 'depth on elements incomplete...'
	    write(6,*) '... using nodal values'
	    nht = 0
	  else
	    write(6,*) 'depth values are incomplete...'
	    stop 'error stop grd_to_basin: depth incomplete'
	  end if
	end if

	write(6,*) 'nht,nhn: ',nht,nhn

	if( nhn == 0 ) then
	  do ie=1,ne
	    hm3v(:,ie) = hhev(ie)
	  end do
	  do il=1,nl
	    ie = ne + il
	    hm3v(:,ie) = hhlv(il)
	  end do
	else
	  do ie=1,nel
	    call basin_get_vertex_nodes(ie,nv,kn)
	    do ii=1,nv
	      k = kn(ii)
	      hm3v(ii,ie) = hhnv(k)
	    end do
	  end do
	end if

	do ie=1,nel
	  if( basin_element_is_1d(ie) ) then
	    hm3v(3,ie) = 1.	!HACK ... set width
	  end if
	end do

	call sp13_set_1d

!--------------------------------------------------------
! debug output
!--------------------------------------------------------

	if( bdebug ) then

	do ie=1,nt,nt/10
	  write(6,*) ie,(nen3v(ii,ie),ii=1,3)
	end do

	do k=1,nk,nk/10
	  write(6,*) k,ippnv(k),xv(k),yv(k)
	end do

	end if

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!***************************************************************

	subroutine basin_to_grd

	use grd
	use basin

	implicit none

	integer nk,ne,nl,nne,nnl,il
	integer k,ie,ii,ibase,n,nv
	integer kn(3)
	integer nhn,nhe
	real flag
	logical bdebug
	logical bconst,bunique
	double precision h,he,hm

!--------------------------------------------------------
! initialize
!--------------------------------------------------------

	flag = -999.
	bdebug = .false.

	ne = nel_2d
	nl = nel - ne
	nne = 3*ne
	nnl = 2*nl

	write(6,*) 'basin_to_grd: ',nkn,ne,nl
	call grd_init(nkn,ne,nl,nne,nnl)

!--------------------------------------------------------
! copy params
!--------------------------------------------------------

	dcor_grd = dcorbas
	dirn_grd = dirnbas
	title_grd = descrr

!--------------------------------------------------------
! copy arrays
!--------------------------------------------------------

	xv(1:nkn) = xgv(1:nkn)
	yv(1:nkn) = ygv(1:nkn)
	ippnv(1:nkn) = ipv(1:nkn)
	ianv(1:nkn) = iarnv(1:nkn)
	ippev(1:ne) = ipev(1:ne)
	iaev(1:ne) = iarv(1:ne)
	ipplv(1:nl) = ipev(ne+1:ne+nl)
	ialv(1:nl) = iarv(ne+1:ne+nl)

!--------------------------------------------------------
! copy element index
!--------------------------------------------------------

	ibase = 0
	ipntev(0) = 0
	do ie=1,ne
	  ibase = ipntev(ie-1)
	  call basin_get_vertex_nodes(ie,nv,kn)
	  do ii=1,nv
	    ibase = ibase + 1
	    inodev(ibase) = kn(ii)
	  end do
	  ipntev(ie) = ibase
	end do

	ibase = 0
	ipntlv(0) = 0
	do il=1,nl
	  ie = ne + il
	  ibase = ipntlv(il-1)
	  call basin_get_vertex_nodes(ie,nv,kn)
	  do ii=1,nv
	    ibase = ibase + 1
	    inodlv(ibase) = kn(ii)
	  end do
	  ipntlv(il) = ibase
	end do

!--------------------------------------------------------
! check depth and copy
!--------------------------------------------------------

	hhnv = flag
	hhev = flag

	bunique = .true.
	bconst = .true.

	do ie=1,nel
	  he = hm3v(1,ie)
	  hm = 0.
	  call basin_get_vertex_nodes(ie,nv,kn)
	  do ii=1,nv
	    h = hm3v(ii,ie)
	    k = kn(ii)
	    if( hhnv(k) == flag ) hhnv(k) = h
	    if( hhnv(k) /= h ) bunique = .false.
	    if( he /= h ) bconst = .false.
	    hm = hm + h
	  end do
	  hm = hm / nv
	  if( nv == 3 ) then
	    hhev(ie) = hm
	  else
	    il = ie - ne
	    hhlv(il) = hm
	  end if
	end do

	if( bconst .and. bunique ) then		!constant depth
	  hhev = hm3v(1,1)
	  hhlv = hm3v(1,1)
	  hhnv = flag
	else if( bunique ) then			!unique depth at nodes
	  hhev = flag
	  hhlv = flag
	else
	  hhnv = flag
	end if
	  
!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!***************************************************************

