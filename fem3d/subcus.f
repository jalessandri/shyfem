c
c $Id: subcus.f,v 1.58 2010-03-08 17:46:45 georg Exp $
c
c custom routines
c
c contents :
c
c subroutine custom( it )	custom routines
c
c subroutine totvol( it )	computes total volume in basin
c subroutine georg( it )	test some nodes and elements
c subroutine fluxw( itact )	tests for mass balance
c subroutine temp( it )		q-flux test
c subroutine bocche		adjust depth at inlets
c
c revision log :
c
c 22.01.1998	ggu	custom routine called at end of time step
c 20.03.1998	ggu	custom routine in own file to avoid compiler warnings
c 08.05.1998	ggu	new custom routine to check for mass balance
c 25.06.1998	ggu	test for q-flux (zerlina)
c 21.08.1998    ggu     xv eliminated
c 03.09.1998    ggu     subroutine bocche to adjust depth at Venice inlets
c 03.10.1998    ggu     subroutine salt to check salt
c 19.04.1999    ggu     subroutine impli to change weighting factor
c 27.05.1999    ggu     use icust to call custom routines
c 05.12.2001    ggu     fixed compiler error with -Wall -pedantic
c 10.08.2003    ggu     use accessor routines for chezy values in anpa()
c 14.08.2003    ggu     new routines test_hakata (icust=26) and node3d
c 14.08.2003    ggu     only 3D in concmass1
c 02.09.2003    ggu     new routine lago (76)
c 04.09.2003    ggu     some fixes in routine anpa
c 05.10.2004    ggu     new routine aldo to set conz in area
c 22.02.2005    ggu     subroutines deleted: salt
c 14.03.2005    ggu     subroutine traccia
c 30.06.2005    ggu     traccia changed, new routines jamal, sedimt
c 01.12.2005    ggu     more changes in traccia
c 23.03.2006    ggu     changed time step to real
c 18.10.2006    ggu     jamal has been updated and comented
c 22.02.2007    ggu     traccie routines updated
c 23.05.2007    ggu     new routines oscillation, kreis, debora
c 03.08.2007    ggu     in jamal reset introduced...
c 12.03.2008    ggu     jamal restructured -> computes projected res time
c 17.03.2008    ggu     new routines zinit, cprint
c 26.06.2008    ggu     routien for testing diffusion
c 09.10.2008    ggu     new call to confop
c 11.10.2008	ggu	pass zfranco into subroutine
c 12.11.2008	ggu	new routine joel
c 19.11.2008	ggu	new routine viscos
c 06.12.2008	ggu	new routine vdiffus()
c 28.01.2009	aac	new routine andreac()
c 24.03.2009	ggu	new zinit, tvd_test
c 23.05.2009	ggu	bug fix in jamal - reset also outer area
c 16.10.2009	ggu	some changes and documentation for traccia
c 18.11.2009	ggu	residence time computation in jamal with correction
c 12.02.2010	ggu	new routines for horizontal diffusion tests (diffus2d)
c 01.03.2010	ggu	new routines jamal_fra() to reset conz
c
c******************************************************************

	subroutine custom( it )

c custom routines

	implicit none

	integer it

	real getpar

	integer icall
	save icall
	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  icall = nint(getpar('icust'))
	  if( icall .le. 0 ) icall = -1
	end if

	if( icall .eq.  6 ) call impli(it)
	if( icall .eq.  7 ) call anpa(it)
	if( icall .eq.  9 ) call channel(it)
	if( icall .eq. 10 ) call sedcust(it)
	if( icall .eq. 11 ) call velchan(it)
	if( icall .eq. 15 ) call concmass(it)
	if( icall .eq. 16 ) call concmass1(it)
	if( icall .eq. 21 ) call hakata(it)
	!if( icall .eq. 25 ) call close_inlets
	if( icall .eq. 25 ) call close_inlets1
	if( icall .eq. 26 ) call test_hakata(it)
	if( icall .eq. 27 ) call traccia
	if( icall .eq. 28 ) call oscillation
	if( icall .eq. 29 ) call kreis
	if( icall .eq. 31 ) call zinit
	if( icall .eq. 32 ) call cprint(it)
	if( icall .eq. 76 ) call lago(it)
	if( icall .eq. 80 ) call aldo(it)
	if( icall .eq. 81 ) call jamal(it)
	if( icall .eq. 811 ) call jamal_fra
	if( icall .eq. 82 ) call sedimt
	if( icall .eq. 83 ) call joel
	if( icall .eq. 90 ) call diffus
	if( icall .eq. 91 ) call viscos
	if( icall .eq. 92 ) call vdiffus(1)
	if( icall .eq. 93 ) call vdiffus(2)
	if( icall .eq. 94 ) call diffus2d
        if( icall .eq. 883 ) call debora(it)
        if( icall .eq. 884 ) call tsinitdebora(it)
        if( icall .eq. 888 ) call uv_bottom
        if( icall .eq. 601 ) call andreac
        if( icall .eq. 602 ) call tvd_test(it)

c	call totvol(it)
c	call georg(it)
c	call fluxw(it)
c	call temp(it)

	end

c*****************************************************************

	subroutine totvol( it )

c computes total volume in basin

	implicit none

	integer it

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer iwegv(1)
	common /iwegv/iwegv

	integer iwet,ie
	real area,voldry,volwet,vol,a,v

	real areaele,volele

	iwet = 0
	area = 0.
	voldry  = 0.
	volwet = 0.

	do ie=1,nel

	  a = areaele(ie)
	  v = volele(ie,+1)

	  if( iwegv(ie) .eq. 0 ) then
	    iwet = iwet + 1
	    area = area + a
	    volwet = volwet + v
	  else
	    voldry = voldry + v
	  end if

	end do

	vol = voldry + volwet

	write(99,'(i10,i5,4e16.8)') 'totvol: '
     +			,it,iwet,area,volwet,voldry,vol

	end

c*****************************************************************

	subroutine georg( it )

c test some nodes and elements

	implicit none

	integer it

	integer ii
	integer ipint,ieint
	integer n,ibase

	integer iwegv(1)
	common /iwegv/iwegv
	real znv(1)
        common /znv/znv
	real zenv(3,1)
	common /zenv/zenv

	integer icall,k1,k2,ie1,ie2,ie3
	save icall,k1,k2,ie1,ie2,ie3
c	data icall,k1,k2,ie1,ie2,ie3 /0,556,555,777,788,788/
	data icall,k1,k2,ie1,ie2,ie3 /0,337,368,513,514,511/

	if( icall .eq. 0 ) then
	  write(97,*) k1,k2,ie1,ie2,ie3
	  k1 = ipint(k1)
	  k2 = ipint(k2)
	  ie1 = ieint(ie1)
	  ie2 = ieint(ie2)
	  ie3 = ieint(ie3)
	  write(97,*) k1,k2,ie1,ie2,ie3
	  icall = 1
	end if

	write(97,*) it,iwegv(ie1),iwegv(ie2),iwegv(ie3)
	write(97,*) znv(k1),znv(k2)
	write(97,*) (zenv(ii,ie1),ii=1,3)
	write(97,*) (zenv(ii,ie2),ii=1,3)
	write(97,*) (zenv(ii,ie3),ii=1,3)

	end

c*****************************************************************

	subroutine fluxw( itact )

c tests for mass balance

	implicit none

	integer itact

	integer ndim,nedim,nscdim
	parameter(ndim=11,nedim=5,nscdim=10)
c	parameter(ndim=13,nedim=7,nscdim=10)
c	parameter(ndim=10,nedim=14,nscdim=10)
c	parameter(ndim=31,nscdim=10)

	logical berror
	integer i,ii,k,ie
	real az,azpar
	real dt
	real dv,dz,fd,ft
	real area

	integer kflxck,ieint
	real flux(nscdim)

	integer nsect
	save nsect
	integer icall
	save icall

	integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer n
	integer kflux(ndim)
	integer iflux(3,ndim)
	save kflux,iflux,n
	integer ielm(nedim)
	save ielm

	real areaele,zetaele

	data icall / 0 /
	data kflux /
     +			 1,81,587,588,589,446,0
     +			,1,447,446,0
     +		/
	data ielm /
     +			 6837,6836,896,902,901
     +		/
c	data kflux /
c     +			 582,534,533,532,568,0
c     +			,582,537,536,535,567,568,0
c     +		/
c	data ielm /
c     +			 890,761,762,760,759,842,843
c     +		/
c	data kflux /
c     +			 1181,1201,1200,1180,0
c     +			,1183,1206,1205,1178,0
c     +		/
c	data ielm /
c     +			 1880,1877,1869,1868,1875,1870,1867
c     +			,1879,1878,1876,1874,1873,1872,1871
c     +		/
c	data kflux /150,173,172,151,0,4311,4315,4314,4309,0
c     +			,1181,1201,1200,1180,0,1182,1204,1203,1202,1179,0
c     +			,1183,1206,1205,1178,0,1185,1208,1207,1177,0
c     +		/

	if( icall .eq. 0 ) then

	   icall = 1

	   do i=1,nedim
	     k = ieint(ielm(i))
	     if( k .le. 0 ) then
		write(6,*) ielm(i),k
		stop 'error stop fluxw(6)'
	     end if
	     ielm(i) = k
	   end do

	   do i=1,ndim
	     if( kflux(i) .ne. 0 ) n = i
	   end do

	   call flxe2i(kflux,n,berror)
	   if( berror ) stop 'error stop fluxw (1)'
	   nsect = kflxck(kflux,n)
	   if( nsect .lt. 0 ) stop 'error stop fluxw (2)'
	   call flxin(kflux,iflux,n)

	   if( nsect .gt. nscdim ) then
	      write(6,*) 'nsect,nscdim ',nsect,nscdim
	      stop 'error stop fluxw (3)'
	   end if

c	   write(97,*) 'initialization of sections : ',nsect,n
c	   write(97,*) (kflux(i),i=1,n)
c	   write(97,*) ((iflux(ii,i),ii=1,3),i=1,n)
	end if

	call getaz(azpar)
	az = azpar

	call flxscs(kflux,iflux,n,az,flux)

	write(97,*) it,(flux(i),i=1,nsect)

	dv = 0.
	do i=1,nedim
	  ie = ielm(i)
	  area = areaele(ie)
	  dz = zetaele(ie,+1) - zetaele(ie,-1)
	  dv = dv + dz * area
	end do

	call get_timestep(dt)
	fd = dt*(flux(1) - flux(2))
	ft = dt*(flux(1) + flux(2))
	write(97,*) dv,fd,ft,(dv-fd)/ft

	end

c*****************************************************************

	subroutine temp( it )

c q-flux test

	implicit none

	integer it

	include 'param.h'

	integer i
	real rit

	integer ndim
	parameter(ndim=5)

	real tempv(nlvdim,1)
	common /tempv/tempv

	integer nodes(ndim)
	save nodes
	data nodes /5,59,113,167,221/

	rit = it/3600.
	write(97,'(6f12.5)') rit,(tempv(1,nodes(i)),i=1,ndim)

	end

c*****************************************************************

	subroutine bocche

c adjust depth at inlets -> must be called from cst routines

	implicit none

	integer ie,ii,itype
	real hm,h,hlido,hmala,hchio
	real getpar

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer iarv(1)
	common /iarv/iarv
	real hev(1)
	common /hev/hev

	write(6,*) 'adjusting depth at inlets...'

	hlido = getpar('hlido')
	hmala = getpar('hmala')
	hchio = getpar('hchio')

	do ie = 1,nel

	  hm = hev(ie)
	  itype = iarv(ie)

	  if( itype .eq. 3 ) then	!chioggia
	    h = min(hm,hchio)
	  else if( itype .eq. 4 ) then	!malamocco
	    h = min(hm,hmala)
	  else if( itype .eq. 5 .or. itype .eq. 9 ) then	!lido
	    h = min(hm,hlido)
	  else
	    h = hm
	  end if

	  if( h .ne. hm ) then
	    write(6,*) 'depth changed: ',ie,itype,hm,h
	  end if

	  call setdepele(ie,h)

	end do

	end

c*****************************************************************

	subroutine impli(it)

c changes implicity

	implicit none

	integer it

	integer it1,it2
	real r
	double precision d

	save d
	data d /1./

	it1 = 30000000
	it2 = 60000000

	if( it .lt. it1 ) then
	  if( d .ne. 1. ) then
	    d = 1.
	    r = d
	    call putpar('azpar',r)
	    call putpar('ampar',r)
	    write(6,*) 'weight changed : ',r
	  end if
	else if( it .gt. it2 ) then
	  if( d .ne. 0.5 ) then
	    d = 0.5
	    r = d
	    call putpar('azpar',r)
	    call putpar('ampar',r)
	    write(6,*) 'weight changed : ',r
	  end if
	else
	  d = it-it1
	  d = d / (it2-it1)
	  d = 1. - 0.5 * d
	  r = d
	  call putpar('azpar',r)
	  call putpar('ampar',r)
	  write(6,*) 'weight changed : ',r
	end if

	end

c*****************************************************************

	subroutine anpa(it)

c changes chezy

	implicit none

	integer it

	integer itmax,itdt
	integer i
	real czmin,czact,czmax
	real czin,czout
	real czp(0:100),czm(0:100)
	save czp,czm
	integer nareas
	save nareas

	integer icall
	save icall
	data icall / 0 /

	itmax = 600000
	itdt = 7200
	czmin = 10.
	czmin = 5.
	czmin = 4.
	czmin = 3.

	if( it .lt. 0 ) return
	if( it .gt. itmax ) return

	if( icall .eq. 0 ) then
	  icall = 1

	  call n_chezy_values(nareas)
	  
	  do i=0,nareas
	    call get_chezy_values(i,czp(i),czm(i))
	  end do

	  write(6,*) 'anpa initialized...  ',nareas
	  write(6,*) (czp(i),i=0,20)
	  write(6,*) (czm(i),i=0,20)
	end if

	if( mod(it,itdt) .ne. 0 ) return

	do i=0,nareas
	  call get_chezy_values(i,czp(i),czm(i))
	end do

	write(15,'(i10,12i5)') it,(nint(czp(i)),i=0,nareas)
	write(15,'(i10,12i5)') it,(nint(czm(i)),i=0,nareas)

	do i=0,nareas
	  if( ( i .ge. 3 .and. i .le. 5 ) .or. i .eq. 9 ) then
	    call get_chezy_values(i,czin,czout)
	    if( czin  .gt. czmin ) czin  = czin  - 1.
	    if( czout .gt. czmin ) czout = czout - 1.
	    call set_chezy_values(i,czin,czout)
	  end if
	end do

	end

c*****************************************************************

	subroutine channel(it)

c channel output

	implicit none

	include 'param.h'

	integer it

        integer ilhkv(1)
        common /ilhkv/ilhkv
        real ulnv(nlvdim,neldim)
        common /ulnv/ulnv
        real vlnv(nlvdim,neldim)
        common /vlnv/vlnv
        real rhov(nlvdim,nkndim)
        common /rhov/rhov
 
        real num(0:nlvdim)
        real nuh(0:nlvdim)
        real tk(0:nlvdim)
        real ep(0:nlvdim)
        real rl(0:nlvdim)

	logical berror
	integer i,k,l,nlev,nlev1
	integer ipext

	integer ndim
	parameter( ndim = 3 )
	integer nodes(ndim)
	save nodes
	integer icall
	save icall

	data nodes / 25 , 28 , 31 /
	data icall / 0 /

	if( icall .eq. 0 ) then
	  icall = 1
	  call n2int(ndim,nodes,berror)
	  if( berror ) stop 'error stop channel: n2int'
	  write(6,*) 'total nodes in channel: ',ndim
	  do i=1,ndim
	    write(6,*) i,nodes(i),ipext(nodes(i))
	  end do
	  open(81,file='gvel.dat',status='unknown',form='formatted')
	  open(82,file='grho.dat',status='unknown',form='formatted')
	  open(89,file='gnum.dat',status='unknown',form='formatted')
	  open(84,file='gnuh.dat',status='unknown',form='formatted')
	  open(85,file='gtke.dat',status='unknown',form='formatted')
	  open(86,file='geps.dat',status='unknown',form='formatted')
	  open(87,file='glen.dat',status='unknown',form='formatted')
	end if
	
c	write(86,*) it,ndim
c	write(88,*) it,ndim
c	do i=1,ndim
c	  k = nodes(i)
c	  nlev = ilhkv(k)
c
c	  write(86,*) i,k,nlev
c	  do l=1,nlev
c	    write(86,*) l,ulnv(l,k),vlnv(l,k)
c	  end do
c
c	  write(88,*) i,k,nlev
c	  do l=1,nlev-1
c	    write(88,*) l,visv(l,k),difv(l,k),tken(l,k),eps(l,k),rls(l,k)
c	  end do
c
c	end do

	k = nodes(2)
	nlev = ilhkv(k)
	nlev1 = nlev - 1

	call gotm_get(k,nlev,num,nuh,tk,ep,rl)

	write(81,*) it,nlev,1.
	write(81,*) (ulnv(l,k),l=1,nlev)

	write(82,*) it,nlev,1.
	write(82,*) (rhov(l,k),l=1,nlev)

	write(89,*) it,nlev1,1.
	write(89,*) (num(l),l=1,nlev1)

	write(84,*) it,nlev1,1.
	write(84,*) (nuh(l),l=1,nlev1)

	write(85,*) it,nlev1,1.
	write(85,*) (tk(l),l=1,nlev1)

	write(86,*) it,nlev1,1.
	write(86,*) (ep(l),l=1,nlev1)

	write(87,*) it,nlev1,1.
	write(87,*) (rl(l),l=1,nlev1)

	end

c*****************************************************************

	subroutine sedcust(it)

c channel output

	implicit none

	include 'param.h'

	integer it

	real unv(1)
	common /unv/unv
	real vnv(1)
	common /vnv/vnv

	integer ie,i
	integer iunit
	real u,v,h
	integer ieint
	real depele

	integer ndim
	parameter(ndim=5)
	integer iesp,iisp
	save iesp,iisp
	integer iespv(ndim),iispv(ndim)
	save iespv,iispv
	integer icall
	save icall

	data iesp / 7106 /
	data iespv / 100 , 6941 , 1989 , 4284 , 6237 /
	data icall / 0 /

	if( icall .eq. 0 ) then
	  iisp = ieint(iesp)
c	  write(78,*) 'element: ',iisp,iesp
	  do i=1,ndim
	    iispv(i) = ieint(iespv(i))
	  end do
	  icall = 1
	end if

	ie = iisp
	u = unv(ie)
	v = vnv(ie)
	h = depele(ie,+1)
	write(78,*) it,u,v,h

	do i=1,ndim
	  ie = iispv(i)
	  u = unv(ie)
	  v = vnv(ie)
	  h = depele(ie,+1)
	  iunit = 80 + i
	  write(iunit,*) it,u,v,h
	end do

	end

c*****************************************************************

	subroutine velchan(it)

c channel velocity output

	implicit none

	include 'param.h'

	integer it

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer ilhkv(1)
        common /ilhkv/ilhkv
        real uprv(nlvdim,1), vprv(nlvdim,1)
	common /uprv/uprv, /vprv/vprv
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real tempv(nlvdim,1)
	common /tempv/tempv
	real znv(1)
	common /znv/znv

	logical berror
	logical b1,b2,b3,b4
	integer k,i,l,nlev
	integer ivert

	integer ndim,ndim1
	parameter(ndim=5,ndim1=1)

	integer k1(ndim),k2(ndim),k3(ndim1),k4(ndim1)
	save k1,k2,k3,k4

	integer icall
	save icall

	data icall / 0 /
	data k1 / 45,34,23,12,1 /
	!data k2 / 11,22,33,44,55 /
	data k2 / 203 , 228 , 253 , 278 , 303 /
	!data k3 / 28 /
	data k3 / 253 /
	data k4 / 558 /

	b1 = .false.
	b2 = .true.
	b3 = .true.
	b4 = .false.

	if( icall .eq. 0 ) then
	  call n2int(ndim,k1,berror)
	  if( berror ) stop 'error stop velchan: k1'
	  call n2int(ndim,k2,berror)
	  if( berror ) stop 'error stop velchan: k2'
	  call n2int(ndim1,k3,berror)
	  if( berror ) stop 'error stop velchan: k3'
	  icall = 1

	  if( b3 ) then
	    k = k3(1)
	    nlev = ilhkv(k)
	    ivert = 1
	    call animhead(94,nlev,ivert,0.,float(nlev),-2.,2.)	!vel
	    call animhead(96,nlev,ivert,0.,float(nlev),-5.,35.)	!temp
	  end if
	end if

	call nantest(nkn*nlvdim,uprv,'uprv')
	call nantest(nkn*nlvdim,vprv,'vprv')
	call nantest(nkn*nlvdim,tempv,'tempv')

	if( b1 ) then
	write(92,*) it,ndim
	do i=1,ndim
	  k = k1(i)
	  nlev = ilhkv(k)
	  write(92,*) k,nlev,znv(k)
	  write(92,*) (uprv(l,k),l=1,nlev)
	  write(92,*) (tempv(l,k),l=1,nlev)
	end do
	end if

	if( b2 ) then
	write(93,*) it,ndim
	do i=1,ndim
	  k = k2(i)
	  nlev = ilhkv(k)
	  write(93,*) k,nlev,znv(k)
	  write(93,*) (uprv(l,k),l=1,nlev)
	  write(93,*) (tempv(l,k),l=1,nlev)
	end do
	end if

	if( b3 ) then
	k = k3(1)
	nlev = ilhkv(k)
	call animdata(94,it,nlev,uprv(1,k))
	call animdata(96,it,nlev,tempv(1,k))
	!write(94,*) it,nlev
	!write(94,*) (uprv(l,k),l=1,nlev)
	!write(96,*) it,nlev,1.
	!write(96,*) (tempv(l,k),l=1,nlev)
	end if

	if( b4 ) then
	k = k4(1)
	nlev = ilhkv(k)
	write(94,*) it,nlev,1.
	write(94,*) (uprv(l,k),l=1,nlev)
	write(96,*) it,nlev,1.
	write(96,*) (tempv(l,k),l=1,nlev)
	end if

	end

c*****************************************************************

	subroutine concmass(it)

c set up zeta and conc for mass conservation test

	implicit none

	include 'param.h'

	integer it

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer ilhkv(1)
        common /ilhkv/ilhkv
        integer nen3v(3,1)
        common /nen3v/nen3v
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real cnv(nlvdim,nkndim)
	common /cnv/cnv
	real znv(1)
	common /znv/znv
	real zenv(3,1)
	common /zenv/zenv

	integer k,ie,ii,n,ip,l
	integer ntot,nlev
	integer ibase,iw
	real y0,dy,dz,y,z
	real conz
	real res

	integer icall
	save icall
	integer kn(5)
	save kn
	data icall / 0 /
	data kn / 3186,2729,2507,2371,2180 /

	if( icall .eq. 0 ) then

	y0 = 28000.
	dy = 25000.
	dz = 0.5
	ntot = 0

	do k=1,nkn
	  y = ygv(k)
	  z = dz*(y-y0)/dy
	  znv(k) = z
	  nlev = ilhkv(k)
	  if( z .gt. 0. ) then
	    conz = 100.
	    ntot = ntot + 1
	  else
	    conz = 0.
	  end if
	  do l=1,nlev
	    cnv(l,k) = conz
	  end do
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

	call setweg(0,iw)

	icall = 1

	write(6,*) 'conz and level initialized... ',ntot,nkn

	end if

	call massconc(1,cnv,nlvdim,res)
	write(6,*) 'total dissolved mass: ',res
	write(6,*) 'values: ',(cnv(1,kn(ii)),ii=1,5)

	end

c*****************************************************************

	subroutine concmass1(it)

c mass conservation test

	implicit none

	integer it

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer iwegv(1)
	common /iwegv/iwegv
	real cnv(nlvdim,nkndim)
	common /cnv/cnv

	real res
	integer i,ie,iweg

	iweg = 0
	do ie=1,nel
	  if(iwegv(ie) .ne. 0 ) iweg = iweg + 1
	end do

	call massconc(1,cnv,nlvdim,res)		!3D

	write(6,*) 'total dissolved mass: ',it,res,iweg
c	call debug_node(4179)
c	call debug_dry

	end

c*****************************************************************

	subroutine hakata(it)

	implicit none

        integer it

        integer ifemopa

        integer iunit
        save iunit
        data iunit / 0 /

        if( iunit .le. 0 ) then
          iunit = ifemopa('Opening Hakata file','.hak','f','new')
        end if

        call hakanode(iunit,it,1608)    !901
c        call hakanode(iunit,it,777)     !178
        call hakanode(iunit,it,585)     !1032 C-10
        call hakanode(iunit,it,383)     !231

	end

c*****************************************************************

	subroutine hakanode(iunit,it,k)

	implicit none

        integer iunit
        integer it
        integer k

        include 'param.h'

        integer ilhkv(1)
        common /ilhkv/ilhkv
        real tempv(nlvdim,1)
        common /tempv/tempv
        real saltv(nlvdim,1)
        common /saltv/saltv
        real uprv(nlvdim,1),vprv(nlvdim,1)
        common /uprv/uprv, /vprv/vprv
        real visv(0:nlvdim,1)
        real difv(0:nlvdim,1)
        common /visv/visv, /difv/difv

        integer nlev

        nlev = ilhkv(k)

        write(iunit,1000) k,'temp'
	call animdata(iunit,it,nlev,tempv(1,k))
        write(iunit,1000) k,'salt'
	call animdata(iunit,it,nlev,saltv(1,k))
        write(iunit,1000) k,'u'
	call animdata(iunit,it,nlev,uprv(1,k))
        write(iunit,1000) k,'v'
	call animdata(iunit,it,nlev,vprv(1,k))
        write(iunit,1000) k,'vis'
	call animdata(iunit,it,nlev+1,visv(0,k))
        write(iunit,1000) k,'dif'
	call animdata(iunit,it,nlev+1,difv(0,k))

        return
 1000   format(i8,2x,a)
	end

c*****************************************************************
c*****************************************************************

	subroutine animhead(iunit,lmax,ivert,xmin,xmax,ymin,ymax)

	implicit none

	integer iunit
	integer lmax,ivert
	real xmin,xmax,ymin,ymax

	integer magic,nvers
	parameter( magic = 46728645 , nvers = 1 )

	write(iunit,*) magic,nvers
	write(iunit,*) lmax,ivert
	write(iunit,*) xmin,xmax,ymin,ymax

	end

c*****************************************************************

	subroutine animdata(iunit,it,n,val)

	implicit none

	integer iunit
	integer it,n
	real val(1)

	integer i

	write(iunit,*) it,n
	write(iunit,*) (val(i),i=1,n)

	end

c*****************************************************************

	subroutine close_inlets

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
c local
        integer ie,ii,k
c	integer ibc,nbc
	integer ibc
	integer ibctot
	real dz,dzcmh
	real zdate
	real dt
	integer nbnds,itybnd
	integer ifemopa

	integer iunit
	save iunit
	data iunit / 0 /

	if( iunit .eq. 0 ) then
	  iunit = ifemopa('','.ccc','form','new')
	end if

c	nbc = nbnds()
	dz = 0.
	dzcmh = 0.
	call get_timestep(dt)

	ibctot=0
	do ibc=1,nbc
	  if( itybnd(ibc) .lt. 0 ) ibctot = ibctot + 1
	end do

	if( ibctot .eq. nbc ) then	!lagoon closed

	    dzcmh = 2.
	    dz = dt * dzcmh / ( 100. * 3600. )

	    call raise_zeta(dz)

	end if

	if( it .eq. 86400 ) call set_zdate(0,0.40)
	call get_zdate(0,zdate)
	write(iunit,*) it,dzcmh,zdate,(itybnd(ibc),ibc=1,nbc)
	
	end

c*****************************************************************

	subroutine close_inlets1

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
	real znv(1)
	common /znv/znv
c local
        integer ie,ii,k
c	integer ibc,nbc
	integer ibc
	integer ibctot
	real dz,dzmmh
	real zdate
	real dt
	integer nbnds,itybnd
	integer ifemopa
	real getpar

	integer ndim
	parameter(ndim=200)
	real vals(ndim)
	integer icl,nbcc,iclose,n,ih,iclass,itcl
	integer itotal,ittot,iclose_old
	real psv,dsv,zrise,zextra,zsalv,zbound,zfranco
	real zmax,wmax,rmax
	real windv,rainv,zlevel
	real t
	character*4 class
	real zvbnds

	save itotal,ittot,iclose_old

	integer iunit,i66
	save iunit,i66
	data iunit / 0 /

        integer icm
        icm(t) = nint(100.*t)

	if( iunit .eq. 0 ) then
	  iunit = ifemopa('','.ccc','form','new')
	  i66 = ifemopa('','.cc1','form','new')
	  itotal = 0
	  ittot = 0
	  iclose_old = 0
	end if

c dsl  pts  btz  cfg  pbo  vgr  tso  fus  pov  ser  tre  gbo  vdg  lsl
c 2750 2449 19   919  533  1385 936  2174 2032 3140 3342 4086 4343 3971

c 1336 1717

	k = 1336
	dsv = znv(k)
	k = 1717
	psv = znv(k)

	zlevel = psv
	call get_timestep(dt)

	zrise   = 0.01*getpar('zrise')
	zfranco = 0.01*getpar('zfranc')
	zsalv   = 0.01*getpar('zsalv')
        zextra  = zsalv - 1.

c iclose=1      => lagoon is completely closed

	call is_closed(0,nbcc,icl)      !nbcc is total number of OBC
	iclose = 0
	if( icl .eq. nbcc ) iclose = 1

          call get_prev(ndim,it,zfranco,n,vals)
          call get_wind_rain(it,windv,rainv,wmax,rmax,dzmmh)

          call class12new(it,n,vals,zsalv,zextra,zrise,wmax,rmax,psv
     +                  ,zmax,zdate,ih,iclass,class)

          call set_zdate(0,zdate)

          if( iclose .eq. 1 ) then      !already closed
            dz = dzmmh*dt/(1000.*3600.)   !convert to [m/timestep]
            call raise_zeta(dz)
          end if

	  zlevel = 0
	  if( ih .gt. 0 ) zlevel = zrise+vals(ih)       !level at time ih
	  zbound = zvbnds(1)                            !level outside
          write(i66,2300) it,iclass,class,icl,ih
     +          ,icm(zmax),icm(zlevel),icm(psv),icm(dsv)
     +		,icm(zbound),icm(zdate),(itybnd(ibc),ibc=1,nbc)

c ittot         number of time steps inlets are completely closed
c itotal        number of total closures

          if( iclose .eq. 1 ) then
            ittot = ittot + 1
	    if( iclose_old .eq. 0 ) itotal = itotal + 1
          end if

	iclose_old = iclose

	write(iunit,*) it,dzmmh,zdate
     +			,(itybnd(ibc),ibc=1,nbc),itotal,ittot
	
 2200       format(20i4)
 2100     format(i10,i2,a5,i2,i3,4f6.2)
 2300     format(i10,i2,a5,i2,i3,3x,10i4)
 2301     format(i10,i2,a5,i2,i3,3x,10i4,2x,3i4)
 2000   format(i10,i4,7f9.2)
 2500   format(a,a4,i10,2f9.2,i10)

	end

c*****************************************************************

	subroutine test_hakata(it)

	implicit none

        integer it

	logical berror
	integer i,k
        integer ifemopa

	integer ndim
	parameter(ndim=4)
	integer nodes(ndim)
        integer iunit
        save iunit
        save nodes
        data iunit / 0 /
        data nodes / 901, 178, 1032, 231 /

        if( iunit .le. 0 ) then
          iunit = ifemopa('Opening out3d','.o3d','f','new')
	  write(6,*) 'Nodes from test_hakata: '
          call n2int(ndim,nodes,berror)
	  do i=1,ndim
	    write(6,*) i,nodes(i)
	  end do
	  if( berror ) stop 'error stop test_hakata'
        end if

	do i=1,ndim
	  k = nodes(i)
	  call node3d(iunit,it,k)
	end do

c        call hakanode(iunit,it,1608)    !901
c        call hakanode(iunit,it,777)     !178
c        call hakanode(iunit,it,585)     !1032 C-10
c        call hakanode(iunit,it,383)     !231

	end

c*****************************************************************

	subroutine node3d(iunit,it,k)

	implicit none

        integer iunit
        integer it
        integer k

        include 'param.h'

        integer ilhkv(1)
        common /ilhkv/ilhkv
        real tempv(nlvdim,1)
        common /tempv/tempv
        real saltv(nlvdim,1)
        common /saltv/saltv
        real uprv(nlvdim,1),vprv(nlvdim,1)
        common /uprv/uprv, /vprv/vprv
	real znv(1)
	common /znv/znv
	real wlnv(0:nlvdim,1)
	common /wlnv/wlnv

        integer nlev
	integer l

        nlev = ilhkv(k)

	write(iunit,1000) 'header (z): ',it,k
	write(iunit,*) 1,znv(k)

	write(iunit,1000) 'header (u): ',it,k
	write(iunit,*) nlev,(uprv(l,k),l=1,nlev)
	write(iunit,1000) 'header (v): ',it,k
	write(iunit,*) nlev,(vprv(l,k),l=1,nlev)
	write(iunit,1000) 'header (w): ',it,k
	write(iunit,*) nlev+1,(wlnv(l,k),l=0,nlev)

	write(iunit,1000) 'header (temp): ',it,k
	write(iunit,*) nlev,(tempv(l,k),l=1,nlev)
	write(iunit,1000) 'header (salt): ',it,k
	write(iunit,*) nlev,(saltv(l,k),l=1,nlev)

        return
 1000   format(a,2i10)
	end

c*****************************************************************

	subroutine lago(it)

	implicit none

	integer it

	include 'param.h'

        real ulnv(nlvdim,neldim)
        common /ulnv/ulnv
        real vlnv(nlvdim,neldim)
        common /vlnv/vlnv
        integer ilhkv(1)
        common /ilhkv/ilhkv
        real uprv(nlvdim,1), vprv(nlvdim,1)
	common /uprv/uprv, /vprv/vprv

	logical berror
	integer i,k
	real u,v
 
	integer ndim
	parameter (ndim=3)
	real vel(ndim)
	integer icall
	integer nodes(ndim)
	save nodes,icall
	data nodes /1020,2021,2527/
	data icall /0/

	if( icall .eq. 0 ) then
	  icall = 1
	  call n2int(ndim,nodes,berror)
	  if( berror ) stop 'error stop'
	end if

	do i=1,ndim
	  k = nodes(i)
	  u=uprv(1,k)
	  v=vprv(1,k)
	  !write(8,*) i,k,u,v
	  vel(i) = sqrt(u*u+v*v)
	end do
	  
	write(8,*) it,(vel(i),i=1,ndim)
	write(9,*) (vel(i),i=1,ndim)

	end

c*****************************************************************

	subroutine aldo(it)

	implicit none

	integer it

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(1)
        common /ilhkv/ilhkv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer l,k,lmax
        real x0,y0,x1,y1,x2,y2,x3,y3
        real xp,yp
        logical is_left

        x0 = 2427.180420
        y0 = 2443.140137
        x1 = 4201.409180
        y1 = 235.025009
        x2 = 4201.409180
        y2 = 235.025009
        x3 = 5664.500000
        y3 = 1452.153198

        if( it .eq. 86400-300 ) then
         do k=1,nkn
          xp = xgv(k)
          yp = ygv(k)
          if( .not. is_left(xp,yp,x0,y0,x1,y1) .and.
     +                  is_left(xp,yp,x2,y2,x3,y3) ) then
            lmax = ilhkv(k)
            do l=1,lmax
              cnv(l,k) = 100.
            end do
          end if
         end do
        end if

        end

c*****************************************************************

        function is_left(xp,yp,x0,y0,x1,y1)

        implicit none

        logical is_left
        real xp,yp,x0,y0,x1,y1
        real dx,dy,dxn,dyn,dxp,dyp
        real scal

        dx = x1-x0
        dy = y1-y0
        dxn = -dy
        dyn = dx

        dxp = xp - x0
        dyp = yp - y0

        scal = dxp*dxn + dyp*dyn

        is_left = scal .gt. 0.

        end

c*****************************************************************
c
c check following routines:
c
c	get_next_traccia  write_traccia  interpolate_traccia
c
c check following logical variables:
c
c usedts	time is in YYMMDD etc.. and must be converted to it
c			-> change in get_next_traccia, write_traccia
c bnearpoint	if true find closest point and use this value
c bintmiss	if true use last good value
c		use (xinit,yinit) if no good value has been found yet
c		in this case also (xinit,yinit) have to be set
c bnodry	treat dry areas as not found elements
c
c if one of bnearpoint or bintmiss is .true. a value will be always found
c if both are .false. points outside of the domain return -999.
c
c fantina:	bnearpoint = .true.   bintmiss = .true.  usedts=.true.
c rachel:	bnearpoint = .false.  bintmiss = .true.  usedts=.false.
c
c for carote di fantina, use usedts=.false.
c
c if there are traccie completely out of domain,
c	use bnearpoint = .false.  bintmiss = .true.
c
c*****************************************************************

	subroutine traccia

	implicit none

	include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        real zeov(3,1)
	common /zeov/zeov
        real zenv(3,1)
	common /zenv/zenv
	integer iwegv(1)
	common /iwegv/iwegv

	integer itnew,itold,ierr
	real zp,zold,znew
	integer ifileo
	integer year,date,time
	character*80 line

	integer iuin,iuout,itp,np,iep
	double precision xdp,ydp
	real xp,yp
	save iuin,iuout,itp,np,iep
	save xp,yp
	save xdp,ydp
	save line

	integer icall
	save icall
	data icall /0/

	if( icall .eq. -1 ) return

c-------------------------------------------------------------------
c initialize
c-------------------------------------------------------------------

	if( icall .eq. 0 ) then
	  iuin = ifileo(55,'traccia.dat','form','old')
	  if( iuin .le. 0 ) goto 99
	  iuout = ifileo(55,'traccia_out.dat','form','new')
	  if( iuout .le. 0 ) goto 99

	  line = ' '
	  write(6,*) 'initializing traccia...'

	  call get_next_traccia(iuin,itp,np,xdp,ydp,line,ierr)
	  xp = xdp
	  yp = ydp
	  if( ierr .ne. 0 ) goto 98

	  iep = 0
	  icall = 1
	end if

c-------------------------------------------------------------------
c prepare time
c-------------------------------------------------------------------

	itnew = it
	itold = it - idt

	if( itp .lt. itold ) goto 89
	ierr = 0

c-------------------------------------------------------------------
c loop on traccie read
c-------------------------------------------------------------------

	do while( ierr .eq. 0 .and. itold .le. itp .and. itp .le. itnew )

	  call interpolate_traccia(iep,itp,itold,itnew,xp,yp,zp)

	  write(6,*) 'traccia: ',itold,itnew,itp,xp,yp,zp,iep
	  write(99,*) 'traccia: ',itold,itnew,itp,xp,yp,zp,iep
     +					,iwegv(abs(iep)),np
	  call write_traccia(iuout,itp,np,xdp,ydp,zp,line)
	  
	  call get_next_traccia(iuin,itp,np,xdp,ydp,line,ierr)
	  xp = xdp
	  yp = ydp

	end do

	if( ierr .gt. 0 ) goto 97
	if( ierr .lt. 0 ) icall = -1

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	return
   99	continue
	stop 'error stop traccia: cannot open file'
   98	continue
	stop 'error stop traccia: first record of traccia'
   97	continue
	write(6,*) itp,np,xp,yp,ierr
	stop 'error stop traccia: read error in traccia'
   89	continue
	write(6,*) itp,itold,itnew
	stop 'error stop traccia: simulation start too late'
        end

c*****************************************************************

	subroutine interpolate_traccia(iep,itp,itold,itnew,xp,yp,zp)

	implicit none

	integer iep
	integer itp,itold,itnew
	real xp,yp
	real zp			!interpolated value

	real xinit,yinit	!lido for rachel
	parameter (xinit=38889.,yinit=32745.2)

        real zeov(3,1)
	common /zeov/zeov
        real zenv(3,1)
	common /zenv/zenv
        real zov(1)
	common /zov/zov
        real znv(1)
	common /znv/znv
	integer iwegv(1)
	common /iwegv/iwegv

	logical bintmiss
	logical bnearpoint
	logical bnodry
	real zold,znew
	integer iepold,k
	real xold,yold
	save xold,yold,iepold

	integer get_nearest_point

	integer icall
	save icall
	data icall /0/

c bnearpoint	if true find closest point and use this value
c bintmiss	if true use last good value
c		use (xinit,yinit) if no good value has been found yet
c		in this case also (xinit,yinit) have to be set
c bnodry	treat dry areas as not found elements
c
c if one of bnearpoint or bintmiss is .true. a value will be always found
c if both are .false. points outside of the domain return -999.
c
c fantina:	bnearpoint = .true.   bintmiss = .true.
c rachel:	bnearpoint = .false.  bintmiss = .true.

	bintmiss = .true.
	bnearpoint = .true.
	bnodry = .true.

	if( icall .eq. 0 ) then		!get first element
	  xold = xinit
	  yold = yinit
	  icall = 1
	  call find_elem_from_old(0,xold,yold,iepold)
	end if

	call find_elem_from_old(iep,xp,yp,iep)

	if( bnodry .and. iep .gt. 0 .and. iwegv(iep) .gt. 0 ) iep = -iep

	if( iep .gt. 0 ) then		!element found -> interpolate
	    call femintp(iep,zeov(1,iep),xp,yp,zold)
	    call femintp(iep,zenv(1,iep),xp,yp,znew)
	    zp = zold + (znew-zold)*(itp-itold)/float(itnew-itold)
	else if( bnearpoint ) then	!interpolate using nearest point
	    k = get_nearest_point(xp,yp)
	    zold = zov(k)
	    znew = znv(k)
	    zp = zold + (znew-zold)*(itp-itold)/float(itnew-itold)
	    write(6,*) 'using nearest point: ',k,zp
	else if( bintmiss ) then	!interpolate using old point
	    iep = iepold
	    if( iep .le. 0 ) goto 99
	    xp = xold
	    yp = yold
	    call femintp(iep,zeov(1,iep),xp,yp,zold)
	    call femintp(iep,zenv(1,iep),xp,yp,znew)
	    zp = zold + (znew-zold)*(itp-itold)/float(itnew-itold)
	else
	    zp = -999.
	end if

	iepold = iep
	xold = xp
	yold = yp

	return
   99	continue
	write(6,*) iepold
	write(6,*) xold,yold,xp,yp
	write(6,*) itp,itold,itnew
	stop 'error stop interpolate_traccia: no element found'
	end

c*****************************************************************

	subroutine get_next_traccia(iuin,itp,np,xdp,ydp,line,ierr)

c 2004 07 18 10 0 41 621 37002.0159999998 39108.9670000002
c 23964600  37778.7 39962.7      60 SA14 2005-10-05::08:50:00

	implicit none

	integer iuin,itp,np
	double precision xdp,ydp
	character*80 line
	integer ierr

	integer year,month,day,hour,min,sec
	integer date,time
	logical dts_initialized

	logical usedts
	parameter (usedts=.true.)

	integer icall
	save icall
	data icall / 0 /

	if( usedts ) then
	  read(iuin,*,iostat=ierr) year,month,day,hour,min,sec,np,xdp,ydp
	else
	  read(iuin,*,iostat=ierr) itp,xdp,ydp,np,line
	  write(6,*) 'new traccia read: ',itp,line
	end if

	if( ierr .ne. 0 ) return
	if( .not. usedts ) return

c	------------------------------------------
c	the next part is only executed if we have to convert
c	the time (given as year month...) to fem_time
c	if we alread read fem_time nothing more is to do
c	------------------------------------------

	if( icall .eq. 0 ) then
	  date = 10000*year + 101
	  time = 0
	  call dtsini(date,time)
	  write(6,*) 'get_next_traccia: initialized dts ',date,time
	  icall = 1
	end if

	call dts2it(itp,year,month,day,hour,min,sec)
	write(6,*) 'new traccia read: ',itp,np

	end

c*****************************************************************

	subroutine write_traccia(iuout,itp,np,xdp,ydp,zp,line)

c writes raw traccia

	implicit none

	integer iuout,itp,np
	double precision xdp,ydp
	real zp
	character*80 line

	integer year,month,day,hour,min,sec
	integer ilen

	logical usedts
	parameter (usedts=.true.)

	if( usedts ) then
	  call dts2dt(itp,year,month,day,hour,min,sec)
	  !write(iuout,*) year,month,day,hour,min,sec,np,xdp,ydp,zp
	  write(iuout,1000) year,month,day,hour,min,sec,np,xdp,ydp,zp
	else
	  call get_string_length(line,ilen)
	  !write(iuout,*) itp,xdp,ydp,np,zp,' ',line(1:ilen)
	  write(iuout,2000) itp,xdp,ydp,np,zp,' ',line(1:ilen)
	end if

	return
 1000	format(i5,5i3,i10,3f14.6)
 2000	format(i10,2f14.6,i10,f14.6,a,a)
	end

c*****************************************************************

	function get_nearest_point(xp,yp)

	implicit none

	integer get_nearest_point
	real xp,yp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(1)
	common /xgv/xgv
	real ygv(1)
	common /ygv/ygv

	integer knear,k
	real dist,d
	real dx,dy

	dist=1.e+30
	knear = 0

	do k=1,nkn
	  dx = xgv(k) - xp
	  dy = ygv(k) - yp
	  d = sqrt( dx*dx + dy*dy )
	  if( d .lt. dist ) then
	    dist = d
	    knear = k
	  end if
	end do

	get_nearest_point = knear
	
	end

c*****************************************************************

	subroutine get_string_length(line,ilen)

	implicit none

	character*80 line
	integer ilen

	integer imax,i

	imax = len(line)
	do i=imax,1,-1
	  if( line(i:i) .ne. ' ' ) goto 1
	end do
    1	continue
	ilen = i
	if( ilen .le. 0 ) ilen = 1

	end
	
c*****************************************************************

	subroutine write_traccia0(iuout,itp,np,xdp,ydp,zp,tz)

c old routine -> do not use anymore

	implicit none

	integer iuout,itp,np
	integer time,ihour
	double precision xdp,ydp
	real zp
	integer tz

	double precision x0,y0
	parameter(x0=2330000. - 50000., y0=5000000.)
	integer itype
	parameter(itype=1)	!0=raw  1=fantina  2=rachel

	integer year,month,day,hour,min,sec
	double precision x,y,z

	call dts2dt(itp,year,month,day,hour,min,sec)
	ihour = hour - tz	!GMT to MET

	if( itype .eq. 0 ) then
	  x = xdp
	  y = ydp
	  z = zp
	  write(iuout,*) year,month,day,hour,min,sec,np,xdp,ydp,itp,zp
	else if( itype .eq. 1 ) then
	  x = xdp + x0
	  y = ydp + y0
	  z = 0.23 - zp
	  time = ihour*10000+min*100+sec
	  write(iuout,'(i8,i8,f14.3,2x,f14.3,2x,f10.3)') time,np,x,y,z
	else if( itype .eq. 2 ) then
c	  (rachel) 2319029.6 5032260.39 38835 27 may 2004
	  x = xdp + x0
	  y = ydp + y0
	  z = zp - 0.23
	  time = ihour*3600+min*60+sec
	  write(iuout,*) x,y,time,day,month,year,z
	else 
	  stop 'error stop write_traccia: invalid itype'
	end if

	write(89,'(6f13.3)') xdp,ydp,x0,y0,x,y

	end

c*****************************************************************

        subroutine jamal(it)

c computes residence time online - one value for whole lagoon

        implicit none

        integer it

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhkv(1)
        common /ilhkv/ilhkv
        integer iarv(1)
        common /iarv/iarv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        real v1v(1)
        common /v1v/v1v

        integer ie,ii,k,lmax,l,ia
        logical bnoret,breset,bstir
        real vol,conz,perc,percmin
        double precision mass,volume
        double precision massaux,volaux
        real volnode
	integer ifemop
	integer iaout,itmin,idtreset
	real dt,c0
	real restime,restime1,restimec
	real remnant,rlast
	real resmed,resstd

	integer iu,it0,ndata
	save iu,it0,ndata
	double precision remint,remlog,remtim
	save remint,remlog,remtim
	double precision rsum,rsumsq
	save rsum,rsumsq
        double precision mass0
        save mass0

        integer icall
        save icall
        data icall / 0 /

c------------------------------------------------------------
c parameters
c------------------------------------------------------------
c
c bnoret	true if no return flow is used (conzentrations outside
c		are explicitly set to 0)
c bstir		simulates completely stirred tank
c		(replaces at every time step conz with average conz)
c percmin	percentage to reach -> after this stop computation
c		use 0 if no premature end is desired
c iaout		area code of elements out of lagoon (used for init and retflow)
c		use -1 to if no outside areas exist
c c0		initial concentration of tracer
c itmin		time from when to compute residence time
c idtreset	time step to reset concentration to c0
c		use 0 if no reset is wanted
c
c--------------------------
c default settings
c--------------------------
        bnoret = .false.
	bstir = .true.
	bstir = .false.
	percmin = 0.
	iaout = -1
	c0 = 1.
	itmin = 0
	idtreset = 0
	idtreset = nint( 3 * 30.5 * 86400 )		!one month is 30.5 days
c--------------------------
c nador
c--------------------------
c	bnoret = .true.
c	percmin = 1.
c	iaout = 0
c--------------------------
c alimini
c--------------------------
c	idtreset = nint(30.5*86400)
c--------------------------
c mar menor
c--------------------------
c 	just use default settings
c--------------------------
c
c------------------------------------------------------------
c do not change anything after this point
c------------------------------------------------------------

c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

        if( it .le. itmin ) return

c------------------------------------------------------------
c initialization -> decide on reset
c------------------------------------------------------------

	breset = .false.		!normally do not reset

        if( icall .eq. 0 ) then
          write(6,*) 'initialization of routine jamal'
	  iu = ifemop('.jam','formatted','new')
	  breset = .true.		!always reset at first call
	  it0 = it
        end if

	if( idtreset .gt. 0 ) then
	  if( it-it0 .ge. idtreset ) breset = .true.
	end if

c------------------------------------------------------------
c flag nodes that are inside lagoon (v1v(k)=1)
c------------------------------------------------------------

        do k=1,nkn
          v1v(k) = 0.
        end do

        do ie=1,nel
          ia = iarv(ie)
          if( ia .ne. iaout ) then
              do ii=1,3
                k = nen3v(ii,ie)
                v1v(k) = 1.
              end do
          end if
        end do

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0
          write(6,*) 'resetting concentrations in jamal at time ',it
          do k=1,nkn
            lmax = ilhkv(k)
	    conz = 0.
            if( v1v(k) .ne. 0. ) conz = c0
            do l=1,lmax
              cnv(l,k) = conz
            end do
          end do
	end if

c------------------------------------------------------------
c total mass (only for nodes inside lagoon)
c------------------------------------------------------------

        mass = 0.
        volume = 0.
        do k=1,nkn
          if( v1v(k) .ne. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              mass = mass + vol*conz
              volume = volume + vol
            end do
          end if
        end do

c------------------------------------------------------------
c stirred tank?
c------------------------------------------------------------

	if( bstir ) then
	  conz = mass / volume
          do k=1,nkn
            if( v1v(k) .ne. 0. ) then
              lmax = ilhkv(k)
              do l=1,lmax
                cnv(l,k) = conz
	      end do
	    end if
	  end do
	end if

        massaux = 0.
        volaux = 0.
        do k=1,nkn
          if( v1v(k) .ne. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              massaux = massaux + vol*conz
              volaux = volaux + vol
            end do
          end if
        end do

	!write(67,*) it,mass,massaux,volume,volaux
	!write(67,*) it,mass,massaux

c------------------------------------------------------------
c------------------------------------------------------------
c reset variables to compute residence time
c------------------------------------------------------------

        if( breset ) then
	  mass0 = mass 
	  remint = 0.
	  remlog = 0.
	  remtim = 0.
	  it0 = it
	  ndata = 0
	  rsum = 0.
	  rsumsq = 0.
	end if

c------------------------------------------------------------
c write to file
c------------------------------------------------------------

	call get_timestep(dt)

	remnant = 0.
        if( mass0 .gt. 0. ) remnant = mass/mass0
        perc = 100.*remnant

	remint = remint + remnant*dt	!integrated remnant function
	restime = remint/86400.		!residence time in days

	rlast = remnant
	if( rlast .ge. 1. ) rlast = 0.
	restimec = restime/(1.-rlast)	!corrected residence time

	remlog = remlog - log(remnant)
	remtim = remtim + (it-it0)
	restime1 = 0.
	if( remlog .gt. 0. ) restime1 = ( remtim / remlog ) / 86400.

	ndata = ndata + 1
	rsum = rsum + restime1
	rsumsq = rsumsq + restime1*restime1
	resmed = rsum / ndata
	resstd = sqrt( rsumsq/ndata - resmed*resmed )

c perc		percentage of mass still in domain
c restime	residence time computed by integrating
c restimec	residence time computed by integrating with correction
c restime1	residence time computed by fitting regression curve
c resmed	average of residence times computed
c resstd	standard deviation of residence time

        write(iu,1000) it,perc,restime,restimec,restime1,resmed,resstd

c------------------------------------------------------------
c finish computation if mass is below threshold
c------------------------------------------------------------

        if( mass0 .ne. 0. .and. perc .lt. percmin ) then
                stop 'finished computing'
        end if

c------------------------------------------------------------
c no return flow -> set outside areas to 0
c------------------------------------------------------------

        if( bnoret ) then
          do k=1,nkn
            if( v1v(k) .eq. 0. ) then
              lmax = ilhkv(k)
              do l=1,lmax
                cnv(l,k) = 0.
              end do
            end if
          end do
        end if

c------------------------------------------------------------
c remember initialization
c------------------------------------------------------------

        icall = icall + 1

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	return
 1000	format(i10,6f10.2)
        end

c*****************************************************************

        subroutine jamal_fra

c reset conz for fra

        implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhkv(1)
        common /ilhkv/ilhkv
        integer iarv(1)
        common /iarv/iarv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        real v1v(1)
        common /v1v/v1v

        integer ie,ii,k,lmax,l,ia
        logical bnoret,breset,bstir
        real vol,conz,perc,percmin
        double precision mass,volume
        double precision massaux,volaux
        real volnode
	integer ifemop
	integer iaout,itmin,idtreset
	real dt,c0
	real restime,restime1,restimec
	real remnant,rlast
	real resmed,resstd

	integer iu,it0,ndata
	save iu,it0,ndata
	double precision remint,remlog,remtim
	save remint,remlog,remtim
	double precision rsum,rsumsq
	save rsum,rsumsq
        double precision mass0
        save mass0

        integer icall
        save icall
        data icall / 0 /

	itmin = 0
	idtreset = 0
	idtreset = nint( 1 * 30.5 * 86400 )		!one month is 30.5 days

c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

        if( it .le. itmin ) return

c------------------------------------------------------------
c initialization -> decide on reset
c------------------------------------------------------------

	breset = .false.		!normally do not reset

        if( icall .eq. 0 ) then
          write(6,*) 'initialization of routine jamal'
	  breset = .true.		!always reset at first call
	  it0 = it
	  it0 = itanf
        end if

	if( idtreset .gt. 0 ) then
	  if( it-it0 .ge. idtreset ) breset = .true.
	end if

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0
          write(6,*) 'resetting concentrations in jamal at time ',it
          write(87,*) 'resetting concentrations in jamal at time ',it
          do k=1,nkn
            lmax = ilhkv(k)
            do l=1,lmax
              cnv(l,k) = 0
            end do
          end do
	  it0 = it
	end if

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c*****************************************************************

        subroutine sedimt

        implicit none

        include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhkv(1)
        common /ilhkv/ilhkv
        integer iarv(1)
        common /iarv/iarv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        real v1v(1)
        common /v1v/v1v

	real conzs(nkndim)
	save conzs
	real conza(nkndim)
	save conza
	real conzh(nkndim)
	save conzh

        integer ie,ii,k,lmax,l,ia
	integer iunit
        logical bnoret
        real vol,conz,perc,wsink,dt,sed,h,r,cnew,rhos
        double precision mass,masss
        real volnode,depnode
	real getpar

	integer iu,id,itmcon,idtcon,itstart
	save iu,id,itmcon,idtcon,itstart

        integer icall
        save icall
        data icall / 0 /

c------------------------------------------------------------
c parameters
c------------------------------------------------------------

        bnoret = .false.
	wsink = 0.
	wsink = 1.e-4
	wsink = 1.e-5
	wsink = 5.e-5
	rhos = 2500.
	call get_timestep(dt)
	call getinfo(iunit)

c------------------------------------------------------------
c initialization
c------------------------------------------------------------

        if( icall .eq. 0 ) then

          write(6,*) 'initialization of routine sedimt: ',wsink

	  do k=1,nkn
	    conzs(k) = 0.
	    conza(k) = 0.
	    conzh(k) = 0.
            lmax = ilhkv(k)
            do l=1,lmax
              cnv(l,k) = 0.
            end do
	  end do

	  itstart = nint(getpar('tcust'))

          iu = 55
          itmcon = nint(getpar('itmcon'))
          idtcon = nint(getpar('idtcon'))
          call confop(iu,itmcon,idtcon,1,3,'set')

          icall = 1

        end if

c------------------------------------------------------------
c is it time ?
c------------------------------------------------------------

        if( it .lt. itstart ) return

c------------------------------------------------------------
c sinking
c------------------------------------------------------------

	if( wsink .gt. 0. ) then
	  l = 1
          do k=1,nkn
              h = depnode(l,k,+1)
              vol = volnode(l,k,+1)
	      r = 0.
	      if( h .gt. 0. ) r = wsink/h
              conz = cnv(l,k)
              conz = max(0.,conz)
	      cnew = conz * exp(-r*dt)
              cnv(l,k) = cnew
              sed = vol * (conz-cnew) 
	      conzs(k) = conzs(k) + sed
              sed = h * (conz-cnew) 
	      conza(k) = conza(k) + sed
              sed = sed / rhos
	      conzh(k) = conzh(k) + sed
          end do
	end if

c------------------------------------------------------------
c total mass
c------------------------------------------------------------

        do k=1,nkn
          v1v(k) = 0.
        end do

        do ie=1,nel
          ia = iarv(ie)
          if( ia .ne. 0 ) then
              do ii=1,3
                k = nen3v(ii,ie)
                v1v(k) = 1.
              end do
          end if
        end do

        mass = 0.
        masss = 0.
        do k=1,nkn
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              mass = mass + vol*conz
            end do
	    masss = masss + conzs(k)
        end do

c------------------------------------------------------------
c write total mass
c------------------------------------------------------------

        write(6,*) 'sedimt: ',it,mass,masss,mass+masss
        write(iunit,*) 'sedimt: ',it,mass,masss,mass+masss

        id = 22       !for sediment -> [kg]
	call confil(iu,itmcon,idtcon,id,1,conzs)
        id = 23       !for sediment -> [kg/m**2]
	call confil(iu,itmcon,idtcon,id,1,conza)
        id = 24       !for sediment -> [m]
	call confil(iu,itmcon,idtcon,id,1,conzh)

c------------------------------------------------------------
c no return flow
c------------------------------------------------------------

        if( bnoret ) then

        do k=1,nkn
          if( v1v(k) .eq. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
                cnv(l,k) = 0.
            end do
          end if
        end do

        end if

c------------------------------------------------------------
c remember initialization
c------------------------------------------------------------

        icall = 1

c------------------------------------------------------------
c end of initialization
c------------------------------------------------------------

        end

c*****************************************************************

	subroutine oscillation

	implicit none

	include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real kenerg,penerg,tenerg

	integer icall
	data icall /0/
	save icall

	if( icall .eq. 0 ) then
	  icall = 1
	end if

	call energ(0,kenerg,penerg)

	tenerg = kenerg + penerg
	write(9,*) it,kenerg,penerg,tenerg

	end

c*****************************************************************

	subroutine kreis

	implicit none

	include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real kenerg,penerg,tenerg

	integer icall
	data icall /0/
	save icall

	if( icall .eq. 0 ) then
	  icall = 1
	  call init_kreis()
	end if

	call energ(0,kenerg,penerg)

	tenerg = kenerg + penerg
	write(9,*) it,kenerg,penerg,tenerg

	end

c*****************************************************************

	subroutine init_kreis

	implicit none

	include 'param.h'

	integer k,ie,ii,l
	real pi,dcor,z0,r0,f,omega,grav
	real aux,r02,r2
	real x,y,z,u,v

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real znv(1)
	common /znv/znv
	real zenv(3,1)
	common /zenv/zenv
	integer nen3v(3,1)
	common /nen3v/nen3v
	real ulnv(nlvdim,neldim)
	common /ulnv/ulnv
	real vlnv(nlvdim,neldim)
	common /vlnv/vlnv

        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real areakv(nlvdim,nkndim)
        common /areakv/areakv

	real v1v(1)
	common /v1v/v1v

	pi = 4.*atan(1.)
	dcor = 45.
	z0 = 0.
	r0 = 0.
	f = 2.0 * 0.729E-4 * sin(dcor*pi/180.)
	omega = 0.5E-5
	grav = 9.81

	aux = (omega*f)/(2.*grav)
	r02 = r0*r0

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  r2 = x*x+y*y
	  z = z0 + aux * ( r2 - r02 )

	  znv(k) = z
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

	call setdepth(nlvdim,hdknv,hdenv,zenv,areakv)

	do ie=1,nel
	  call baric(ie,x,y)
	  u = -omega * y
	  v = +omega * x
	  do l=1,nlv
	    ulnv(l,ie) = u
	    vlnv(l,ie) = v
	  end do
	end do

	call vtot
	call uvint
	call uvtop0(v1v)
	call uvtopr(v1v)

	end

c*****************************************************************

c*****************************************************************
        subroutine bclevvar

        implicit none 

        include 'param.h' 

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        integer edim
        parameter(edim = 20) !20!number of element near the open boundary
        integer elese(edim)!elements near the boundary
c  !      data elese / 290,321,326,335,308,2112,324,336,334,325,
c !    +     2118,2113,328,295,2119,2116,333,2121,2120,2114
c!     +      /                                                !grid chao0_primapubb.grd
        data elese / 3002,3001,2902,2901,2802,2801,2702,2701,
     +  2602,2601,3004,3003,2904,2903,2804,2803,2704,2703,2604,2603
     +      /    !grid chao0_ordinato.grd
c     !   data elese /382,381,380,379,378,377,376,375,374,373,372,371,
c     !+  369,368,367,366,365,364,363,362,361,360,359,358,357,356,
c     !+  355,354,353
c     !+         /!grid bati_gradino.grd
        integer eleint(edim)
        integer ie,l,i
        logical berror          !true on return if error
        
        real ulnv(nlvdim,neldim),vlnv(nlvdim,neldim)
        common /ulnv/ulnv, /vlnv/vlnv
        real utlnv(nlvdim,neldim), vtlnv(nlvdim,neldim)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real hdeov(nlvdim,neldim)
        common /hdeov/hdeov

        real u,v
        real ut,vt
        real h

        integer nlev
        parameter (nlev = 5)
      
        real upresc(nlev),vpresc(nlev)
        real umax
        real alpha

        umax = 0.1 !unity m/s
        upresc(1) = umax !vertical profile u
        upresc(2) = umax/2
        upresc(3) = -umax/2!-umax/2
        upresc(4) = -umax/2!-umax/2
        upresc(5) = -umax/2!-umax/2!deb220307
        
        vpresc(1) = 0. !vertical profile v
        vpresc(2) = 0.
        vpresc(3) = 0.
        vpresc(4) = 0.
        vpresc(5) = 0.

        do i=1,edim 
            eleint(i)=elese(i)
        enddo
                call e2int(edim,eleint,berror)
            
        !print*,eleint 
        if (it .le. 43200)then
                alpha=it/43200.
                !write(6,*)'funziona',it,alpha
        else
                alpha=1.
        endif

        do i=1,edim
                ie=eleint(i)
                do l=1,nlev
                       h = hdeov(l,ie)
                       ulnv(l,ie) = upresc(l)*alpha
                       vlnv(l,ie) = vpresc(l)*alpha
                !   write(6,*)'funzionaDEB',ie,l,ulnv(l,ie),vlnv(l,ie)
                       utlnv(l,ie) = ulnv(l,ie)*h
                       vtlnv(l,ie) = vlnv(l,ie)*h
         !  write(6,*) l,' h ',h,' ulnv ',ulnv(l,ie),' vlnv ',vlnv(l,ie)
                   !write(6,*) l,'utlnv',utlnv(l,ie),'vtlnv',vtlnv(l,ie)
                enddo
        enddo

        end
c*****************************************************************
        subroutine bclevvar_ini !deb

        implicit none 

        include 'param.h' 

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer edim
        parameter(edim = 20) !number of element near the open boundary
        integer elese(edim)!elements near the boundary
c!        data elese / 290,321,326,335,308,2112,324,336,334,325,
c!     +     2118,2113,328,295,2119,2116,333,2121,2120,2114
c!     +      /                                                !grid chao0_primapubb.grd
          data elese / 3002,3001,2902,2901,2802,2801,2702,2701,
     +  2602,2601,3004,3003,2904,2903,2804,2803,2704,2703,2604,2603 
     +  /    !grid chao0_ordinato.grd
c    !   data elese /382,381,380,379,378,377,376,375,374,373,372,371,
c    ! +  369,368,367,366,365,364,363,362,361,360,359,358,357,356,
c    ! +  355,354,353
c    ! +     /!grid bati_gradino.grd

        integer eleint(edim)
        integer ie,l,i,ii
        logical berror          !true on return if error
       
        real ulov(nlvdim,neldim),vlov(nlvdim,neldim)
        common /ulov/ulov, /vlov/vlov
        real utlov(nlvdim,neldim), vtlov(nlvdim,neldim)
        common /utlov/utlov, /vtlov/vtlov

        real ulnv(nlvdim,neldim),vlnv(nlvdim,neldim)
        common /ulnv/ulnv, /vlnv/vlnv
        real utlnv(nlvdim,neldim), vtlnv(nlvdim,neldim)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real hdeov(nlvdim,neldim)
        common /hdeov/hdeov
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv


        integer iuvfix(neldim)
        common /iuvfix/iuvfix
        save /iuvfix/

        real u,v
        real ut,vt
        real h

        integer nlev
        parameter (nlev = 5)
      
        real upresc(nlev),vpresc(nlev)
        real umax

        do ii=1,neldim
                iuvfix(ii)= 0
        enddo

        umax = 0.!0.1 !unity m/s
        upresc(1) = umax !vertical profile u
        upresc(2) = umax/2
        upresc(3) = -umax/2
        upresc(4) = -umax/2
        upresc(5) = -umax/2
        
        vpresc(1) = 0. !vertical profile v
        vpresc(2) = 0.
        vpresc(3) = 0.
        vpresc(4) = 0.
        vpresc(5) = 0.

        do i=1,edim 
            eleint(i)=elese(i)
        enddo
                call e2int(edim,eleint,berror)
            
        !print*,eleint 
        do i=1,edim
                ie=eleint(i)
                iuvfix(ie)=1
                do l=1,nlev
                       h = hdeov(l,ie)
                       ulov(l,ie) = upresc(l)!deb 12feb2007
                       vlov(l,ie) = vpresc(l)!deb 12feb2007
                       ulnv(l,ie) = upresc(l)
                       vlnv(l,ie) = vpresc(l)
                       utlov(l,ie) = ulov(l,ie)*h
                       utlnv(l,ie) = ulnv(l,ie)*h
                       vtlov(l,ie) = vlov(l,ie)*h
                       vtlnv(l,ie) = vlnv(l,ie)*h
           !write(6,*)'inizializzo velocita bclevvar_ini'
           !write(6,*) l,' h ',h,' ulnv ',ulov(l,ie),' vlnv ',vlov(l,ie)
           !write(6,*) l,' h ',h,' ulnv ',ulnv(l,ie),' vlnv ',vlnv(l,ie)
                   !write(6,*) l,'utlnv',utlnv(l,ie),'vtlnv',vtlnv(l,ie)
                enddo
        enddo

        end


c****************************************************************************


c*****************************************************************

	subroutine debora(it)

	implicit none

	integer it

	integer icall
	save icall
	data icall /0/

	if( icall .gt. 0 ) return
	if( it .lt. 258000 ) return
	!if( it .lt. 3000 ) return

	icall = 1

	write(6,*) 'changed in debora...'
	call putpar('ilin',0.)
	call putpar('itsplt',2.)

	end

c*****************************************************************

        subroutine tsinitdebora(it)

	implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlv,nlvdi
        common /level/ nlvdi,nlv
	integer ie,k,i,l,it
        integer itype
        integer iarv(1)
        common /iarv/iarv
        integer nen3v(3,1)
        common /nen3v/nen3v
        real tempv(nlvdim,1)
        common /tempv/tempv
        real saltv(nlvdim,1)
        common /saltv/saltv

	integer icall
	save icall
	data icall /0/

        if(icall .gt. 0) return

        do l=1,nlv
             do ie=1,nel
                itype = iarv(ie)        
               do i=1,3
                        k=nen3v(i,ie)
                        if( itype .lt. 90) then
                          saltv(l,k)=30.
                          tempv(l,k)=14.
                        else
                          saltv(l,k)=35.
                          tempv(l,k)=14.
                        endif
                end do
             end do
        end do

        icall = 1

        end

c*****************************************************************

	subroutine zinit

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real znv(1)
	common /znv/znv
	real zenv(3,1)
	common /zenv/zenv
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer k,ie,ii
	real xmin,xmax,ymin,ymax
	real zav,dz,pi
	real x,y,r,z

	integer icall
	save icall
	data icall / 0 /

	if( icall .gt. 0 ) return

	icall = 1

	zav = 0.
	dz = .20

	pi = 4.*atan(1.)

	xmax = xgv(1)
	xmin = xgv(1)
	ymax = ygv(1)
	ymin = ygv(1)

	do k=1,nkn
	  xmax = max(xmax,xgv(k))
	  xmin = min(xmin,xgv(k))
	  ymax = max(ymax,ygv(k))
	  ymin = min(ymin,ygv(k))
	end do

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  r = (y - ymin)/(ymax-ymin)	! r is in [0:+1]
	  z = zav + (2*r-1.)*dz		! linear
	  z = zav - dz * cos(r*pi)	! sinusoidal
	  znv(k) = z
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

	end

c****************************************************************

	subroutine cprint(it)

	implicit none

	include 'param.h'

	integer it

	integer i,k,lmax,l
	logical berror

	real cnv(nlvdim,nkndim)
	common /cnv/cnv
	real saltv(nlvdim,1)
	common /saltv/saltv
	real tempv(nlvdim,1)
	common /tempv/tempv
        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer icall
	save icall
	data icall / 0 /

	integer ndim
	parameter (ndim=5)
	integer nodes(ndim)
	save nodes
	!data nodes / 984, 4860, 4636, 4585 /
	data nodes / 984, 4860, 4636, 4585 , 3452 /

	if( icall .eq. 0 ) then
	  icall = 1
	  call n2int(ndim,nodes,berror)
	  if( berror ) stop 'error stop cprint'
	end if

	write(87,*)
	write(87,*) 'time = ',it
	write(87,*)

	do i=1,ndim
	  k = nodes(i)
	  lmax = ilhkv(k)
	  write(87,*) i,k,lmax,(cnv(l,k),l=1,lmax)
	  write(81,*) i,k,lmax,(saltv(l,k),l=1,lmax)
	  write(82,*) i,k,lmax,(tempv(l,k),l=1,lmax)
	end do

	write(86,*) it,(cnv(1,nodes(i)),i=1,ndim)

	end

c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************

	subroutine diffus2d

c test for 2D horizontal diffusion algorithm

c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )		(n=n)
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )		(n=1)
c
c C(x,t) =  1/(4*pi*a*t) * exp( -|x|**2/(4*a*t) )		(n=2)

	implicit none

	include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real cnv(nlvdim,nkndim)
        common /cnv/cnv

	integer ndim
	parameter (ndim=101)

	integer k,k0,i,kin,n
	integer kmin,kmax,kstep
	integer itt,idtfrq,it0
	real c0,ctot,ct0,c
	real dx,rk
	real depth
	real caux(ndim)
	real xv(2*ndim)
	real yv(2*ndim)
	character*80 file

	integer ipint
	real getpar

	save ct0

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------
c initialization
c----------------------------------------------------------

	idtfrq = 500
	it0 = 1000		!initialize with this time

	c0 = 100.
	k0 = 51051
	dx = 100.
	depth = 10.

	if( icall .eq. 0 ) then
	  rk = getpar('chpar')
	  kin = ipint(k0)
	  do k=kin,kin
	    cnv(1,k) = c0
	  end do
	  call tsmass(cnv,+1,nlvdim,ctot)
	  write(6,*) 'diffus2d initialized (1): ',c0,ctot
	  call cv_init(it0,k0,rk,c0,ct0,cnv)
	  call tsmass(cnv,+1,nlvdim,ctot)
	  write(6,*) 'diffus2d initialized: ',it0,it,itanf,idt
	  write(6,*) 'diffus2d initialized: ',c0,ctot,ct0
	  ct0 = ctot/depth
	end if

	icall = icall + 1
	itt = it + it0 - idt

c----------------------------------------------------------
c compute average over section
c----------------------------------------------------------

	if( mod(itt,idtfrq) .ne. 0 ) return

	rk = getpar('chpar')
	call tsmass(cnv,+1,nlvdim,ctot)
	kin = ipint(k0)
	c = cnv(1,kin)
	write(6,*) 'diffus2d: ',c0,c,ctot,ct0

	kmin = 51001
	kmax = 51101
	kstep = 1
	call extract_conz(cnv,ndim,xv,yv,k0,kmin,kmax,kstep)
        call make_filename(itt,file,'hor','.txt')
	call write_file(file,ndim,xv,yv)

	kmin =   1051
	kmax = 101051
	kstep = 1000
	call extract_conz(cnv,ndim,xv,yv,k0,kmin,kmax,kstep)
        call make_filename(itt,file,'ver','.txt')
	call write_file(file,ndim,xv,yv)

	kmin =   1001
	kmax = 101101
	kstep = 1001
	call extract_conz(cnv,ndim,xv,yv,k0,kmin,kmax,kstep)
        call make_filename(itt,file,'dia','.txt')
	call write_file(file,ndim,xv,yv)

	call make_anal(itt,ndim,xv,yv,ct0,rk,dx)
        call make_filename(itt,file,'ana','.txt')
	call write_file(file,ndim,xv,yv)

	call extract_irreg(cnv,2*ndim,n,xv,yv,k0,dx,60.)
        call make_filename(itt,file,'p60','.txt')
	call write_file(file,n,xv,yv)

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c****************************************************************

	subroutine write_file(file,ndim,xv,yv)

	implicit none

	character*(*) file
	integer ndim
	real xv(ndim)
	real yv(ndim)

	integer i,iu

	iu = 88

	open(iu,file=file,status='unknown',form='formatted')
	do i=1,ndim
	  write(iu,*) xv(i),yv(i)
	end do
	close(iu)

	end

c****************************************************************

	subroutine make_anal(it,ndim,xv,yv,c0,rk,dx)

c C(x,t) =  1/(4*pi*a*t) * exp( -|x|**2/(4*a*t) )		(n=2)

	implicit none

	include 'param.h'

	integer it
	integer ndim
	real xv(ndim),yv(ndim)
	real c0,rk,dx

	integer i,nc
	real pi,aux,x,y,r
	real ctot,dxx,daux,fact
	double precision dtot

	pi = 4. * atan(1.)
	aux = 4. * rk * it
	nc = 1 + ndim / 2

	do i=1,ndim
	  x = (i-nc)*dx
	  y = c0 * exp(-x*x/aux) / (pi*aux)
	  xv(i) = x
	  yv(i) = y
	end do

	fact = 0.9
	fact = 0.95
	r = nc * dx
	dtot = 0.
	do i=1,1000
	  x = r
	  y = c0 * exp(-x*x/aux) / (pi*aux)
	  dxx = (1.-fact) * r
	  daux = 2.*pi*r*y*dxx
	  dtot = dtot + daux
	  if( dtot .gt. 100 .and. daux .lt. 1.e-10 ) goto 1
	  !write(6,*) i,x,y,daux,real(dtot)
	  r = fact * r
	end do
	write(6,*) i,daux
	stop 'error stop make_anal: error too high'

    1	continue
	ctot = dtot
	write(6,*) 'analytic solution: ',i,y,daux,ctot

	end

c****************************************************************

	subroutine cv_init(it,k0,rk,c0,ct0,cv)

c initializes cnv with analytic solution (2D)

	implicit none

	include 'param.h'

	integer it,k0
	real rk,c0,ct0
        real cv(nlvdim,nkndim)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer i,k,kin
	real x,y,dx,dy,r
	real x0,y0
	real pi,aux,ctot

	integer ipint

	pi = 4. * atan(1.)
	aux = 4. * rk * it
	ctot = c0 * pi * aux
	ct0 = ctot

	kin = ipint(k0)
	x0 = xgv(kin)
	y0 = ygv(kin)

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  dx = x-x0
	  dy = y-y0
	  r = sqrt( dx*dx + dy*dy )
	  y = ctot * exp(-r*r/aux) / (pi*aux)
	  cv(1,k) = y
	end do

	end

c****************************************************************

	subroutine extract_irreg(cv,ndim,n,xv,yv,k0,dx,phi)

c extracts conz from array

	implicit none

	include 'param.h'

        real cv(nlvdim,nkndim)
	integer ndim
	integer n
	real xv(ndim)
	real yv(ndim)
	integer k0
	real dx,phi

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer i,kin0,ie,ip
	real x0,y0
	real pi,rad,phi0
	real r,x,y
	real zp

	integer ipint

	write(6,*) 'starting to extract... phi = ',phi

	kin0 = ipint(k0)
	x0 = xgv(kin0)
	y0 = ygv(kin0)
	write(6,*) 'center... ',k0,kin0,x0,y0

c------------------------------------------------------------
c prepare parameters
c------------------------------------------------------------

	xv(1) = 0.
	yv(1) = cv(1,kin0)

	pi = 4.*atan(1.)
	rad = pi/180.
	phi0 = phi * rad

c------------------------------------------------------------
c go backward
c------------------------------------------------------------

	ie = 0
	do i=2,ndim
	  r = -(i-1)*dx
	  x = x0 + r * cos(phi0)
	  y = y0 + r * sin(phi0)
	  call intp_on_node(ie,x,y,cv,zp)
	write(6,*) ie,r,x,y,zp
	  if( ie .le. 0 ) goto 1
	  xv(i) = r
	  yv(i) = zp
	end do
	stop 'error stop extract_irreg: too many nodes (1)'
    1	continue
	n = i-1
	write(6,*) 'backward... ',n
	
c------------------------------------------------------------
c invert arrays
c------------------------------------------------------------

	do i=1,n/2
	  ip = n+1-i
	  call swap(xv(i),xv(ip))
	  call swap(yv(i),yv(ip))
	end do

c------------------------------------------------------------
c go foreward
c------------------------------------------------------------

	ie = 0
	do i=n+1,ndim
	  r = (i-n)*dx
	  x = x0 + r * cos(phi0)
	  y = y0 + r * sin(phi0)
	  call intp_on_node(ie,x,y,cv,zp)
	  if( ie .le. 0 ) goto 2
	  xv(i) = r
	  yv(i) = zp
	end do
	stop 'error stop extract_irreg: too many nodes (2)'
    2	continue
	n = i-1
	write(6,*) 'foreward... ',n
	
	write(6,*) 'finished to extract... n = ',n

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c****************************************************************

	subroutine swap(a,b)

	implicit none

	real a,b

	real aux

	aux = a
	a = b
	b = aux

	end

c****************************************************************

	subroutine intp_on_node(ie,x,y,cv,zp)

	implicit none

	include 'param.h'

	integer ie
	real x,y
	real cv(nlvdim,nkndim)
	real zp

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ii,k
	real z(3)

	call find_elem_from_old(ie,x,y,ie)

	if( ie .le. 0 ) return

	do ii=1,3
	  k = nen3v(ii,ie)
	  z(ii) = cv(1,k)
	end do

	call femintp(ie,z,x,y,zp)

	end

c****************************************************************

	subroutine extract_conz(cv,ndim,xv,yv,k0,kmin,kmax,kstep)

c extracts conz from array

	implicit none

	include 'param.h'

        real cv(nlvdim,nkndim)
	integer ndim
	real xv(ndim)
	real yv(ndim)
	integer k0,kmin,kmax,kstep

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer i,k,kin
	real x,y,dx,dy,dxy
	real x0,y0
	real fact

	integer ipint

	kin = ipint(k0)
	x0 = xgv(kin)
	y0 = ygv(kin)

	i = 0
	fact = -1
	do k=kmin,kmax,kstep
	  kin = ipint(k)
	  if( k .eq. k0 ) fact = +1
	  x = xgv(kin)
	  y = ygv(kin)
	  dx = x-x0
	  dy = y-y0
	  dxy = sqrt( dx*dx + dy*dy )
	  i = i + 1
	  if( i .gt. ndim ) stop 'error stop extract_conz: ndim'
	  xv(i) = fact*dxy
	  yv(i) = cv(1,kin)
	end do

	end

c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************

	subroutine diffus

c test for horizontal diffusion algorithm

	implicit none

	include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real cnv(nlvdim,nkndim)
        common /cnv/cnv

	integer k,k0,k1,i
	real c0,cmed
	real caux(25)

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------
c initialization
c----------------------------------------------------------

	c0 = 100.

	if( icall .eq. 0 ) then
	  do k=109,117
	    cnv(1,k) = c0
	  end do
	end if

	icall = icall + 1

c----------------------------------------------------------
c compute average over section
c----------------------------------------------------------

	do i=1,25
	  k1 = i*9
	  k0 = k1 - 8
	  cmed = 0.
	  do k=k0,k1
	    cmed = cmed + cnv(1,k)
	  end do
	  cmed = cmed / 9
	  caux(i) = cmed
	end do

	write(88,*) it,caux

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c****************************************************************

	subroutine viscos

c viscosity algorithm

	implicit none

	include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real vprv(nlvdim,1)
	common /vprv/vprv

	integer k,i,l

	write(88,*) it
	do k=109,117
	  write(88,'(12i5)') k,nlv,(nint(1000*vprv(l,k)),l=1,nlv)
	end do

	write(89,*) it
	do l=1,nlv
	  write(89,'(12i5)') l,(nint(1000*vprv(l,k)),k=109,117)
	end do

	end

c****************************************************************

	subroutine vdiffus(mode)

c diffusion algorithm
c
c solution of purely diffusional part :
c
c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
c
c for n-dimensions and
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
c
c for 1 dimension
c
c the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area

	implicit none

	integer mode

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer k,i,l,lc,lmax
	integer it0
	real c0,cmed,pi,a,t,aux,x2,cc,dc
	real caux(nlvdim)
	save it0

	real getpar

	integer icall
	save icall
	data icall / 0 /

c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )

	lc = (nlv+1)/2
	c0 = 100.
	pi = 4. * atan(1.)
	a = getpar('diftur')
	dc = -0.5

	if( icall .eq. 0 ) then
	  if( mode .eq. 1 ) then
	    it0 = it
	    do k=1,nkn
	      cnv(lc,k) = c0
	    end do
	  else
	    do k=1,nkn
	      lmax = ilhkv(k)
	      do l=1,lmax
	        cnv(l,k) = c0 + l * dc
	      end do
	    end do
	  end if
	  icall = 1
	  write(68,*) 'vdiffus: ',mode,lc,c0,pi,a,it0  
	end if

	t = it - it0
	aux = 4. * a * t
	do l=1,nlv
	  x2 = (l-lc)**2
	  cc = 0.
	  if( aux .gt. 0. ) cc = c0 * exp( -x2 / aux ) / sqrt( pi * aux )
	  caux(l) = cc
	end do
	if( aux .eq. 0 ) caux(lc) = c0

	write(68,*) it,nlv
	do l=1,nlv
	  !k = i*9 - 4
	  write(68,'(i4,6f12.4)') l,(cnv(l,i*9-4),i=1,25,6),caux(l)
	end do

	end

c****************************************************************

        subroutine uv_bottom

        implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(1)
        common /ilhkv/ilhkv

        integer k,l,m
        real u,v,uz,cdir

        real rdebug(nkndim)
        integer iudeb
        save rdebug,iudeb
        data iudeb /0/

        do k = 1,nkn

          l = ilhkv(k)                          !bottom layer
          m = l

          call getuv(m,k,u,v)                   !gets velocities u/v
          call getmd(u,v,UZ,CDIR)               !gets UZ, CDIR

          rdebug(k) = uz
        end do

        write(6,*) 'debug value written... ',iudeb
        call conwrite(iudeb,'.ggu',1,888,1,rdebug)

        end

c********************************************************************

	subroutine joel

c writes node list for Malta Coastal Model

	implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhkv(1)
        common /ilhkv/ilhkv
        real hlv(1)
        common /hlv/hlv
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real hev(1)
	common /hev/hev
	real v1v(1)
	common /v1v/v1v

        integer k,l,ie,ii,i
        real x,y,h

	real x0,y0,phi0
	real pi,rad
	real xfact,yfact

	integer ibc,nnodes,kext
	integer nkbnds,kbnds,ipext

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------
c do only at first call
c----------------------------------------------------------

	if( icall .ne. 0 ) return
	icall = 1

c----------------------------------------------------------
c initialize geographic projection
c----------------------------------------------------------

	x0 = 14.047
	y0 = 35.672
	phi0 = 36.

	pi = 4. * atan(1.)
	rad = pi/180.

	yfact = 60*1852
	xfact = yfact * cos(phi0*rad)

c----------------------------------------------------------
c get deepest depth on node
c----------------------------------------------------------

        do k = 1,nkn
	  v1v(k) = -999.
	end do

        do ie = 1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    v1v(k) = max(v1v(k),hev(ie))
	  end do
	end do

c----------------------------------------------------------
c write out values for total domain
c----------------------------------------------------------

	open(1,file='joel_total.dat')

	write(1,*) nkn
        do k = 1,nkn
	  kext = ipext(k)
	  l = ilhkv(k)
	  h = v1v(k)
	  x = xgv(k) / xfact + x0
	  y = ygv(k) / yfact + y0
	  write(1,*) k,kext,l,x,y,h
	end do

	close(1)

c----------------------------------------------------------
c write out values for open boundary
c----------------------------------------------------------

	if( nbc .gt. 0 ) then
	open(1,file='joel_bound.dat')

	ibc = 1
	nnodes = nkbnds(ibc)

	write(1,*) nnodes
        do i = 1,nnodes
	  k = kbnds(ibc,i)			!get boundary nodes
	  kext = ipext(k)
	  l = ilhkv(k)
	  h = v1v(k)
	  x = xgv(k) / xfact + x0
	  y = ygv(k) / yfact + y0
	  write(1,*) k,kext,l,x,y,h
	end do

	close(1)
	end if

c----------------------------------------------------------
c write out values z-levels
c----------------------------------------------------------

	open(1,file='joel_zlevels.dat')

	write(1,*) nlv
        do l = 1,nlv
	  h = hlv(l)
	  write(1,*) l,h
	end do

	close(1)

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c********************************************************************

        subroutine andreac

c introduce le condizioni iniziali per la concentrazione

        implicit none

        include 'param.h'

        real cnv(nlvdim,nkndim)
        common /cnv/cnv

        integer icall
        save icall
        data icall / 0 /

        if(icall.eq.0)then
          cnv(1,4311)=100.
          cnv(1,4312)=100.
       	  icall=1
        end if

	end

c********************************************************************

        subroutine tvd_test(it)

c introduce le condizioni iniziali per la concentrazione

        implicit none

	integer it

        include 'param.h'

	integer ndim
	parameter (ndim=2000)

	integer nelems
	integer iel(ndim)
	real conz(ndim)
	real yco(ndim)
	save nelems,iel,conz,yco

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        integer nen3v(3,1)
        common /nen3v/nen3v

	integer i,ii,ie,k,ien
	integer igrad
	real grad,gradaux,gradmax
	real dy,x0,y,ymin,ymax
	real cc(3)
	real ccc,cmin,cmax

        integer icall
        save icall
        data icall / 0 /

	dy = 100.
	x0 = 1100.
	call mima(ygv,nkn,ymin,ymax)

        if( icall .eq. 0 ) then
       	  icall = 1

	  y = ymin
	  i = 0
	  ie = 0
	  do while( y+dy .le. ymax )
	    i = i + 1
	    if( i .gt. ndim ) stop 'error stop tvd_test: ndim'
	    y = ymin + (i-1) * dy + dy/2.
	    call find_elem_from_old(ie,x0,y,ien)
	    !write(6,*) ymin,ymax,x0,y,ie,ien
	    if( ien .le. 0 ) stop 'error stop tvd_test (1): iel'
	    iel(i) = ien
	    yco(i) = y
	    ie = ien
	  end do

	  nelems = i
	  write(6,*) 'set up of tvd test: ',ymin,ymax,nelems

	  write(68,*) 46728645,1
	  write(68,*) nelems,0
	  write(68,*) 0,nelems+1,-10,110.

        end if

	do i=1,nelems

	  y = yco(i)
	  ie = iel(i)

	  do ii=1,3
	    k = nen3v(ii,ie)
	    cc(ii) = cnv(1,k)
	  end do

	  call femintp(ie,cc,x0,y,ccc)

	  conz(i) = ccc
	end do

	cmin = +1.e+30
	cmax = -1.e+30
	do k=1,nkn
	  cmin = min(cmin,cnv(1,k))
	  cmax = max(cmax,cnv(1,k))
	end do
	write(71,*) it,cmin,cmax	!only for check
	cmin = min(cmin,0.)
	cmax = max(cmax,100.)
	write(70,*) it,cmin,cmax-100.

	write(68,*) it,nelems
	write(68,*) (conz(i),i=1,nelems)

	grad = 0.
	igrad = 0
	do i=4,nelems-1			!really start from 2 -> bug
	  gradaux = abs(conz(i+1)-conz(i-1))
	  if( gradaux .gt. grad ) then
	    grad = gradaux
	    igrad = i
	  end if
	end do

	gradmax = 100./(2.*dy)
	grad = grad / (2.*dy)
	grad = 100. * grad / gradmax	!is fraction of max grad in percentage
	write(69,*) it,grad,igrad

        end

c********************************************************************

        subroutine make_filename(it,file,pre,ext)

        implicit none

        integer it
        character*(*) file,pre,ext

        integer i,n
        integer ichanm

        write(file,'(a,i10,a)') pre,it,ext
        n = ichanm(file)
        do i=n,1,-1
          if( file(i:i) .eq. ' ' ) file(i:i) = '0'
        end do

        end

c**********************************************************************
