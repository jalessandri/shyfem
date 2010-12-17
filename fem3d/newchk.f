c
c $Id: newchk.f,v 1.40 2010-02-26 17:35:06 georg Exp $
c
c routines for various checks
c
c contents :
c
c subroutine test3d(iunit,nn)           test output for new variables
c subroutine check_all			checks arrays for sanity (shell)
c subroutine check_fem			checks arrays for sanity
c subroutine check_values		checks important variables
c subroutine tsmass(ts,z,nlvdim,tstot)  computes mass of T/S or any conc. ts
c subroutine debug_dry			writes debug information on dry areas
c subroutine debug_node(k)		writes debug information on node k
c subroutine mimafem(string)		writes some min/max values to stdout
c subroutine mass_conserve(vf,va)	checks mass conservation
c
c subroutine check_node(k)		debug info on node k
c subroutine check_elem(ie)		debug info on element ie
c subroutine check_nodes_in_elem(ie)	debug info on nodes in element ie
c subroutine check_elems_around_node(k) debug info on elements around node k
c
c revision log :
c
c 24.08.1998    ggu     levdbg used for debug
c 26.08.1998    ggu     subroutine tsmass transferred from newbcl0
c 26.08.1998    ggu     subroutine convol, tstvol transferred from newcon1
c 22.10.1999    ggu     igrv, ngrv eliminated (subst by ilinkv, lenkv)
c 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
c 05.12.2001    ggu     always execute tstvol, more debug info
c 11.10.2002    ggu     call to setaix deleted
c 09.01.2003    ggu     some variables saved in contst
c 27.03.2003    ggu     new routine value_check
c 31.07.2003    ggu     eliminate compiler warnings
c 31.07.2003    ggu     eliminate useless variables
c 10.08.2003    ggu     new routine check_fem
c 03.09.2003    ggu     routines check and sanity_check deleted
c 03.09.2003    ggu     renamed value_check to check_values, new check_all
c 13.03.2004    ggu     write total volume to inf file
c 15.03.2005    ggu     call to check_austausch() eliminated
c 23.03.2006    ggu     changed time step to real
c 23.08.2007    ggu     test for boundary nodes using routines in testbndo.h
c 27.09.2007    ggu     deleted tstvol,tstvol1,contst, new mass_conserve
c 24.06.2008    ggu     bpresv deleted
c 06.12.2008    ggu     read vreps from STR file
c 21.01.2009    ggu     new var vrerr to stop if mass error is too high
c 23.03.2009    ggu     more debug for vrerr, new routine check_node()
c 02.04.2009    ggu     new routine check_elem()
c 06.04.2009    ggu     new check_elems_around_node, check_nodes_in_elem
c 26.02.2010    ggu     in test3d() write also meteo data
c 08.04.2010    ggu     more info in checks (depth and area)
c
c*************************************************************

	subroutine test3d(iunit,nn)

c test output for new variables
c
c nn	number of first array elements to be printed

	implicit none

c argument
	integer iunit
	integer nn
c parameter
	include 'param.h'
c common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/nlvdi,nlv

	integer ilhv(1)
	real fcorv(1)
	real zov(1),znv(1)
	real ulov(nlvdim,1)
	real ulnv(nlvdim,1)
	real vlov(nlvdim,1)
	real vlnv(nlvdim,1)
	real wlov(0:nlvdim,1)
	real wlnv(0:nlvdim,1)
	real utlov(nlvdim,1)
	real utlnv(nlvdim,1)
	real vtlov(nlvdim,1)
	real vtlnv(nlvdim,1)
	real hlv(1),hldv(1)
	real hlhv(1)
	real hev(1)
	common /ilhv/ilhv, /fcorv/fcorv
	common /zov/zov, /znv/znv
	common /ulov/ulov
	common /vlov/vlov
	common /wlov/wlov
	common /ulnv/ulnv
	common /vlnv/vlnv
	common /wlnv/wlnv
	common /utlov/utlov
	common /vtlov/vtlov
	common /utlnv/utlnv
	common /vtlnv/vtlnv
	common /hlv/hlv, /hldv/hldv
	common /hlhv/hlhv
	common /hev/hev

	real visv(0:nlvdim,1)
	common /visv/visv
	real difv(0:nlvdim,1)
	common /difv/difv
        real saltv(nlvdim,1)
        real tempv(nlvdim,1)
        common /saltv/saltv
        common /tempv/tempv
        real zeov(3,1), zenv(3,1)
        common /zeov/zeov, /zenv/zenv
        real difhv(nlvdim,1)
        common /difhv/difhv
        integer iwegv(1)
        common /iwegv/iwegv

        real tauxnv(nkndim), tauynv(nkndim)
        real wxv(nkndim), wyv(nkndim)
        real ppv(nkndim)
        common /tauxnv/tauxnv, /tauynv/tauynv
        common /wxv/wxv, /wyv/wyv
        common /ppv/ppv

c local
	logical bmeteo
	integer i,l,nk,ne
	integer iu,ii

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in test3d'

	bmeteo = .false.

	iu = iunit
	if( iu .le. 0 ) iu = 6

	write(iu,*) 'time:',it
	write(iu,*) 'nn  :',nn
	write(iu,*) 'nkn :',nkn
	write(iu,*) 'nel :',nel
	write(iu,*) 'nlvdi,nlv :',nlvdi,nlv	

	write(iu,*) 'hlv :'
	write(iu,*) (hlv(l),l=1,nlv)
	write(iu,*) 'hldv :'
	write(iu,*) (hldv(l),l=1,nlv)

	if(nn.eq.0) then
		nk=nkn
		ne=nel
	else
		nk=min(nn,nkn)		!$$cmplerr
		ne=min(nn,nel)
	end if

	write(iu,*) 'ilhv :'
	write(iu,*) (ilhv(i),i=1,ne)
	write(iu,*) 'fcorv :'
	write(iu,*) (fcorv(i),i=1,ne)
	write(iu,*) 'hlhv :'
	write(iu,*) (hlhv(i),i=1,ne)
	write(iu,*) 'hev :'
	write(iu,*) (hev(i),i=1,ne)
	write(iu,*) 'iwegv :'
	write(iu,*) (iwegv(i),i=1,ne)

	write(iu,*) 'zov :'
	write(iu,*) (zov(i),i=1,nk)
	write(iu,*) 'znv :'
	write(iu,*) (znv(i),i=1,nk)
	write(iu,*) 'zeov :'
	write(iu,*) ((zeov(ii,i),ii=1,3),i=1,ne)
	write(iu,*) 'zenv :'
	write(iu,*) ((zenv(ii,i),ii=1,3),i=1,ne)

	if( bmeteo ) then
	write(iu,*) 'ppv :'
	write(iu,*) (ppv(i),i=1,nk)
	write(iu,*) 'wxv :'
	write(iu,*) (wxv(i),i=1,nk)
	write(iu,*) 'wyv :'
	write(iu,*) (wyv(i),i=1,nk)
	write(iu,*) 'tauxnv :'
	write(iu,*) (tauxnv(i),i=1,nk)
	write(iu,*) 'tauynv :'
	write(iu,*) (tauynv(i),i=1,nk)
	end if

	do l=1,nlv
	write(iu,*)
	write(iu,*) 'level :',l
	write(iu,*) 'ulov :'
	write(iu,*) (ulov(l,i),i=1,ne)
	write(iu,*) 'vlov :'
	write(iu,*) (vlov(l,i),i=1,ne)
	write(iu,*) 'wlov :'
	write(iu,*) (wlov(l-1,i),i=1,nk)
	write(iu,*) 'ulnv :'
	write(iu,*) (ulnv(l,i),i=1,ne)
	write(iu,*) 'vlnv :'
	write(iu,*) (vlnv(l,i),i=1,ne)
	write(iu,*) 'wlnv :'
	write(iu,*) (wlnv(l-1,i),i=1,nk)
	write(iu,*) 'utlov :'
	write(iu,*) (utlov(l,i),i=1,ne)
	write(iu,*) 'vtlov :'
	write(iu,*) (vtlov(l,i),i=1,ne)
	write(iu,*) 'utlnv :'
	write(iu,*) (utlnv(l,i),i=1,ne)
	write(iu,*) 'vtlnv :'
	write(iu,*) (vtlnv(l,i),i=1,ne)
	write(iu,*) 'visv :'
	write(iu,*) (visv(l,i),i=1,nk)
	write(iu,*) 'difv :'
	write(iu,*) (difv(l,i),i=1,nk)
	write(iu,*) 'tempv :'
	write(iu,*) (tempv(l,i),i=1,nk)
	write(iu,*) 'saltv :'
	write(iu,*) (saltv(l,i),i=1,nk)
	write(iu,*) 'difhv :'
	write(iu,*) (difhv(l,i),i=1,ne)
	end do

	end

c******************************************************************

	subroutine check_all

c checks arrays for sanity

	implicit none

        integer levdbg
        real getpar

        levdbg = nint(getpar('levdbg'))

        if( levdbg .ge. 5 ) call check_fem
        if( levdbg .ge. 2 ) call check_values

	!call mimafem('panic')

	end

c******************************************************************

	subroutine check_fem

c checks arrays for sanity

	implicit none

c-------------------------------------------------------
c check geom arrays
c-------------------------------------------------------

	call check_ev
	call check_geom

c-------------------------------------------------------
c check vertical structure
c-------------------------------------------------------

	call check_vertical

c-------------------------------------------------------
c check various arrays
c-------------------------------------------------------

	call check_eddy
	call check_coriolis
	call check_chezy

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c******************************************************************

	subroutine check_values

c checks important variables

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real hdknv(nlvdim,1)
        common /hdknv/hdknv
        real hdkov(nlvdim,1)
        common /hdkov/hdkov

        real hdenv(nlvdim,1)
        common /hdenv/hdenv
        real hdeov(nlvdim,1)
        common /hdeov/hdeov

        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv

        real zov(1),znv(1)
        common /zov/zov, /znv/znv
	real unv(1),vnv(1)
        common /unv/unv, /vnv/vnv

        real saltv(nlvdim,1)
        common /saltv/saltv
        real tempv(nlvdim,1)
        common /tempv/tempv

        real utlov(nlvdim,1),vtlov(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	common /utlnv/utlnv, /vtlnv/vtlnv           

        real ulov(nlvdim,1),vlov(nlvdim,1)
        common /ulov/ulov, /vlov/vlov
        real ulnv(nlvdim,1),vlnv(nlvdim,1)
	common /ulnv/ulnv, /vlnv/vlnv           

	character*16 text

	text = '*** value check'

	call check1Dr(nkn,zov,-10.,+10.,text,'zov')
	call check1Dr(nkn,znv,-10.,+10.,text,'znv')

	call check2Dr(3,3,nel,zeov,-10.,+10.,text,'zeov')
	call check2Dr(3,3,nel,zenv,-10.,+10.,text,'zenv')

	call check1Dr(nkn,unv,-10000.,+10000.,text,'unv')
	call check1Dr(nkn,vnv,-10000.,+10000.,text,'vnv')

	call check2Dr(nlvdim,nlv,nkn,utlnv,-10000.,+10000.,text,'utlnv')
	call check2Dr(nlvdim,nlv,nkn,vtlnv,-10000.,+10000.,text,'vtlnv')

	call check2Dr(nlvdim,nlv,nkn,ulnv,-10.,+10.,text,'ulnv')
	call check2Dr(nlvdim,nlv,nkn,vlnv,-10.,+10.,text,'vlnv')

	call check2Dr(nlvdim,nlv,nkn,tempv,-10.,+50.,text,'tempv')
	call check2Dr(nlvdim,nlv,nkn,saltv,-1.,+50.,text,'saltv')

	call check2Dr(nlvdim,nlv,nkn,hdknv,0.,+10000.,text,'hdknv')
	call check2Dr(nlvdim,nlv,nkn,hdkov,0.,+10000.,text,'hdkov')

	call check2Dr(nlvdim,nlv,nel,hdenv,0.,+10000.,text,'hdenv')
	call check2Dr(nlvdim,nlv,nel,hdeov,0.,+10000.,text,'hdeov')

	end

c**********************************************************************

        subroutine tsmass(ts,mode,nlvdi,tstot)

c computes mass of T/S or any concentration ts

        implicit none

        integer nlvdi          !dimension of levels
        real ts(nlvdi,1)       !concentration on nodes
c        real z(3,1)             !water level
	integer mode
        real tstot              !total computed mass of ts
c
	include 'param.h'

	double precision scalcont

	if( nlvdi .ne. nlvdim ) then
	  write(6,*) nlvdi,nlvdim
	  stop 'error stop tsmass: nlvdim'
	end if
	if( mode .ne. 1 .and. mode .ne. -1 ) then
	  write(6,*) 'mode = ',mode
	  stop 'error stop tsmass: wrong value for mode'
	end if

	tstot = scalcont(mode,ts)

        end

c************************************************************

	subroutine debug_dry

c writes debug information on dry areas

	implicit none

	include 'param.h'

c common
	integer itanf,itend,idt,nits,niter,it
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer iwegv(1)
	common /iwegv/iwegv

	integer ie,iweg

	iweg = 0
	do ie=1,nel
	  if( iwegv(ie) .ne. 0 ) iweg = iweg + 1
	end do

	write(6,*) 'drydry... ',it,iweg

	end

c*************************************************************

	subroutine debug_node(k)

c writes debug information on final volume around node k (internal)

	implicit none

	integer k

	include 'param.h'

c common
	integer itanf,itend,idt,nits,niter,it
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	real eps1,eps2,pi,flag,high
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	common /mkonst/ eps1,eps2,pi,flag,high
	integer nen3v(3,1)
	integer ilhv(1)
	real hev(1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	real utlov(nlvdim,1),vtlov(nlvdim,1)
	real wlov(0:nlvdim,1),wlnv(0:nlvdim,1)
	real uprv(nlvdim,1),vprv(nlvdim,1),wprv(0:nlvdim,1)
	common /nen3v/nen3v
	common /hev/hev
	common /utlnv/utlnv, /vtlnv/vtlnv, /wlnv/wlnv
	common /utlov/utlov, /vtlov/vtlov
	common /wlov/wlov
	common /uprv/uprv, /vprv/vprv, /wprv/wprv
	common /ilhv/ilhv
	real zov(1), znv(1)
	common /zov/zov, /znv/znv
        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv
	real unv(1),uov(1)
	common /unv/unv, /uov/uov
	real vnv(1),vov(1)
	common /vnv/vnv, /vov/vov
        integer iwegv(1)
        common /iwegv/iwegv
	include 'ev.h'

	integer ie,ii,kk,l,i
	integer ilevel
	integer iweg
	real flux,dzvol,avvol
	real diff,rdiff
	real aj,uv0,uv1
	real b,c
	real dt,az,azt,azpar

	real getpar
	integer ipext,ieext

	integer netot
	integer kinf(2,ngrdim)
	save netot,kinf

	integer kmem
	save kmem
	data kmem / 0 /

	if( k .ne. kmem ) then
	  netot = 0
          do ie=1,nel
            do ii=1,3
	      kk = nen3v(ii,ie)
	      if( kk .eq. k ) then
	        netot = netot + 1
	        if( netot .gt. ngrdim ) then
		  stop 'error stop debug_node: ngrdim'
	        end if
	        kinf(1,netot) = ie
	        kinf(2,netot) = ii
	      end if
	    end do
	  end do
	  kmem = k
	  write(6,*) 'new node for debug...'
	  write(6,*) k,ipext(k),netot
	  do i=1,netot
	    ie = kinf(1,i)
	    ii = kinf(2,i)
	    write(6,*) ie,ieext(ie),ii
	  end do
	end if

c compute inflow into column and volume of column

	call get_timestep(dt)
	call getaz(azpar)
	az = azpar
	azt = 1. - az

	flux = 0.	! flux into water column
	dzvol = 0.	! volume change due to water level change
	avvol = 0.	! average volume of water column
	iweg = 0

	do i=1,netot
	  ie = kinf(1,i)
	  ii = kinf(2,i)
          aj=ev(10,ie)
          ilevel=ilhv(ie)
          kk=nen3v(ii,ie)
	  if( kk .ne. k ) stop 'error stop debug_node: internal error'
          b=ev(ii+3,ie)
          c=ev(ii+6,ie)
          uv0=0.
          uv1=0.
          do l=ilevel,1,-1
            uv1=uv1+utlnv(l,ie)*b+vtlnv(l,ie)*c
            uv0=uv0+utlov(l,ie)*b+vtlov(l,ie)*c
          end do
          uv1=unv(ie)*b+vnv(ie)*c
          uv0=uov(ie)*b+vov(ie)*c
          flux  = flux  + dt*12.*aj*(uv0*azt+uv1*az)
          dzvol = dzvol + 4.*aj*(zenv(ii,ie)-zeov(ii,ie))
          avvol = avvol + 4.*aj*hev(ie)
	  if( iwegv(ie) .ne. 0 ) iweg = iweg + 1
        end do

	diff = abs(flux-dzvol)
	rdiff = 0.
	if( avvol .gt. 0. ) rdiff = diff/avvol

	write(6,*) 'debug... ',it,diff,rdiff,iweg

	end

c*************************************************************

        subroutine mimafem(string)

c writes some min/max values to stdout

        implicit none

	include 'param.h'

        character*(*) string
c common
        integer itanf,itend,idt,nits,niter,it
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
        common /femtim/ itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv

        real znv(1)
        real ulnv(nlvdim,1),vlnv(nlvdim,1)
        real wlnv(0:nlvdim,1)
        real up0v(1),vp0v(1)
        real saltv(nlvdim,1),tempv(nlvdim,1)
        common /znv/znv
        common /ulnv/ulnv, /vlnv/vlnv
        common /wlnv/wlnv
        common /up0v/up0v, /vp0v/vp0v
        common /saltv/saltv, /tempv/tempv

        real hdknv(nlvdim,1)
        common /hdknv/hdknv
        real hdkov(nlvdim,1)
        common /hdkov/hdkov

        real hdenv(nlvdim,1)
        common /hdenv/hdenv
        real hdeov(nlvdim,1)
        common /hdeov/hdeov

        real hlhv(1)
        common /hlhv/hlhv

        real utlov(nlvdim,1),vtlov(nlvdim,1)
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov, /utlnv/utlnv, /vtlnv/vtlnv           

c local
	integer ie,l,k
	real u,v,w,z,s,t,c
        real high
        real zmin,zmax,umin,umax,vmin,vmax
        real hknmax,hkomax,henmax,heomax
        real utomax,utnmax,vtomax,vtnmax
        real hlvmax,h1vmax
        real bprmax
c functions
	integer ipext,ieext

c-----------------------------------------------------
c initial check and write
c-----------------------------------------------------

	if(nlvdim.ne.nlvdi) stop 'error stop mimafem: level dimension'

        !return  !FIXME
        write(6,*) '------------------ ',string,' ',it

c-----------------------------------------------------
c check water levels and barotropic velocities
c-----------------------------------------------------

        high = 1.e+30

        zmin =  high
        zmax = -high
        umin =  high
        umax = -high
        vmin =  high
        vmax = -high

	do k=1,nkn
	  z = znv(k)
	  u = up0v(k)
	  v = vp0v(k)
          zmin = min(zmin,z)
          zmax = max(zmax,z)
          umin = min(umin,u)
          umax = max(umax,u)
          vmin = min(vmin,v)
          vmax = max(vmax,v)
	end do

        write(6,*) zmin,zmax,umin,umax,vmin,vmax

c-----------------------------------------------------
c check of layer thickness
c-----------------------------------------------------

        hknmax = -high
        hkomax = -high
        do k=1,nkndim
          do l=1,nlvdim
            hknmax = max(hknmax,hdknv(l,k))
            hkomax = max(hkomax,hdkov(l,k))
          end do
        end do

        henmax = -high
        heomax = -high
        do ie=1,neldim
          do l=1,nlvdim
            henmax = max(henmax,hdenv(l,ie))
            heomax = max(heomax,hdeov(l,ie))
          end do
        end do

        write(6,*) hknmax,hkomax,henmax,heomax

        hlvmax = -high
        do ie=1,neldim
            hlvmax = max(hlvmax,hlhv(ie))
        end do

        write(6,*) hlvmax

c-----------------------------------------------------
c check of transports
c-----------------------------------------------------

        utomax = -high
        utnmax = -high
        vtomax = -high
        vtnmax = -high
        do ie=1,neldim
          do l=1,nlvdim
            utomax = max(utomax,abs(utlov(l,ie)))
            utnmax = max(utnmax,abs(utlnv(l,ie)))
            vtomax = max(vtomax,abs(vtlov(l,ie)))
            vtnmax = max(vtnmax,abs(vtlnv(l,ie)))
          end do
        end do

        write(6,*) utomax,utnmax,vtomax,vtnmax

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

        write(6,*) '------------------'

        end

c*************************************************************

	subroutine mass_conserve(vf,va)

c checks mass conservation of single boxes (finite volumes)

	implicit none

	include 'param.h'

	real vf(nlvdim,1)
	real va(nlvdim,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        integer ilhv(1)
        integer ilhkv(1)
        real utlov(nlvdim,1),vtlov(nlvdim,1)
        real utlnv(nlvdim,1),vtlnv(nlvdim,1),wlnv(0:nlvdim,1)
        common /nen3v/nen3v
        common /utlov/utlov, /vtlov/vtlov
        common /utlnv/utlnv, /vtlnv/vtlnv, /wlnv/wlnv
        common /ilhv/ilhv
        common /ilhkv/ilhkv
        real mfluxv(nlvdim,1)
        common /mfluxv/mfluxv
	real rqv(1)
	common /rqv/rqv
	include 'ev.h'

	logical berror
	integer ie,l,ii,k,lmax,mode
	real am,az,azt,dt,azpar,ampar
	real areafv,b,c
	real ffn,ffo,ff
	real vmax,vrmax,vdiv,vdiff,vrdiff
	real abot,atop
	real volo,voln
	real ubar,vbar
	real vbmax,vlmax,vrbmax,vrlmax
	real vrwarn,vrerr

	real volnode,areanode,getpar
	include 'testbndo.h'

	integer ninfo
	save ninfo
	data ninfo /0/

c----------------------------------------------------------------
c initialize
c----------------------------------------------------------------

	if( ninfo .eq. 0 ) call getinfo(ninfo)

	vrwarn = getpar('vreps')
	vrerr = getpar('vrerr')

	mode = +1
        call getazam(azpar,ampar)
	az = azpar
	am = ampar
        azt = 1. - az
	call get_timestep(dt)

	do k=1,nkn
          lmax = ilhkv(k)
	  do l=1,lmax
	    vf(l,k) = 0.
	    va(l,k) = 0.
	  end do
	end do

c----------------------------------------------------------------
c compute horizontal divergence
c----------------------------------------------------------------

        do ie=1,nel
          areafv = 4. * ev(10,ie)               !area of triangle / 3
          lmax = ilhv(ie)
          do l=1,lmax
            do ii=1,3
                k=nen3v(ii,ie)
                b = ev(ii+3,ie)
                c = ev(ii+6,ie)
                ffn = utlnv(l,ie)*b + vtlnv(l,ie)*c
                ffo = utlov(l,ie)*b + vtlov(l,ie)*c
                ff = ffn * az + ffo * azt
                vf(l,k) = vf(l,k) + 3. * areafv * ff
                va(l,k) = va(l,k) + areafv
            end do
          end do
        end do

c----------------------------------------------------------------
c include vertical divergence
c----------------------------------------------------------------

	do k=1,nkn
          lmax = ilhkv(k)
	  abot = 0.
	  do l=lmax,1,-1
	    atop = va(l,k)
	    vdiv = wlnv(l,k)*abot - wlnv(l-1,k)*atop
	    vf(l,k) = vf(l,k) + vdiv + mfluxv(l,k)
	    abot = atop
	  end do
	end do

c----------------------------------------------------------------
c check mass balance in boxes
c----------------------------------------------------------------

	berror = .false.
	vrmax = 0.
	vmax = 0.
	do k=1,nkn
	 !if( is_inner(k) ) then
	 if( .not. is_external_boundary(k) ) then
	  berror = .false.
          lmax = ilhkv(k)
	  do l=1,lmax
	    voln = volnode(l,k,+1)
	    volo = volnode(l,k,-1)
	    vdiv = vf(l,k)
	    vdiff = voln - volo - vdiv * dt
	    vdiff = abs(vdiff)
	    vrdiff = vdiff / volo
	    vmax = max(vmax,vdiff)
	    vrmax = max(vrmax,vrdiff)
	    if( vrdiff .gt. vrerr ) then
		berror = .true.
	        write(6,*) l,k,volo,voln,vdiff,vrdiff
	    end if
	  end do
	  if( berror ) call check_node(k)
	 end if
	end do

	vlmax = vmax		!absolute error for each box
	vrlmax = vrmax		!relative error for each box

c----------------------------------------------------------------
c barotropic
c----------------------------------------------------------------

	do k=1,nkn
	    vf(1,k) = 0.
	    va(1,k) = 0.
	end do

        do ie=1,nel
          areafv = 4. * ev(10,ie)               !area of triangle / 3

	  ubar = 0.
	  vbar = 0.
          lmax = ilhv(ie)
          do l=1,lmax
	    ubar = ubar + az * utlnv(l,ie) + azt * utlov(l,ie)
	    vbar = vbar + az * vtlnv(l,ie) + azt * vtlov(l,ie)
	  end do

          do ii=1,3
                k=nen3v(ii,ie)
                b = ev(ii+3,ie)
                c = ev(ii+6,ie)
		ff = ubar * b + vbar * c
                vf(1,k) = vf(1,k) + 3. * areafv * ff
                va(1,k) = va(1,k) + areafv * ff
          end do
        end do

	vrmax = 0.
	vmax = 0.
	do k=1,nkn
	 !if( is_inner(k) ) then
	 if( .not. is_external_boundary(k) ) then
	    voln = volnode(1,k,+1)
	    volo = volnode(1,k,-1)
	    vdiv = vf(1,k) + rqv(k)
	    vdiff = voln - volo - vdiv * dt
	    vdiff = abs(vdiff)
	    vrdiff = vdiff / volo
	    vmax = max(vmax,vdiff)
	    vrmax = max(vrmax,vrdiff)
	    !write(78,*) l,k,volo,voln,vdiff,vrdiff
	 end if
	end do

	vbmax = vmax		!absolute error for water column
	vrbmax = vrmax		!relative error for water column

c----------------------------------------------------------------
c write diagnostic output
c----------------------------------------------------------------

c	vbmax 		!absolute error for water column
c	vrbmax 		!relative error for water column
c	vlmax 		!absolute error for each box
c	vrlmax 		!relative error for each box

	if( vrlmax .gt. vrwarn ) then
	  write(6,*) 'mass error: ',vbmax,vlmax,vrbmax,vrlmax
	  if( vrlmax .gt. vrerr ) then
	    write(6,*) 'mass error of matrix solution is very high'
	    write(6,*) 'the relative mass error is = ',vrlmax
	    write(6,*) 'the limit of the mass error is vrerr = ',vrerr
	    write(6,*) 'Probably there is some problem with the solution'
	    write(6,*) 'of the system matrix. However, if you think'
	    write(6,*) 'you can live with this mass error, then please'
	    write(6,*) 'increase the value of vrerr in the STR file.'
	    stop 'error stop mass_conserve: mass error too high'
	  end if
	end if

	write(ninfo,*) 'mass_balance: ',vbmax,vlmax,vrbmax,vrlmax

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*************************************************************
c*************************************************************
c*************************************************************
c*************************************************************
c*************************************************************

	subroutine check_node(k)

c writes debug information on node k

	implicit none

	integer k

	include 'param.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhv(1)
        common /ilhv/ilhv
        integer ilhkv(1)
        common /ilhkv/ilhkv
	integer inodv(1)
	common /inodv/inodv

        real wlnv(0:nlvdim,1)
        common /wlnv/wlnv
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
        real mfluxv(nlvdim,1)
        common /mfluxv/mfluxv
	real znv(1), zov(1)
	common /znv/znv, /zov/zov
        real saltv(nlvdim,1)
        common /saltv/saltv
        real tempv(nlvdim,1)
        common /tempv/tempv
	real visv(0:nlvdim,1)
	common /visv/visv
	real difv(0:nlvdim,1)
	common /difv/difv

        real hdknv(nlvdim,1)
        common /hdknv/hdknv
        real hdkov(nlvdim,1)
        common /hdkov/hdkov
        real areakv(nlvdim,1)
        common /areakv/areakv

	integer iu
	integer l,lmax,kk

	integer ipext
	real volnode

	iu = 16
	iu = 6
	lmax = ilhkv(k)

	write(iu,*) '-------------------------------- check_node'
	write(iu,*) 'it,k,kext:  ',it,k,ipext(k)
	write(iu,*) 'lmax,inodv: ',lmax,inodv(k)
	write(iu,*) 'xgv,ygv:    ',xgv(k),ygv(k)
	write(iu,*) 'zov,znv:    ',zov(k),znv(k)
	write(iu,*) 'hdkov:      ',(hdkov(l,k),l=1,lmax)
	write(iu,*) 'hdknv:      ',(hdknv(l,k),l=1,lmax)
	write(iu,*) 'areakv:     ',(areakv(l,k),l=1,lmax)
	write(iu,*) 'volold:     ',(volnode(l,k,-1),l=1,lmax)
	write(iu,*) 'volnew:     ',(volnode(l,k,+1),l=1,lmax)
	write(iu,*) 'wlnv:       ',(wlnv(l,k),l=0,lmax)
	write(iu,*) 'mfluxv:     ',(mfluxv(l,k),l=1,lmax)
	write(iu,*) 'tempv:      ',(tempv(l,k),l=1,lmax)
	write(iu,*) 'saltv:      ',(saltv(l,k),l=1,lmax)
	write(iu,*) 'visv:       ',(visv(l,k),l=1,lmax)
	write(iu,*) 'difv:       ',(difv(l,k),l=1,lmax)
	write(iu,*) '-------------------------------------------'

	end

c*************************************************************

	subroutine check_elem(ie)

c writes debug information on element ie

	implicit none

	integer iunit
	integer ie

	include 'param.h'
	include 'ev.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer nen3v(3,1)
	common /nen3v/nen3v
	integer ilhv(1)
	common /ilhv/ilhv
	integer ilhkv(1)
	common /ilhkv/ilhkv
	integer iwegv(1)
	common /iwegv/iwegv

	real hev(1)
	common /hev/hev
	real hlhv(1)
	common /hlhv/hlhv

	real zenv(3,1)
	common /zenv/zenv
	real zeov(3,1)
	common /zeov/zeov
	real znv(1), zov(1)
	common /znv/znv, /zov/zov

        real utlov(nlvdim,1),vtlov(nlvdim,1)
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov
        common /utlnv/utlnv, /vtlnv/vtlnv

        real ulnv(nlvdim,1)
        common /ulnv/ulnv
        real vlnv(nlvdim,1)
        common /vlnv/vlnv
        real wlnv(0:nlvdim,1)
        common /wlnv/wlnv

        real hdenv(nlvdim,1)
        common /hdenv/hdenv
        real hdeov(nlvdim,1)
        common /hdeov/hdeov

	integer iu
	integer l,lmax,ii

	integer ieext

	iu = 16
	iu = 6
	lmax = ilhv(ie)

	write(iu,*) '-------------------------------- check_elem'
	write(iu,*) 'it,ie,ieext: ',it,ie,ieext(ie)
	write(iu,*) 'lmax,iwegv:  ',lmax,iwegv(ie)
	write(iu,*) 'area:        ',ev(10,ie)*12.
	write(iu,*) 'nen3v  :     ',(nen3v(ii,ie),ii=1,3)
	write(iu,*) 'hev,hlhv:    ',hev(ie),hlhv(ie)
	write(iu,*) 'zeov:        ',(zeov(ii,ie),ii=1,3)
	write(iu,*) 'zenv:        ',(zenv(ii,ie),ii=1,3)
	write(iu,*) 'hdeov:       ',(hdeov(l,ie),l=1,lmax)
	write(iu,*) 'hdenv:       ',(hdenv(l,ie),l=1,lmax)
	write(iu,*) 'utlov:       ',(utlov(l,ie),l=1,lmax)
	write(iu,*) 'vtlov:       ',(vtlov(l,ie),l=1,lmax)
	write(iu,*) 'utlnv:       ',(utlnv(l,ie),l=1,lmax)
	write(iu,*) 'vtlnv:       ',(vtlnv(l,ie),l=1,lmax)
	write(iu,*) 'ulnv:        ',(ulnv(l,ie),l=1,lmax)
	write(iu,*) 'vlnv:        ',(vlnv(l,ie),l=1,lmax)
	write(iu,*) '-------------------------------------------'

	end

c*************************************************************

	subroutine check_nodes_in_elem(ie)

c writes debug information on nodes in element ie

	implicit none

	integer ie

        integer nen3v(3,1)
        common /nen3v/nen3v

	integer ii,k,iu
	integer ieext

	iu = 16
	iu = 6
	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking nodes in element: ',ie,ieext(ie)
	write(iu,*) '-------------------------------------------'

	do ii=1,3
	  k = nen3v(ii,ie)
	  call check_node(k)
	end do

	end

c*************************************************************

	subroutine check_elems_around_node(k)

c writes debug information on elements around node k

	implicit none

	integer k

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v

	integer ie,ii,kk,iu
	logical bdebug

	integer ipext

	iu = 16
	iu = 6
	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking elements around node: ',k,ipext(k)
	write(iu,*) '-------------------------------------------'

	do ie=1,nel
	  bdebug = .false.
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk .eq. k ) bdebug = .true.
	  end do
	  if( bdebug ) call check_elem(ie)
	end do

	end

c*************************************************************

