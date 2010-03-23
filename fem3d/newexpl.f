c
c $Id: newexpl.f,v 1.10 2010-03-08 17:46:45 georg Exp $
c
c explicit term routines
c
c contents :
c
c subroutine set_explicit
c subroutine viscous_stability(ahpar,ahstab)	computes stability for viscosity
c subroutine set_diff_horizontal_new1
c subroutine set_diff_horizontal_new
c subroutine set_advective
c subroutine set_semi_lagrange
c subroutine set_barocl
c subroutine set_barocl_new
c subroutine set_barocl_new1
c
c revision log :
c
c 01.05.2007	ggu	new file -> all explicit terms here
c 28.09.2007	ggu	semi-lagrangian part introduced
c 16.04.2008	ggu	bugfix in set_barocl (real do indices!!)
c 14.07.2008	ggu&ccf	ahpar is real in set_diff_horizontal_new
c 03.11.2008    ggu&dbf nudging implemented (call to bclnudge)
c 09.11.2008	ggu	set_barocl_new (cleaned version of set_barocl)
c 19.11.2008	ggu	new set_diff_horizontal_new1(), viscous_stability()
c 19.02.2010	ggu	in viscous_stability() for dt=1
c 26.02.2010	ggu	new call to momentum_viscous_stability()
c 26.02.2010	ggu	set_advective() cleaned up
c 26.02.2010	ggu	new momentum_advective_stability()
c 08.03.2010	ggu	run only down to avail layers (bug fix)
c
c notes :
c
c explicit term is on left side
c
c******************************************************************

	subroutine set_explicit

	implicit none
        
        include 'param.h'

	integer ie,l
        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real fxv(nlvdim,1)      !new HYDRO deb
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
        save /fxv/,/fyv/
        
        logical bbarcl
        integer ilin,itlin,ibarcl
        real getpar

c-------------------------------------------
c parameters
c-------------------------------------------

        ilin = nint(getpar('ilin'))
        itlin = nint(getpar('itlin'))
        ibarcl = nint(getpar('ibarcl'))
        bbarcl = ibarcl .eq. 1 .or. ibarcl .eq. 2

c-------------------------------------------
c initialization
c-------------------------------------------

	do ie=1,nel
	  do l=1,nlv
	    fxv(l,ie) = 0.
	    fyv(l,ie) = 0.
	  end do
	end do

c-------------------------------------------
c nudging
c-------------------------------------------

	call bclnudge

c-------------------------------------------
c horizontal diffusion
c-------------------------------------------

	call set_diff_horizontal_new1

c-------------------------------------------
c advective (non-linear) terms
c-------------------------------------------

        if( ilin .eq. 0 ) then
          if( itlin .eq. 0 ) then
	    call set_advective
	  else
	    call set_semi_lagrange
	  end if
	end if

c-------------------------------------------
c baroclinic contribution
c-------------------------------------------

        !if( bbarcl ) call set_barocl
        if( bbarcl ) call set_barocl_new

c-------------------------------------------
c end of routine
c-------------------------------------------

	end

c******************************************************************

	subroutine momentum_viscous_stability(ahpar,rindex,dstab)

c computes stability for viscosity

	implicit none

        include 'param.h'
	include 'ev.h'

	real ahpar
	real rindex
	real dstab(nlvdim,neldim)

        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ieltv(3,1)
	common /ieltv/ieltv
	integer ilhv(1)
	common /ilhv/ilhv
	real difhv(nlvdim,1)
	common /difhv/difhv
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,iei,l,lmax,k
	integer noslip
	real u,v,ui,vi
	real anu,ax,ay
	real area,areai
	real dt
	real a,ai,amax,afact

	rindex = 0.
	if( ahpar .le. 0 ) return

	amax = 0.

	do ie=1,nel

	  lmax = ilhv(ie)
	  area = 12. * ev(10,ie)

	  do l=1,lmax

	    a = 0.
	    do ii=1,3
              iei = ieltv(ii,ie)
              if( iei .le. 0 ) iei = ie

              areai = 12. * ev(10,iei)

	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              a = a + ai
	    end do

	    amax = max(amax,a)
	    dstab(l,ie) = a

          end do

	end do

	rindex = amax

	end

c******************************************************************

	subroutine set_diff_horizontal_new1

	implicit none

        include 'param.h'
	include 'ev.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlv,nlvdi
        common /level/ nlvdi,nlv

        real fxv(nlvdim,1)
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
        save /fxv/,/fyv/

	integer ieltv(3,1)
	common /ieltv/ieltv
	integer ilhv(1)
	common /ilhv/ilhv
        real utlov(nlvdim,1),vtlov(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov
	real difhv(nlvdim,1)
	common /difhv/difhv

	integer ie,ii,iei,l,lmax
	integer noslip
	real u,v,ui,vi
	real anu,ahpar,ax,ay
	real area,areai
	real dt
	real a,ai,amax,afact
	logical bnoslip

	real getpar

	call get_timestep(dt)

        ahpar = getpar('ahpar')
	if( ahpar .le. 0 ) return

        noslip = nint(getpar('noslip'))
	bnoslip = noslip .ne. 0

	amax = 0.

	do ie=1,nel

	  lmax = ilhv(ie)
	  area = 12. * ev(10,ie)

	  do l=1,lmax
	    u  = utlov(l,ie)
	    v  = vtlov(l,ie)

	    a = 0.
	    do ii=1,3
              iei = ieltv(ii,ie)
              afact = 1.
              if( bnoslip .and. iei .eq. 0 ) afact = -1.
              if( iei .le. 0 ) iei = ie

              areai = 12. * ev(10,iei)

	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              a = a + ai

	      ui = afact * utlov(l,iei)
	      vi = afact * vtlov(l,iei)

	      ax = ai * ( ui - u )
	      ay = ai * ( vi - v )

	      fxv(l,ie) = fxv(l,ie) - ax	!- because f is on left side
	      fyv(l,ie) = fyv(l,ie) - ay
	    end do
	    amax = max(amax,a)
          end do

	end do

	amax = amax * dt
	!write(99,*) 'stability viscosity: ',amax

	end

c******************************************************************

	subroutine set_diff_horizontal_new

	implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlv,nlvdi
        common /level/ nlvdi,nlv

        real saux1(nlvdim,1),saux2(nlvdim,1)
        real saux3(nlvdim,1),saux4(nlvdim,1)
        common /saux1/saux1, /saux2/saux2
        common /saux3/saux3, /saux4/saux4
        real fxv(nlvdim,1)      !new HYDRO deb
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
        save /fxv/,/fyv/

        real cb,cd 
        real ahpar,khpar 
        real vismol                     !constant vertical molecular viscosity
        integer noslip                  !no slip BC
        common /noslip/noslip
        common /rpara1/cb,cd,ahpar,khpar,vismol
        integer inodv(1)
        common /inodv/inodv
        integer ibound(1)
        common /ibound/ibound
        save /inodv/,/ibound/

        real v1v(1)
        common /v1v/v1v
        real v2v(1)
        common /v2v/v2v
        real v3v(1)
        common /v3v/v3v
	include 'ev.h'
        integer ilhv(1), ilhkv(1)
        common /ilhv/ilhv, /ilhkv/ilhkv
        integer nen3v(3,1), iarv(1)
        common /nen3v/nen3v, /iarv/iarv
        real utlov(nlvdim,1),vtlov(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov

	integer ie,k,ii,l,lmax
	real ao,u,v
	real anu,rv,acux,acuy
	real area,w,um,vm
	real dt
	real wm,wmax
        real getpar

	logical is_material_boundary_node

	is_material_boundary_node(k) = inodv(k) .ne. 0 
     +					.and. ibound(k) .eq. 0

	call get_timestep(dt)

        ahpar = getpar('ahpar')
	anu = ahpar
        !write(6,*)ahpar
	if( anu .le. 0 ) return
	!write(6,*) 'anu: ',anu

	do k=1,nkn
	  v1v(k) = 0.		!inverse area
	  v2v(k) = 0.		!grade
	  do l=1,nlv
	    saux1(l,k) = 0.		!average u
	    saux2(l,k) = 0.		!average v
	  end do
	end do

	do ie=1,nel
	  ao = ev(10,ie)
	  area = 12. * ao
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    v1v(k) = v1v(k) + anu / area
	    v2v(k) = v2v(k) + 1.
	    do l=1,lmax
	      u = utlov(l,ie)
	      v = vtlov(l,ie)
	      saux1(l,k) = saux1(l,k) + u
	      saux2(l,k) = saux2(l,k) + v
	    end do
	  end do
	end do

	do k=1,nkn
	  rv = 1. / v2v(k)
	  !if( noslip .ne. 0 .and. is_material_boundary_node(k) ) rv = 0.
	  lmax = ilhkv(k)
	  do l=1,lmax
	    saux1(l,k) = saux1(l,k) * rv
	    saux2(l,k) = saux2(l,k) * rv
	  end do
	end do

	wmax = 0.
	do ie=1,nel

	  wm = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    wm = wm + v1v(k)
	  end do
	  wmax = max(wmax,wm)

	  lmax = ilhv(ie)
	  do l=1,lmax
	    acux = 0.
	    acuy = 0.
	    u = utlov(l,ie)
	    v = vtlov(l,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      w = v1v(k)
	      um = saux1(l,k)
	      vm = saux2(l,k)
	      acux = acux - w * u + w * um
	      acuy = acuy - w * v + w * vm
	    end do
	    fxv(l,ie) = fxv(l,ie) - acux	!- because f is on left side
	    fyv(l,ie) = fyv(l,ie) - acuy
	  end do
	end do

	wmax = wmax * dt

	!write(99,*) 'stability diffusion: ',wmax

	end

c******************************************************************

        subroutine set_advective

        implicit none

        include 'param.h'

        integer nlv,nlvdi
        !common /nlv/nlv
        common /level/ nlvdi,nlv

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real saux1(nlvdim,1),saux2(nlvdim,1)
        real saux3(nlvdim,1),saux4(nlvdim,1)
        common /saux1/saux1, /saux2/saux2
        common /saux3/saux3, /saux4/saux4

        real fxv(nlvdim,1)      !new HYDRO deb
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
        save /fxv/,/fyv/
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real uprv(nlvdim,1),vprv(nlvdim,1)
        common /uprv/uprv, /vprv/vprv
        real hdenv(nlvdim,1)
        common /hdenv/hdenv
        real hdknv(nlvdim,1)
        common /hdknv/hdknv
	include 'ev.h'
        integer nen3v(3,1), iarv(1)
        common /nen3v/nen3v, /iarv/iarv
        integer ilhv(1)
        common /ilhv/ilhv
        integer ilhkv(1)
        common /ilhkv/ilhkv
        real radv
        common /radv/radv

        integer ii,ie,k,l,lmax
        real b,c
        real ut,vt
        real uc,vc
        real up,vp
        real um,vm
        real f,h
	real xadv,yadv
	real area,vol

	!write(6,*) 'set_advective called...'

c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
            saux1(l,k) = 0.
	    saux2(l,k) = 0.
	    saux3(l,k) = 0.
	  end do
	end do

c---------------------------------------------------------------
c accumulate momentum that flows into nodes (weighted by flux)
c---------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
            h = hdenv(l,ie)
	    ut = utlnv(l,ie)
	    vt = vtlnv(l,ie)
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
                f = ut * b + vt * c	! f>0 => flux into node
                if( f .gt. 0. ) then
		  saux1(l,k) = saux1(l,k) + f
		  saux2(l,k) = saux2(l,k) + f * ut
		  saux3(l,k) = saux3(l,k) + f * vt
                end if
	    end do
          end do
	end do

c---------------------------------------------------------------
c compute average momentum for every node
c---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
            h = hdknv(l,k)
	    if( saux1(l,k) .gt. 0 ) then	!flux into node
	      saux2(l,k) = saux2(l,k) / saux1(l,k)
	      saux3(l,k) = saux3(l,k) / saux1(l,k)
	    else				!only flux out of node
	      saux2(l,k) = uprv(l,k) * h
	      saux3(l,k) = vprv(l,k) * h
	    end if
	  end do
	end do

c---------------------------------------------------------------
c compute advective contribution
c---------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax

	    area = 12. * ev(10,ie)
            h = hdenv(l,ie)
	    vol = area * h
  	    ut = utlnv(l,ie)
  	    vt = vtlnv(l,ie)
            uc = ut / h
            vc = vt / h

	    xadv = 0.
	    yadv = 0.
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
                up = saux2(l,k) / h		!NEW
                vp = saux3(l,k) / h
                f = ut * b + vt * c
                if( f .lt. 0. ) then	!flux out of node => into element
                  xadv = xadv + f * ( up - uc )
                  yadv = yadv + f * ( vp - vc )
                end if
            end do

	    fxv(l,ie) = fxv(l,ie) + xadv
	    fyv(l,ie) = fyv(l,ie) + yadv

	  end do
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

        end

c******************************************************************

	subroutine momentum_advective_stability(rindex,astab)

c computes courant number of advective terms in momentum equation

	implicit none

        include 'param.h'

	real rindex
	real astab(nlvdim,neldim)

        integer nlv,nlvdi
        !common /nlv/nlv
        common /level/ nlvdi,nlv

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'ev.h'
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real hdenv(nlvdim,1)
        common /hdenv/hdenv
        integer ilhv(1)
        common /ilhv/ilhv
        integer nen3v(3,1)
        common /nen3v/nen3v

	integer ie,l,ii,k,lmax
	real cc,cmax
	real ut,vt
	real area,h,vol
	real b,c,f,ftot

	cmax = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax

	    area = 12. * ev(10,ie)
            h = hdenv(l,ie)
	    vol = area * h

  	    ut = utlnv(l,ie)
  	    vt = vtlnv(l,ie)

	    ftot = 0.
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
                f = ut * b + vt * c
                if( f .lt. 0. ) ftot = ftot - f
            end do

	    cc = area*ftot/vol
	    astab(l,ie) = cc
	    cmax = max(cmax,cc)

	  end do
	end do

	rindex = cmax

	end

c******************************************************************

	subroutine set_semi_lagrange

        implicit none
         
        include 'param.h'
        
	integer ie,l
	real xadv,yadv,dt
        real uadv(neldim),vadv(neldim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real fxv(nlvdim,1)
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv

	call get_timestep(dt)

        call back_trace(uadv,vadv)

	l = 1			!only for one layer

	do ie=1,nel
          xadv = uadv(ie) / dt
          yadv = vadv(ie) / dt

	  fxv(l,ie) = fxv(l,ie) + xadv
	  fyv(l,ie) = fyv(l,ie) + yadv
	end do

	end

c******************************************************************

        subroutine set_barocl

        implicit none
         
        include 'param.h'
	include 'ev.h'
        
        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
        real fxv(nlvdim,1)      !new HYDRO deb
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
        save /fxv/,/fyv/
        real rhov(nlvdim,1)
        common /rhov/rhov
        !integer itanf,itend,idt,nits,niter,it
        !integer itemp,isalt
        !real k,l,ie,ii				!BUG
        integer k,l,ie,ii			!BUG
        real dt
        real rrho0
        real salref,temref,sstrat,tstrat
        real hlayer
        real hhi
        real hdeov(nlvdim,1)
        common /hdeov/hdeov
        real hdkov(nlvdim,1)
        common /hdkov/hdkov

        real xbcl,ybcl
        real bpresxv(nlvdim,1),bpresyv(nlvdim,1)!deb 100407
        common /bpresxv/bpresxv, /bpresyv/bpresyv!deb 100407
        integer lmax
        integer ilhv(1)
        common /ilhv/ilhv
        integer nen3v(3,1)
        common /nen3v/nen3v

        real rhop,presbt,presbcx,presbcy,dprescx,dprescy,br,cr!deb
        real b,c
c        integer icall
c        save icall
c        data icall /0/

        if(nlvdim.ne.nlvdi)stop 'error stop : level dimension in barocl'

c        if(icall.eq.-1) return


        call get_timestep(dt)

        do ie=1,nel
                do l=1,nlv
                        bpresxv(l,ie) = 0.
                        bpresyv(l,ie) = 0.
                enddo
        enddo

        do ie=1,nel
            presbcx = 0.
            presbcy = 0.
            lmax=ilhv(ie)
            !print*,lmax,' lmax'
            do l=1,lmax            
                hlayer = 0.5 * hdeov(l,ie)
                
		br = 0.
		cr = 0.                 
                do ii=1,3                 
                        k = nen3v(ii,ie)
                        rhop = rhov(l,k) ! rho^prime for each node of element 
                        !print*,'rhov ', l,k,rhov(l,k)
                        b = ev(3+ii,ie)!gradient in x della funz di forma
                        c = ev(6+ii,ie)!gradient in y della funz di forma
                        br = br + (b*rhop) 
                        cr = cr + (c*rhop)
                end do
                presbcx = presbcx + br*hlayer
                presbcy = presbcy + cr*hlayer
                bpresxv(l,ie) = presbcx
                bpresyv(l,ie) = presbcy
                presbcx = presbcx + br*hlayer
                presbcy = presbcy + cr*hlayer
           end do
        end do
        
        rrho0=1./rowass
        !print*,'rowass ', rowass, grav,hhi 
        do ie=1,nel
            lmax=ilhv(ie)
                do l=1,lmax
                     hhi = hdeov(l,ie)
                     !print*, 'hhi ',hhi,bpresxv(l,ie),l,ie
                     xbcl =  rrho0*grav*hhi*bpresxv(l,ie)
                     ybcl =  rrho0*grav*hhi*bpresyv(l,ie)
                     fxv(l,ie) = fxv(l,ie) +xbcl
                     fyv(l,ie) = fyv(l,ie) +ybcl
                     !write(6,*)'fxv ',fxv,ie
                enddo
        enddo

        end

c**********************************************************************

        subroutine set_barocl_new

        implicit none
         
        include 'param.h'
	include 'ev.h'
        
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real fxv(nlvdim,1)      !new HYDRO deb
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
        real rhov(nlvdim,1)
        common /rhov/rhov
	real bpresv(nlvdim,1)
	common /bpresv/bpresv
        real hdeov(nlvdim,1)
        common /hdeov/hdeov
	real hldv(1)
	common /hldv/hldv
        integer ilhv(1)
        common /ilhv/ilhv
        integer ilmv(1)
        common /ilmv/ilmv
        integer nen3v(3,1)
        common /nen3v/nen3v

        integer k,l,ie,ii,lmax,lmin
        double precision hlayer,hhi
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr

        if(nlvdim.ne.nlvdi) stop 'error stop set_barocl_new: nlvdi'

        raux=grav/rowass

        do ie=1,nel
          presbcx = 0.
          presbcy = 0.
	  lmin = ilmv(ie)
          lmax = ilhv(ie)
          do l=1,lmax
            hhi = hdeov(l,ie)
            hhi = hldv(l)
            hlayer = 0.5 * hhi
                
	    br = 0.
	    cr = 0.                 
            do ii=1,3                 
              k = nen3v(ii,ie)
              rhop = rhov(l,k)		!rho^prime for each node of element 
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y
              br = br + b * rhop
              cr = cr + c * rhop
            end do

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer

            xbcl =  raux * hhi * presbcx
            ybcl =  raux * hhi * presbcy
            fxv(l,ie) = fxv(l,ie) + xbcl
            fyv(l,ie) = fyv(l,ie) + ybcl

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer
          end do
        end do
        
        end

c**********************************************************************

        subroutine set_barocl_new1

        implicit none
         
        include 'param.h'
	include 'ev.h'
        
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real fxv(nlvdim,1)      !new HYDRO deb
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
        real rhov(nlvdim,1)
        common /rhov/rhov
	real bpresv(nlvdim,1)
	common /bpresv/bpresv
        real hdeov(nlvdim,1)
        common /hdeov/hdeov
	real hldv(1)
	common /hldv/hldv
        integer ilhv(1)
        common /ilhv/ilhv
        integer ilmv(1)
        common /ilmv/ilmv
        integer nen3v(3,1)
        common /nen3v/nen3v

        integer k,l,ie,ii,lmax,lmin
        double precision hlayer,hhi
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr

	double precision px(0:nlvdim)
	double precision py(0:nlvdim)

        if(nlvdim.ne.nlvdi) stop 'error stop set_barocl_new: nlvdi'

        raux=grav/rowass

        do ie=1,nel
          presbcx = 0.
          presbcy = 0.
	  lmin = ilmv(ie)
          lmax = ilhv(ie)

	  px(0) = presbcx
	  py(0) = presbcy
          do l=1,lmax
            hhi = hdeov(l,ie)
            hhi = hldv(l)
            hlayer = 0.5 * hhi
            hlayer = hhi
                
	    br = 0.
	    cr = 0.                 
            do ii=1,3                 
              k = nen3v(ii,ie)
              rhop = rhov(l,k)		!rho^prime for each node of element 
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y
              br = br + b * rhop
              cr = cr + c * rhop
            end do

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer
	    px(l) = presbcx
	    py(l) = presbcy

          end do

          do l=1,lmax
	    presbcx = 0.5*(px(l) + px(l-1))
	    presbcy = 0.5*(py(l) + py(l-1))
            xbcl =  raux * hhi * presbcx
            ybcl =  raux * hhi * presbcy
            fxv(l,ie) = fxv(l,ie) + xbcl
            fyv(l,ie) = fyv(l,ie) + ybcl
	  end do
        end do
        
        end

c**********************************************************************
