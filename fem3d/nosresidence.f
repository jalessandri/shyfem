c
c $Id: nosresidence.f,v 1.9 2010-03-22 15:29:31 georg Exp $
c
c analyzes NOS file for residence time
c
c revision log :
c
c 05.03.2004    ggu     copied from nosextr
c 18.10.2006    ggu     heavily commented
c 13.02.2009    ggu     compute min/max of residence time
c 29.10.2009    ggu     compute and write 2D fields, use log to compute restime
c 12.01.2010    ggu     lots of changes
c
c****************************************************************

	program nosresidence

c analyzes CON file for residence time

	implicit none

	include 'param.h'
	include 'evmain.h'
	include 'basin.h'

c--------------------------------------------------
        character*80 descrr
        common /descrr/descrr
c--------------------------------------------------

	character*80 title
	real cv3(nlvdim,nkndim)
	real cv3d(nlvdim,nkndim)
	real cv2(nkndim)
	real cvres3(nlvdim,nkndim)
	real cvres2(nkndim)
	real cvrestot(nlvdim,nkndim)
        double precision cvacu(nlvdim,nkndim)

	integer ilhkv(nkndim)
	integer ilhv(neldim)
	real hev(neldim)
	real hlv(nlvdim)

	common /ilhkv/ilhkv
	common /ilhv/ilhv
	common /hev/hev
	common /hlv/hlv

	logical bminmax,balways,breset
	logical blog,badj
        integer l,lmax,k
	integer itstart,itend
	integer nin,nvers,nlv,nvar
	integer it,ivar,ierr,it0
	integer nread,nused,nrepl
	integer minused
        real conz,c0,conze,conzold
	real dt,ddt,secs_in_day,eps
	real xmin,xmax
	real rl,clog,cl
	real ctop,ccut
	real res
        double precision tacu

	integer iapini,ideffi

c---------------------------------------------------------------
c parameters to be changed
c---------------------------------------------------------------

c dt		time step used (output)
c c0		initial concentration
c it0		start of simulation (where conz = c0)
c bminmax	true  -> compute only between itstart and itend
c		false -> compute always
c breset	true  -> concentration gets reset during simulation
c itstart	start of computation if bminmax is true
c itend		end of computation if bminmax is true
c
c blog		use logarithmic regression to compute residence time
c badj		adjust residence time for tail of distribution
c
c ctop		maximum to be used for frequency curve
c ccut		cut residence time at this level (cor res time computation)

c        dt = 600
c        dt = 3600
c	c0 = 100.

        dt = 3600
        dt = 86400
	c0 = 1.
	it0 = 0
	ctop = 150.
	ccut = 150.
	ctop = 7000.
	ccut = 7000.
	ctop = 400.
	ccut = 500.
	minused = 5

	bminmax = .true.
	bminmax = .false.
	itstart = 86400
	itend = 30000000

	blog = .false.
	badj = .false.
	badj = .true.

c---------------------------------------------------------------
c do not change anything beyond this point
c---------------------------------------------------------------

	secs_in_day = 86400.
	eps = 0.3
	conzold = 1.e+30

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call set_ev

c---------------------------------------------------------------
c open input file and read headers
c---------------------------------------------------------------

	nin=ideffi('datdir','runnam','.con','unform','old')
	if(nin.le.0) goto 100

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

	call set_ilhv_post(nlv,nkn,nel,hlv,hev,nen3v,ilhkv,ilhv)

c---------------------------------------------------------------
c initialize variables and arrays
c---------------------------------------------------------------

        it = 0
        ivar = 99
        open(3,file='nosres.nos',status='unknown',form='unformatted')
        call wfnos(3,3,nkn,nel,nlv,1,title,ierr)
        call wsnos(3,ilhkv,hlv,hev,ierr)
        open(2,file='nosres2d.nos',status='unknown',form='unformatted')
        call wfnos(2,3,nkn,nel,1,1,title,ierr)
        call wsnos(2,ilhkv,hlv,hev,ierr)

	balways = .not. bminmax

	nread = 0
	nused = 0
	nrepl = 0

	call acu_reset(tacu,cvacu,ilhkv)

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cv3d(l,k) = 0.
            cvrestot(l,k) = 0.
          end do
        end do

        open(76,file='resi.txt',status='unknown',form='formatted')

c---------------------------------------------------------------
c time loop -> accumulate concentrations
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	!write(6,*) 'time : ',it,ivar

	if( balways .or. it .ge. itstart .and. it .le. itend ) then

	  nused = nused + 1

	  call acu_aver(conz,cv3,ilhkv,hlv)
          call mima3d(cv3d,nlvdim,nkn,ilhkv,xmin,xmax)
	  write(68,*) it,conz,xmin,xmax
	  if( conz - conzold .gt. eps ) then	!reset
	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing res time: ',it,conzold,conz,nused
	    write(6,*) '-------------------------------------------'
	    if( nused .gt. minused ) then
	      call acu_comp(blog,badj,it,dt,c0,ccut,ilhkv,hlv,tacu,cvacu
     +				,cv3d,cvres3,cvres2)
	      call acu_freq(it,ctop,ilhkv,cvres3)
	      nrepl = nrepl + 1
              do k=1,nkn
                lmax = ilhkv(k)
                do l=1,lmax
                  cvrestot(l,k) = cvrestot(l,k) + cvres3(l,k)
	        end do
	      end do
	    end if
	    call acu_reset(tacu,cvacu,ilhkv)
	    it0 = it
	    nused = 0
	  end if
	  !write(68,*) it,conz,conzold
	  conzold = conz

	  tacu = tacu + (it-it0)

          do k=1,nkn
            lmax = ilhkv(k)
            do l=1,lmax
              conz = cv3(l,k)
	      if ( conz .gt. 0. .and. conz .le. 1. ) then
		if( blog ) then
	          rl = - log(conz/c0)
                  cvacu(l,k) = cvacu(l,k) + rl
		else
                  cvacu(l,k) = cvacu(l,k) + conz
		end if
              else if( conz .lt. 0. ) then
		write(140,*) it,k,conz
		conz = 0.
      	      endif  
	      cv3d(l,k) = conz
            end do
          end do

	end if

	goto 300

c---------------------------------------------------------------
c end of time loop
c---------------------------------------------------------------

  100	continue

c---------------------------------------------------------------
c compute last residence time
c---------------------------------------------------------------

	    conz = 0.

	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing res time: ',it,conzold,conz,nused
	    write(6,*) '-------------------------------------------'

	    if( nused .gt. minused ) then
	      nrepl = nrepl + 1
	      call acu_comp(blog,badj,it,dt,c0,ccut,ilhkv,hlv,tacu,cvacu
     +				,cv3d,cvres3,cvres2)
	      call acu_freq(it,ctop,ilhkv,cvres3)
              do k=1,nkn
                lmax = ilhkv(k)
                do l=1,lmax
                  cvrestot(l,k) = cvrestot(l,k) + cvres3(l,k)
	        end do
	      end do
	    end if

	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing final res time: ',nrepl
	    write(6,*) '-------------------------------------------'

	    if( nrepl .gt. 0 ) then
	      it = 0
	      call acu_final(it,ilhkv,hlv,nrepl,cvrestot,cv2)
	      call acu_freq(it,ctop,ilhkv,cvrestot)
	    end if

c---------------------------------------------------------------
c finish up
c---------------------------------------------------------------

	write(6,*)
	write(6,*) 'parameters   dt = ',dt,'     c0 = ',c0
	write(6,*) nread,' records read'
	write(6,*) nused,' records used'
        write(6,*) 'levels    : ',lmax
	if( bminmax ) write(6,*) 'min/max time used : ',itstart,itend
	write(6,*)
	write(6,*) 'output written to nosres.nos and nosres2d.nos'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c***************************************************************

	subroutine average_residence_time(nlvdi,nkn,ilhkv,hlv,rl,res)

c averages residence time for whole basin

	implicit none

	integer nlvdi,nkn
	integer ilhkv(1)
	real hlv(1)
	real rl(nlvdi,1)
	real res

	include 'ev.h'

	integer k,lmax,l
	real htop,hbot,htot,h
	real ctot,cg,r
	double precision rtot

	rtot = 0.
        do k=1,nkn
          lmax = min(nlvdi,ilhkv(k))
          htop = 0.
          htot = 0.
          ctot = 0.
          do l=1,lmax
            hbot = hlv(l)
            h = hbot - htop
            cg = rl(l,k)
            ctot = ctot + h * cg
            htot = htot + h
            htop = hbot
          end do
          r = ctot / htot
	  rtot = rtot + r
        end do

	res = rtot / nkn

	end


c***************************************************************

	subroutine make_2D(nlvdi,nkn,ilhkv,hlv,cv3,cv2)

c convert 3D field to 2D

	implicit none

	integer nlvdi,nkn
	integer ilhkv(1)
	real hlv(1)
	real cv3(nlvdi,1)
	real cv2(1)

	integer k,lmax,l
	real htop,hbot,htot,h
	real ctot,cg

        do k=1,nkn
          lmax=ilhkv(k)
          htop = 0.
          htot = 0.
          ctot = 0.
          do l=1,lmax
            hbot = hlv(l)
            h = hbot - htop
            cg = cv3(l,k)
            ctot = ctot + h * cg
            htot = htot + h
            htop = hbot
          end do
          cv2(k) = ctot / htot
        end do

	end

c***************************************************************

        subroutine mima3d(xx,nlvdim,nkn,ilhkv,xmin,xmax)

c computes min/max of 3d matrix
c
c xx            matrix
c nlvdim	dimension of levels
c nkn           total number of nodes
c ilhkv		number of layers in node
c xmin,xmax     min/max value in matrix

	implicit none

	integer nlvdim,nkn
	real xx(nlvdim,nkn)
	integer ilhkv(nkn)
	real xmin,xmax

	integer k,l,lmax
	real c

	xmin = 1.e+30
	xmax = -xmin

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    c = xx(l,k)
	    xmin = min(xmin,c)
	    xmax = max(xmax,c)
	  end do
	end do

	end

c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

	do i=1,n
	  if(xx(i).ne.rnull) goto 1
	end do
    1	continue

	if(i.le.n) then
	  xmax=xx(i)
	  xmin=xx(i)
	else
	  xmax=rnull
	  xmin=rnull
	end if

	nmin=i+1

        do i=nmin,n
          x=xx(i)
	  if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
	  end if
        end do

        end

c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************

	subroutine acu_aver_elem(conz,cv)

	implicit none

	include 'param.h'
	include 'ev.h'

	real conz
	real cv(nlvdim,nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer nen3v(3,neldim)
	common /nen3v/nen3v
	integer ilhv(neldim)
	common /ilhv/ilhv

	logical bdebug
	integer k,lmax,l,ie,ii
	real h,area
	real hd(nlvdim)
	double precision cc,hh

	cc = 0.
	hh = 0.

        do ie=1,nel
	  bdebug = ie .eq. -1
          lmax = ilhv(ie)
	  area = 4. * ev(10,ie)
	  call get_elem_depth(ie,hd)
	if( bdebug ) write(6,*) ie,lmax,area
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      h = hd(l)
	      cc = cc + cv(l,k) * h * area
	      hh = hh + h * area
	if( bdebug ) write(6,*) l,k,h,cc,hh
	    end do
	  end do
	end do

	conz = cc / hh

	end

c***************************************************************

	subroutine acu_aver(conz,cv,ilhkv,hlv)

	implicit none

	include 'param.h'

	real conz
	real cv(nlvdim,nkndim)
	integer ilhkv(nkndim)
	real hlv(nlvdim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	logical busedep
	integer k,lmax,l,n
	real hup,hbot,h
	double precision cc,hh

	busedep = .false.
	busedep = .true.

	cc = 0.
	hh = 0.
	h = 1.
	n = 0

        do k=1,nkn
          lmax = ilhkv(k)
	  hup = 0.
          do l=1,lmax
	    hbot = hlv(l)
	    if( busedep ) h = hbot - hup
	    cc = cc + cv(l,k) * h
	    hh = hh + h
	    n = n + 1
	    hup = hbot
	  end do
	end do

	conz = cc / hh
	!conz = cc / n

	end

c***************************************************************

	subroutine acu_reset(tacu,cvacu,ilhkv)

	implicit none

	include 'param.h'

	double precision tacu
	double precision cvacu(nlvdim,nkndim)
	integer ilhkv(nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,lmax,l

	tacu = 0.
        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cvacu(l,k) = 0.
          end do
        end do

	end

c***************************************************************

	subroutine acu_comp(blog,badj,it,dt,c0,ccut,ilhkv,hlv,tacu,cvacu
     +				,cv3d,cv3,cv2)

c compute residence time

	implicit none

	include 'param.h'

	logical blog,badj
	integer it
	real dt
	real c0
	real ccut  !cut residence time at this level (cor res time computation)
	integer ilhkv(nkndim)
	real hlv(nlvdim)
	double precision tacu
	double precision cvacu(nlvdim,nkndim)
	real cv3d(nlvdim,nkndim)
	real cv3(nlvdim,nkndim)
	real cv2(nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,lmax,l,ivar,ierr
	real conz,conze,res,rese
	real xmin,xmax,cmax
	real secs_in_day,ddt

c---------------------------------------------------------------
c set parameters
c---------------------------------------------------------------

	secs_in_day = 86400.

	tacu = tacu / secs_in_day
	ddt = dt / secs_in_day

c---------------------------------------------------------------
c compute residence times
c---------------------------------------------------------------

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            conz = cvacu(l,k)
	    if( blog ) then
	      conz = tacu / conz
	    else
              conz = ddt * conz / c0		!convert to res time
	      if( badj ) then
	        conze = cv3d(l,k) / c0
	        if( conze .ge. 1 ) conze = 0.
	        if( conze .le. 0 ) conze = 0.
                conz = conz / ( 1. - conze )	!adjusted res time
	      end if
	    end if
	    if ( ccut .gt. 0. .and. conz .gt. ccut ) conz = ccut
            cv3(l,k) = conz
          end do
        end do

c---------------------------------------------------------------
c compute average residence times
c---------------------------------------------------------------

	write(6,*) 'basin wide residence times:'

	call acu_aver_elem(rese,cv3)
	call average_residence_time(nlvdim,nkn,ilhkv,hlv,cv3,res)
        call mima3d(cv3,nlvdim,nkn,ilhkv,xmin,xmax)
	write(6,*) ' (aver/min/max): ',res,rese,xmin,xmax
	cmax = xmax

	call make_2D(nlvdim,nkn,ilhkv,hlv,cv3,cv2)

	call average_residence_time(1,nkn,ilhkv,hlv,cv2,res)
        call mima3d(cv2,1,nkn,ilhkv,xmin,xmax)
	write(6,*) '  original 2D  (aver/min/max): ',res,xmin,xmax

c---------------------------------------------------------------
c write to file and terminal
c---------------------------------------------------------------

	ivar = 99

        call wrnos(3,it,ivar,nlvdim,ilhkv,cv3,ierr)
        call wrnos(2,it,ivar,1,ilhkv,cv2,ierr)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine acu_final(it,ilhkv,hlv,nrepl,cv3,cv2)

	implicit none

	include 'param.h'

	integer it
	integer ilhkv(nkndim)
	real hlv(nlvdim)
	integer nrepl
	real cv3(nlvdim,nkndim)
	real cv2(nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,lmax,l,ivar,ierr
	real conz,conze,res,rese
	real xmin,xmax,cmax

c---------------------------------------------------------------
c compute residence times
c---------------------------------------------------------------

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
	    cv3(l,k) = cv3(l,k) / nrepl
	  end do
	end do

c---------------------------------------------------------------
c compute average residence times
c---------------------------------------------------------------

	write(6,*) 'basin wide residence times (average):'

	call acu_aver_elem(rese,cv3)
	call average_residence_time(nlvdim,nkn,ilhkv,hlv,cv3,res)
        call mima3d(cv3,nlvdim,nkn,ilhkv,xmin,xmax)
	write(6,*) ' (aver/min/max): ',res,rese,xmin,xmax
	cmax = xmax

	call make_2D(nlvdim,nkn,ilhkv,hlv,cv3,cv2)

	call average_residence_time(1,nkn,ilhkv,hlv,cv2,res)
        call mima3d(cv2,1,nkn,ilhkv,xmin,xmax)
	write(6,*) '  original 2D  (aver/min/max): ',res,xmin,xmax

c---------------------------------------------------------------
c write to file and terminal
c---------------------------------------------------------------

	ivar = 99

        call wrnos(3,it,ivar,nlvdim,ilhkv,cv3,ierr)
        call wrnos(2,it,ivar,1,ilhkv,cv2,ierr)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine acu_freq(it,ctop,ilhkv,cv3)

	implicit none

	include 'param.h'

	integer it
	real ctop			!cut at this value of residence time
	integer ilhkv(nkndim)
	real cv3(nlvdim,nkndim)

	integer ndim
	parameter (ndim=100)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,lmax,l,i,ic
	integer icount(0:ndim)
	real conz,c,amax,dc,val,tot,cmax
	character*50 file

c---------------------------------------------------------------
c compute maximum
c---------------------------------------------------------------

	cmax = 0.
        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
	    cmax = max(cmax,cv3(l,k))
	  end do
	end do

	amax = cmax
	if( ctop .gt. 0. .and. amax .gt. ctop ) amax = ctop

c---------------------------------------------------------------
c compute frequency curve
c---------------------------------------------------------------

	ic = 0
	do i=0,ndim
	  icount(i) = 0
	end do

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
	    conz = cv3(l,k)
	    i = nint(ndim*conz/amax)
	    if (i .lt. 0) i = 0
	    if (i .gt. ndim) i = ndim
	    icount(i) = icount(i) + 1
	    ic = ic + 1
	  end do
	end do

c---------------------------------------------------------------
c write frequency curve to file
c---------------------------------------------------------------

	dc = 1.
	tot = 0.
	call make_name(it,file,'freq_by_bin_','.his')
	write(6,*) 'writing to file: ',file
	open(11,file=file,status='unknown',form='formatted')
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/float(ic)
	  tot = tot + val*dc
	  write(11,*) i,val,icount(i)
	end do
	close(11)
	write(6,*) 'writing finished: ',tot

	dc = amax/100.
	tot = 0.
	call make_name(it,file,'freq_by_res_','.his')
	write(6,*) 'writing to file: ',file
	open(11,file=file,status='unknown',form='formatted')
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/(float(ic)*dc)
	  tot = tot + val*dc
	  write(11,*) c,val,icount(i)
	end do
	close(11)
	write(6,*) 'writing finished: ',tot

c---------------------------------------------------------------
c write out all data to file (for debug and median)
c---------------------------------------------------------------

	write(76,*) nkn,it
        do k=1,nkn
	  conz = cv3(1,k)
	  write(76,*) conz
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine make_name(it,file,pre,ext)

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
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine set_ilhv_post(nlv,nkn,nel,hlv,hev,nen3v,ilhkv,ilhv)

c simplicistic

	implicit none

	integer nlv,nkn,nel
	real hlv(1)
	real hev(1)
	integer nen3v(3,1)
	integer ilhkv(1)
	integer ilhv(1)

	integer ie,ii,k,l
	integer lmin,lmax

	do ie=1,nel
	  lmin = nlv
	  do ii=1,3
	    k = nen3v(ii,ie)
	    lmin = min(lmin,ilhkv(k))
	  end do
	  l = 1
	  do while( hlv(l) .lt. hev(ie) )
	    l = l + 1
	  end do
	  lmax = l
	  if( lmax .gt. lmin ) lmax = lmin
	  ilhv(ie) = lmax
	end do

	end

c**********************************************************************

	subroutine get_elem_depth(ie,hd)

	implicit none

	include 'param.h'

	integer ie
	real hd(1)

	integer l,lmax
	real hmax,htop,hbot

	integer ilhkv(nkndim)
	integer ilhv(neldim)
	real hev(neldim)
	real hlv(nlvdim)

	common /ilhkv/ilhkv
	common /ilhv/ilhv
	common /hev/hev
	common /hlv/hlv

	lmax = ilhv(ie)
	hmax = hev(ie)

	htop = 0.
	do l=1,lmax
	  hbot = hlv(l)
	  if( l .eq. lmax ) hbot = hmax
	  hd(l) = hbot - htop
	  htop = hbot
	end do

	end

c**********************************************************************


