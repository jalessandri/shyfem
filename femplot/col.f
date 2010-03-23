c
c $Id: col.f,v 1.8 2009-05-26 15:58:20 georg Exp $
c
c include file for colors
c
c revision log :
c
c 05.08.2003    ggu     changed adjcoltab, new cfast, init_cfast
c 20.08.2003    ggu     some more comments
c 21.08.2003    ggu     colini is called in colsetup
c 04.10.2004    ggu     adjust amin and decimal point in colauto
c 21.04.2009    ggu     allow for constant values (make 2 isolines)
c
c=====================================================================
c
c---------------------------------------------------------------------
c
c	integer isodim
c	parameter(isodim=256)
c
c        integer isopar,isoanz,niso,ncol
c        real fnull
c        real fiso(isodim)
c        real ciso(isodim+1)
c
c        common /isolin/ isopar,isoanz,niso,ncol
c        common /fsolin/ fnull
c        common /isofol/ fiso
c        common /isocol/ ciso
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c typical usage:
c
c $color
c         x0col = 0.1  y0col = 0.65
c         x1col = 0.45 y1col = 0.7
c  
c         legcol = 'Bathymetry [m]'
c         ndccol = -1  faccol = 1
c         icolor = 1
c  
c         colmin = 0   colmax = 0.666
c         valmin = 0   valmax = 300
c         niso = 50   nisomx = 70
c         nctick = 5
c $end
c
c other variables:
c
c	 isoval =    1.   3.   5.
c	 color  = 0.3  0.4  0.5  0.6
c
c**********************************************************

	subroutine colsetup

c sets up color table after read from $color section

	implicit none

	include 'color.h'

	logical bdebug
	integer niso,nisomx,nisodf
	integer icolor
	real colmin,colmax,valmin,valmax,dval
	logical bisord,bcolrd

	real getpar

	bdebug = .true.
	bdebug = .false.

c-------------------------------------------------
c initialize value and check for error
c-------------------------------------------------

	icauto = 0

	call colini

	colmin = getpar('colmin')	!min color [0..1]
	colmax = getpar('colmax')	!max color [0..1]
	valmin = getpar('valmin')	!min isovalue
	valmax = getpar('valmax')	!max isovalue
	dval   = getpar('dval')		!distance of isovalues
	niso   = nint(getpar('niso'))	!total number of isovalues to use
	nisomx = nint(getpar('nisomx'))	!maximum number of isovalues allowed
	nisodf = nint(getpar('nisodf'))	!default number of isovalues to use
	icolor = nint(getpar('icolor'))	!color table

	if( bdebug ) then
	  write(6,*) 'colsetup'
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	  write(6,*) 'nisord,ncolrd: ',nisord,ncolrd
	end if

	if( colmin .eq. colmax ) then
	  write(6,*) 'colmin/max : ',colmin,colmax
	  write(6,*) 'values must be different...'
	  stop 'error stop colsetup: colmin/max'
	end if

	bisord = nisord .gt. 0		!isovalues have been read
	bcolrd = ncolrd .gt. 0		!colors have been read

c-----------------------------------------------------
c see what has been read -> compute niso
c-----------------------------------------------------

	if( bisord .or. bcolrd ) then

c	  --------------------------------------------
c	  at least one type of values has been read
c	  --------------------------------------------

	  if( bisord .and. bcolrd ) then
	    if( nisord + 1 .ne. ncolrd ) then
		write(6,*) 'Number of isolines and colors read'
		write(6,*) '  are not compatible. There must be'
		write(6,*) '  one color more than isovalues.'
		write(6,*) 'nisord, ncolrd : ',nisord,ncolrd
		stop 'error stop colsetup: nisord/ncolrd'
	    end if
	  end if

c	  --------------------------------------------
c	  set niso
c	  --------------------------------------------

	  if( bisord ) niso = nisord
	  if( bcolrd ) niso = ncolrd - 1

c	  --------------------------------------------
c	  check computed number of isovalues
c	  --------------------------------------------

	  if( niso .le. 0 ) then
	    write(6,*) 'Cannot use only one color.'
	    stop 'error stop colsetup: niso'
	  else if( niso .eq. 1 ) then
	    if( .not. bisord .and. valmax .ne. valmin ) then
	      write(6,*) 'One isovalue only possible for valmin=valmax'
	      stop 'error stop colsetup: niso'
	    end if
	  end if

c	  --------------------------------------------
c	  if we cannot determine isovalues -> automatic
c	  --------------------------------------------

	  if( .not. bisord .and. valmax .eq. valmin ) icauto = 1

	else if( valmin .ne. valmax ) then

c	  --------------------------------------------
c	  isovalues can be determined
c	  --------------------------------------------

	  if( dval .ne. 0. ) then
	    niso = nint((valmax-valmin)/dval) + 1
	    if( niso .le. 0 ) then
	      write(6,*) 'Error computing niso.'
	      write(6,*) 'If dval < 0 then valmax < valmin.'
	      stop 'error stop colsetup: niso'
	    end if
	  else
	    if( niso .le. 0 ) niso = nisodf
	  end if

	else ! if( valmin .eq. valmax ) then

c	  --------------------------------------------
c	  here we can do only automatic evaluation
c	  --------------------------------------------

	  if( dval .ne. 0. ) then
	    niso = 0
c	  else
c	    use given value or 0
	  end if

	  icauto = 1

	end if

c-----------------------------------------------------
c check computed niso value and store
c-----------------------------------------------------

	if( niso .gt. nisomx ) then
	  write(6,*) 'Value of niso > nisomx.'
	  write(6,*) 'niso,nisomx : ',niso,nisomx
	  write(6,*) 'Either reduce iso-values or increase nisomx.'
	  stop 'error stop colsetup: nisomx'
	end if

	call putpar('niso',float(niso))
	isoanz = niso

c-----------------------------------------------------
c compute all other values that can be set up
c-----------------------------------------------------

	if( icauto .eq. 0 ) then

c	  --------------------------------------------
c	  no automatic evaluation -> set up everything
c	  --------------------------------------------

	  if( niso .le. 0 ) stop 'internal error colsetup (1)'

	  if( .not. bcolrd ) then
	    call mkval(niso+1,ciso,colmin,colmax,0.)
	    call adjcoltab(icolor,niso+1,ciso)
	  end if
	  if( .not. bisord ) call mkval(niso,fiso,valmin,valmax,dval)

	else ! if( icauto .ne. 0 ) then

c	  --------------------------------------------
c	  automatic evaluation -> set up as much as possible
c	  --------------------------------------------

	  if( niso .gt. 0 ) then

	    if( .not. bcolrd ) then
	      call mkval(niso+1,ciso,colmin,colmax,0.)
	      call adjcoltab(icolor,niso+1,ciso)
	    end if
	    if( bisord ) stop 'internal error colsetup (2)'
	    if( valmin .ne. valmax ) stop 'internal error colsetup (3)'

	  end if

	end if

c-----------------------------------------------------
c debug output
c-----------------------------------------------------

	if( bdebug ) then
	  write(6,*) 'colsetup: ',icauto
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	end if

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	end

c**********************************************************

	subroutine colauto(valmin,valmax)

c sets up color table in automatic mode

	implicit none

	real valmin,valmax		!min/max for array of values

	include 'color.h'

	logical bdebug
	integer i
	integer niso,nisomx,nisodf
	integer icolor
	real colmin,colmax,dval
	real amin,amax
        real fact,val,fdec
        integer ndec,idec

	real getpar
	real rnext,rnexta

	bdebug = .true.

	amin = 0.
	amax = 0.

c-------------------------------------------------
c check if we are in automatic mode
c-------------------------------------------------

	if( bdebug ) then
	  write(6,*) 'colauto: ',icauto
	  if( icauto .eq. 0 ) then
	    write(6,*) 'colors: ',isoanz+1
	    write(6,*) (ciso(i),i=1,isoanz+1)
	    write(6,*) 'isovals: ',isoanz
	    write(6,*) (fiso(i),i=1,isoanz)
	  end if
	end if

	if( icauto .lt. 0 ) stop 'error stop colauto: no initialization'
	if( icauto .eq. 0 ) return

c-------------------------------------------------
c initialize values
c-------------------------------------------------

	colmin = getpar('colmin')	!min color [0..1]
	colmax = getpar('colmax')	!max color [0..1]
	dval   = getpar('dval')		!distance of isovalues
	niso   = nint(getpar('niso'))	!total number of isovalues
	nisomx = nint(getpar('nisomx'))	!maximum number of isovalues allowed
	nisodf = nint(getpar('nisodf'))	!default number of isovalues to use
	icolor = nint(getpar('icolor'))	!color table

	if( bdebug ) then
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	end if

c-------------------------------------------------
c set up data structures
c-------------------------------------------------

	if( niso .gt. 0 ) then

c	  ----------------------------------------
c	  we already know how many isolines to draw (and color table)
c	  ----------------------------------------

	  dval = 0.
	  call mkval(niso,fiso,valmin,valmax,dval)

	else

c	  ----------------------------------------
c	  we must compute everything -> niso and dval first
c	  ----------------------------------------

	  if( dval .ne. 0. ) then
	    niso = nint((valmax-valmin)/dval) + 1
	  else			!find a good dval
	    niso = nisodf
	    dval = (valmax-valmin)/niso
	    if( dval .ne. 0. ) then
	      dval = rnext(dval,+3)	!round dval upwards
	      niso = nint((valmax-valmin)/dval) + 1
	    else
	      dval = 1.
	      niso = 1
	    end if
	  end if

c	  ----------------------------------------
c	  check niso if in range
c	  ----------------------------------------

	  if( niso .le. 0 ) then
	    write(6,*) 'Error computing niso.'
	    write(6,*) dval,valmin,valmax,niso
	    write(6,*) 'If dval < 0 then valmax < valmin.'
	    stop 'error stop colauto: niso'
	  else if( niso .gt. nisomx ) then
	    write(6,*) 'Value of niso > nisomx.'
	    write(6,*) 'niso,nisomx : ',niso,nisomx
	    write(6,*) 'Either increase dval or increase nisomx.'
	    stop 'error stop colauto: nisomx'
	  end if

c	  ----------------------------------------
c	  now we have a good niso and dval -> get amin
c	  ----------------------------------------

	  amin = valmin / dval
	  amin = rnexta(amin,+4)
	  amin = amin * dval
          do while( amin .gt. valmin )          !find starting point
            amin = amin - dval
          end do
          if( amin + niso*dval .lt. amax ) then
            amin = amin + dval
          end if
          amax = amin + (niso-1)*dval   !this is maximum isoline
	  amax = valmax

c	  ----------------------------------------
c	  special case - constant array
c	  ----------------------------------------

	  if( valmin .eq. valmax ) then
	    niso = 2
	    dval = valmax / 5.
	    if( dval .eq. 0 ) dval = 1.
	    amin = valmax - 0.5*dval
	    amax = amin + dval
	  end if

c	  ----------------------------------------
c	  adjust decimal point
c	  ----------------------------------------

          ndec = getpar('ndccol')
          fact = getpar('faccol')
          val = fact * dval
          fdec = log10(val)
          if( fdec .lt. 0. ) then
            idec = 1 + ifix(-fdec)      !1 for 0.1-0.9
            call putpar('ndccol',float(idec))
            write(6,*) 'decimal point adjusted: ',fact,dval,ndec,idec
          end if

c	  ----------------------------------------
c	  set up arrays
c	  ----------------------------------------

	  call mkval(niso,fiso,amin,amax,dval)
	  call mkval(niso+1,ciso,colmin,colmax,0.)
	  call adjcoltab(icolor,niso+1,ciso)

	end if

	isoanz = niso

	if( bdebug ) then
	  write(6,*) 'colauto final...'
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	  write(6,*) 'amin,amax: ',amin,amax
	  write(6,*) 'colors: ',niso+1
	  write(6,*) (ciso(i),i=1,niso+1)
	  write(6,*) 'isovals: ',niso
	  write(6,*) (fiso(i),i=1,niso)
	  write(6,*) 'isoanz: ',isoanz
	end if

	end

c**********************************************************

	subroutine mkval(n,array,amin,amax,da)

c makes array of values

	implicit none

	integer n
	real array(1)
	real amin,amax
	real da

	integer i
	real dval

	if( n .le. 0 ) then
	  return
	else if( n .eq. 1 ) then
	  dval = 0.
	else
	  dval = da
	  if( dval .eq. 0. ) dval = (amax-amin) / (n-1)
	end if

	do i=1,n
	  array(i) = amin + (i-1) * dval
	end do

	end

c**********************************************************

	subroutine adjcoltab(icolor,ncol,ciso)

c adjusts color table (only hsb and gray)

	implicit none

	integer icolor
	integer ncol
	real ciso(1)

	integer i
        integer n
        real x1,x2,y1,y2
        real rm,aux,x,y
        real rs
        real a(3,10)

        real cfast

	if( icolor .eq. 0 ) then
c	  gray original
	else if( icolor .eq. 1 ) then
c	  hsb original
	else if( icolor .eq. 2 ) then
	  do i=1,ncol
	    ciso(i) = ciso(i)**0.5
	  end do
	else if( icolor .eq. 3 ) then
          x1 = 0.2
          x2 = 0.4
          rm = 0.8
          y1 = rm*x1
          y2 = 1. - rm*(1.-x2)
          aux = (y2-y1)/(x2-x1)
	  do i=1,ncol
            x = ciso(i)
            if( x .le. x1 ) then
              y = rm * x
            else if( x .ge. x2 ) then
              y = 1. - rm*(1.-x)
            else
              y = y1 + (x-x1) * aux
            end if
	    ciso(i) = y
	  end do
	else if( icolor .eq. 4 ) then
          n = 2
          a(1,1) = 0.2
          a(2,1) = 0.4
          a(3,1) = 2
          a(1,2) = 0.6
          a(2,2) = 0.75
          a(3,2) = 2
          call init_cfast(n,a,rs)
          write(12,*) rs
	  do i=1,ncol
            x = ciso(i)
            y = cfast(x,n,a,rs)
	    ciso(i) = y
	  end do
	  do i=0,100
            x = i/100.
            y = cfast(x,n,a,rs)
            write(9,*) x,x,y
	  end do
        else
          write(6,*) 'icolor = ',icolor
          stop 'error stop adjcoltab: icolor'
	end if

	end

c**********************************************************

        function cfast(x,n,a,rs)

c computes rm from array a

        implicit none

        real cfast      !y value for x
        real x          !x value where y has to be computed
        integer n       !entries into a
        real a(3,1)     !x1,x2,rm for every region
        real rs         !derivative for regions not in a

        integer i
        real x0,y0
        real x1,x2,rm

        x0 = 0.
        y0 = 0.

        do i=1,n
          x1 = a(1,i)
          x2 = a(2,i)
          rm = a(3,i)
          if( x .le. x1 ) then
            cfast = y0 + (x-x0) * rs
            write(11,*) 1,x,x0,y0,x1,x2,rm
            return
          else if( x .le. x2 ) then
            y0 = y0 + (x1-x0) * rs
            cfast = y0 + (x-x1) * rm
            write(11,*) 2,x,x0,y0,x1,x2,rm
            return
          else
            y0 = y0 + (x1-x0) * rs
            y0 = y0 + (x2-x1) * rm
            x0 = x2
          end if
        end do

        cfast = y0 + (x-x0) * rs
        write(11,*) 3,x,x0,y0,x1,x2,rm

        end

c**********************************************************

        subroutine init_cfast(n,a,rs)

c computes rm from array a

        implicit none

        integer n       !entries into a
        real a(3,1)     !x1,x2,rm for every region
        real rs         !computed derivative for regions not in a

        integer i
        real x,y

        x = 0.
        y = 0.

        do i=1,n
          x = x + a(2,i) - a(1,i)
          y = y + a(3,i) * ( a(2,i) - a(1,i) )      !rm*(x2-x1)
        end do

        if( x .gt. 1. .or. y .gt. 1. ) then
          stop 'error stop init_cfast: cumulative value too high'
        end if

        rs = (1.-y) / (1.-x)

        end

c**********************************************************