c
c $Id: lapwin.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 20.08.2003	ggu	new laplacian interpolation
c 02.09.2003	ggu	some comments, write to .dat file
c 30.10.2003	ggu	subroutine prepare_bc included in this file
c 31.10.2003	ggu	wind interpolation is finished
c
c notes :
c
c please prepare file like this:
c
c----------------- start
c
c time  n
c k1	valx1   valy1   valp1
c k2	valx2   valy2   valp2
c ...
c kn	valxn   valyn   valpn
c time  n
c k1	valx1   valy1   valp1
c k2	valx2   valy2   valp2
c ...
c kn	valxn   valyn   valpn
c----------------- end
c
c first line of file must be empty !!!
c in every record header time and number of nodes to read must be given
c in every record wx,wy and p (pressure) must be given (set p=0 if not needed)
c
c run memory and set the basin ( memory -b venlag62 )
c run lapwin with input file ( lapwin < input.dat )
c
c****************************************************************

        program lapwin

c laplacian interpolation for wind field

	implicit none

	include 'param.h'

	integer matdim
	parameter (matdim = nkndim*100)

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

	include 'evmain.h'

	real rmat(matdim)
	real zv(nkndim,3)
	real rzv(nkndim,3)

	integer k,i
        integer itime,ierr
        integer nvars
        integer iunit,nout
	real flag
	real zmin,zmax

	integer iapini

	flag = 1.23456e+23
        iunit = 5
        nvars = 3
        nout = 1

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn = ',nkn,'  nel = ',nel
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*)
        write(6,*) ' dcor = ',dcor,'  dirn = ',dirn
        write(6,*)

c-----------------------------------------------------------------
c check dimension
c-----------------------------------------------------------------

	if( nkn*(mbw+1) .gt. matdim ) then
	  write(6,*) nkn*(mbw+1),matdim
	  stop 'error stop laplap: matdim too small'
	end if

c-----------------------------------------------------------------
c open output file
c-----------------------------------------------------------------

        open(nout,file='wind.win',form='unformatted',status='unknown')

c-----------------------------------------------------------------
c set up ev
c-----------------------------------------------------------------

	call set_ev
	call check_ev

c-----------------------------------------------------------------
c read BC and interpolate
c-----------------------------------------------------------------

    1   continue
	  call prepare_bc(iunit,nvars,nkndim,nkn,rzv,flag,itime,ierr)
          if( ierr .ne. 0 ) goto 2
          do i=1,3
	    call lapint(rmat,zv(1,i),rzv(1,i),flag)
            call mima(zv(1,i),nkn,zmin,zmax)
            write(6,*) 'min/max: ',i,zmin,zmax
          end do
          call wrwin(nout,itime,nkn,zv(1,1),zv(1,2),zv(1,3))
          goto 1
    2   continue

        if( ierr .gt. 0 ) then
          write(6,*) 'ierr = ',ierr
          stop 'error stop lapwin: error in prepare_bc'
        end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c****************************************************************

	subroutine prepare_bc(iunit,ncols,nkndim,nkn,rzv,flag,itime,ierr)

c reads boundary conditions from file and sets up array
c
c file must be made like this:
c
c----------------- start
c
c time  nodes
c k1	valx1   valy1   valp1
c k2	valx2   valy2   valp2
c ...
c kn	valxn   valyn   valpn
c----------------- end

	implicit none

	integer iunit
	integer ncols
	integer nkndim
	integer nkn
	real rzv(nkndim,ncols)
	real flag
        integer itime
        integer ierr

	integer k,kn,i,n
        integer ndata
	real val(10)

	integer ipint

        if( ncols .gt. 10 ) then
          stop 'error stop prepare_bc: ncols'
        end if

	do i=1,ncols
	  do k=1,nkn
	    rzv(k,i) = flag
          end do
	end do

	read(iunit,*,end=100,err=91) itime,ndata
	write(6,*) itime,ndata

        do n=1,ndata
	  read(iunit,*) k,(val(i),i=1,ncols)
	  kn = ipint(k)
	  if( kn .le. 0 ) goto 99
	  if( kn .gt. nkn ) goto 98
          do i=1,ncols
	    rzv(kn,i) = val(i)
          end do
	  write(6,*) k,kn,ncols,(val(i),i=1,ncols)
        end do

	return
  100   continue
        ierr = -1
        return
   91   continue
        ierr = +1
        return
   98	continue
	write(6,*) k,kn,nkn,val
	stop 'error stop prepare_bc: error in internal node number'
   99	continue
	write(6,*) k,kn,nkn,val
	stop 'error stop prepare_bc: no such node'
	end

c******************************************************************
