c
c $Id: subdef.f,v 1.11 2008-04-01 11:58:12 georg Exp $
c
c create default names
c
c contents :
c
c function idefop(type,form,status)		  opens default file (outdated)
c subroutine mkname(dir,name,ext,file)            makes file name
c subroutine defmak(defdir,defnam,ext,file)	  makes file with defaults
c subroutine deffile(ext,file)			  makes file with fem defaults
c function ideffi(defdir,defnam,ext,form,status)  opens file in default dir
c function ifemop(ext,form,status)		  opens file with default name
c function ifemopa(text,ext,form,status)	  opens file with default name
c function ifem_open_file(ext,status)		  opens unformated file
c
c revision log :
c
c 23.05.1997	ggu     $$EXTENS - default extension may be overwritten
c 18.06.1997	ggu     restructured - idefna,idefop,idefts,idefun deleted
c 16.01.1998	ggu     idefop reintroduced -> to avoid link error
c 21.01.1998	ggu     in mkname: give extension with or without dot
c 08.08.2000	ggu     new routine ifemop
c 27.11.2001	ggu     error message rewritten
c 11.10.2002	ggu	new subroutine deffile
c 07.03.2007	ggu	new routine ifem_open_file
c
c****************************************************************
c
c	function idefna(type,name)
c
c basin files
c
c	call getfnm('basdir',dir) call getfnm('basnam',nam) ext='.bas'
c	call getfnm('basdir',dir) call getfnm('basnam',nam) ext='.geo'
c	call getfnm('basdir',dir) call getfnm('basnam',nam) ext='.dep'
c	call getfnm('basdir',dir) call getfnm('basnam',nam) ext='.bnd'
c	call getfnm('basdir',dir) call getfnm('basnam',nam) ext='.isl'
c
c run output files
c
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.out'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.ous'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.ext'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.flt'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.fdf'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.rst'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.res'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.elv'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.nov'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.con'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.rms'
c	call getfnm('datdir',dir) call getfnm('basnam',nam) ext='.el0'
c
c temporary files
c
c	call getfnm('tmpdir',dir) call getfnm('runnam',nam) ext='.cpu'
c	call getfnm('tmpdir',dir) call getfnm('runnam',nam) ext='.prg'
c	call getfnm('tmpdir',dir) call getfnm('runnam',nam) ext='.err'
c	call getfnm('datdir',dir) call getfnm('runnam',nam) ext='.inf'
c	call getfnm('runnam',nam) ext='.66'
c
c**************************************************************

	function idefop(type,form,status)

c opens default file -> should not be used any more

	implicit none

	integer idefop
	character*(*) type,form,status

	character*4 ext
	character*6 defdir,defnam

	integer ideffi

	defdir = 'datdir'
	defnam = 'runnam'

	if( type .eq. 'bas' ) then
	  defdir = 'basdir'
	  defnam = 'runnam'
	end if

	ext(1:1) = '.'
	ext(2:) = type

        idefop = ideffi(defdir,defnam,ext,form,status)

	end

c**************************************************************

        subroutine mkname(dir,name,ext,file)

c makes file name given its constituents
c
c dir   directory
c name  name
c ext   extension (with or without dot)
c file  created file name (return)

        implicit none

c arguments
        character*(*) dir,name,ext,file
c local
        integer nall,nstart,nend,naux
c function
        integer ichafs,ichanm

        nall=1
        file=' '

        nstart=ichafs(dir)
        nend=ichanm(dir)
        if(nend.gt.0) then
		file(nall:)=dir(nstart:nend)
        	nall=nall+nend-nstart+1
	end if

        nstart=ichafs(name)
        nend=ichanm(name)
        if(nend.gt.0) then
		file(nall:)=name(nstart:nend)
       		nall=nall+nend-nstart+1
	end if

c new code: if extension given in name do not add default extension
c		extension MUST have 3 chars	$$EXTENS
c FIXME

	naux = nall - 4
	if( naux .gt. 0 .and. file(naux:naux) .eq. '.' ) return

        nstart=ichafs(ext)
        nend=ichanm(ext)

	if( nend .gt. 0 ) then
	   if( ext(nstart:nstart) .ne. '.' ) then !add dot if not in ext
		file(nall:nall) = '.'
		nall = nall + 1
	   end if
	   file(nall:)=ext(nstart:nend)
       	   nall=nall+nend-nstart+1
	end if

	end

c**************************************************************

        subroutine defmak(defdir,defnam,ext,file)

c makes file with defaults supplied
c
c defdir   directory
c defnam  name
c ext   extension (with dot)
c file  created file name (return)

        implicit none

        character*(*) defdir,defnam,ext,file
	character*80 dir,nam

	dir = ' '
        if( defdir .ne. ' ' ) call getfnm(defdir,dir)
        call getfnm(defnam,nam)

	call mkname(dir,nam,ext,file)

	end

c**************************************************************

        subroutine deffile(ext,file)

c makes file with defaults for fem model supplied
c
c ext   extension (with dot)
c file  created file name (return)

        implicit none

        character*(*) ext,file

	call defmak('datdir','runnam',ext,file)

	end

c**************************************************************

        function ideffi(defdir,defnam,ext,form,status)

c opens file in default dir

c defdir   directory
c defnam  name
c ext   extension (with dot)
c form  formatted ?
c status open status

        implicit none

	integer ideffi
        character*(*) defdir,defnam,ext,status,form
	character*80 file
	integer ifileo

	call defmak(defdir,defnam,ext,file)
	ideffi=ifileo(0,file,form,status)

	end

c**************************************************************

        function ifemop(ext,form,status)

c opens file with default name (run) and extension given for fem model

c ext   extension (with dot)
c form  formatted ?
c status open status

        implicit none

	integer ifemop
        character*(*) ext,status,form

	character*80 file,defdir,defnam
	integer ifileo

	defdir = ' '
	defnam = 'runnam'
	call defmak(defdir,defnam,ext,file)
	ifemop=ifileo(0,file,form,status)

	end

c**************************************************************

        function ifemopa(text,ext,form,status)

c opens file with default name (run) and extension given for fem model
c in case of error exits with error message text

c text  error message
c ext   extension (with dot)
c form  formatted ?
c status open status

        implicit none

	integer ifemopa
        character*(*) text,ext,status,form

	character*80 file,defdir,defnam
	integer ifileo

	defdir = ' '
	defnam = 'runnam'
	call defmak(defdir,defnam,ext,file)
	ifemopa=ifileo(0,file,form,status)

	if( ifemopa .le. 0 ) then
	  write(6,*) 'error opening file ',file
	  write(6,*) text
	  stop 'error stop ifemopa'
	end if

	end

c**************************************************************

        function ifem_open_file(ext,status)

c opens unformated file with default name (run) and extension given
c in case of error exits 

c ext   extension (with dot)
c status open status

        implicit none

	integer ifem_open_file
        character*(*) ext,status

	character*80 file,defdir,defnam
	character*80 form
	integer ifileo

	form = 'unform'
	defdir = ' '
	defnam = 'runnam'
	call defmak(defdir,defnam,ext,file)
	ifem_open_file = ifileo(0,file,form,status)

	if( ifem_open_file .le. 0 ) then
	  write(6,*) 'error opening file ',file
	  stop 'error stop ifem_open_file'
	end if

	end

c**************************************************************
