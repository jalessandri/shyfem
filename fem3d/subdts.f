c
c $Id: subdts.f,v 1.9 2009-05-29 15:40:42 georg Exp $
c
c routines for date and time
c
c revision log :
c
c 01.12.2003    ggu     new routine month_from_name(), make_lower_case()
c 09.03.2004    ggu     new routine month_name()
c 22.09.2004    ggu     rewritten, integrated file newdat.f
c 05.10.2004    ggu     dtsadj got some FIXME
c 11.03.2005    ggu     new test added
c 13.11.2005    ggu     small bug fix (UNPACK)
c 01.12.2005    ggu     new routine dts_initialized and block data
c 07.05.2009    ggu     new routine dtsyear() and date_compute()
c
c notes :
c
c important routines:
c
c subroutine dtsini(date,time)		    initializes dts routines
c subroutine dtsyear(year)		    initializes dts routines (for year)
c subroutine dts2dt(it,year,month,day,hour,min,sec) converts it to date and time
c subroutine dts2it(it,year,month,day,hour,min,sec) converts date and time to it

c************************************************************************
c************************************************************************
c************************************************************************
c
c packing routines
c
c************************************************************************
c************************************************************************
c************************************************************************

	subroutine packdate(date,year,month,day)

c packs year, month, day into date [YYYYMMDD]

	implicit none

	integer date
        integer year,month,day

	date = day + 100*month + 10000*year

	end

c************************************************************************

	subroutine unpackdate(date,year,month,day)

c splits date [YYYYMMDD] into year, month, day

	implicit none

	integer date
        integer year,month,day

        month = date / 100
        day = date - month * 100

        year = month / 100
        month = month - year * 100

	end

c************************************************************************
 
	subroutine packtime(time,hour,min,sec)

c packs hour, min, sec into time [HHMMSS]

	implicit none

	integer time
        integer hour,min,sec

	time = sec + 100*min + 10000*hour

	end

c************************************************************************

	subroutine unpacktime(time,hour,min,sec)

c splits time [HHMMSS] into hour, min, sec

	implicit none

	integer time
        integer hour,min,sec

        min = time / 100
        sec = time - min * 100

        hour = min / 100
        min = min - hour * 100

	end

c************************************************************************
c************************************************************************
c************************************************************************
c
c dts routines
c
c************************************************************************
c************************************************************************
c************************************************************************

	subroutine dtsgf(it,line)

c formats date and time given it (secs)

	implicit none

	integer it
	character*(*) line

	integer year,month,day,hour,min,sec

	call dts2dt(it,year,month,day,hour,min,sec)
	call dtsform(year,month,day,hour,min,sec,line)

	end

c************************************************************************

	subroutine dtsform(year,month,day,hour,min,sec,line)

c formats date and time

	implicit none

	integer year,month,day,hour,min,sec
	character*(*) line

	character*20 lineaux
	integer i

        lineaux = ' '

	if( year .gt. 0 ) then
	  write(lineaux,1000) year,month,day,hour,min,sec
	  do i=1,20
	    if( lineaux(i:i) .eq. ' ' ) lineaux(i:i) = '0'
	  end do
	  line = lineaux
	else
	  write(lineaux,2000) day,hour,min,sec
	  do i=11,20
	    if( lineaux(i:i) .eq. ' ' ) lineaux(i:i) = '0'
	  end do
	  line = lineaux
	end if

	return
 1000	format(i4,1h-,i2,1h-,i2,2h::,i2,1h:,i2,1h:,i2)
 2000	format(i10,2h::,i2,1h:,i2,1h:,i2)
	end

c************************************************************************
c************************************************************************
c************************************************************************

	subroutine adjmod(low,high,modulo)

c re-distributes values (internal routine)

	implicit none

	integer low,high,modulo

	integer ia

	if( low .ge. 0 ) then
	  ia = low/modulo
	  low = low - ia*modulo
	  high = high + ia
	else
	  ia = low/modulo
	  if( mod(low,modulo) .ne. 0 ) ia = ia - 1
	  low = low - ia*modulo
	  high = high + ia
	end if

	end

c************************************************************************

	subroutine dtsadj(year,month,day,hour,min,sec)

c adjusts date and time (internal routine)

	implicit none

        integer year,month,day
        integer hour,min,sec

	integer ia,jd
	logical bdebug

	integer days_in_year

	bdebug = .true.
	bdebug = .false.

        if( bdebug ) then
          write(6,*) 'enter: ',year,month,day,hour,min,sec
        end if

c---------------------------------------------------------
c adjusts seconds
c---------------------------------------------------------

	if( sec .lt. 0 .or. sec .ge. 60 ) then
	  if( bdebug ) write(6,*) day,hour,min,sec
	  call adjmod(sec,min,60)
	  if( bdebug ) write(6,*) day,hour,min,sec,ia
	end if

c---------------------------------------------------------
c adjusts minutes
c---------------------------------------------------------

	if( min .lt. 0 .or. min .ge. 60 ) then
	  if( bdebug ) write(6,*) day,hour,min,sec
	  call adjmod(min,hour,60)
	  if( bdebug ) write(6,*) day,hour,min,sec,ia
	end if

c---------------------------------------------------------
c adjusts hours
c---------------------------------------------------------

	if( hour .lt. 0 .or. hour .ge. 24 ) then
	  if( bdebug ) write(6,*) day,hour,min,sec
	  call adjmod(hour,day,24)
	  if( bdebug ) write(6,*) day,hour,min,sec,ia
	end if

        if( bdebug ) then
          write(6,*) 'hours: ',year,month,day,hour,min,sec
        end if

	if( year .eq. 0 ) return

c---------------------------------------------------------
c adjusts negative days
c---------------------------------------------------------

	if( day .le. 0 ) then
          call date2j(year,month,1,jd)    !FIXME
	  day = day + jd - 1
	  month = 1
	  if( bdebug ) write(6,*) 'negative: ',jd,day,month,year
	end if

	do while ( day .le. 0 )
	  year = year - 1
	  day = day + days_in_year(year)
	  if( bdebug ) write(6,*) 'negative: ',jd,day,month,year
	end do

        if( bdebug ) then
          write(6,*) 'days: ',year,month,day,hour,min,sec
        end if

c---------------------------------------------------------
c go back to beginning of year
c---------------------------------------------------------

        call date2j(year,month,1,jd)    !FIXME
	day = day + jd - 1
	month = 1
	if( bdebug ) write(6,*) 'positive: ',jd,day,month,year

	do while( day .gt. days_in_year(year) )
	  day = day - days_in_year(year)
	  year = year + 1
	  if( bdebug ) write(6,*) 'positive: ',jd,day,month,year
	end do

        if( bdebug ) then
          write(6,*) 'jdays: ',year,month,day,hour,min,sec
        end if

c---------------------------------------------------------
c adjust days finally
c---------------------------------------------------------

	jd = day
        call j2date(jd,year,month,day)
	if( bdebug ) write(6,*) 'final: ',jd,day,month,year

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c************************************************************************

	subroutine dts2dt(it,year,month,day,hour,min,sec)

c converts it to date and time

	implicit none

	integer it
        integer year,month,day
        integer hour,min,sec

        integer year0,month0,day0
        integer hour0,min0,sec0
	integer last
	common /dtsdts/ year0,month0,day0,hour0,min0,sec0,last
	save /dtsdts/

	last = it

	year=year0
	month=month0
	day=day0
	hour=hour0
	min=min0
	sec=sec0

	sec = sec + it

	call dtsadj(year,month,day,hour,min,sec)

	end

c************************************************************************

	subroutine dts2it(it,year,month,day,hour,min,sec)

c converts date and time to it

	implicit none

	integer it
        integer year,month,day
        integer hour,min,sec

	integer days,days0,secs,secs0

        integer year0,month0,day0
        integer hour0,min0,sec0
	integer last
	common /dtsdts/ year0,month0,day0,hour0,min0,sec0,last
	save /dtsdts/

	call date2days(days,year,month,day)
	call date2days(days0,year0,month0,day0)

	!write(6,*) days,days0

	secs = 3600 * hour + 60 * min + sec
	secs0 = 3600 * hour0 + 60 * min0 + sec0

	it = ( days - days0 ) * 86400 + secs - secs0

	end

c************************************************************************

	subroutine dtsini(date,time)

c sets date and time for 0 fem time

	implicit none

	integer date,time
	logical bdebug

        integer year0,month0,day0
        integer hour0,min0,sec0
	integer last
	common /dtsdts/ year0,month0,day0,hour0,min0,sec0,last
	save /dtsdts/

	bdebug = .true.
	bdebug = .false.

	call unpackdate(date,year0,month0,day0)
	!call unpackdate(time,hour0,min0,sec0)	!UNPACK
	call unpacktime(time,hour0,min0,sec0)

	if( bdebug ) then
	  write(6,*) '-------- dtsini -------------'
	  write(6,*) year0,month0,day0,hour0,min0,sec0
	  write(6,*) '-----------------------------'
	end if

	last = 0

	end

c************************************************************************

	subroutine dtsyear(year)

c initialization by only giving year

        implicit none

        integer year

        integer date,time

        date = 10000*year + 101
        time = 0

	call dtsini(date,time)

        end

c************************************************************************

	function dts_initialized()

c checks if dts routines are already initialized

	implicit none

	logical dts_initialized

        integer year0,month0,day0
        integer hour0,min0,sec0
	integer last
	common /dtsdts/ year0,month0,day0,hour0,min0,sec0,last
	save /dtsdts/

	if( year0 .eq. 0 ) then
	  dts_initialized = .false.
	else
	  dts_initialized = .true.
	end if

	end

c************************************************************************
c
c common block /dtsdts/ is initialized in block data statement
c
c************************************************************************
c************************************************************************
c************************************************************************
c
c convert days and seconds
c
c************************************************************************
c************************************************************************
c************************************************************************

        subroutine hms2secs(secs,hour,min,sec)

c packs hour, min, sec into secs

        implicit none

        integer secs
        integer hour,min,sec

        secs = sec + 60 * ( min + 60 * hour )

        end

c************************************************************************

        subroutine secs2hms(secs,hour,min,sec)

c splits secs into hour, min, sec

        implicit none

        integer secs
        integer hour,min,sec

        min = secs / 60
        sec = secs - min * 60

        hour = min / 60
        min = min - hour * 60

        end

c************************************************************************
c************************************************************************
c************************************************************************
 
        subroutine date2days(days,year,month,day)
 
c given date returns total days from 1/1/1
 
	implicit none

	integer days
        integer year,month,day

	integer iy,jd
 
        iy=year-1
        days=365*iy+iy/4-iy/100+iy/400

c other possibility to sum up years (saver but slower)
c       days=0
c       do i=1,iyear-1
c         days=days+days_in_year(i)
c	end do

	call date2j(year,month,day,jd)

	days = days + jd
 
        end
 
c************************************************************************

        subroutine days2date(days,year,month,day)

c splits days from 1/1/1 into year, month, day

        implicit none

        integer year,month,day
        integer days

        integer iy,auxdays,days1
        integer idmon

        days1 = days - 1

c       estimate year (using average days in 400 years)

        iy = (days*400) / (365*400+97) + 2

c       correct year

        auxdays = days1 + 1
        do while( auxdays .gt. days1 )
          iy = iy - 1
          auxdays = 365*iy + iy/4 - iy/100 + iy/400
        end do

c       year found, auxdays days left

        year = iy + 1
        auxdays = days1 - auxdays

c       find month

        month = 0
        do while( auxdays .ge. 0 )
          month = month + 1
          auxdays = auxdays - idmon(year,month)
        end do

c       find day

        day = auxdays + idmon(year,month) + 1

        end

c************************************************************************
c************************************************************************
c************************************************************************
c
c julian date routines
c
c************************************************************************
c************************************************************************
c************************************************************************

        subroutine date2j(year,month,day,jd)
 
c converts date to julian days

        implicit none
 
        integer year,month,day
	integer jd

        integer jdmon,idmon
 
	if( month .lt. 1 .or. month .gt. 12 ) then
	  write(6,*) year,month,day
	  stop 'error stop date2j: month'
        else if( day .lt. 1 .or. day .gt. idmon(year,month) ) then
	  write(6,*) year,month,day
	  stop 'error stop date2j: day'
	end if

	jd = day + jdmon(year,month-1)

	end

c************************************************************************

        subroutine j2date(jd,year,month,day)
 
c given julian days and year returns month and day

        implicit none
 
	integer jd
        integer year,month,day
 
        integer jdmon,days_in_year
 
        if( jd .lt. 1 .or. jd .gt. days_in_year(year) ) then
	  write(6,*) year,jd,days_in_year(year)
	  stop 'error stop j2date: jd'
        end if

	month = jd / 30
	if( jd .le. jdmon(year,month) ) month = month - 1
	day = jd - jdmon(year,month)
	month = month + 1
 
        end
 
c************************************************************************

        function jdmon(year,month)

c returns total days from January to end of month

        implicit none

        integer jdmon
        integer year,month

        logical bises

        integer jmm(0:12)
        save jmm
        data jmm /0,31,59,90,120,151,181,212,243,273,304,334,365/

        jdmon=jmm(month)

        if( month .ge. 2 .and. bises(year) ) jdmon = jdmon + 1

        end

c************************************************************************

        function idmon(year,month)

c returns days in month

        implicit none

        integer idmon
        integer year,month

        logical bises

        integer imm(0:12)
        save imm
        data imm /0,31,28,31,30,31,30,31,31,30,31,30,31/

        idmon=imm(month)

        if( month .eq. 2 .and. bises(year) ) idmon = idmon + 1

        end

c************************************************************************
 
        function bises(year)
 
c true if year is bisestial
 
        logical bises
        integer year
 
        if(  ( mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0 )
     +                  .or. mod(year,400) .eq. 0 ) then
                bises = .true.
        else
                bises = .false.
        end if
 
        end
 
c************************************************************************
 
        function days_in_year(year)
 
c returns total days in year
 
	integer days_in_year
        integer year

	logical bises
 
        if( bises(year) ) then
	  days_in_year = 366
        else
	  days_in_year = 365
        end if
 
        end
 
c************************************************************************
c************************************************************************
c************************************************************************
c
c month name routines
c
c************************************************************************
c************************************************************************
c************************************************************************

        subroutine month_name(im,mname)
 
c returns month name with im given (1-12)
 
	implicit none

        integer im              ![1-12]
        character*(*) mname

        integer lang            !1=it  2=de  else=en
        parameter (lang=0)

        character*3 months_en(12),months_it(12),months_de(12)
        save months_en,months_it,months_de
        data months_en
     +          /'Jan','Feb','Mar','Apr','May','Jun'
     +          ,'Jul','Aug','Sep','Oct','Nov','Dec'/
        data months_it
     +          /'Gen','Feb','Mar','Apr','Mag','Giu'
     +          ,'Lug','Ago','Set','Ott','Nov','Dic'/
        data months_de
     +          /'Jan','Feb','Mar','Apr','Mai','Jun'
     +          ,'Jul','Aug','Sep','Okt','Nov','Dez'/

        if( im .lt. 1 .or. im .gt. 12 ) then
          write(6,*) 'error in month : ',im
          stop 'error stop monthname'
        end if

        if( lang .eq. 1 ) then
          mname = months_it(im)
        else if( lang .eq. 2 ) then
          mname = months_de(im)
        else
          mname = months_en(im)
        end if

        end

c************************************************************************

        function month_from_name(month_name)
 
c returns month number from name [1-12] , 0 for error
 
	implicit none

        integer month_from_name
        character*(*) month_name

        character*3 name
        integer i

        character*3 months(12)
        save months
        data months
     +          /'jan','feb','mar','apr','may','jun'
     +          ,'jul','aug','sep','oct','nov','dec'/

        name = month_name
        call make_lower_case(name)

        do i=1,12
          if( months(i) .eq. name ) then
            month_from_name = i
            return
          end if
        end do

        month_from_name = 0

        return
        end

c************************************************************************

        subroutine make_lower_case(name)

        implicit none

        character*(*) name

        integer upper_start,upper_end,lower_start,lower_end
        integer n,i,ic

        upper_start = ichar('A')
        upper_end = ichar('Z')
        lower_start = ichar('a')
        lower_end = ichar('z')

        n = len(name)
        do i=1,n
          ic = ichar(name(i:i))
          if( ic .ge. upper_start .and. ic .le. upper_end ) then
            ic = ic + lower_start - upper_start
          end if
          name(i:i) = char(ic)
        end do

        end

c************************************************************************
c************************************************************************
c************************************************************************

        block data dts_blockdata

        implicit none

        integer year0,month0,day0
        integer hour0,min0,sec0
	integer last

	common /dtsdts/ year0,month0,day0,hour0,min0,sec0,last
	save /dtsdts/

        data year0,month0,day0,hour0,min0,sec0,last /0,0,0,0,0,0,0/

        end

c************************************************************************
c************************************************************************
c************************************************************************
c
c test routines
c
c************************************************************************
c************************************************************************
c************************************************************************

	subroutine datetest

c tests date routines

	implicit none

	integer date,time
        integer year,month,day
	integer hour,min,sec

	character*40 line
	integer jd,days,it
	integer days_in_year

	date = 19970823
	time = 101030

	call unpackdate(date,year,month,day)
	write(6,*) date
	write(6,*) year,month,day
	date = 0
	call packdate(date,year,month,day)
	write(6,*) date

	call unpacktime(time,hour,min,sec)
	write(6,*) time
	write(6,*) hour,min,sec
	time = 0
	call packtime(time,hour,min,sec)
	write(6,*) time

	write(6,*)
	write(6,*) '*** ',year,month,day,hour,min,sec
	sec = sec + 10805
	call dtsadj(year,month,day,hour,min,sec)
	write(6,*) '*** ',year,month,day,hour,min,sec

	write(6,*)
	write(6,*) '*** ',year,month,day,hour,min,sec
	sec = sec - 3715
	call dtsadj(year,month,day,hour,min,sec)
	write(6,*) '*** ',year,month,day,hour,min,sec

	write(6,*)
	write(6,*) '*** ',year,month,day,hour,min,sec
	sec = sec - 82800
	call dtsadj(year,month,day,hour,min,sec)
	write(6,*) '*** ',year,month,day,hour,min,sec

	year=1996
	days = days_in_year(year)
	days = 0
	do jd=1,days
	  call j2date(jd,year,month,day)
	  write(6,*) jd,year,month,day
	end do

	date = 19980201
	date = 19980801
	time = 0
	it = 0
	call dtsini(date,time)
	call dts2dt(it,year,month,day,hour,min,sec)
	call dtsform(year,month,day,hour,min,sec,line)
	write(6,'(a)') line
	write(6,*)

	it = -2678400

	do while( it .le. 86400 )
	  call dts2dt(it,year,month,day,hour,min,sec)
	  call dtsform(year,month,day,hour,min,sec,line)
	  write(6,'(a)') line
	  it = it + 3600
	end do

	date = 0
	call dtsini(date,time)

	it = -86400
	do while( it .le. 86400 )
	  call dts2dt(it,year,month,day,hour,min,sec)
	  call dtsform(year,month,day,hour,min,sec,line)
	  write(6,'(a)') line
	  it = it + 3600
	end do

	end

c************************************************************************

        subroutine check_j2date(jd,year,iprint)
        integer jdd,jd,year,month,day,iprint
        call j2date(jd,year,month,day)
        call date2j(year,month,day,jdd)
        if( iprint .ne. 0 ) write(6,*) year,month,day,jd,jdd
        if( jd .ne. jdd ) write(6,*) '*** error : ',jd,jdd
        end

c************************************************************************

        subroutine test_j2date

        implicit none
        integer year,j,i

        year = 2001
        write(6,*) 'checking j2date for complete year ',year
        do j=1,365
          call check_j2date(j,year,0)
        end do
        year = 2004
        write(6,*) 'checking j2date for complete year ',year
        do j=1,366
          call check_j2date(j,year,0)
        end do

        do i=3,4
          year = 2000 + i
          write(6,*) 'checking j2date for single days in year ',year
          call check_j2date(1,year,1)
          call check_j2date(29,year,1)
          call check_j2date(30,year,1)
          call check_j2date(31,year,1)
          call check_j2date(59,year,1)
          call check_j2date(60,year,1)
          call check_j2date(61,year,1)
          call check_j2date(89,year,1)
          call check_j2date(90,year,1)
          call check_j2date(91,year,1)
          call check_j2date(365,year,1)
        end do

        end

c************************************************************************

        subroutine check_secs_pack(secs,iprint)
        integer secs,secs0,iprint
        integer hour,min,sec
        call secs2hms(secs,hour,min,sec)
        call hms2secs(secs0,hour,min,sec)
        if( secs .ne. secs0 ) then
          write(6,*) '*** error : ',hour,min,sec,secs,secs0
        end if
        if( iprint .ne. 0 ) write(6,*) secs,hour,min,sec
        end

c************************************************************************

        subroutine check_days_pack(days,iprint)
        integer days,days0,iprint
        integer year,month,day
        call days2date(days,year,month,day)
        call date2days(days0,year,month,day)
        if( days .ne. days0 ) then
          write(6,*) '*** error : ',year,month,day,days,days0
        end if
        if( iprint .ne. 0 ) write(6,*) days,year,month,day
        end

c************************************************************************

        subroutine test_sd_pack

        implicit none

        integer j,year,month,day
        integer iprint

        write(6,*) 'checking packing of seconds of whole day'
        do j=1,86400
          iprint = 0
          if( mod(j,7200) .eq. 0 ) iprint = 1
          call check_secs_pack(j,iprint)
        end do

        write(6,*) 'checking packing of days from year 1-2004'
        call check_days_pack(1,1)
        do j=1,365*2005
          iprint = 0
          if( mod(j,30000) .eq. 0 ) iprint = 1
          call check_days_pack(j,iprint)
        end do
        call check_days_pack(j,1)

        end

c************************************************************************

        subroutine check_dpack(year,month,day,iprint)
        integer year,month,day,iprint
        integer year0,month0,day0
        integer date
        call packdate(date,year,month,day)
        call unpackdate(date,year0,month0,day0)
        if( year .ne. year0 .or. month .ne. month0 
     +          .or. day .ne. day0 ) then
          write(6,*) '*** error : ',year,month,day,year0,month0,day0
        end if
        if( iprint .ne. 0 ) write(6,*) year,month,day,date
        end
        
c************************************************************************

        subroutine check_tpack(year,month,day,iprint)
        integer year,month,day,iprint
        integer year0,month0,day0
        integer date
        call packdate(date,year,month,day)
        call unpackdate(date,year0,month0,day0)
        if( year .ne. year0 .or. month .ne. month0 
     +          .or. day .ne. day0 ) then
          write(6,*) '*** error : ',year,month,day,year0,month0,day0
        end if
        if( iprint .ne. 0 ) write(6,*) year,month,day,date
        end
        
c************************************************************************

        subroutine test_dt_pack

        implicit none

        integer j,year,month,day

        year = 2001
        write(6,*) 'checking packing of complete year ',year
        do j=1,365
          call j2date(j,year,month,day)
          call check_dpack(year,month,day,0)
        end do

        year = 2004
        write(6,*) 'checking complete year ',year
        do j=1,366
          call j2date(j,year,month,day)
          call check_dpack(year,month,day,0)
        end do

        end

c************************************************************************

        subroutine test_var

        implicit none

        integer year,month,day
        integer hour,min,sec
	integer date,time,it
	character*30 line

	date = 20050209
	time = 0

	call dtsini(date,time)

	it = 0
	call dts2dt(it,year,month,day,hour,min,sec)
        !call dts2it(it,year,month,day,hour,min,sec)
	call dtsgf(it,line)

	write(6,*) 'test_var:'
	write(6,*) date,time
	write(6,*) it,year,month,day,hour,min,sec
	write(6,*) line

	end

c************************************************************************

        subroutine date_compute

c interactively check date/time and seconds

        implicit none

        integer year,month,day,hour,min,sec
        integer date,time
        integer mode,it

    3   continue

        write(6,*) 'Enter reference year (CR to exit): '
        read(5,'(i10)') year
        if( year .eq. 0 ) goto 2

        call dtsyear(year)

    4   continue

        write(6,*) 'date to secs: -1   secs to date: +1     end: CR'
        read(5,'(i10)') mode
        if( mode .eq. 0 ) goto 3

    5   continue

        if( mode .gt. 0 ) then
          write(6,*) 'Enter seconds: (CR to end):'
          read(5,'(i10)') it
          if( it .eq. 0 ) goto 4
	  call dts2dt(it,year,month,day,hour,min,sec)
        else
          write(6,*) 'Enter date (YYYYMMDD): (CR to end):'
          read(5,'(i10)') date
          if( date .eq. 0 ) goto 4
          time = 0
	  call unpackdate(date,year,month,day)
	  call unpacktime(time,hour,min,sec)
	  call dts2it(it,year,month,day,hour,min,sec)
        end if

        write(6,*) it,year,month,day,hour,min,sec

        goto 5
        
    2   continue

        end

c************************************************************************

        subroutine test_date_all

	call test_j2date
        call test_dt_pack
        call test_sd_pack
        call test_var

        end

c************************************************************************
c	program datet
c	!call datetest
c        call test_date_all
c        call date_compute
c	end
c************************************************************************