c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c Program to post-process output files from bootstrap iterations
c
c Created by JAS - June 2012
c
      include 'commons.h'
c
c external
      integer lnblnk
c
c local arrays
      integer m
      parameter ( m = 1000 )
      real*8 afitval( m, mtot ), aparams( m, md ), lfrac( 3, m )
      real*8 eparams( md ), params( md )
c
c local variables
      integer i, ii, j, k
      real*8 a, adev, ave, b, curt, sdev, skew, var
      character*100 bootfile, line
c
      print *, 'DiskFit (v1.2)  Copyright (C) 2015,' //
     +         ' Jerry Sellwood and Kristine Spekkens'
      print *
      print *, 'This program comes with ABSOLUTELY NO WARRANTY.' //
     +         '  It is free software'
      print *, 'and you are welcome to' //
     +         ' redistribute it under certain conditions.'
      print *, 'See <http://www.gnu.org/licenses/> for details.'
      print *
c
c read input values
      print *,'Enter input parameter file name'
      read( *, * )infile
      call readinp
c
      if( .not. luncert )then
        print *, 'bootlace run when uncertainties are not requested'
        call crash( 'bootlace' )
      end if
c
c initialize error vectors
      do i = 1, nd
        eparams( i ) = 0
      end do
      do i = 1, ntot
        efitval( i ) = 0
      end do
c
      if( nunc .gt. m )then
        print *, 'too many bootstrap iterations requested', nunc
        print *, 'decrease nunc or contact code administrator'
        call crash( 'bootlace' )
      end if
c
c open bootstrap file
      k = lnblnk( outroot ) + 1
      bootfile = outroot
      write( bootfile( k:k+4 ), '( a5 )' )'bstrp'
      open( 4, file = bootfile, status = 'old' )
c check first line of header
      read( 4, '( a )' )line
      if( line .ne. infile )then
        print *, 'first line of header file appears wrong'
        print '( ''    read: '', a )', line( 1:50 )
        print '( ''expected: '', a )', infile( 1:50 )
c        call crash( 'bootlace' )
      end if
c skip remainder of header lines
      do i = 1, 4
        read( 4, * )
      end do
c bootstrap loop
      j = 0
      read( 4, '( a )', iostat = k )line
c read best-fit parameters only once
      read( line( 19:100 ), * )ii
      if( ii .eq. 0 )then
        print *, 'Reading best fit parameters'
        read( 4, * )( params( i ), i = 1, nd )
        read( 4, * )( fitval( i ), i = 1, ntot )
        if( lphot .and. ldisk .and. ( lnax .or. lbulge ) )then
          read( 4, * )dskfrac, barfrac, blgfrac
        end if
        read( 4, * )
        read( 4, '( a )', iostat = k )line
      else
        print *, 'Best fitting parameter values not found'
        call crash( 'bootlace' )
      end if
      do while ( k .eq. 0 )
c read next fitted parameters from bootstrap file
        read( line( 19:100 ), * )ii
        print *, 'Reading bootstrap number: ', ii
        j = j + 1
        if( j .gt. m )then
          print *, 'arrays too small - too many iterations'
          call crash( 'bootlace' )
        end if
        read( 4, * )( aparams( j, i ), i = 1, nd )
        read( 4, * )( afitval( j, i ), i = 1, ntot )
        if( lphot .and. ldisk .and. ( lnax .or. lbulge ) )then
          read( 4, * )lfrac( 1, j ), lfrac( 2, j ), lfrac( 3, j )
        end if
        read( 4, * )
c do not want best-fit parameters flagged with ii = 0, so over-write them
        if( ii .eq. 0 )j = j - 1
        read( 4, '( a )', iostat = k )line
c skip any header info from concatenating files
        if( ( k .eq. 0 ) .and. ( line .eq. infile ) )then
          do i = 1, 4
            read( 4, * )
          end do
          read( 4, '( a )', iostat = k )line
        end if
c end bootstrap loop here
      end do
c close input file
      close( 4 )
      print *, 'Actual number of bootstraps read', j
c process data
      nunc = j
      call geterrs( aparams, afitval, lfrac, eparams, m )
      call bestfit( params, eparams )
c revise errors in optimal output file
      call revisout
      end

