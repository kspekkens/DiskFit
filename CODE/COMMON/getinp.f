      subroutine getinp( fname )
c Copyright (C) 2020, Jerry Sellwood
      implicit none
c Plagiarized from subroutine readinp of DiskFit package
c
c This subroutine reads the named input file. The structure of the input file
c   is given in the documentation - when tat is ready!
c
c   lines of characters in the input file can be either in single quotes
c   or plain text.  If the first character of a line is ' [char(39)], the
c   string is presumed to start at the second character
c
c calling argument
      character*(*) fname
c
c common block
      include 'comcube.h'
c
c externals
      integer lnblnk
c
c local variables
      character line*120, string*10
      integer i, ic, istat, j, k, mGH
c
c set version number
      version = 1.0
c try to open specified file
      open( unit = 13, file = fname, status = 'old', err = 2 )
      infile = fname
c skip the first line
      read( 13, * )
c initialize line counter
      ic = 1
c line 2: name of data cube file
      call getline( 13, line, istat, ic )
      if( line( 1:1 ) .eq. char( 39 ) )then
        read( line, *, err = 1, end = 1 )datfile
      else
        datfile = line
      end if
      i = lnblnk( datfile )
      print *, 'Name of data cube file ' // datfile( 1:i )
c line 3: names of mask file
      call getline( 13, line, istat, ic )
      if( line( 1:1 ) .eq. char( 39 ) )then
        read( line, *, err = 1, end = 1 )maskfile
      else
        maskfile = line
      end if
      i = lnblnk( maskfile )
      print *, 'Name of mask file ' // maskfile( 1:i )
c line 4: root names of maps to be created
      call getline( 13, line, istat, ic )
      if( line( 1:1 ) .eq. char( 39 ) )then
        read( line, *, err = 1, end = 1 )outpfile
      else
        outpfile = line
      end if
c parse name of output file for subdirectories and root name
      j = lnblnk( outpfile )
      i = index( outpfile, '/' )
      k = i + 1
      do while ( i .gt. 0 )
        i = index( outpfile( k:j ), '/' )
        k = k + i
      end do
      i = index( outpfile( k:j ), '.' )
      k = k + i - 1
      if( k .gt. 0 )then
        outroot = outpfile( 1:k )
      else
        outroot = ' '
      end if
      print *, 'Names of map files ' // outpfile( 1:j )
c line 5: flag for vels in m/s
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )VELmps
      if( VELmps )then
        print *, 'Velocity channel values in m/s'
      else
        print *, 'Velocity channel values in km/s'
      end if
c line 6: noise threshold parameters - suggest 3 and 2
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )imin, thresh
      write( string, '( f4.1 )' ), imin
      print *, 'Line fits not attempted if peak intensity <' //
     +         string( 1:4 ) // ' times noise'
      write( string, '( f4.1 )' ), thresh
      print *, 'Line fits discarded if max fitted intensity <' //
     +         string( 1:4 ) // ' times noise'
c line 7: order of Gauss-Hermite fit
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )mGH
      if( mGH .eq. 0 )then
        mfit = 4
      else
        mfit = mGH + 2
      end if
      if( mfit .eq. 4 )then
        print *, 'Simple Gaussian fits to each line profile'
      else if( mfit .eq. 5 )then
        print *, 'h3 only fits to each line profile'
      else if( mfit .eq. 6 )then
        print *, 'h3 and h4 fits to each line profile'
      else
        print *, 'mGH not one of 0, 3 or 4'
        go to 1
      end if
c line 8: limits on line width and GH coefficients - suggest 20 and 1
      call getline( 13, line, istat, ic )
      if( mfit .eq. 4 )then
        read( line, *, err = 1, end = 1 )wmax
        hmax = 0
      else
        read( line, *, err = 1, end = 1 )wmax, hmax
      end if
      write( string, '( f4.1 )' ), wmax
      print *, 'Max allowed line width ' //
     +                    string( 1:4 ) // ' velocity channels'
      if( mfit .gt. 4 )then
        write( string, '( f4.2 )' ), hmax
        print *, 'Bound on abs values of GH coefficients ' // string
      end if
c line 9: number of bootstraps for each line fit
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )nboot
      if( nboot .eq. 0 )then
        print *, 'No uncertainties will be estimated'
      else
        print *, 'Uncertainties estimated from', nboot, ' iterations'
      end if
c line 10: list of toggles for creating FITS files
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )lfout
      extn( 1, 1 ) = 'mnV'
      extn( 1, 2 ) = 'emV'
      extn( 2, 1 ) = 'fit'
      extn( 2, 2 ) = '   '
      print *, 'The following FITS files will be created:'
      j = lnblnk( outpfile )
      if( nboot .eq. 0 )then
        print *, outpfile( 1:j ) // '.' // extn( 1, 1 ) // '.fits'
      else
        print *, outpfile( 1:j ) // '.' // extn( 1, 1 ) // '.fits' //
     + '  and ' //  outpfile( 1:j ) // '.' // extn( 1, 2 ) // '.fits'
      end if
      if( lfout )then
        print *, outpfile( 1:j ) // '.' // extn( 2, 1 ) // '.fits'
      end if
      close( 13 )
      print *, 'Successfully read input file'
      print *
      return
c read failed
    1 print '( 2x, a )', line
      print *, 'Problem with line', ic
    2 print '( ''Error reading file '', a )', fname
      call crash( 'getinp' )
      end
