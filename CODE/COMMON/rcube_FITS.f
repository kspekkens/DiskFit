      subroutine rcube_FITS
c Copyright (C) 2020, Jerry Sellwood
c
c reads velocity data cube
c
c   Created by JAS January 2020
c
      use aacube
      implicit none
      include 'comcube.h'
c
c local FITS variables, from FITS cookbook
      character*4 comment
      integer bitpix, blocksize, firstpix, group, nfound
      integer readwrite, status, unit
      integer naxes( 3 )
      logical anynull
      real dv, nullval, v1
      real, allocatable :: buffer( : )
      real*8, allocatable :: buf( : )
c
c local variables
      character string*8
      integer i, j, k, lnblnk, ioffs, nchan, nline
      logical vdecr, lfreq
      real dmax, dmin, frest, vfac, y
c
c the status parameter must always be initialized
      status = 0
c
c get unused logical unit numbers for the FITS file
      call ftgiou( unit, status )
c
      readwrite = 0
      call ftopen( unit, datfile, readwrite, blocksize, status )
c this is fatal - could not open the file named as the data cube
      if( status .ne. 0 )then
        i = lnblnk( datfile )
        print *, 'failed to open FITS cube file ' // datfile( 1:i )
        call crash( 'rcube_FITS' )
      end if
c determine the size of the image
      call ftgknj( unit, 'NAXIS', 1, 3, naxes, nfound, status )
c check that it found all NAXIS keywords
      if( nfound .ne. 3 )then
        print *, 'Failed to read the NAXISn keyword'
        call crash( 'rcube_FITS' )
      end if
c      print *, 'cube size', naxes
c determine the type of data
      call ftgidt( unit, bitpix, status )
c      print *, 'bitpix =', bitpix
c get pixel scale and convert to arcsec
      call ftgkye( unit, 'CDELT2', v1, comment, status )
      pixscale = 3600. * abs( v1 )
c get nature of 3rd axis
      call ftgkys( unit, 'CTYPE3', string, comment, status )
      if( string .eq. 'FELO-HEL' .or. string .eq. 'VELO-HEL' )then
        lfreq = .false.
      else if( string( 1:4 ) .eq. 'FREQ' )then
        lfreq = .true.
        frest = 1420405751.7667
      else
        print *, 'Unrecognized keyword string CTYPE3 = ' // string
        call crash( 'rcube_FITS' )
      end if
c get velocity channel parameters
      call ftgkye( unit, 'CRPIX3', v1, comment, status )
      ioffs = -nint( v1 ) + 1
      call ftgkye( unit, 'CRVAL3', v1, comment, status )
      call ftgkye( unit, 'CDELT3', dv, comment, status )
c relativistic formula good for all v, speed of light in m/s
      if( lfreq )then
        y = ( frest - v1 ) / frest
        v1 = 2.998e8 * ( 2. * y - y**2 ) / ( 2. - 2 * y + y**2 )
        y = dv / frest
        dv = -2.998e8 * ( 2. * y - y**2 ) / ( 2. - 2 * y + y**2 )
      end if
      vdecr = dv .lt. 0.
      vfac = 1
      if( VELmps )vfac = 1.e-3
      if( .not. vdecr )then
        print *, 'Warning, positive dv is untested!'
        print *, 'Velocity range',
     +        vfac * ( v1 + dv * real( ioffs ) ),
     +        vfac * ( v1 + dv * real( naxes( 3 ) + ioffs - 1 ) )
        call gtintg( 'Enter 1 to continue', i )
        if( i .ne. 1 )stop
c      else
c        print *, 'Velocity range',
c     +        vfac * ( v1 + dv * real( naxes( 3 ) + ioffs - 1 ) ),
c     +        vfac * ( v1 + dv * real( ioffs ) )
      end if
c
      nsizex = naxes( 1 )
      nsizey = naxes( 2 )
      nsizev = naxes( 3 )
      group = 1
      firstpix = 1
      nullval = 0 ! no checks for undefined pixels - done elsewhere
c input buffers
      allocate ( buffer( nsizex ) )
      if( bitpix .eq. -64 )allocate ( buf( nsizex ) )
c raw data arrays
      allocate ( dcube( nsizex, nsizey, nsizev ) )
      allocate ( xval( nsizex ) )
      allocate ( yval( nsizey ) )
      allocate ( vval( nsizev ) )
c
      j = 0
      dmin = 1.e6
      dmax = -dmin
c loop through velocity channels
      nchan = 1
      do while ( nchan .le. nsizev )
        k = nchan
        if( vdecr )k = nsizev + 1 - nchan
        vval( k ) = v1 + dv * real( nchan + ioffs - 1 )
c loop through y
        nline = 1
        do while ( nline .le. nsizey )
          if( bitpix .eq. -32 )then
            call ftgpve( unit, group, firstpix, nsizex, nullval,
     +                   buffer, anynull, status )
          else if( bitpix .eq. -64 )then
            call ftgpvd( unit, group, firstpix, nsizex, nullval,
     +                   buf, anynull, status )
            do i = 1, nsizex
              buffer( i ) = buf( i )
            end do
          else
            print *, 'Unrecognized value of BITPIX:', bitpix
            call crash( 'rcube_FITS' )
          end if
          yval( nline ) = nline
c place data in array - arranging planes in ascending order of velocity
          i = 1
          do while ( i .le. nsizex )
            if( nline .eq. 1 )xval( i ) = i
            dmin = min( dmin, buffer( i ) )
            dmax = max( dmax, buffer( i ) )
            dcube( i, nline, k ) = buffer( i )
            j = j + 1
            i = i + 1
          end do
c end loop over y
          firstpix = firstpix + nsizex
          nline = nline + 1
        end do
c end loop over velocity channels
        nchan = nchan + 1
      end do
      print *, j, ' intensity values, min and max:', dmin, dmax
      print *, nsizev, ' velocity channels from',
     +  vfac * vval( 1 ), ' to', vfac * vval( nsizev ), ' km/s'
c the FITS files must always be closed before exiting the program
c   any unit numbers allocated with FTGIOU are freed with FTFIOU
      call ftclos( unit, status )
      call ftfiou( unit, status )
c print out error messages, if any error was flagged
      if( status .gt. 0 )call printerror( status )
      return
      end
