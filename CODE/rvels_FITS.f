      subroutine rvels_FITS
c reads velocity map, and possibly velocity error map, as FITS files
c
c   Created by RZS March 2009
c   Polished by JAS October 2009
c   Revised by JAS March 2011
c   Allow for 64-bit reals JAS April 2013
c
      include 'commons.h'
c
c local FITS variables, from FITS cookbook
      integer bitpix, blocksize, firstpix, group, nbuffer, nfound
      integer readwrite, status, unit, unit2
      integer naxes( 2 ), naxes2( 2 )
      logical anynull
      real evelnow, eveltmp, nullval
      real buffer( mpx ), buffer2( mpx )
      real*8 buf( mpx )
c
c local variables
      integer i, istepin, j, nline
      real cospa, mpsfac, r, sinpa, u, v, velmax, velmin, x, y
c      real nanum
      parameter ( istepin = 1 )
c
      l2D = .true.
c process inner region variables which are used only here
      b_over_a = 1 - regeps
      regpa = ( regpa + 90. ) * pi / 180.
c
      mpsfac = 1
      if( VELmps )mpsfac = 1000
c
c read in velocities ( vsys - 500.0 ) < V < ( vsys + 500.0 ) only
      velmin = ( vsys - 500.0 ) * mpsfac
      velmax = ( vsys + 500.0 ) * mpsfac
c use eISM to control the fake error value. Adopt 10 if none
c   is given. evelnow is used only when no error file exists
      eveltmp = errtol / 2.0 ! important in IF below
      evelnow = 0.0
      if( eISM .le. 0.0 )evelnow = 10.0
c the status parameter must always be initialized
      status = 0
c
c get unused logical unit numbers for the FITS file
      call ftgiou( unit, status )
      call ftgiou( unit2, status )
c
      readwrite = 0
      call ftopen( unit, invfile, readwrite, blocksize, status )
c this is fatal - could not open the file named as the velocity map
      if( status .ne. 0 )then
        print *, 'failed to open FITS velocity map'
        call crash( 'rvels_FITS' )
      end if
c attempt to read error file. If it doesn't exist then set all errors to 0
c   and rely on turbulence term to set a contant error for all input values
      call ftopen( unit2, inevfile, readwrite, blocksize, status )
      if( status .eq. 0 )then
c file was opened - OK to proceed
        evelfile = .true.
      else if( status .eq. 104 )then
c file doesn't exist, so will assume constant 0 error
        evelfile = .false.
        call ftfiou( unit2, status )
c        print *, 'Adopting error ', evelnow, ', eISM ', eISM, 'km/s'
        status = 0
      else
c there was some other error opening the file - warn and carry on
        evelfile = .false.
        call ftfiou( unit2, status )
c        print *, 'No error file found.'
c        print *, 'Adopting error ', evelnow, ', eISM ', eISM, 'km/s'
        status = 0
      end if
c
c determine the size of the image
      call ftgknj( unit, 'NAXIS', 1, 2, naxes, nfound, status )
c check that it found both NAXIS1 and NAXIS2 keywords
      if( nfound .ne. 2 )then
        print *, 'Failed to read the NAXISn keywords from vels file'
        call crash( 'rvels_FITS' )
      end if
c determine the type of data
      call ftgidt( unit, bitpix, status )
c
      if( evelfile )then
        call ftgknj( unit2, 'NAXIS', 1, 2, naxes2, nfound, status )
c check that it found both NAXIS1 and NAXIS2 keywords
        if( nfound .ne. 2 )then
          print *, 'Failed to read the NAXISn keywords from error file'
          call crash( 'rvels_FITS' )
        end if
        if( ( naxes( 1 ) .ne. naxes2( 1 ) ) .and.
     +      ( naxes( 2 ) .ne. naxes2( 2 ) ) )then
          print *, naxes, naxes2
          print *, 'Your vels and evels images differ in size. Fix it.'
          call crash( 'rvels_FITS' )
        end if
c check the type of data
        call ftgidt( unit2, i, status )
        if( i .ne. bitpix )then
          print *, i, bitpix
          print *,
     +         'Your vels and evels images differ in data type. Fix it.'
          call crash( 'rvels_FITS' )
        end if
      end if
c
      nsizex = naxes( 1 )
      nsizey = naxes( 2 )
      if( ( nsizex .gt. mpx ) .or. ( nsizey .gt. mpy ) )then
        print *, 'Too many pixels needed', nsizex, nsizey
        call crash( 'rvels_FITS' )
      end if
      group = 1
      firstpix = 1
      nullval = 0 ! no checks for undefined pixels - done elsewhere
c
      nbuffer = min( naxes( 1 ), MPX )
      if( naxes( 1 ) .gt. MPX )then
        print *, '****'
        print *, 'Your image is too large. Increase MPX'
        print *, 'Reading only box ', nbuffer, ' x', nbuffer
      end if
c
      cospa = cos( regpa )
      sinpa = sin( regpa )
c
      j = 0
c loop though array twice - first get pixels within regrad
      nline = 1
      do while ( nline .le. naxes( 2 ) )
        if( bitpix .eq. -32 )then
          call ftgpve( unit, group, firstpix, nbuffer, nullval,
     +                 buffer, anynull, status )
        else if( bitpix .eq. -64 )then
          call ftgpvd( unit, group, firstpix, nbuffer, nullval,
     +                 buf, anynull, status )
          do i = 1, nbuffer
            buffer( i ) = buf( i )
          end do
        else
          print *, 'Unrecognized value of BITPIX:', bitpix
          call crash( 'rvels_FITS' )
        end if
c read uncertainties only if the file exists
        if( evelfile )then
          if( bitpix .eq. -32 )then
            call ftgpve( unit2, group, firstpix, nbuffer, nullval,
     +                   buffer2, anynull, status )
          else if( bitpix .eq. -64 )then
            call ftgpvd( unit2, group, firstpix, nbuffer, nullval,
     +                   buf, anynull, status )
            do i = 1, nbuffer
              buffer2( i ) = buf( i )
            end do
          end if
        end if
        yval( nline ) = nline
c place data in array
        i = 1
        do while ( i .le. nbuffer )
          xval( i ) = i
          x = real( i ) - xcen
          y = real( nline ) - ycen
          u = x * cospa + y * sinpa
          v = -x * sinpa + y * cospa
          r = sqrt( u**2 + ( v / b_over_a )**2 )
          if( evelfile )then
            evelnow = buffer2( i ) / mpsfac
            eveltmp = evelnow  ! take the error seriously
          end if
c set initial dummy values for all pixels
          ldat( i, nline ) = -10000
          ldate( i, nline ) = 0
          lgpix( i, nline ) = .false.
          if( r .le. regrad )then
c use pixel only if value is valid
            if( ( buffer( i ) .gt. velmin ) .and.
     +          ( buffer( i ) .lt. velmax ) .and.
     +          ( eveltmp .gt. 0.01 ) .and.
     +          ( eveltmp .lt. errtol ) )then
              ldat( i, nline ) = buffer( i ) / mpsfac
              ldate( i, nline ) = evelnow
              lgpix( i, nline ) = .true.
              j = j + 1
            end if
          end if
          i = i + istepin
        end do
c increment pointers and loop back to read the next group of pixels
        firstpix = firstpix + naxes( 1 ) * istepin
        nline = nline + istepin
      end do
c
c now get pixels outside regrad
      firstpix = 1
      nline = 1
      do while ( nline .le. naxes( 2 ) )
        call ftgpve( unit, group, firstpix, nbuffer, nullval,
     +               buffer, anynull, status )
        if( evelfile )then
          call ftgpve( unit2, group, firstpix, nbuffer, nullval,
     +                 buffer2, anynull, status )
        end if
c
c place data in array
        i = 1
        do while ( i .le. nbuffer )
          x = real( i ) - xcen
          y = real( nline ) - ycen
          u = x * cospa + y * sinpa
          v = -x * sinpa + y * cospa
          r = sqrt( u**2 + ( v / b_over_a )**2 )
          if( evelfile )then
            evelnow = buffer2( i ) / mpsfac
            eveltmp = evelnow   ! take error seriously
          end if
          if( r .gt. regrad )then
c use pixel only if value is valid
            if( ( buffer( i ) .gt. velmin ) .and.
     +          ( buffer( i ) .lt. velmax ) .and.
     +          ( eveltmp .gt. 0.01 ) .and.
     +          ( eveltmp .lt. errtol ) )then
              ldat( i, nline ) = buffer( i ) / mpsfac
              ldate( i, nline ) = evelnow
              lgpix( i, nline ) = .true.
              j = j + 1
            end if
          end if
          i = i + istepout
        end do
c increment pointers and loop back to read the next group of pixels
        firstpix = firstpix + naxes( 1 ) * istepout
        nline = nline + istepout
      end do
c total number of points read
      inp_pts = j - 1
c the FITS files must always be closed before exiting the program
c   any unit numbers allocated with FTGIOU are freed with FTFIOU
      call ftclos( unit, status )
      call ftfiou( unit, status )
      if( evelfile )then
        call ftclos( unit2, status )
        call ftfiou( unit2, status )
      end if
c print out error messages, if any error was flagged
      if( status .gt. 0 )call printerror( status )
      return
      end
