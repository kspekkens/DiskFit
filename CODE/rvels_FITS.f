      subroutine rvels_FITS
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c reads velocity map, and possibly velocity error map, as FITS files
c
c   Created by RZS March 2009
c   Polished by JAS October 2009
c   Revised by JAS March 2011
c   Allow for 64-bit reals JAS April 2013
c   Updated to f90 - JAS Jan 2015
c
      include 'commons.h'
c
c local FITS variables, from FITS cookbook
      integer bitpix, blocksize, firstpix, group, nfound
      integer readwrite, status, unit, unit2
      integer naxes( 2 ), naxes2( 2 )
      logical anynull
      real evelnow, eveltmp, nullval
      real, allocatable :: buffer( : ), buffer2( : )
      real*8, allocatable :: buf( : )
c
c local variables
      integer i, istepin, j, lnblnk, nline
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
      call ftopen( unit, datfile, readwrite, blocksize, status )
c this is fatal - could not open the file named as the velocity map
      if( status .ne. 0 )then
        i = lnblnk( datfile )
        print *,
     +        'failed to open FITS velocity map file ' // datfile( 1:i )
        call crash( 'rvels_FITS' )
      end if
c attempt to read error file - if one is expected
      if( lerrfile )then
        call ftopen( unit2, errfile, readwrite, blocksize, status )
        if( status .ne. 0 )then
c named file was not found
          i = lnblnk( errfile )
          print *,
     +        'Unable to open the uncertainties file ' // errfile( 1:i )
          call crash( 'rvels_FITS' )
        end if
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
      if( lerrfile )then
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
      group = 1
      firstpix = 1
      nullval = 0 ! no checks for undefined pixels - done elsewhere
c input buffers
      allocate ( buffer( nsizex ) )
      if( lerrfile )allocate ( buffer2( nsizex ) )
      if( bitpix .eq. -64 )allocate ( buf( nsizex ) )
c raw data arrays
      allocate ( ldat( nsizex, nsizey ) )
      allocate ( ldate( nsizex, nsizey ) )
      allocate ( lgpix( nsizex, nsizey ) )
      allocate ( xval( nsizex ) )
      allocate ( yval( nsizey ) )
c
      cospa = cos( regpa )
      sinpa = sin( regpa )
c
      j = 0
c loop though array twice - first get pixels within regrad
      nline = 1
      do while ( nline .le. nsizey )
        if( bitpix .eq. -32 )then
          call ftgpve( unit, group, firstpix, nsizex, nullval,
     +                 buffer, anynull, status )
        else if( bitpix .eq. -64 )then
          call ftgpvd( unit, group, firstpix, nsizex, nullval,
     +                 buf, anynull, status )
          do i = 1, nsizex
            buffer( i ) = buf( i )
          end do
        else
          print *, 'Unrecognized value of BITPIX:', bitpix
          call crash( 'rvels_FITS' )
        end if
c read uncertainties only if the file exists
        if( lerrfile )then
          if( bitpix .eq. -32 )then
            call ftgpve( unit2, group, firstpix, nsizex, nullval,
     +                   buffer2, anynull, status )
          else if( bitpix .eq. -64 )then
            call ftgpvd( unit2, group, firstpix, nsizex, nullval,
     +                   buf, anynull, status )
            do i = 1, nsizex
              buffer2( i ) = buf( i )
            end do
          end if
        end if
        yval( nline ) = nline
c place data in array
        i = 1
        do while ( i .le. nsizex )
          xval( i ) = i
          x = real( i ) - xcen
          y = real( nline ) - ycen
          u = x * cospa + y * sinpa
          v = -x * sinpa + y * cospa
          r = sqrt( u**2 + ( v / b_over_a )**2 )
          if( lerrfile )then
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
        firstpix = firstpix + nsizex * istepin
        nline = nline + istepin
      end do
c
c now get pixels outside regrad
      firstpix = 1
      nline = 1
      do while ( nline .le. nsizey )
        if( bitpix .eq. -32 )then
          call ftgpve( unit, group, firstpix, nsizex, nullval,
     +                 buffer, anynull, status )
        else if( bitpix .eq. -64 )then
          call ftgpvd( unit, group, firstpix, nsizex, nullval,
     +                 buf, anynull, status )
          do i = 1, nsizex
            buffer( i ) = buf( i )
          end do
        end if
c read uncertainties only if the file exists
        if( lerrfile )then
          if( bitpix .eq. -32 )then
            call ftgpve( unit2, group, firstpix, nsizex, nullval,
     +                   buffer2, anynull, status )
          else if( bitpix .eq. -64 )then
            call ftgpvd( unit2, group, firstpix, nsizex, nullval,
     +                   buf, anynull, status )
            do i = 1, nsizex
              buffer2( i ) = buf( i )
            end do
          end if
        end if
c
c place data in array
        i = 1
        do while ( i .le. nsizex )
          x = real( i ) - xcen
          y = real( nline ) - ycen
          u = x * cospa + y * sinpa
          v = -x * sinpa + y * cospa
          r = sqrt( u**2 + ( v / b_over_a )**2 )
          if( lerrfile )then
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
        firstpix = firstpix + nsizex * istepout
        nline = nline + istepout
      end do
c total number of points read
      inp_pts = j - 1
c the FITS files must always be closed before exiting the program
c   any unit numbers allocated with FTGIOU are freed with FTFIOU
      call ftclos( unit, status )
      call ftfiou( unit, status )
      if( lerrfile )then
        call ftclos( unit2, status )
        call ftfiou( unit2, status )
      end if
c print out error messages, if any error was flagged
      if( status .gt. 0 )call printerror( status )
      return
      end
