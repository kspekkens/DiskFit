      subroutine readimg
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c reads intensity data as a FITS image
c
c   Created by EIB
c   Polished by JAS for ASR
c   Adapted for inclusion in velfit by JAS Mar 2011
c   Added option to read an image mask: JAS Jun 13
c   Updated to f90 - JAS Jan 2015
c   Allow for user-supplied uncertainties - JAS Sep 15
c
      include 'commons.h'
c
c local array
      integer*2, allocatable :: timg16( :, : )
c
c local variables
      integer bitpix, blocksize, group, i, j, lnblnk, nfound, nullval
      integer readwrite, status, unit, naxes( 2 )
      logical anyf
c
      l2D = .true.
c get an unused logical unit number for the FITS file
      status = 0
      call ftgiou( unit, status )
c open the FITS file
      readwrite = 0
      blocksize = 1
      call ftopen( unit, datfile, readwrite, blocksize, status )
c this is fatal
      if( status .ne. 0 )then
        print *, 'failed to open FITS intensity map'
        call crash( 'readimg' )
      end if
c determine the size of the image
      call ftgknj( unit, 'NAXIS', 1, 2, naxes, nfound, status )
c check that it found both keywords
      if( nfound .ne. 2 )then
	print *, 'failed to read the NAXIS keywords'
        call crash( 'readimg' )
      end if
c check local arrays are large enough for this image
      nsizex = naxes( 1 )
      nsizey = naxes( 2 )
c determine the type of data in the file
      bitpix = 0
      call ftgidt( unit, bitpix, status )
c initialize variables
      group = 1
      nullval = 0
      anyf = .false.
c allocate space
      allocate ( ldat( nsizex, nsizey ) )
      if( bitpix .eq. 16 )then
        allocate ( timg16( nsizex, nsizey ) )
        call ftg2dj( unit, group, nullval, nsizex, nsizex,
     +               nsizey, timg16, anyf, status )
        do j = 1, nsizey
          do i = 1, nsizex
            ldat( i, j ) = timg16( i, j )
          end do
        end do
        deallocate ( timg16 )
      else if( bitpix .eq. -32 )then
        call ftg2de( unit, group, nullval, nsizex, nsizex,
     +               nsizey, ldat, anyf, status )
      else
        print *, 'Unrecognized value for bitpix', bitpix
        call crash( 'readimg' )
      end if
c close the FITS file and release the unit number
      call ftclos( unit, status )
      call ftfiou( unit, status )
c check for any error, and if so print out error messages
      if( status .gt. 0 )call printerror( status )
c      print *, 'read in FITS image'
c set xval & yval arrays
      allocate ( xval( nsizex ) )
      allocate ( yval( nsizey ) )
      do i = 1, nsizex
        xval( i ) = i
      end do
      do i = 1, nsizey
        yval( i ) = i
      end do
c
c read uncertainties file if one is supplied
c
      if( lerrfile )then
c get an unused logical unit number for the FITS file
        status = 0
        call ftgiou( unit, status )
c open the FITS file
        readwrite = 0
        blocksize = 1
c attempt to read error file
        call ftopen( unit, errfile, readwrite, blocksize, status )
        if( status .ne. 0 )then
c named file was not found
          i = lnblnk( errfile )
          print *,
     +        'Unable to open the uncertainties file ' // errfile( 1:i )
          call crash( 'readimg' )
        end if
c determine the size of the image
        call ftgknj( unit, 'NAXIS', 1, 2, naxes, nfound, status )
c check that it found both keywords
        if( nfound .ne. 2 )then
          print *,
     +         'failed to read the NAXIS keywords of uncertainties file'
          call crash( 'readimg' )
        end if
c check local arrays are large enough for this image
        if( ( naxes( 1 ) .ne. nsizex ) .or.
     +      ( naxes( 2 ) .ne. nsizey ) )then
          print *, 'uncertanties size differs from image size - fix it'
          call crash( 'readimg' )
        end if
c determine the type of data in the file
        bitpix = 0
        call ftgidt( unit, bitpix, status )
c initialize variables
        group = 1
        nullval = 0
        anyf = .false.
c allocate space
        allocate ( ldate( nsizex, nsizey ) )
c read data
        if( bitpix .eq. 16 )then
          call ftg2dj( unit, group, nullval, nsizex, nsizex,
     +                 nsizey, timg16, anyf, status )
          allocate ( timg16( nsizex, nsizey ) )
          do j = 1, nsizey
            do i = 1, nsizex
              ldate( i, j ) = timg16( i, j )
            end do
          end do
          deallocate ( timg16 )
        else if( bitpix .eq. -32 )then
          call ftg2de( unit, group, nullval, nsizex, nsizex,
     +                 nsizey, ldate, anyf, status )
        else
          print *,
     +     'Unrecognized value for bitpix in uncertainties file', bitpix
          call crash( 'readimg' )
        end if
c close the FITS file and release the unit number
        call ftclos( unit, status )
        call ftfiou( unit, status )
c check for any error, and if so print out error messages
        if( status .gt. 0 )call printerror( status )
c        print *, 'read uncertainties file'
      end if
c
c read mask file if one is supplied
c
      if( lmask )then
c get an unused logical unit number for the FITS file
        status = 0
        call ftgiou( unit, status )
c open the FITS file
        readwrite = 0
        blocksize = 1
        call ftopen( unit, mskfile, readwrite, blocksize, status )
        if( status .ne. 0 )then
c named file was not found
          i = lnblnk( mskfile )
          print *, 'Unable to open the mask file ' // mskfile( 1:i )
          call crash( 'readimg' )
        end if
c determine the size of the image
        call ftgknj( unit, 'NAXIS', 1, 2, naxes, nfound, status )
c check that it found both keywords
        if( nfound .ne. 2 )then
          print *, 'failed to read the NAXIS keywords of mask file'
          call crash( 'readimg' )
        end if
c check local arrays are large enough for this image
        if( ( naxes( 1 ) .ne. nsizex ) .or.
     +      ( naxes( 2 ) .ne. nsizey ) )then
          print *, 'mask size differs from image size - fix it'
          call crash( 'readimg' )
        end if
c determine the type of data in the file
        bitpix = 0
        call ftgidt( unit, bitpix, status )
c initialize variables
        group = 1
        nullval = 0
        anyf = .false.
c allocate space
        allocate ( dmask( nsizex, nsizey ) )
c read data
        if( bitpix .eq. 16 )then
          call ftg2dj( unit, group, nullval, nsizex, nsizex,
     +                 nsizey, timg16, anyf, status )
          allocate ( timg16( nsizex, nsizey ) )
          do j = 1, nsizey
            do i = 1, nsizex
              dmask( i, j ) = timg16( i, j )
            end do
          end do
          deallocate ( timg16 )
        else if( bitpix .eq. -32 )then
          call ftg2de( unit, group, nullval, nsizex, nsizex,
     +                 nsizey, dmask, anyf, status )
        else
          print *, 'Unrecognized value for bitpix in mask', bitpix
          call crash( 'readimg' )
        end if
c close the FITS file and release the unit number
        call ftclos( unit, status )
        call ftfiou( unit, status )
c check for any error, and if so print out error messages
        if( status .gt. 0 )call printerror( status )
c        print *, 'read mask file'
      end if
      return
      end
