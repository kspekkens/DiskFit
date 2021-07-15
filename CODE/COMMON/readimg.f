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
c   Use FITS routines for data type conversions - JAS Apr 16
c
      include 'commons.h'
c
c local variables
      integer blocksize, group, i, j, lnblnk, nfound, nullval
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
	print *, 'failed to read the NAXIS keywords', nfound
        call crash( 'readimg' )
      end if
c set array size for this image
      nsizex = naxes( 1 )
      nsizey = naxes( 2 )
c initialize variables
      group = 1
      nullval = 0
      anyf = .false.
c allocate space and red the image
      allocate ( ldat( nsizex, nsizey ) )
      call ftg2de( unit, group, nullval, nsizex, nsizex,
     +             nsizey, ldat, anyf, status )
c close the FITS file and release the unit number
      call ftclos( unit, status )
      call ftfiou( unit, status )
c check for any error, and if so print out error messages
      if( status .gt. 0 )call printerror( status )
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
c check the size of the image
        call ftgknj( unit, 'NAXIS', 1, 2, naxes, nfound, status )
        if( nfound .ne. 2 )then
          print *,
     +            'value of NAXIS keyword in uncertainties file', nfound
          call crash( 'readimg' )
        end if
        if( ( naxes( 1 ) .ne. nsizex ) .or.
     +      ( naxes( 2 ) .ne. nsizey ) )then
          print *, 'uncertanties size differs from image size - fix it'
          call crash( 'readimg' )
        end if
c initialize variables
        group = 1
        nullval = 0
        anyf = .false.
c allocate space and read data
        allocate ( ldate( nsizex, nsizey ) )
        call ftg2de( unit, group, nullval, nsizex, nsizex,
     +               nsizey, ldate, anyf, status )
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
c check the size of the image
        call ftgknj( unit, 'NAXIS', 1, 2, naxes, nfound, status )
        if( nfound .ne. 2 )then
          print *, 'value of NAXIS keyword in mask file', nfound
          call crash( 'readimg' )
        end if
        if( ( naxes( 1 ) .ne. nsizex ) .or.
     +      ( naxes( 2 ) .ne. nsizey ) )then
          print *, 'mask size differs from image size - fix it'
          call crash( 'readimg' )
        end if
c initialize variables
        group = 1
        nullval = 0
        anyf = .false.
c allocate space and read data
        allocate ( dmask( nsizex, nsizey ) )
        call ftg2de( unit, group, nullval, nsizex, nsizex,
     +               nsizey, dmask, anyf, status )
c close the FITS file and release the unit number
        call ftclos( unit, status )
        call ftfiou( unit, status )
c check for any error, and if so print out error messages
        if( status .gt. 0 )call printerror( status )
c        print *, 'read mask file'
      end if
      return
      end
