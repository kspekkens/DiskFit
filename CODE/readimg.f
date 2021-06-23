      subroutine readimg
c reads intensity data as a FITS image
c
c   Created by EIB
c   Polished by JAS for ASR
c   Adapted for inclusion in velfit by JAS Mar 2011
c
      include 'commons.h'
c
c local array
      integer*2 timg16( mpx, mpy )
c
c local variables
      integer bitpix, blocksize, group, i, j, nfound, nullval
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
      call ftopen( unit, inpfile, readwrite, blocksize, status )
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
      if( ( naxes( 1 ) .gt. mpx ) .or. ( naxes( 2 ) .gt. mpy ) )then
        print *, 'local parameters mpx or mpy too small'
        print *, 'increase and recompile'
        call crash( 'readimg' )
      end if
      nsizex = naxes( 1 )
      nsizey = naxes( 2 )
c determine the type of data in the file
      bitpix = 0
      call ftgidt( unit, bitpix, status )
c initialize variables
      group = 1
      nullval = 0
      anyf = .false.
      if( bitpix .eq. 16 )then
        call ftg2dj( unit, group, nullval, mpx, naxes( 1 ),
     +               naxes( 2 ), timg16, anyf, status )
        do j = 1, naxes( 2 )
          do i = 1, naxes( 1 )
            ldat( i, j ) = timg16( i, j )
          end do
        end do
      else if( bitpix .eq. -32 )then
        call ftg2de( unit, group, nullval, mpx, naxes( 1 ),
     +               naxes( 2 ), ldat, anyf, status )
      else
        print *, 'Unrecognized value for bitpix', bitpix
        call crash( 'readimg' )
      end if
c close the FITS file and release the unit number
      call ftclos( unit, status )
      call ftfiou( unit, status )
c check for any error, and if so print out error messages
      if( status .gt. 0 )call printerror( status )
c set xval & yval arrays
      do i = 1, nsizex
        xval( i ) = i
      end do
      do i = 1, nsizey
        yval( i ) = i
      end do
      return
      end
