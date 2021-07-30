      subroutine get_mask
c Copyright (C) 2020, Jerry Sellwood
c
c Routine to read a previously-created mask file and to set values in the
c   logical array llmask
      use aacube
      implicit none
c
      include 'comcube.h'
c
c local FITS variables, from FITS cookbook
      integer bitpix, blocksize, firstpix, group, nfound
      integer readwrite, status, unit
      integer naxes( 2 )
      logical anynull
      real nullval
      real, allocatable :: buffer( : )
      real*8, allocatable :: buf( : )
c
c local variables
      integer i, lnblnk, nline
c
c the status parameter must always be initialized
      status = 0
c
c get unused logical unit number for the FITS file
      call ftgiou( unit, status )
c
      readwrite = 0
      call ftopen( unit, maskfile, readwrite, blocksize, status )
c this is fatal - could not open the file named as the velocity map
      if( status .ne. 0 )then
        i = lnblnk( maskfile )
        print *, 'failed to open mask file ' // maskfile( 1:i )
        call crash( 'GET_MASK' )
      end if
c
c determine the size of the image
      call ftgknj( unit, 'NAXIS', 1, 2, naxes, nfound, status )
c check that it found both NAXIS1 and NAXIS2 keywords
      if( nfound .ne. 2 )then
        print *, 'Failed to read the NAXISn keywords from mask file'
        call crash( 'GET_MASK' )
      end if
c determine the type of data
      call ftgidt( unit, bitpix, status )
c
      nsizex = naxes( 1 )
      nsizey = naxes( 2 )
      if( ( naxes( 1 ) .gt. nsizex ) .or. ( naxes( 2 ) .gt. nsizey )
     +                      )call crash( 'Array too small in GET_MASK' )
      group = 1
      firstpix = 1
      nullval = 0 ! no checks for undefined pixels - done elsewhere
c input buffers
      allocate ( buffer( nsizex ) )
      if( bitpix .eq. -64 )allocate ( buf( nsizex ) )
c
c loop though array
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
          call crash( 'GET_MASK' )
        end if
c set up mask
        do i = 1, nsizex
          llmask( i, nline ) = buffer( i ) .gt. 0.
        end do
c increment pointers and loop back to read the next group of pixels
        firstpix = firstpix + nsizex
        nline = nline + 1
      end do
c the FITS files must always be closed before exiting the program
c   any unit numbers allocated with FTGIOU are freed with FTFIOU
      call ftclos( unit, status )
      call ftfiou( unit, status )
c print out error messages, if any error was flagged
      if( status .gt. 0 )call printerror( status )
      return
      end
