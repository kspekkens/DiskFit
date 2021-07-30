      subroutine write_mask( work, nx, ny )
c Copyright (C) 2020, Jerry Sellwood
      use aacube
      implicit none
c
c creates a FITS image of the blanking mask
c
c calling arguments
      integer nx, ny
      real work( nx, ny )
c
      include 'comcube.h'
c
c externals
      integer lnblnk
c
c local FITS variables
      integer bitpix, blocksize, fpixel, naxis, readwrite, status
      integer naxes( 2 )
      logical anynull
      real nullval
c
c local arrays
      integer*2, allocatable :: outl16( : )
      real, allocatable :: outline( : )
      real*8, allocatable :: outlin2( : )
c
c local variables
      character file*100, outfits*100
      integer i, iout, j, jj, k
      integer unitD, unitM, unitO
      real diff, x
c
      if( nx .ne. nsizex .or. ny .ne. nsizey )then
        print *, 'Expected image size', nsizex, nsizey
        print *, 'Output array dimens', nx, ny
        call crash( 'Muddle in WRITE_MASK' )
      end if
c create fits files
      nullval = -600000
      status = 0
c reopen the input data file as readonly
      call ftgiou( unitD, status )
      readwrite = 0
      call ftopen( unitD, datfile, readwrite, blocksize, status )
c determine the type of data
      call ftgidt( unitD, bitpix, status )
      allocate ( outline( nsizex ) )
      if( bitpix .eq. -64 )then
        allocate ( outlin2( nsizex ) )
      else if( bitpix .eq. 16 )then
        allocate ( outl16( nsizex ) )
      end if
c construct file name
      k = lnblnk( maskfile )
      if( k .eq. 0 )then
        k = lnblnk( outroot )
        outfits = outroot( 1:k )
        outfits( k+1:k+5 ) = '.' // 'mask'
        k = k + 5
      else
        outfits = maskfile
      end if
      if( outfits( k-4:k ) .ne. '.fits' )then
        outfits( k+1:k+6 ) = '.fits'
        k = k + 6
      end if
c delete the file if it already exists, so we can then recreate it
      call deletefile( outfits, status )
c start the fits writing process
      call ftgiou( unitO, status )
      blocksize = 1             ! needed for historical reasons.
      call ftinit( unitO, outfits, blocksize, status )
c create minimal required header
      naxis = 2
      naxes( 1 ) = nsizex
      naxes( 2 ) = nsizey
      call ftphps( unitO, bitpix, naxis, naxes, status )
c this is the main loop - select fitted Vp (and its uncertainty) only
      fpixel = 1
      do j = 1, naxes( 2 )
        if( bitpix .eq. -32 )then
          do i = 1, naxes( 1 )
            outline( i ) = work( i, j )
          end do
          call ftppne( unitO, 0, fpixel, naxes( 1 ), outline,
     +                  nullval, status )
        else if( bitpix .eq. -64 )then
          do i = 1, naxes( 1 )
            outlin2( i ) = work( i, j )
          end do
          call ftppnd( unitO, 0, fpixel, naxes( 1 ), outlin2,
     +                 nullval, status )
        else if( bitpix .eq. 16 )then
          do i = 1, naxes( 1 )
            outl16( i ) = work( i, j )
          end do
          call ftppnj( unitO, 0, fpixel, naxes( 1 ), outl16,
     +                 nullval, status )
        end if
        fpixel = fpixel + naxes( 1 )
      end do
c add some useful header information
      call ftpcom( unitO, '******************************', status )
      write( file, '( f5.3 )' )version
      call ftpcom( unitO, ' MAP CREATED BY BLANK v' //
     +             file( 1:5 ), status )
      call ftpcom( unitO, '******************************', status )
c record date and time
      call timestamp( file )
      k = lnblnk( file )
      call ftpcom( unitO,
     +               'Date map created: ' // file( 1:k ), status )
c file names
      file = infile
      call ftpkys( unitO, 'INFILE', file,
     +               'Input text file', status )
      file = outroot
      call ftpkys( unitO, 'OUTFILE', file,
     +               'Root name for output files', status )
c put some useful text in the string file
      file = ' '
      call ftpkye( unitO, 'NOISE', noise, 4,
     +      'estimated from edge planes of smoothed data cube', status )
      call ftpkye( unitO, 'Smooth', smsig, 3,
     +                    'Gaussian smoothing sigma in pixels', status )
c close file
      call ftclos( unitO, status )
      call ftfiou( unitO, status )
      if( status .gt. 0 )call printerror( status )
      i = lnblnk( outfits )
      print *, 'Wrote:  ', outfits( 1:i )
      return
      end
