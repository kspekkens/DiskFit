      subroutine writemap( work, nx, ny, l )
c Copyright (C) 2020, Jerry Sellwood
      use aacube
      implicit none
c
c creates FITS images of the maps and optionally uncertainties of the same
c   x-y dimensions as the input data cube
c plagiarized from writemod from the DiskFit package
c
c calling arguments
      integer l, nx, ny
      real work( nx, ny, l )
c
      include 'comcube.h'
c
c externals
      integer lnblnk
c
c local FITS variables
      integer bitpix, blocksize, fpixel, naxis, readwrite, status
      integer naxes( 3 )
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
        call crash( 'Muddle in WRITEMAP' )
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
c loop over output files
      jj = 1
      if( nboot .gt. 0 )jj = 2
      do iout = 1, jj
c construct file names
        k = lnblnk( outpfile )
        if( k .gt. 0 )outfits = outpfile( 1:k )
        outfits( k+1:k+4 ) = '.' // extn( 1, iout )
        outfits( k+5:k+9 ) = '.fits'
c delete the file if it already exists, so we can then recreate it
        call deletefile( outfits, status )
c start the fits writing process
        call ftgiou( unitO, status )
        blocksize = 1           ! needed for historical reasons.
        call ftinit( unitO, outfits, blocksize, status )
c create minimal required header
        naxis = 2
        naxes( 1 ) = nsizex
        naxes( 2 ) = nsizey
        call ftphps( unitO, bitpix, naxis, naxes, status )
c this is the main loop - select fitted Vp (and its uncertainty) only
        k = mfit + 1
        if( iout .eq. 2 )k = 2 * mfit + 3
        fpixel = 1
        do j = 1, naxes( 2 )
          if( bitpix .eq. -32 )then
            do i = 1, naxes( 1 )
              outline( i ) = work( i, j, k )
            end do
            call ftppne( unitO, 0, fpixel, naxes( 1 ), outline,
     +                   nullval, status )
          else if( bitpix .eq. -64 )then
            do i = 1, naxes( 1 )
              outlin2( i ) = work( i, j, k )
            end do
            call ftppnd( unitO, 0, fpixel, naxes( 1 ), outlin2,
     +                   nullval, status )
          else if( bitpix .eq. 16 )then
            do i = 1, naxes( 1 )
              outl16( i ) = work( i, j, k )
            end do
            call ftppnj( unitO, 0, fpixel, naxes( 1 ), outl16,
     +                   nullval, status )
          end if
          fpixel = fpixel + naxes( 1 )
        end do
c add some useful header information
        call ftpcom( unitO, '******************************', status )
        write( file, '( f5.3 )' )version
        call ftpcom( unitO, ' MAP CREATED BY MAKEMAP v' //
     +               file( 1:5 ), status )
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
        file = outpfile
        call ftpkys( unitO, 'OUTFILE', file,
     +               'Root name for output files', status )
c put some useful text in the string file
        file = ' '
        k = lnblnk( file )
        if( k .gt. 0 )call ftpkys( unitO, 'COMPS', file,
     +                 'Fit parameters', status )
        k = 0
        if( mfit .gt. 4 )k = mfit - 2
        call ftpkyj( unitO, 'FIT ', k,
     +               'Order of GH expansion', status )
        call ftpkyj( unitO, 'BOOTS', nboot,
     +               'number of bootstrap iterations per line', status )
        call ftpkye( unitO, 'NOISE', noise, 4,
     +               'estimated from edge planes of data cube', status )
        call ftpkye( unitO, 'MIN_I', Imin, 1,
     +         'times noise: min intensity required for a fit', status )
        call ftpkye( unitO, 'StoN', thresh, 1,
     +               'Min S/N of fitted lines', status )
        x = wmax * abs( vval( 1 ) - vval( 2 ) )
        if( VELmps )x = 1.e-3 * x
        call ftpkye( unitO, 'VWIDTH', x, 2,
     +             'Max velocity width (km/s) of fitted lines', status )
        if( mfit .gt. 4 )then
          x = hmax
          call ftpkye( unitO, 'MAXH', x, 2,
     +                 'Max abs values of h3 and h4', status )
        end if
c close files
        call ftclos( unitO, status )
        call ftfiou( unitO, status )
        if( status .gt. 0 )call printerror( status )
        i = lnblnk( outfits )
        print *, 'Wrote:  ', outfits( 1:i )
      end do
c
      if( lfout )then
c construct file names
        k = lnblnk( outpfile )
        if( k .gt. 0 )outfits = outpfile( 1:k )
        outfits( k+1:k+4 ) = '.' // extn( 2, 1 )
        outfits( k+5:k+9 ) = '.fits'
c delete the file if it already exists, so we can then recreate it
        call deletefile( outfits, status )
c start the fits writing process
        call ftgiou( unitM, status )
        blocksize = 1           ! needed for historical reasons.
        call ftinit( unitM, outfits, blocksize, status )
c create minimal required header
        naxis = 3
        naxes( 1 ) = nsizex
        naxes( 2 ) = nsizey
c save fit parameters and chisq
        jj = mfit + 1
c and uncertainties in fit parameters, if computed
        if( nboot .gt. 0 )jj = jj + mfit
        naxes( 3 ) = jj
        call ftphps( unitM, bitpix, naxis, naxes, status )
c this is the main loop
        fpixel = 1
        jj = 0
        do k = 1, naxes( 3 )
          jj = jj + 1
c skip Vp
          if( k .eq. mfit + 1 )jj = jj + 1
          do j = 1, naxes( 2 )
            if( bitpix .eq. -32 )then
              do i = 1, naxes( 1 )
                outline( i ) = work( i, j, jj )
              end do
              call ftppne( unitM, 0, fpixel, naxes( 1 ), outline,
     +                     nullval, status )
            else if( bitpix .eq. -64 )then
              do i = 1, naxes( 1 )
                outlin2( i ) = work( i, j, jj )
              end do
              call ftppnd( unitM, 0, fpixel, naxes( 1 ), outlin2,
     +                     nullval, status )
            else if( bitpix .eq. 16 )then
              do i = 1, naxes( 1 )
                outl16( i ) = work( i, j, jj )
              end do
              call ftppnj( unitM, 0, fpixel, naxes( 1 ), outl16,
     +                     nullval, status )
            end if
            fpixel = fpixel + naxes( 1 )
          end do
        end do
c add some useful header information
        call ftpcom( unitM, '*********************************',
     +               status )
        write( file, '( f5.3 )' )version
        call ftpcom( unitM, ' FITTED CUBE CREATED BY MAKEMAP v' //
     +               file( 1:5 ), status )
        call ftpcom( unitM, '*********************************',
     +               status )
c record date and time
        call timestamp( file )
        k = lnblnk( file )
        call ftpcom( unitM,
     +               'Date cube created: ' // file( 1:k ), status )
c file names
        file = infile
        call ftpkys( unitM, 'INFILE', file,
     +               'Input text file', status )
        file = outpfile
        call ftpkys( unitM, 'OUTFILE', file,
     +               'Root name for output files', status )
c put some useful text in the string file
        file = ' '
        k = lnblnk( file )
        if( k .gt. 0 )call ftpkys( unitM, 'COMPS', file,
     +                 'Fit parameters', status )
        k = 0
        if( mfit .gt. 4 )k = mfit - 2
        call ftpkyj( unitM, 'FIT ', k,
     +               'Order of GH expansion', status )
        call ftpkyj( unitM, 'BOOTS', nboot,
     +               'number of bootstrap iterations per line', status )
        call ftpkye( unitM, 'NOISE', noise, 4,
     +               'estimated from edge planes of data cube', status )
        call ftpkye( unitM, 'MIN_I', imin, 1,
     +         'times noise: min intensity required for a fit', status )
        call ftpkye( unitM, 'StoN', thresh, 1,
     +               'Min S/N of fitted lines', status )
        x = wmax * abs( vval( 1 ) - vval( 2 ) )
        if( VELmps )x = 1.e-3 * x
        call ftpkye( unitM, 'VWIDTH', x, 2,
     +             'Max velocity width (km/s) of fitted lines', status )
        if( mfit .gt. 4 )then
          x = hmax
          call ftpkye( unitM, 'MAXH', x, 2,
     +                 'Max abs values of h3 and h4', status )
        end if
        i = lnblnk( outfits )
        print *, 'Wrote:  ', outfits( 1:i )
        call ftclos( unitM, status )
        call ftfiou( unitM, status )
      end if
      call ftclos( unitD, status )
      call ftfiou( unitD, status )
      if( status .gt. 0 )call printerror( status )
      return
      end
