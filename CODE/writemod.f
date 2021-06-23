      subroutine writemod
c creates FITS images of the model and residuals of the same size as the
c   input image, with values for all pixels, not just the sampled ones.
c   Uses same formulae as in func.f
c
c    Written by RZS summer 2009
c    Polished by - JAS October 2009
c    Improved by - JAS March 2011
c    Adapted for warp models - JAS July 2011
c    Added output components of the photometry - JAS March 2012
c    Output predicted velocity values for ALL pixels, set nullval to NaN,
c        and set residuals only for pixels used in the fit - JAS March 2012
c    Polished by KS March 2012
c    Allow for 64-bit reals JAS April 2013
c
      include 'commons.h'
c
c externals
      integer lnblnk
c      real nanum
c
c local FITS variables
      integer bitpix, blocksize, fpixel, naxis, readwrite, status
      integer naxes( 2 )
      integer*2 datl16( mpx ), outl16( mpx )
      logical anynull
      real nullval
      real datline( mpx ), outline( mpx )
      real*8 datlin2( mpx ), outlin2( mpx )
c
c local variables
      character*100 datfile, outfits, file
      integer i, ii, iout, jj, k, kk, kout
      integer unitD, unitO
      logical lout
      real diff, val, x, y
c
c create fits files
c      nullval = nanum( 0. )
      nullval = -600000
      status = 0
c reopen the input data file as readonly
      call ftgiou( unitD, status )
      readwrite = 0
      if( lvels )datfile = invfile
      if( lphot )datfile = inpfile
      call ftopen( unitD, datfile, readwrite, blocksize, status )
c determine the type of data
      call ftgidt( unitD, bitpix, status )
c start a big loop over output files
      kout = 2
      if( lphot )kout = 5
      do iout = 1, kout
c construct file names
        k = lnblnk( outroot )
        outfits = outroot( 1:k )
        if( iout .eq. 1 )outfits( k+1:k+3 ) = 'mod'
        if( iout .eq. 2 )outfits( k+1:k+3 ) = 'res'
        lout = iout .lt. 3
        if( iout .eq. 3 .and. ( lnax .or. lbulge ) )then
          outfits( k+1:k+3 ) = 'dsk'
          lout = .true.
        end if
        if( iout .eq. 4 .and. lnax )then
          outfits( k+1:k+3 ) = 'bar'
          lout = .true.
        end if
        if( iout .eq. 5 .and. lbulge )then
          outfits( k+1:k+3 ) = 'blg'
          lout = .true.
        end if
        if( lout )then
          outfits( k+4:k+8 ) = '.fits'
c set output file names to report in output parameter file
          if( iout .eq. 1)outmfile = outfits
          if( iout .eq. 2)outrfile = outfits
c delete the file if it already exists, so we can then recreate it
          call deletefile( outfits, status )
c start the fits writing process
          call ftgiou( unitO, status )
          blocksize = 1         ! needed for historical reasons.
          call ftinit( unitO, outfits, blocksize, status )
c create minimal required header.
          naxis = 2
          naxes( 1 ) = nsizex
          naxes( 2 ) = nsizey
          call ftphps( unitO, bitpix, naxis, naxes, status )
c this is the main loop
          do ii = 1, naxes( 2 )
            jj = ii - loy + 1
            fpixel = 1 + ( ii - 1 ) * naxes( 1 )
c get next row
            if( bitpix .eq. -32 )then
              call ftgpve( unitD, 1, fpixel, naxes( 1 ), nullval,
     +                     datline, anynull, status )
            else if( bitpix .eq. -64 )then
              call ftgpvd( unitD, 1, fpixel, naxes( 1 ), nullval,
     +                     datlin2, anynull, status )
              do i = 1, naxes( 1 )
                datline( i ) = datlin2( i )
              end do
            else if( bitpix .eq. 16 )then
              call ftgpvj( unitD, 1, fpixel, naxes( 1 ), nullval,
     +                     datl16, anynull, status )
              do i = 1, naxes( 1 )
                datline( i ) = datl16( i )
              end do
            else
              print *, 'Unrecognized value of BITPIX:', bitpix
              call crash( 'writemod' )
            end if
            do i = 1, naxes( 1 )
              kk = i - lox + 1
c set default value for pixels outside the fitted region etc
              outline( i ) = nullval
c skip if this pixel is outside fitted region
              if( inmask( kk, jj ) )then
c use precalculated model values
                if( lvels )then
c save model value and get data value
                  val = model( kk, jj )
                  if( .not. lsystemic )val = model( kk, jj ) + vsys
                  if( VELmps )val = val * 1000.0
                  if( iout .eq. 1 )outline( i ) = val
c                if( iout .eq. 1 )outline( i ) = datline( i )
c compute residual only if data value is reasonable
                  if( iout .eq. 2 )then
                    if( datline( i ) .ne. nullval) then
                      if( Velmps )then
                        diff = .001 * ( datline( i ) - vsys * 1000.0 )
                      else
                        diff = datline( i ) - vsys
                      end if
                      if( abs( diff ) .lt. 500. )then
                        if( bitpix .eq. -32 )then
                          outline( i ) = datline( i ) - val
                        else if( bitpix .eq. -64 )then
                          outlin2( i ) = datlin2( i ) - val
                        end if
                      end if
                    end if
                  end if
                else if( lphot )then
                  if( iout .eq. 1 )outline( i ) = model( kk, jj )
c bad pixels are flagged with a large negative value
                  if( ( iout .eq. 2 ) .and.
     +                ( datline( i ) .gt. -500.0 ) )outline( i ) =
     +                                    datline( i ) - model( kk, jj )
                  if( iout .eq. 3 )outline( i ) = diskint( kk, jj )
                  if( iout .eq. 4 )outline( i ) = barint( kk, jj )
                  if( iout .eq. 5 )outline( i ) = bulgeint( kk, jj )
                  if( bitpix .eq. -64 )outlin2( i ) = outline( i )
                  if( bitpix .eq. 16 )outl16( i ) = outline( i )
                end if
              end if
            end do
c end of this row, write it out using group=0
            if( bitpix .eq. -32 )then
              call ftppne( unitO, 0, fpixel, naxes( 1 ), outline,
     +                     nullval, status )
            else if( bitpix .eq. -64 )then
              call ftppnd( unitO, 0, fpixel, naxes( 1 ), outlin2,
     +                     nullval, status )
            else if( bitpix .eq. 16 )then
              call ftppnj( unitO, 0, fpixel, naxes( 1 ), outl16,
     +                     nullval, status )
            end if
          end do
c add some useful header information
          call ftpcom( unitO, '***************************', status )
          call ftpcom( unitO, ' FILE GENERATED BY DISKFIT ', status )
          call ftpcom( unitO, '***************************', status )
          file = infile
          call ftpkys( unitO, 'INFILE', file,
     +                 'Input text file', status )
          file = outpfile
          call ftpkys( unitO, 'OUTPFILE', file,
     +                 'Output parameter file', status )
c Components included in model:
          file = ' '
          if( lvels )then
            k = 1
            if( ldisk )then
              file(1:4) = 'disk'
              k = 5
            end if
            if( lradial )then
              file(k:k+8) = ' + m=0 flow'
              k = k + 10
            end if
            if( lnax )then
              if( order .eq. 1 )file(k:k+10) = ' + m=1 flow'
              if( order .eq. 2 )file(k:k+10) = ' + m=2 flow'
              k = k + 10
            end if
            call ftpkys( unitO, 'COMPS', file,
     +                   'Model components', status )
          end if
          if( lphot )then
            k = 1
            if( ldisk )then
              file( 1:4 ) = 'disk'
              k = 5
            end if
            if( lbulge )then
              file( k:k+8 ) = ' + bulge'
              k = k + 8
            end if
            if( lnax )then
              file(k:k+5) = ' + bar'
              k = k + 5
            end if
            call ftpkys( unitO, 'COMPS', file,
     +                   'Model components', status )
          end if
c Basic model properties
          x = smaf
          call ftpkyf( unitO, 'SMAF', x, 2,
     +                 'Maximum semimajor axis', status )
          x = pixscale
          call ftpkyf( unitO, 'PIXSCALE', x, 2,
     +                 'Pixel scale (asec/pix)', status )
          if( lcentre )then
            x = newxcen
            y = newycen
            call ftpkyf( unitO, 'XCENOUT', x, 2,
     +                   'best fitting X center', status )
            call ftpkyf( unitO, 'YCENOUT', y, 2,
     +                  'best fitting Y center', status )
          end if
          if( lpa )then
            x = ( newpa * 180. / pi ) - 90.
            call ftpkyf( unitO, 'PA_OUT', x, 2,
     +                   'best fitting PA', status )
          end if
          if( leps )then
            x = neweps
            call ftpkyf( unitO, 'EPS_OUT', x, 2,
     +                   'Best fitting ellipticity', status )
          end if
c Close files
          call ftclos( unitO, status )
          call ftfiou( unitO, status )
          if( status .gt. 0 )call printerror( status )
          i = lnblnk( outfits )
          print *, 'Wrote:  ', outfits( 1:i )
        end if
      end do
c
      call ftclos( unitD, status )
      call ftfiou( unitD, status )
      if( status .gt. 0 )call printerror( status )
      return
      end
