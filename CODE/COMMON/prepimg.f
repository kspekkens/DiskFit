      subroutine prepimg( gtdst )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c puts raw image into a more convenient format and sets uncertainties
c
c   Originally written by EIB (named skysub!)
c   Polished for ASR by JAS
c   Renamed and included in velfit by JAS March 2011
c   Find median pixel value - JAS July 2011
c   Added calling argument - JAS Sep 2011
c   Added sparse sampling in outer disk - JAS Jan 2012
c   Used preset mask - JAS March 2012
c   Modified to use an image mask file - JAS June 2013
c   Updated to f90 - JAS Jan 2015
c   Allow for user-supplied uncertainties - JAS Sep 15
c
      include 'commons.h'
c
c calling argument - returns sma of largest active data point
      real gtdst
c
c external - from Numerical Recipes, single precision
      real select
c
c local variables
      integer i, j, k
      real ai, cospa, maxim
      real r, sinpa, xel, xmloc, xr, yel, ymloc, yr
c
c the compiler needs this piece of code for some reason
      if( allocated( dmask ) )then
        print *, 'dmask is allocated'
        print *, size( dmask, 1 )
      else
        if( lmask )print *, 'lmask is .T. but dmask is not allocated'
      end if
c copy the image, or excise selected piece, and find the brightest pixel
      maxim = 0
      k = 0
      do i = 1, xrange
        xval( i ) = xval( i + lox - 1 )
      end do
      allocate ( sdat( xrange, yrange ) )
      allocate ( sigma( xrange, yrange ) )
      do j = 1, yrange
        yval( j ) = yval( j + loy - 1 )
        do i = 1, xrange
          sdat( i, j ) = ldat( i + lox - 1, j + loy - 1 )
          if( lerrfile )sigma( i, j ) =
     +                  ldate( i + lox - 1, j + loy - 1 )
          if( lmask )then
            if( dmask( i + lox - 1, j + loy - 1 ) .ne. 0. )then
              sdat( i, j ) = -9999
              k = k + 1
            end if
          end if
          if( sdat( i, j ) .gt. maxim )then
            maxim = sdat( i, j )
            xmloc = xval( i )
            ymloc = yval( j )
          end if
        end do
      end do
      if( lmask )print *, k, ' pixels masked in prepimg'
      inp_pts = xrange * yrange
c check that brightest pixel is close to the selected image centre
      i = xrange / 2
      j = yrange / 2
      xr = min( xrange, yrange ) / 10
      if( ( abs( xmloc - xval( i ) ) .gt. xr ) .or.
     +    ( abs( ymloc - yval( j ) ) .gt. xr ) )then
        print *, 'Selected region not centered on the brightest pixel?'
c        print *, 'check lox, loy, xrange, yrange in input data'
c        call crash( 'prepimg' )
      end if
c select pixels to be used
      call setmsk( gtdst )
c parameters for sparse sampling
      b_over_a = 1 - regeps
      regpa = ( regpa + 90. ) * pi / 180.
      cospa = cos( regpa )
      sinpa = sin( regpa )
c allocate space
      allocate ( lgpix( xrange, yrange ) )
      allocate ( x2res( xrange * yrange, 1 ) )
c pixct = # pixels in fit.  pix_ind = # of indep. points
      pixct = 0
      pix_ind = 0
      do j = 1, yrange
        yr = yval( j ) - ycen
        do i = 1, xrange
          lgpix( i, j ) = .false.
c work over all but bad pixels
          if( sdat( i, j ) .ge. -4.0 * skysig )then
            xr = xval( i ) - xcen
c select only pixels within the mask
            if( inmask( i, j ) )then
              if( istepout .gt. 1 )then
c sparse sampling in outer image
                xel = xr * cospa + yr * sinpa
                yel = -xr * sinpa + yr * cospa
                r = sqrt( xel**2 + ( yel / b_over_a )**2 )
                if( r .lt. regrad )then
                  lgpix( i, j ) = .true.
                else
                  if( ( mod( i - 1, istepout ) .eq. 0 ) .and.
     +                ( mod( j - 1, istepout ) .eq. 0 ) )then
                    lgpix( i, j ) = .true.
                  else
c blank out data value for unused pixels
                    sdat( i, j ) = -10000
                  end if
                end if
              else
                lgpix( i, j ) = .true.
              end if
              if( lgpix( i, j ) )then
                pixct = pixct + 1
                pix_ind = pix_ind + 1
c save input intensity to find the median
                x2res( pixct, 1 ) = sdat( i, j )
              end if
            end if
          end if
        end do
      end do
c find the median using routine select from Numerical Recipes
      typcl = select( pixct / 2, pixct, x2res )
c if no error file supplied, define the measurement uncertainties for each pixel
      if( .not. lerrfile )then
        do j = 1, yrange
          do i = 1, xrange
            if( lgpix( i, j ) )then
c assume photon counting statistics dominate
              sigma( i, j ) = sqrt( ( sdat( i, j ) + sky ) / gain )
            end if
          end do
        end do
      end if
c increase uncertainties of pixels outside largest ellipse
c   to reduce their weight in the fit
      cospa = cos( pa )
      sinpa = sin( pa )
      do j = 1, yrange
        yr = yval( j ) - ycen
        do i = 1, xrange
          if( lgpix( i, j ) )then
            xr = xval( i ) - xcen
            xel = xr * cospa + yr * sinpa
            yel = -xr * sinpa + yr * cospa
            ai = sqrt( xel**2 + ( yel / ( 1.0 - eps ) )**2 )
            if( ai .gt. sma( nellip ) )then
c scale sigmai by (1+x), where x rises linearly from 0 to 1
c    over the range 1 to 1.05 times smaf
              ai = 20. * ( ai / sma( nellip ) - 1. )
              sigma( i, j ) = sigma( i, j ) * ( 1. + ai )
            end if
          end if
        end do
      end do
c return space that is no longer needed
      deallocate ( ldat )
      if( lerrfile )deallocate ( ldate )
      if( lmask )deallocate ( dmask )
      return
      end
