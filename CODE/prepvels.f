      subroutine prepvels( gtdst )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c puts raw data in more convenient format
c
c   Originally written by KS
c   Modified by RZS summer 2009
c   Polished by JAS October 2009
c   Renamed to prepvels by JAS March 2011
c   Calling argument added, select subsection of vely map, and set typcl
c     for smoothing - JAS Sep 2011
c   Implement preset mask - JAS March 2012
c   Updated to f90 - JAS Jan 2015
c
      include 'commons.h'
c
c calling argument - returns sma of largest active data point
      real gtdst
c
c local variables
      integer i, j
c
c 2D image case
      if( l2D )then
c copy the velocity map or extract selected piece
        if( lFITS )then
          allocate ( sdat( xrange, yrange ) )
          allocate ( sdate( xrange, yrange ) )
          inp_pts = 0
          do j = 1, yrange
            yval( j ) = yval( j + loy - 1 )
            do i = 1, xrange
              sdat( i, j )  =  ldat( i + lox - 1, j + loy - 1 )
              sdate( i, j ) = ldate( i + lox - 1, j + loy - 1 )
              lgpix( i, j ) = lgpix( i + lox - 1, j + loy - 1 )
              if( lgpix( i, j ) )inp_pts = inp_pts + 1
            end do
          end do
          do i = 1, xrange
            xval( i ) = xval( i + lox - 1 )
          end do
        end if
c set the mask now xval and yval arrays are right
        call setmsk( gtdst )
c pixct = # pixels in fit. pix_ind = # of indep. points
        allocate ( x2res( xrange, yrange ) )
        allocate ( sigma( xrange, yrange ) )
        pixct = 0
        pix_ind = 0
        do j = 1, yrange
          do i = 1, xrange
c skip bad pixels and points outside mask
            if( lgpix( i, j ) .and. inmask( i, j ) )then
              pixct = pixct + 1
              pix_ind = pix_ind + 1
              x2res( pixct, 1 ) = sdat( i, j )
            else
c flag pixels outside the mask as bad
              lgpix( i, j ) = .false.
            end if
          end do
        end do
c get the measurement uncertainties for each pixel
        do j = 1, yrange
          do i = 1, xrange
            if( lgpix( i, j ) )then
              sigma( i, j ) = sqrt( sdate( i, j )**2 + eISM**2 )
c penalizing data away from major axis - code not used
c              xr = xval( i ) - xcen
c              yr = yval( j ) - ycen
c              xel = xr * cos( maskp ) + yr * sin( maskp )
c              yel = -xr * sin( maskp ) + yr * cos( maskp )
c              tmpdst = sqrt( xel**2 + ( yel / ( 1. - eps ) )**2 )
c              penalty = 5.0 / abs( xel / tmpdst )
c              sigma( i, j ) = sqrt( sdate( i, j )**2 + eISM**2 + penalty**2 )
              if( sigma( i, j ) .eq. 0. )then
                print *, 'At least one sigma_ij in dataset is'
                print *, 'exactly equal to 0: minimization will crash'
                print *, 'Change data or set deltaISM>0 in input file'
                call crash( 'prepvels' )
              end if
            end if
          end do
        end do
      else
c 1D pixel list - xval and yval preset
        call setmsk( gtdst )
        allocate ( lgpix( inp_pts, 1 ) )
        allocate ( sigma( inp_pts, 1 ) )
        allocate ( x2res( inp_pts, 1 ) )
c count pixels
        pixct = 0
        pix_ind = 0
        do i = 1, inp_pts
          lgpix( i, 1 ) = .false.
c theta = acos( abs( xel / tmpdst ) )
          if( inmask( i, 1 ) )then
            lgpix( i, 1 ) = .true.
            pixct = pixct + 1
            pix_ind = pix_ind + 1
            x2res( pixct, 1 ) = sdat( i, 1 )
c set uncertainty
            sigma( i, 1 ) = sqrt( sdate( i, 1 )**2 + eISM**2 )
            if( sigma( i, 1 ) .eq. 0. )then
              print *, 'At least one sigma_i in dataset is'
              print *, 'exactly equal to 0: minimization will crash'
              print *, 'Change data or set deltaISM>0 in input file'
              call crash( 'prepvels' )
            end if
          end if
        end do
      end if
c set typical velocity for smoothing - use Heapsort from Numerical Recipes
      if( ( lambda1 .gt. 0. ) .or. ( lambda2 .gt. 0. ) )then
        call hpsort( pixct, x2res )
        i = pixct / 10
        j = pixct - i
        typcl = .5 * abs( x2res( j, 1 ) - x2res( i, 1 ) )
        if( VELmps )typcl = .001 * typcl
c round down
        if( typcl .gt. 100. )then
          i = 0.1 * typcl
          typcl = 10 * i
        else if( typcl .gt. 10. )then
          i = typcl
          typcl = i
        else
          i = 10. * typcl
          if( i .gt. 0 )typcl =  .1 * real( i )
        end if
      end if
      return
      end
