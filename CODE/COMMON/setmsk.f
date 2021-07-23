      subroutine setmsk( gtdst )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c set mask to select the region of image that will be used for the fit
c
c created by JAS - March 2012
c   Updated to f90 - JAS Jan 2015
c   Change to in mask for vels and linter0 - Jul 2015
c
      include 'commons.h'
c
c calling argument
      real gtdst
c
c local variables
      integer i, j
      real ai, cmp, maskein, maskeout, maskp, masksmout, masksmin, smp
      real xel, xr, yel, yr
c
c outer mask boundary larger and a little rounder than the outermost ellipse
      masksmout = 1.1 * smaf
      maskeout = 0.9 * eps
c position angle is the initial guess
      maskp = pa
      cmp = sngl( cos( dble( maskp ) ) )
      smp = sngl( sin( dble( maskp ) ) )
c ellipticities of inner mask boundary - needed for vels only
      if( lvels .and. ( .not. linter0 ) )then
        masksmin = sma0
        maskein = min( 1.05 * eps, 0.9 )
      end if
      gtdst = 0
c 2D image case
      if( l2D )then
        allocate ( inmask( xrange, yrange ) )
        do j = 1, yrange
          yr = yval( j ) - ycen
          do i = 1, xrange
            inmask( i, j ) = .false.
            xr = xval( i ) - xcen
            xel = xr * cmp + yr * smp
            yel = -xr * smp + yr * cmp
            ai = sqrt( xel**2 + ( yel / ( 1. - maskeout ) )**2 )
            if( ai .lt. masksmout )then
              gtdst = max( ai, gtdst )
c mask out inner pixels, if interpolation to zero for vels is not selected
              if( lvels .and. ( .not. linter0 ) )then
                ai = sqrt( xel**2 + ( yel / ( 1. - maskein ) )**2 )
                if( ai .gt. masksmin )inmask( i, j ) = .true.
              else
                inmask( i, j ) = .true.
              end if
            end if
          end do
        end do
c if user supplied an image mask - flag these pixels too
        if( lmask )then
          do j = 1, yrange
            do i = 1, xrange
              if( sdat( i, j ) .lt. -9998. )inmask( i, j ) = .false.
            end do
          end do
        end if
      else
c 1D pixel list - used only for vels
        allocate ( inmask( inp_pts, 1 ) )
        do i = 1, inp_pts
          xr = xval( i ) - xcen
          yr = yval( i ) - ycen
          xel = xr * cmp + yr * smp
          yel = -xr * smp + yr * cmp
          ai = sqrt( xel**2 + ( yel / ( 1. - maskeout ) )**2 )
          inmask( i, 1 ) = .false.
          if( ai .lt. masksmout )then
            gtdst = max( ai, gtdst )
            if( linter0 )then
              inmask( i, 1 ) = .true.
            else
c skip points inside inner radius, if interpolation to zero is not selected
              ai = sqrt( xel**2 + ( yel / ( 1. - maskein ) )**2 )
              if( ai .gt. masksmin )inmask( i, 1 ) = .true.
            end if
          end if
        end do
      end if
      return
      end
