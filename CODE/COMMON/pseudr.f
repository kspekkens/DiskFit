      subroutine pseudr( jdum )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c Uses a bootstrap technique to generate a new realization of the velocity
c   field to be fitted - an alternative to pseudp written by KS
c Creates a new velocity map by adding errors to the model velocity
c   distribution from a random rotation and radial scaling of the original
c   x2res map
c
c   Created by RZS August 2009
c   Polished by JAS October 2009
c   Updated to f90 - JAS Jan 2015
c
      include 'commons.h'
c
c calling argument
      integer jdum
c
c externals
      integer neari
      real*8 ran1_dbl
c
c local arrays
      real annwidth( 3 ), cosrpa( 3 ), sinrpa( 3 ), rini( 3 ), rout( 3 )
      real, allocatable :: newx2res( :, : )
c
c local variables
      integer i, ix, iy, j, jind, jx, jy, kx, ky, nreg
      logical lp
      real cospa, fring, meps, newr, rpa, rrad, rtmp, sinpa, sysv
      real x, xr, xel, y, yr, yel
c
c this routine assumes a 2D map
      if( .not. l2D )then
        print *,
     +    'This boostrap not available for text format: set junc >0'
        call crash( 'PSEUDR' )
      end if
c
      sysv = 0
      if( lvels )sysv = vsys
      if( lsystemic )sysv = 0
      fring = sma0
      if( linter0 )fring = 0
c determine regions. If bar 1, 2 or 3. If no bar, just one.
      nreg = 1
      rini( 1 ) = fring
      if( lnax )then
        if( ( nminr .eq. 1 ) .neqv. ( nmaxr .eq. nellip ) )then
          nreg = 2
          if( nminr .eq. 1 ) then
            rout( 1 ) = sma( nmaxr )
            rini( 2 ) = rout( 1 )
          else
            rout( 1 ) = sma( nminr )
            rini( 2 ) = rout( 1 )
          end if
        else if(
     +      .not. ( ( nminr .eq. 1 ) .and. ( nmaxr .eq. nellip ) ) )then
          nreg = 3
          rout( 1 ) = sma( nminr )
          rini( 2 ) = rout( 1 )
          rout( 2 ) = sma( nmaxr )
          rini( 3 ) = rout( 2 )
        end if
      end if
      rout( nreg ) = sma( nellip )
c width of each annulus is the total size of radial range
      do i = 1, nreg
        annwidth( i ) = rout( i ) - rini( i )
      end do
c write out region boundaries
c      print *, 'regions, boundaries & random angle:'
      do i = 1, nreg
c create random PA angles for each annulus/region
        rpa = ran1_dbl( jdum ) * 2.d0 * pi
c        print '( i3, 3f10.2)', i, rini( i ), rout( i ), 180. * rpa / pi
        cosrpa( i ) = cos( rpa )
        sinrpa( i ) = sin( rpa )
      end do
c local disk geometric parameters
      meps = 0.9 * eps
      cospa = cos( pa )
      sinpa = sin( pa )
c find new random beginning of annulus
      rrad = ran1_dbl( jdum )
c      print '( a, f10.3 )', 'new random radius', rrad
c initialize newx2res map
      allocate ( newx2res( xrange, yrange ) )
      do iy = 1, yrange
        do ix = 1, xrange
          newx2res( ix, iy ) = x2res( ix, iy )
        end do
      end do
c rotate each annulus in x2res by a random angle
      do iy = 1, yrange
        do ix = 1, xrange
          if( lgpix( ix, iy ) )then
c deproject current pixel
            yr = yval( iy ) - ycen
            xr = xval( ix ) - xcen
            xel =    xr * cospa + yr * sinpa
            yel = ( -xr * sinpa + yr * cospa ) / ( 1.d0 - meps )
c find which annulus and region the pixel is in
            rtmp = sqrt( xel**2.0 + yel**2.0 )
c region 1 is the default
            j = 1
            do jind = 1, nreg
              if( rtmp .gt. rini( jind ) )j = jind
            end do
c rotate pixel through a random angle
            x =  xel * cosrpa( j ) + yel * sinrpa( j )
            y = -xel * sinrpa( j ) + yel * cosrpa( j )
c shift pixel in radius
            newr = rtmp + rrad * ( rout( j ) - rini( j ) )
            if( newr .gt. rout( j ) )newr = newr + rini( j ) - rout( j )
c compute projected position (recycle xr, yr variables)
            rtmp = max( rtmp, 0.01 )
            x = ( newr / rtmp ) * x
            y = ( newr / rtmp ) * y * ( 1.d0 - meps )
            xr = x * cospa - y * sinpa + xcen
            yr = y * cospa + x * sinpa + ycen
c the position of the closest pixel in the original map
            jx = neari( xval, xrange, xr )
            jy = neari( yval, yrange, yr )
            lp = .true.
            j = 2
            do while ( lp )
c the ibl & jbl arrays are preset to cycle through ever more distant points
              kx = max( jx + ibl( j ), 1 )
              kx = min( kx, xrange )
              ky = max( jy + jbl( j ), 1 )
              ky = min( ky, yrange )
              if( lgpix( kx, ky ) )then
                newx2res( ix, iy ) = x2res( kx, ky )
                lp = .false.
              end if
              j = j + 1
c blur list ran out before a dependent point was found
c   but the previously assigned value of newx2res can suffice
              if( j .gt. mblur )lp = .false.
            end do
          end if
        end do
      end do
c generate realization of velocity field + errors:
      do iy = 1, yrange
        do ix = 1, xrange
          if( lgpix( ix, iy ) )then
            sdat( ix, iy ) = model( ix, iy ) + sysv
     *                            + newx2res( ix, iy ) * sigma( ix, iy )
          end if
        end do
      end do
      return
      end
