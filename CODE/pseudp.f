      subroutine pseudp( jdum )
c This subroutine uses a bootstrap technique to generate a new realization
c   of the velocity field to be fitted
c Creates areas with constant fractional residuals as described in SS07
c
c   Originally written by KS 2006
c   Debugged & polished by JAS May 2009 & October 2009
c   Alternative option pseudr created by RZS September 2009
c   Converted to 2-D pixel map and largely rewritten by JAS March 2011
c   Debugged again by JAS July 2011
c
      include 'commons.h'
c
c calling argument
      integer jdum
c
c externals
c      integer nearpt
      real*8 ran1_dbl
c
c local arrays
      integer newx( mapx, mapy ), newy( mapx, mapy )
      logical dflag( mapx, mapy )
      real newx2res( mapx, mapy )
c
c local variables
      integer ix, iy, j, jind, jx, jy, k, kx, ky
      logical lp, in
      real frac, sysv, val
c
      sysv = 0
      if( lvels )sysv = vsys
      if( lsystemic )sysv = 0
c swap points at random
      print*, xrange, yrange
      do iy = 1, yrange
        do ix = 1, xrange
          lp = .false.
c ensure pixel selected to swap is good
          do while ( .not. lp )
            jx = int( ran1_dbl( jdum ) * xrange ) + 1
            jy = int( ran1_dbl( jdum ) * yrange ) + 1
            lp = lgpix( jx, jy )
          end do
c save these, set new residual, and mark all pixels as independent
          newx( ix, iy ) = jx
          newy( ix, iy ) = jy
          newx2res( ix, iy ) = x2res( jx, jy )
          dflag( ix, iy ) = .false.
        end do
      end do
c set jind - junc is real!!
      jind = nint( junc )
c re-create quasi-coherent residuals if jind > 1
      if( jind .gt. 1 )then
        frac = 1. / real( jind )
c raise flag for randomly chosen independent points
        do iy = 1, yrange
          do ix = 1, xrange
            if( lgpix( ix, iy ) )then
              val = ran1_dbl( jdum )
              if( val .gt. frac )dflag( ix, iy ) = .true.
            end if
          end do
        end do
c copy residuals of dependent points from around source of independent point
        do iy = 1, yrange
          do ix = 1, xrange
            if( lgpix( ix, iy ) .and. dflag( ix, iy ) )then
c the ibl & jbl arrays are preset to cycle through ever more distant points
              lp = .true.
              j = 2
              k = 2
c find jind-1 neighboring pixels surrounding the flagged point
              do while ( lp )
c select the pixel neighboring that which was the source of the residual
                jx = newx( ix, iy ) + ibl( j ) * istepout
                in = ( jx .gt. 0 ) .and. ( jx .le. xrange )
                jy = newy( ix, iy ) + jbl( j ) * istepout
                in = in .and. ( jy .gt. 0 ) .and. ( jy .le. yrange )
                if( lgpix( jx, jy ) .and. in )then
c identify the neighboring pixel to receive the correlated residual
                  kx = max( ix + ibl( k ) * istepout, 1 )
                  kx = min( kx, xrange )
                  ky = max( iy + jbl( k ) * istepout, 1 )
                  ky = min( ky, yrange )
c assign the correlated residual to the appropriate pixel
                  if( lgpix( kx, ky ) .and. .not. dflag( kx, ky ) )
     +                              newx2res( kx, ky ) = x2res( jx, jy )
c keep going until the patch is complete
                  k = k + 1
                  lp = k .le. jind
                end if
                j = j + 1
c blur list may run out before the patch is complete
c   but previously assigned values of newx2res can suffice
                lp = j .le. mblur
              end do
            end if
          end do
        end do
      end if
c generate realization of velocity field with new errors
      do iy = 1, yrange
        do ix = 1, xrange
          if( lgpix( ix, iy ) )then
            sdat( ix, iy ) = model( ix, iy ) + sysv
     +                            + newx2res( ix, iy ) * sigma( ix, iy )
          end if
        end do
      end do
      return
      end
