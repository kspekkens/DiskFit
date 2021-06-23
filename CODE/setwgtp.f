      subroutine setwgtp( xpos, ypos, km, iwt, wt )
c Sets the weight factors for a point at xpos, ypos relative to the
c   current disk center.  Called both from func & from writemod
c
c   Created by JAS by adapting ASR's code for photometry - March 2011
c
c common block
      include 'commons.h'
c
c calling arguments
      integer km, iwt( mk )
      real xpos, ypos, wt( mk )
c
c local variables
      integer i, j, jj, k, kk
      real ai, arg, bn, expon, frac1, frac2, xel, yel
c
c clear the wgt array
      do j = 1, mk
        iwt( j ) = 0
        wt( j ) = 0
      end do
      jj = 0
      kk = 0
      if( ldisk )then
c compute semi-major axis of disk ellipse through this pixel
        xel = xpos * cphi + ypos * sphi
        yel = -xpos * sphi + ypos * cphi
        ai = sqrt( xel**2 + ( yel / b_over_a )**2 )
c find nearest interior ellipse by bi-section - assumed in ascending order
        k = 1
        j = nellip
        do while ( k + 1 .lt. j )
          i = ( k + j ) / 2
          if( sma( i ) .lt. ai )then
            k = i
          else
            j = i
          end if
        end do
c linear interpolation between rings and linear extrapolation to values
c   inside the smallest or outside the largest rings
        frac2 = ( ai - sma( k ) ) / ( sma( k + 1 ) - sma( k ) )
        frac1 = 1. - frac2
c calculate the weights for disk component
        iwt( kk + 1 ) = k
        iwt( kk + 2 ) = k + 1
        wt( kk + 1 ) = frac1
        wt( kk + 2 ) = frac2
        kk = kk + 2
        jj = nellip
      end if
c
      if( lnax )then
c compute semi-major axis of bar ellipse through this pixel
        xel = xpos * cosphib + ypos * sinphib
        yel = -xpos * sinphib + ypos * cosphib
        ai = sqrt( xel**2 + ( yel / b_over_a2 )**2 )
c find nearest interior ellipse by bi-section - assumed in ascending order
        k = nminr
        j = nminr - 1 + nbar
        if( ( ai .gt. sma( k ) ) .and. ( ai .lt. sma( j ) ) )then
          do while ( k + 1 .lt. j )
            i = ( k + j ) / 2
            if( sma( i ) .lt. ai )then
              k = i
            else
              j = i
            end if
          end do
c linear interpolation
          frac2 = ( ai - sma( k ) ) / ( sma( k + 1 ) - sma( k ) )
          frac1 = 1. - frac2
          k = k + 1 - nminr
c calculate the weights for the bar component
          iwt( kk + 1 ) = jj + k
          iwt( kk + 2 ) = jj + k + 1
          wt( kk + 1 ) = frac1
          wt( kk + 2 ) = frac2
          kk = kk + 2
        end if
        jj = jj + nbar
      end if
c bulge weights
      if( lbulge )then
c for 1 < n < 10 Graham 01 gives k = 1.9992n - 0.3271
        bn = 1.9992 * bulgen - 0.3271
        expon = 1. / bulgen
c calculate the bulge weight array
        xel = xpos * cphi + ypos * sphi
        yel = -xpos * sphi + ypos * cphi
        ai = sqrt( xel**2 + ( yel / ( 1. - bulgel ) )**2 )
c calculate weights for the bulge component
        do k = 1, nbulge
          arg = -bn * ( ( ai / rbulge )**expon - 1. )
          arg = min( arg, 20. )
          iwt( kk + 1 ) = jj + k
          wt( kk + 1 ) = exp( arg )
          kk = kk + 1
        end do
      end if
      km = max( kk, km )
      return
      end
