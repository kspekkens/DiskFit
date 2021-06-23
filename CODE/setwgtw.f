      subroutine setwgtw( xpos, ypos, km, iwt, wt )
c Sets the weight factors for a point at xpos, ypos relative to the
c   current disk center.  Called both from func & from writemod
c   This version differs from setwgtv because it allows for a warp
c
c   Created for a simple warp model - JAS July 2011
c
c common block
      include 'commons.h'
c
c calling arguments
      integer km, iwt( mk )
      real xpos, ypos, wt( mk )
c
c local arrays
      real r2( mellip ), cst( mellip )
c
c local variables
      integer j, jg, jl, k
      real ct, frac0, frac1, rg2, rl2, r2max, r2min, rp2, st, theta
      real tiny, x, y
      parameter ( tiny = 1.e-2 )
c
c clear the wgt array
      do j = 1, mk
        iwt( j ) = 0
        wt( j ) = 0
      end do
c angle of this pixel from the center
      theta = atan2( ypos, xpos )
      rp2 = xpos**2 + ypos**2
      ct = cos( theta )
      st = sin( theta )
c tabulate projected radii of rings along this radius vector
      r2max = 0
      r2min = 2 * sma( nellip )**2
      do j = 1, nellip
        x = ct * wcp( j ) + st * wsp( j )
        cst( j ) = x
        y = ( st * wcp( j ) - ct * wsp( j ) ) / wba( j )
        r2( j ) = sma( j )**2 / ( x**2 + y**2 )
        r2max = max( r2max, r2( j ) )
        r2min = min( r2min, r2( j ) )
      end do
c point on or outside largest ellipse
      if( rp2 .gt. r2max - tiny )then
        rg2 = r2max
        rl2 = 0
        do j = 1, nellip
          if( abs( r2( j ) - r2max ) .lt. tiny )then
            jg = j
          else
            if( r2( j ) .gt. rl2 )then
              rl2 = r2( j )
              jl = j
            end if
          end if
        end do
c point on or inside innermost ellipse
      else if( rp2 .lt. r2min + tiny )then
        if( linter0 )then
c interpolate to zero
          jl = 0
          rl2 = 0
          do j = 1, nellip
            if( r2( j ) .eq. r2min )jg = j
          end do
          rg2 = r2min
        else
c find two smallest ellipses
          rg2 = r2max
          rl2 = r2min
          do j = 1, nellip
            if( r2( j ) .eq. r2min )then
              jl = j
            else
              if( r2( j ) .lt. rg2 )then
                rg2 = r2( j )
                jg = j
              end if
            end if
          end do
        end if
      else
c find the closest ellipses on either side of this point
        jg = nellip
        jl = 1
        rl2 = 0
        rg2 = 2 * sma( nellip )**2
        do j = 1, nellip
          if( r2( j ) .gt. rp2 )then
            if( r2( j ) .lt. rg2 )then
              rg2 = r2( j )
              jg = j
            end if
          else
            if( r2( j ) .gt. rl2 )then
              rl2 = r2( j )
              jl = j
            end if
          end if
        end do
      end if

      k = abs( jg - jl )
      if( k .eq. 0 )then
        print *, 'selected rings', jl, jg, rp2
        jl = max( jl - 5, 1 )
        jg = min( jg + 5, nellip )
        if( jl .gt. jg )then
          j = jl
          jl = jg
          jg = j
        end if
        do j = jl, jg
          print *, j, r2( j )
        end do
        call crash( 'setwgtw' )
      end if

      rl2 = sqrt( rl2 )
      rp2 = sqrt( rp2 )
      rg2 = sqrt( rg2 )
c weights for interpolation or extrapolation
      if( jl .gt. 0 )then
        frac0 = ( rp2 - rl2 ) / ( rg2 - rl2 )
        frac1 = 1. - frac0
      else
c interpolate to zero at center inside innermost ring
        frac0 = rp2 / rg2
        frac1 = 0
      end if
      iwt( 1 ) = jl
      iwt( 2 ) = jg
      if( fronly )then
        wt( 1 ) = frac1
        wt( 2 ) = frac0
      else
c velocity weights sin(i) * cos(theta), with theta the in-plane angle
        if( jl .gt. 0 )wt( 1 ) =
     +            frac1 * wsi( jl ) * cst( jl ) * rl2 / sma( jl )
        wt( 2 ) = frac0 * wsi( jg ) * cst( jg ) * rg2 / sma( jg )
      end if
      j = 2
c systemic velocity weights - they are all 1
      if( lsystemic )then
        j = j + 1
        iwt( j ) = ntot
        wt( j ) = 1
      end if
      km = max( j, km )
      return
      end
