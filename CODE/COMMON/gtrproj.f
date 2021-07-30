      subroutine gtrproj( xpos, ypos, rproj )
c Copyright (C) 2020, Jerry Sellwood and Kristine Spekkens
c
c finds the deprojected radius in a warped disk of the pixel at xpos, ypos
c   relative to the current disk center
c Most of the code here is copied from setwgtw
c
c   Created - JAS Aug 2020
c
c common block
      include 'commons.h'
c
c calling arguments
      real rproj, xpos, ypos
c
c local array!s
      real r2( mellip )
c
c local variables
      integer j, jg, jl, k
      real ct, frac0, frac1, rg2, rl2, r2max, r2min, rp2, st, theta
      real tiny, x, y
      parameter ( tiny = 1.e-2 )
c
c angle of this pixel from the center
      theta = atan2( ypos, xpos )
      rp2 = xpos**2 + ypos**2
      ct = cos( theta )
      st = sin( theta )
c tabulate radii of rings along this radius vector
      r2max = 0
      r2min = 2 * sma( nellip )**2
      do j = 1, nellip
        x = ct * wcp( j ) + st * wsp( j )
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
        call crash( 'GTRPROJ' )
      end if
c
      rl2 = sqrt( rl2 )
      rp2 = sqrt( rp2 )
      rg2 = sqrt( rg2 )
c weights for interpolation or extrapolation
      if( jl .gt. 0 )then
        if( abs( rg2 - rl2 ) .gt. tiny )then
          frac0 = ( rp2 - rl2 ) / ( rg2 - rl2 )
        else
          frac0 = 1
        end if
        frac1 = 1. - frac0
      else
c interpolate to zero at center inside innermost ring
        frac0 = rp2 / rg2
        frac1 = 0
      end if
c get projected radius
      rproj = frac1 * sma( jl ) + frac0 * sma( jg )
      return
      end
