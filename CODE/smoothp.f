      subroutine smoothp( wa, mt )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c Purpose is to add the penalty for large second derivatives between
c   the values of the mean orbital speed and non-axisymmetric velocities
c   WARNING - assumes sma values are equally spaced, but need not start
c   from the center
c
c Routine written by JAS May 2009
c
c
c common blocks
      include 'commons.h'
c
c calling argument is the matrix
      integer mt
      real*8 wa( mt, mt )
c
c local variables
      integer j, jj, k, kk, iuse
      real da, fac
      save iuse, jj, da
      data iuse / 0 /
c
c skip if no smoothing
      if( ( lambda1 .le. 0. ) .and. ( lambda2 .le. 0. ) )return
c checks
      if( iuse .eq. 0 )then
c determine ring spacing da and check they are equally spaced
        da = ( sma( nellip ) - sma( 1 ) ) / dble( nellip - 1 )
        do j = 2, nellip
          fac = sma( 1 ) + real( j - 1 ) * da
          if( abs( sma( j ) - fac ) .gt. .01 * da )then
            print *, 'Ellipses not equally spaced'
            print *, 'expected sma for ellipse', j, ' is', fac,
     +                              ' actual value is', sma( j )
            call crash( 'SMOOTHP' )
          end if
        end do
c sma of first ellipse determines the smoothing rule at the inner boundary
        if( ( abs( da - sma( 1 ) ) .lt. .01 * da ) .and. lvels )then
c v_t constrained to start from zero
          jj = 1
        else if( ( sma( 1 ) .gt. 1.01 * da ) .or. lphot )then
c allow first value to float freely
          jj = 3
        else
          print *,
     +           'First ellipse is closer to zero than the mean spacing'
          print *, 'Need more code to handle this case'
          call crash( 'SMOOTHP' )
        end if
c guess a typical velocity in km/s or intensity in ADU only if none is preset
        if( typcl .le. 0. )typcl = 100
        if( lvels )then
          print *, 'typical velocity adopted is', typcl,  ' km/s'
        else
          print *, 'typical intensity is', typcl,  ' ADU'
        end if
        iuse = 1
      end if
c axisymmetric terms
      if( lambda1 .gt. 0. )then
c scale out ring spacing
        fac = 2. * lambda1 * real( nellip**4 ) / typcl**2
        if( jj .eq. 3 )then
c inner two ellipses
          j = 1
          wa( j, j     ) = wa( j, j     ) +      fac
          wa( j, j + 1 ) = wa( j, j + 1 ) - 2. * fac
          wa( j, j + 2 ) = wa( j, j + 2 ) +      fac
          j = 2
          wa( j, j - 1 ) = wa( j, j - 1 ) - 2. * fac
          wa( j, j     ) = wa( j, j     ) + 5. * fac
          wa( j, j + 1 ) = wa( j, j + 1 ) - 4. * fac
          wa( j, j + 2 ) = wa( j, j + 2 ) +      fac
        end if
        do j = jj, nellip - 2
c constrain curve to be antisymmetric about the origin
          if( j .eq. 1 )then
            wa( 1, 1 ) = wa( 1, 1 ) + 5. * fac
          else if( j .eq. 2 )then
            wa( 2, 1 ) = wa( 2, 1 ) - 4. * fac
            wa( 2, 2 ) = wa( 2, 2 ) + 6. * fac
          else
c normal case
            wa( j, j - 2 ) = wa( j, j - 2 ) +      fac
            wa( j, j - 1 ) = wa( j, j - 1 ) - 4. * fac
            wa( j, j     ) = wa( j, j     ) + 6. * fac
          end if
          wa( j, j + 1 ) = wa( j, j + 1 ) - 4. * fac
          wa( j, j + 2 ) = wa( j, j + 2 ) +      fac
        end do
c outer two ellipses
        j = nellip - 1
        wa( j, j - 2 ) = wa( j, j - 2 ) +      fac
        wa( j, j - 1 ) = wa( j, j - 1 ) - 4. * fac
        wa( j, j     ) = wa( j, j     ) + 5. * fac
        wa( j, j + 1 ) = wa( j, j + 1 ) - 2. * fac
        j = nellip
        wa( j, j - 2 ) = wa( j, j - 2 ) +      fac
        wa( j, j - 1 ) = wa( j, j - 1 ) - 2. * fac
        wa( j, j     ) = wa( j, j     ) +      fac
      end if
      if( lambda2 .gt. 0. )then
c smooth axisymmetric radial velocity terms
c   without the constraint that these curves should pass through the origin
        if( lradial .and. ( nradial .gt. 4 ) )then
c allow for a different lambda
          fac = 2. * lambda2 * real( nellip**4 ) / typcl**2
c inner two ellipses
          k = nellip + 1
          wa( k, k     ) = wa( k, k     ) +      fac
          wa( k, k + 1 ) = wa( k, k + 1 ) - 2. * fac
          wa( k, k + 2 ) = wa( k, k + 2 ) +      fac
          k = nellip + 2
          wa( k, k - 1 ) = wa( k, k - 1 ) - 2. * fac
          wa( k, k     ) = wa( k, k     ) + 5. * fac
          wa( k, k + 1 ) = wa( k, k + 1 ) - 4. * fac
          wa( k, k + 2 ) = wa( k, k + 2 ) +      fac
c generic ellipses
          do j = 3, nradial - 2
            k = nellip + j
            wa( k, k - 2 ) = wa( k, k - 2 ) +      fac
            wa( k, k - 1 ) = wa( k, k - 1 ) - 4. * fac
            wa( k, k     ) = wa( k, k     ) + 6. * fac
            wa( k, k + 1 ) = wa( k, k + 1 ) - 4. * fac
            wa( k, k + 2 ) = wa( k, k + 2 ) +      fac
          end do
c outer two ellipses
          k = nellip + nradial - 1
          wa( k, k - 2 ) = wa( k, k - 2 ) +      fac
          wa( k, k - 1 ) = wa( k, k - 1 ) - 4. * fac
          wa( k, k     ) = wa( k, k     ) + 5. * fac
          wa( k, k + 1 ) = wa( k, k + 1 ) - 2. * fac
          k = nellip + nradial
          wa( k, k - 2 ) = wa( k, k - 2 ) +      fac
          wa( k, k - 1 ) = wa( k, k - 1 ) - 2. * fac
          wa( k, k     ) = wa( k, k     ) +      fac
        end if
c smooth non-axisymmetric velocity terms, both radial and tangential
c   without the constraint that these curves should pass through the origin
        if( lnax .and. ( nnasymm .gt. 4 ) )then
c allow for a different lambda
          fac = 2. * lambda2 / da**2
          do kk = 1, 2
c inner two ellipses
            k = nellip + ( kk - 1 ) * nnasymm + 1
            wa( k, k     ) = wa( k, k     ) +      fac
            wa( k, k + 1 ) = wa( k, k + 1 ) - 2. * fac
            wa( k, k + 2 ) = wa( k, k + 2 ) +      fac
            k = nellip + ( kk - 1 ) * nnasymm + 2
            wa( k, k - 1 ) = wa( k, k - 1 ) - 2. * fac
            wa( k, k     ) = wa( k, k     ) + 5. * fac
            wa( k, k + 1 ) = wa( k, k + 1 ) - 4. * fac
            wa( k, k + 2 ) = wa( k, k + 2 ) +      fac
c generic ellipses
            do j = 3, nnasymm - 2
              k = nellip + ( kk - 1 ) * nnasymm + j
              wa( k, k - 2 ) = wa( k, k - 2 ) +      fac
              wa( k, k - 1 ) = wa( k, k - 1 ) - 4. * fac
              wa( k, k     ) = wa( k, k     ) + 6. * fac
              wa( k, k + 1 ) = wa( k, k + 1 ) - 4. * fac
              wa( k, k + 2 ) = wa( k, k + 2 ) +      fac
            end do
c outer two ellipses
            k = nellip + kk * nnasymm - 1
            wa( k, k - 2 ) = wa( k, k - 2 ) +      fac
            wa( k, k - 1 ) = wa( k, k - 1 ) - 4. * fac
            wa( k, k     ) = wa( k, k     ) + 5. * fac
            wa( k, k + 1 ) = wa( k, k + 1 ) - 2. * fac
            k = nellip + kk * nnasymm
            wa( k, k - 2 ) = wa( k, k - 2 ) +      fac
            wa( k, k - 1 ) = wa( k, k - 1 ) - 2. * fac
            wa( k, k     ) = wa( k, k     ) +      fac
          end do
        end if
      end if
      return
      end
