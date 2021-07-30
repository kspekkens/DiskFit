      subroutine groots( nroot, vm )
c  Copyright (C) 2020, Jerry Sellwood
      use aacube
      implicit none
c
c calling arguments
      integer nroot
      real vm
c
c common blocks
c
      include 'comcube.h'
c
      real amp, h3, h4, mean, sigma, cont
      common / fitGH / amp, mean, sigma, h3, h4, cont
c
c externals
      real dLbydv, dL2dv2, GHfit
c
c local variables
      integer ifail, ind, ir, j, k, l, nmin
      logical keep( 5 )
      real a( 5 ), al, am, ar, vd, vmn( 5 ), vmx( 5 )
      real*8 dv, err, v, v1, v2, tol, w( 4 )
c
      nroot = 0
      nmin = 0
c
      j = 0
      am = -10
      v2 = max( mean - 5 * sigma, vval( 1 ) )
      dv = .11 * sigma
      do while
     +      ( sngl( v2 ) .lt. min( mean + 5. * sigma, vval( nsizev ) ) )
        v = v2
        v2 = v2 + dv
        v2 = min( v2, dble( mean + 5. * sigma ) )
        v1 = v2
c set precision
        tol = 1.e-3
        ir = 1
        ind = 1
        ifail = 1
c find zero if one exists in this range
        do while ( ind .gt. 0 )
          call fndzro( v1, v, err, tol, ir, w, ind, ifail )
          err = dLbydv( sngl( v1 ) )
        end do
c check whether a root was found (ifail=4, a pole, is sometimes a root)
        if( ifail .eq. 0 .or. ifail .eq. 4 )then
c select maxima
          if( dL2dv2( sngl( v1 ) ) .lt. 0. )then
            nroot = nroot + 1
            if( nroot .gt. 5 )call crash( 'Too many maxs in GROOTS' )
            vmx( nroot ) = v1
            keep( nroot ) = .true.
c find the greatest max
            a( nroot ) = GHfit( vmx( nroot ) )
            if( a( nroot ) .gt. am )then
              am = a( nroot )
              vm = vmx( nroot )
              j = nroot
            end if
c remember minima
          else
            nmin = nmin + 1
            if( nmin .gt. 5 )call crash( 'Too many mins in GROOTS' )
            vmn( nmin ) = v1
          end if
        end if
      end do
      if( nroot .gt. 1 )then
c check for closely-spaced pairs
        do k = 1, nroot - 1
c find close velocity pairs - assumes they are consecutive roots
          vd = abs( vmx( k ) - vmx( k + 1 ) )
          if( vd .lt. .5 * sigma )then
c discard lesser max of close pairs - that will never be the abs largest
            if( a( k ) .gt. a( k + 1 ) )then
              keep( k + 1 ) = .false.
            else
              keep( k ) = .false.
            end if
          end if
c discard lesser maximum when the intervening minimum is not deep enough
          do ind = 1, nmin
            if( ( vmn( ind ) .gt. vmx( k ) ) .and. 
     +          ( vmn( ind ) .lt. vmx( k + 1 ) ) )then
              ar = GHfit( vmn( ind ) )
c limiting depth is 80% of lower max or 50% of average, whichever is the less
              al = min( a( k ), a( k + 1 ) )
              al = min( 0.5 * ( a( k ) + a( k + 1 ) ), 0.8 * al )
              if( ar .gt. al )then
                if( a( k ) .lt. a( k + 1 ) )then
                  keep( k ) = .false.
                else
                  keep( k + 1 ) = .false.
                end if
              end if
            end if
          end do
        end do
c discard any maxima that are less than 20% of the largest
        do ir = 1, nroot
          if( a( ir ) .lt. 0.2 * am )then
            keep( ir ) = .false.
          end if
        end do
c return velocity of greatest maximum
        vm = vmx( j )
c count surviving roots
        k = 0
        do ir = 1, nroot
          if( keep( ir ) )k = k + 1
        end do
        nroot = max( k, 1 )
c single maximum
      else
        vm = vmx( 1 )
      end if
      if( nroot .eq. 0 )vm = vdef
      return
      end

      real function dLbydv( v )
c  Copyright (C) 2020, Jerry Sellwood
      implicit none
c
c calling argument
      real v
c
c common block
c
      real amp, h3, h4, mean, sigma, cont
      common / fitGH / amp, mean, sigma, h3, h4, cont
c
c local variables
      real alpha, da, dG, G, H3f, H4f, w
      real*8 pi
      parameter ( pi = 3.141592653589793238462643d0 )
c
      w = ( v - mean ) / sigma
      alpha = exp( -.5 * w**2 ) / sqrt( 2. * pi )
      da = -w * alpha
      H3f = ( 2. * w**3 - 3. * w ) / sqrt( 3. )
      H4f = ( 4. * w**4 - 12. * w**2 + 3. ) / sqrt( 24. )
      G = 1 + h3 * H3f + h4 * H4f
      dG = h3 * ( 6 * w**2 - 3 ) / sqrt( 3. ) +
     +     h4 * ( 8 * w**3 - 12 * w ) / sqrt( 6. )
      dLbydv = amp * ( G * da + alpha * dg ) / sigma**2
      return
      end

      real function dL2dv2( v )
c  Copyright (C) 2020, Jerry Sellwood
      implicit none
c
c calling argument
      real v
c
c common block
c
      real amp, h3, h4, mean, sigma, cont
      common / fitGH / amp, mean, sigma, h3, h4, cont
c
c local variables
      real alpha, da, da2, dG, dG2, G, H3f, H4f, w
      real*8 pi
      parameter ( pi = 3.141592653589793238462643d0 )
c
      w = ( v - mean ) / sigma
      alpha = exp( -.5 * w**2 ) / sqrt( 2. * pi )
      da = -w * alpha
      da2 = alpha * ( w**2 - 1 )
      H3f = ( 2. * w**3 - 3. * w ) / sqrt( 3. )
      H4f = ( 4. * w**4 - 12. * w**2 + 3. ) / sqrt( 24. )
      G = 1 + h3 * H3f + h4 * H4f
      dG = h3 * ( 6 * w**2 - 3 ) / sqrt( 3. ) +
     +     h4 * ( 8 * w**3 - 12 * w ) / sqrt( 6. )
      dG2 = h3 * 12 * w / sqrt( 3. ) +
     +      h4 * ( 24 * w**2 - 12 ) / sqrt( 6. )
      dL2dv2 = amp * ( G * da2 + 2 * da * dg + alpha * dG2 ) / sigma**3
      return
      end
