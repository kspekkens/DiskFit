      real*8 function func( p )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c Function minimized by the call to Powell
c   Routine originally written by EIB
c   Modified extensively for non-axisymmetric photometry by JAS November 2004
c   Tweaked by ASR
c   Adapted for non-axisymmetric velocities by KS 2006
c   Revised again by JAS in consultation with RZS June 2009
c   Much of the work delegated to subroutine setwgt - JAS October 2009
c   Photometry fitting option resurrected - JAS March 2011
c   Simple warp model added - JAS July 2011
c   Calls to setwgtx moved to setwgt - JAS March 2012
c   Updated to f90 - JAS Jan 2015
c   Set an upper bound on bulge r_e of 10x initial guess - JAS Mar 2017
c   Increased Avalue from 400 to 4000 - JAS Aug 2017
c   Made chi^2 walls cumulative, and improved verbose option - JAS May 2020
c   Moved chi^2 wall on rwarp in to half axis of first ellipse - JAS Aug 2020
c
c real*8 precision needed for matrix arithmetic
c
c common block
      include 'commons.h'
c
c calling argument
      real*8 p( md )
c
c local arrays
      integer, allocatable :: ipiv( : )
      logical, allocatable :: nskip( : )
      real*8, allocatable :: rhs( : )
      real*8, allocatable :: wa( :, : )
      real*8, allocatable :: work( : )
c
c local variables
      character*120 outh, outv
      integer i, ii, ij, ik, iuse, ix, iy, j, jj, k, kk, km, rank !, info
      logical creath, flag
      real best
      real*8 Avalue, rcond, re_ini, resd, sysv, tmod
      save best, iuse, outh, re_ini
c
      parameter ( Avalue = 400 )
      data iuse / -1 /, best / 1.e6 /
c
c flag first call
      creath = iuse .eq. -1
c initialize function value
      func = best
      flag = .false.
c build output string
      k = 0
c output calling parameter list and check for reasonalble values
      i = 0
      if( lpa )then
        i = 1
        if( creath )write( outh( k+1:k+10 ), '( 3x, a7 )' )'diskpa:'
        write( outv( k+1:k+10 ), '( f10.4 )' )
     +                                        p( i ) * 180.d0 / pi - 90.
        k = k + 10
c prevent fit wandering outside 0 < disk PA < 2pi
        if( abs( p( i ) - pi ) .gt. pi )then
          func = func + Avalue * ( 1. + ( abs( p( i ) - pi ) - pi )**2 )
        end if
      end if
      if( leps )then
        i = i + 1
        if( creath )write( outh( k+1:k+10 ), '( 3x, a7 )' )'diskel:'
        write( outv( k+1:k+10 ), '( f10.6 )' )p( i )
        k = k + 10
c prevent fit wandering outside 0.05 < disk ellipticity < 0.95
        if( p( i ) .lt. 0.05d0 )then
          func = func + Avalue * ( 1. + ( p( i ) - .05 )**2 )
        else if( p( i ) .gt. 0.95d0 )then
          func = func + Avalue * ( 1. + ( p( i ) - .95 )**2 )
        end if
      end if
      if( lcentre )then
        if( creath )write( outh( k+1:k+20 ), '( 5x, a15 )' )
     +                                                 'xcen:     ycen:'
        write( outv( k+1:k+20 ), '( 2f10.4 )' )p( i + 1 ), p( i + 2 )
        k = k + 20
c prevent fit wandering more than 50 pixels from the estimated center
        resd = ( xcen - p( i + 1 ) )**2 + ( ycen - p( i + 2 ) )**2
        if( resd .gt. 2500.d0 )func = func +
     +                     Avalue * ( 1. + ( sqrt( resd ) - 50.d0 )**2 )
        i = i + 2
      end if
      if( lvels )then
        if( lnax )then
          if( lphib )then
            i = i + 1
            if( creath )write( outh( k+1:k+10 ), '( 4x, a6 )' )'phi_b:'
            write( outv( k+1:k+10 ), '( f10.4 )' )p( i ) * 180.d0 / pi
            k = k + 10
          end if
        end if
c set velocity offset depending on whether we fit for it
        sysv = vsys
        if( lsystemic )sysv = 0
      end if
c options for photometry
      if( lphot )then
        if( lnax )then
          if( lphib )then
            i = i + 1
            if( creath )write( outh( k+1:k+10 ), '( 3x, a7 )' )'bar-pa:'
            write( outv( k+1:k+10 ),
     +                           '( f10.4 )' )p( i ) * 180.d0 / pi - 90.
            k = k + 10
          end if
          if( lepsb )then
            i = i + 1
            if( creath )write( outh( k+1:k+10 ), '( 3x, a7 )' )'bar-el:'
            write( outv( k+1:k+10 ), '( f10.4 )' )p( i )
            k = k + 10
c prevent fit wandering outside 0.05 < bar ellipticity < 0.95
            if( p( i ) .lt. 0.05d0 )then
              func = func + Avalue * ( 1. + ( p( i ) - .05 )**2 )
            else if( p( i ) .gt. 0.95d0 )then
              func = func + Avalue * ( 1. + ( p( i ) - .95 )**2 )
            end if
          end if
        end if
        if( lbulge )then
          if( lbleps )then
            i = i + 1
            if( creath )write( outh( k+1:k+10 ), '( 3x, a7 )' )'blg-el:'
            write( outv( k+1:k+10 ), '( f10.4 )' )p( i )
            k = k + 10
c KS: prevent fit wandering outside 0. < bulge ellipticity < 0.95
            if( p( i ) .lt. 0.d0 )then
              func = func + Avalue * ( 1. + p( i )**2 )
            else if( p( i ) .gt. 0.95d0 )then
              func = func + Avalue * ( 1. + ( p( i ) - .95 )**2 )
            end if
          end if
          if( lsersn )then
            i = i + 1
            if( creath )write( outh( k+1:k+10 ), '( 4x, a6 )' )'blg-n:'
            write( outv( k+1:k+10 ), '( f10.4 )' )p( i )
            k = k + 10
c prevent fit from wandering to absurdly small Sersic indices
            if( p( i ) .lt. 0.17 )then
              func = func + Avalue * ( 1. + ( p( i ) - .17 )**2 )
            end if
c prevent fit from wandering to absurdly large Sersic indices
            if( p( i ) .gt. 10. )then
              func = func + Avalue * ( 1. + ( p( i ) - 10. )**2 )
            end if
          end if
          if( lr_e )then
            i = i + 1
            if( creath )write( outh( k+1:k+10 ), '( 4x, a6 )' )'r-blg:'
            write( outv( k+1:k+10 ), '( f10.4 )' )p( i )
            k = k + 10
c prevent fit wandering to negative r_eff
            if( p( i ) .lt. 0.05 )then
              func = func + Avalue * ( 1. + ( p( i ) - .05 )**2 )
            end if
c prevent fit from wandering to very large values
            if( creath )re_ini = p( i )            
            flag = p( i ) .gt. 9. * re_ini
            if( p( i ) .gt. 10. * re_ini )then
              func = func + Avalue * ( 1. + ( p( i ) - re_ini )**2 )
            end if
          end if
        end if
      end if
c simple warp model
      if( lwarp )then
c inner radius of warp
        if( lrwarp )then
          i = i + 1
          if( creath )write( outh( k+1:k+10 ), '( 4x, a6 )' )'rwarp:'
          write( outv( k+1:k+10 ), '( f10.4 )' )p( i )
          k = k + 10
c prevent rwarp from wandering outside allowed range
          if( p( i ) .lt. .5 * sma( 1 ) )then
           func = func + Avalue * ( 1. + ( p( i ) - .5 * sma( 1 ) )**2 )
          else if( p( i ) .gt. sma( nellip ) )then
            func = func +
     +               Avalue * ( 1. + ( p( i ) - sma( nellip ) )**2 )
          end if
        end if
c maximum eccentricity of warp
        if( lwepsm )then
          i = i + 1
          if( creath )write( outh( k+1:k+10 ), '( 4x, a6 )' )'wepsm:'
          write( outv( k+1:k+10 ), '( f10.4 )' )p( i )
          k = k + 10
c prevent warp parameters from wandering outside reasonable ranges
          if( abs( p( i ) ) .gt. .5 )then
            func = func + Avalue * ( 1. + ( abs( p( i ) ) - .5 )**2 )
          end if
        end if
c max position angle of warp
        if( lwpm )then
          i = i + 1
          if( creath )write( outh( k+1:k+10 ), '( 4x, a6 )' )'wphim:'
          write( outv( k+1:k+10 ), '( f10.4 )' )p( i )
          k = k + 10
c prevent warp parameters from wandering outside reasonable ranges
          if( abs( p( i ) ) .gt. 1. )then
            func = func + Avalue * ( 1. + ( abs( p( i ) ) - 1. )**2 )
          end if
        end if
      end if
c check that the number of parameters set agrees with the number expected
      if( i .ne. nd )then
        print *, 'logical error in number of parameters'
        call crash( 'func' )
      end if
c output header if this is a new iteration
      if( iuse .ne. iter )then
        if( iter .eq. 0 )then
          best = 1.e6
        else
          print '( a, i5, 1x, a, 1pe10.3 )', 'starting iteration', iter,
     +             'best reduced chisq so far:', best
        end if
        print '( a )', outh( 1:k )
        if( .not. verbose )print '( a )', outv( 1:k )
c output warning if needed
        if( flag )then
          print *, 'Warning: bulge effective radius near its maximum' //
     +                         ' allowed value of', 10. * sngl( re_ini )
          print *, '   the minimum could be artificial'
        end if
        iuse = iter
      end if
c print out calling arguments only if verbose
      if( verbose )print '( a )', outv( 1:k )
c there's no point in evaluating func if we have set it already
      if( sngl( func ) - best .gt. 1.e-5 * best )then
        if( verbose )then
          print *, 'skipping function evaluation', sngl( func ), best
        end if
      else
c allocate space
        allocate ( ipiv( ntot ) )
        allocate ( nskip( ntot ) )
        allocate ( rhs( ntot ) )
        allocate ( wa( ntot, ntot ) )
        allocate ( work( ntot ) )
c set the weights
        fronly = .false.
        call setwgt( p, lseeing, km )
c start by assuming all intensities are positive - relevant for photometry only
        do i = 1, ntot
          nskip( i ) = .true.
        end do
c initialize the wa matrix
    1   k = 0
        do kk = 1, ntot
          if( nskip( kk ) )then
            j = k
            k = k + 1
            do jj = kk, ntot
              j = j + 1
              wa( j, k ) = 0
            end do
          end if
        end do
c make the wa matrix
        if( l2D )then
          do iy = 1, yrange
            do ix = 1, xrange
              if( lgpix( ix, iy ) )then
                k = 0
                do kk = 1, ntot
                  if( nskip( kk ) )then
                    ik = 0
                    do i = 1, km
                      if( iwgt( i, ix, iy ) .eq. kk )ik = i
                    end do
                    j = k
                    k = k + 1
                    do jj = kk, ntot
                      if( nskip( jj ) )then
                        j = j + 1
                        if( ik .gt. 0 )then
                          ij = 0
                          do i = 1, km
                            if( iwgt( i, ix, iy ) .eq. jj )ij = i
                          end do
                          if( ij .gt. 0 )then
                            wa( j, k ) = wa( j, k ) +
     +        wgt( ij, ix, iy ) * wgt( ik, ix, iy ) / sigma( ix, iy )**2
                          end if
                        end if
                      end if
                    end do
                  end if
                end do
              end if
            end do
          end do
        else
          do i = 1, inp_pts
            if( lgpix( i, 1 ) )then
              k = 0
              do kk = 1, ntot
                if( nskip( kk ) )then
                  ik = 0
                  do ii = 1, km
                    if( iwgt( ii, i, 1 ) .eq. kk )ik = ii
                  end do
                  j = k
                  k = k + 1
                  do jj = kk, ntot
                    if( nskip( jj ) )then
                      j = j + 1
                      if( ik .gt. 0 )then
                        ij = 0
                        do ii = 1, km
                          if( iwgt( ii, i, 1 ) .eq. jj )ij = ii
                        end do
                        if( ij .gt. 0 )then
                          wa( j, k ) = wa( j, k ) +
     +              wgt( ij, i, 1 ) * wgt( ik, i, 1 ) / sigma( i, 1 )**2
                        end if
                      end if
                    end if
                  end do
                end if
              end do
            end if
          end do
        end if
c fill out symmetric part
        k = 0
        do kk = 1, ntot
          if( nskip( kk ) )then
            j = k
            k = k + 1
            do jj = kk, ntot
              j = j + 1
              wa( k, j ) = wa( j, k )
            end do
          end if
        end do
        rank = k
c set the rhs vector
        j = 0
        do jj = 1, ntot
          if( nskip( jj ) )then
            j = j + 1
            rhs( j ) = 0
          end if
        end do
        if( l2D )then
          do iy = 1, yrange
            do ix = 1, xrange
              if( lgpix( ix, iy ) )then
                j = 0
                do jj = 1, ntot
                  if( nskip( jj ) )then
                    ij = 0
                    do i = 1, km
                      if( iwgt( i, ix, iy ) .eq. jj )ij = i
                    end do
                    j = j + 1
                    if( ij .gt. 0 )then
                      if( lvels )then
                        rhs( j ) = rhs( j ) + wgt( ij, ix, iy ) *
     +                    ( sdat( ix, iy ) - sysv ) / sigma( ix, iy )**2
                      else
                        rhs( j ) = rhs( j ) +
     +           wgt( ij, ix, iy ) * sdat( ix, iy ) / sigma( ix, iy )**2
                      end if
                    end if
                  end if
                end do
              end if
            end do
          end do
        else
          do i = 1, inp_pts
            if( lgpix( i, 1 ) )then
              j = 0
              do jj = 1, ntot
                if( nskip( jj ) )then
                  ij = 0
                  do ii = 1, km
                    if( iwgt( ii, i, 1 ) .eq. jj )ij = ii
                  end do
                  j = j + 1
                  if( ij .gt. 0 )then
                    if( lvels )then
                      rhs( j ) = rhs( j ) + wgt( ij, i, 1 ) *
     +                        ( sdat( i, 1 ) - sysv ) / sigma( i, 1 )**2
                    else
                      rhs( j ) = rhs( j ) +
     +                  wgt( ij, i, 1 ) * sdat( i, 1 ) / sigma( i,1 )**2
                    end if
                  end if
                end if
              end do
            end if
          end do
        end if
        j = 0
        do jj = 1, ntot
          if( nskip( jj ) )then
            j = j + 1
            if( rhs( j ) .eq. 0 )then
              print '( a2, $ )', '-'
            end if
          end if
        end do
c add smoothing penalty
        call smoothp( wa, ntot )
c solve the linear system
CCCC LAPACK:
c      info = 0
c      call f07adf( rank, rank, wa, ntot, ipiv, info )
c      if( info .eq. 0 )then
c        call f07aef( 'N', rank, 1, wa, ntot, ipiv, rhs, ntot, info )
c        if( info .ne. 0 )print *, 'f07aef failed with info =', info
c      else
c        print *, 'the factor U is singular'
c        print '( a2, $ )', '.'
c      end if
CCCC LINPACK:
        rcond = 0
        call dgeco( wa, ntot, rank, ipiv, rcond, work )
        if( abs( rcond ) .gt. 1.d-10 )then
          call dgesl( wa, ntot, rank, ipiv, rhs, 0 )
        else
c          print *, 'the matrix may be singular'
c          print *,
c     +        'this problem could be caused by too closely spaced rings'
c          print *, ' try adjusting the ring spacing'
c          call crash( 'func' )
          print '( a2, $ )', 's '
        end if
c insert zeros for skipped rows
        j = 0
        flag = .false.
        do i = 1, ntot
          if( nskip( i ) )then
            j = j + 1
c nskip is always true for velocities
            fitval( i ) = rhs( j )
c flag new negative values
            nskip( i ) = lvels .or. rhs( j ) .gt. 0.d0
            flag = flag .or. ( .not. nskip( i ) )
          else
            fitval( i ) = 0
          end if
        end do
c loop back if some intensities are still negative
        if( flag )go to 1
c define chisq
        func = 0
        if( l2D )then
          do iy = 1, yrange
            do ix = 1, xrange
              if( lgpix( ix, iy ) )then
c sum fitted ellipse values to find total model value
                tmod = 0
                do j = 1, ntot
                  ij = 0
                  do i = 1, km
                    if( iwgt( i, ix, iy ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )tmod = tmod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
                end do
c difference from observed value
                if( lvels )then
                  resd = ( sdat( ix, iy ) - sysv ) - tmod
                else
c could still have negative intensities in extrapolated region
                  tmod = max( tmod, 0.d0 )
                  resd = sdat( ix, iy ) - tmod
                end if
c
c add this term only for pixels that made some contribution to the weights.
c   In the code above, this is all pixels in the list.  If that code is
c   changed to ignore data anywhere in the fit for the circular velocity,
c   then the ignored pixels must not contribute to the chi^2, since the
c   predicted velocity will be small and the residuals very large
c
                func = func + ( resd / sigma( ix, iy ) )**2
              end if
            end do
          end do
        else
          do i = 1, inp_pts
            if( lgpix( i, 1 ) )then
c sum fitted ellipse values to find total model value
              tmod = 0
              do j = 1, ntot
                ij = 0
                do ii = 1, km
                  if( iwgt( ii, i, 1 ) .eq. j )ij = ii
                end do
                if(
     +            ij .gt. 0 )tmod = tmod + wgt( ij, i, 1 ) * fitval( j )
              end do
c difference from observed value
              if( lvels )then
                resd = ( sdat( i, 1 ) - sysv ) - tmod
              else
                tmod = max( tmod, 0.d0 )
                resd = sdat( i, 1 ) - tmod
              end if
              func = func + ( resd / sigma( i, 1 ) )**2
            end if
          end do
        end if
c find reduced chi^2
        dof = pix_ind - ( ntot + nd )
        func = func / dof
        if( verbose )then
          if( sngl( func ) .lt. best )then
            print *, 'new best value of func', sngl( func ), best
          else
            print *, 'no improvement', sngl( func )
          end if
        end if
      end if
      best = min( best, sngl( func ) )
      return
      end
