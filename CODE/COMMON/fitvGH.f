      subroutine fitvGH( dat, n, pars, ifail )
c  Copyright (C) 2020, Jerry Sellwood
      use aacube
      implicit none
c
c calling arguments
      integer n, ifail( 7 ) 
      real dat( n ), pars( 8, 2 )
c
c common blocks
c
      include 'comcube.h'
c
      real amp, h3, h4, mean, sigma, cont
      common / fitGH / amp, mean, sigma, h3, h4, cont
c
c externals
      external GHfit, GHfun2, GHgrd2, GHparm
      real GHfit
c
c local arrays
      integer ia( nsize )
      real bstfit( nsize ), resids( nsize ), smth( nsize )
      real*8 d( 6 ), h34( 6 )
      integer, allocatable, save :: iw(:)
      real*8, allocatable, save :: w(:)
c
c local variables
      integer i, iboot, j, k, m, liw, lw
      integer uiparm
      logical bimod, fail, firstc, nset
      real amn, amx, chisq, dv, f, sg, smin, smx, sn, tgt, v
      real vl, vm, vu
      real*8 a, s2, urparm
      save firstc, smin, liw, lw
c
      data firstc / .true. /
c
      fail = .false.
      bimod = .false.
c find largest value
      amx = -100
      do i = 1, n
        amx = max( dat( i ), amx )
      end do
c skip if maximum is within the noise
      amx = amx / noise
      fail = amx .lt. Imin
      if( fail )then
        ifail( 1 ) = ifail( 1 ) + 1
      else
c select data
        nv = 0
        do i = 1, n
          if( dat( i ) .gt. -100. )then
            nv = nv + 1
            vals( 1, nv ) = vval( i )
            vals( 2, nv ) = dat( i ) / noise
            smth( nv ) = vals( 2, nv )
          end if
        end do
        if( nv .lt. 10 )then
          fail = .true.
          ifail( 1 ) = ifail( 1 ) + 1
        else
c estimate inital parameter values from smoothed data
          call bcsmth( smth, nv, 10 )
          smx = -1000
          do i = 1, nv
            if( smth( i ) .gt. smx )then
              k = i
              smx = smth( i )
            end if
          end do
          dv = abs( vval( 1 ) - vval( 2 ) )
c get full width at half max
          tgt = .5 * smx
          if( k .lt. nv )then
            nset = .true.
            do i = k, nv
              if( nset )then
c increasing i is working away from the max
                if( ( smth( i ) .gt. tgt ) .and.
     +              ( smth( i + 1 ) .lt. tgt ) )then
c linear interpolation
                  f = ( smth( i ) - tgt ) /
     +                                     ( smth( i ) - smth( i + 1 ) )
                  vu = ( 1. - f ) * vals( 1, i ) +
     +                        f   * vals( 1, i + 1 )
                  nset = .false.
                end if
              end if
            end do
c no values below target
            if( nset )vu = vals( 1, nv )
          else
            vu = vals( 1, nv ) + .5 * dv
          end if
          if( k .gt. 1 )then
            nset = .true.
c decreasing i is working away from the max
            do i = k - 1, 1, -1
              if( nset )then
                if( ( smth( i ) .lt. tgt ) .and.
     +              ( smth( i + 1 ) .gt. tgt ) )then
c linear interpolation
                  f = ( smth( i + 1 ) - tgt ) /
     +                                    ( smth( i + 1 ) - smth( i ) )
                  vl = ( 1. - f ) * vals( 1, i + 1 ) +
     +                        f   * vals( 1, i )
                  nset = .false.
                end if
              end if
            end do
c no values below target
            if( nset )vl = vals( 1, 1 )
          else
            vl = vals( 1, 1 ) - .5 * dv
          end if
c set first guess - sigma ~ half-width at half max
          sg = .5 * max( dv, vu - vl )
c fit profile
          if( firstc )then
c set minimum acceptable sigma - 2 channel spacings for a min threshold line
            smin = 2. * dv * thresh
c allocate space for minimizer
            liw = 60
            lw = 71 + mfit * ( mfit + 15 ) / 2
c            lw = 77 + mfit * ( mfit + 17 ) / 2
            allocate ( iw( liw ) )
            allocate ( w( lw ) )
c allocate space for bootstraps if needed
            if( nboot .gt. 0 )then
              allocate ( parms( nboot, mfit + 1 ) )
            end if
            firstc = .false.
          end if
c initialize
          call deflt( 2, iw, liw, lw, w )
c parameters to suppress all output
          do i = 19, 24
            iw( i ) = 0
          end do
          iw( 23 ) = -1
c allow more iterations
          iw( 17 ) = 5000
          iw( 18 ) = mfit * 5000
c set absolute and relative tolerances
          w( 31 ) = 1.d-7
          w( 32 ) = 1.d-5
c insert first guess
          h34( 1 ) = 0
          h34( 2 ) = smx
          h34( 3 ) = .5 * ( vl + vu )
          amn = h34( 3 )
          h34( 4 ) = sg
          h34( 5 ) = 0
          h34( 6 ) = 0
          do i = 1, mfit
            d( i ) = 1
          end do
          iboot = 0
          do while ( iboot .le. nboot )
c flag this as a fresh start
            iw( 1 ) = 12
c least squares fitting routine
            call sumsl( mfit, d, h34, GHfun2, GHgrd2, iw, liw, lw, w,
     +                  uiparm, urparm, GHparm )
            if( iboot .eq. 0 )then
c get reduced chi-square value at minimum
              call GHfun2( mfit, h34, k, s2, uiparm, urparm, GHparm )
              chisq = s2 / real( nv - mfit )
c flag from minimizer
              fail = iw( 1 ) .gt. 7
              if( fail )ifail( 2 ) = ifail( 2 ) + 1
c fitted amplitude is within noise
              if( .not. fail )then
                s2 = h34( 4 )
                if( s2 .eq. 0.d0 )s2 = sg
c intensity values are already scaled by the noise
                sn = abs( h34( 2 ) / ( s2 * sqrt( 2. * pi ) ) )
                fail = sn .lt. thresh
                if( fail )ifail( 3 ) = ifail( 3 ) + 1
              end if
c fitted mean V is outside the vel range
              if( .not. fail )then
                f = h34( 3 )
                fail = ( f .lt. vval( 1 ) ) .or.
     +                 ( f .gt. vval( nsizev ) )
                if( fail )ifail( 4 ) = ifail( 4 ) + 1
              end if
c fitted sigma is too large or too small for a low S/N line
              if( .not. fail )then
                fail = ( sngl( h34( 4 ) ) .gt. wmax * dv ) .or.
     +                 ( sngl( h34( 4 ) ) .lt. smin / amx )
                if( fail )ifail( 5 ) = ifail( 5 ) + 1
              end if
c nonsense values of h3 or h4
              if( .not. fail )then
                fail = ( abs( h34( 5 ) ) .gt. hmax ) .or.
     +                 ( abs( h34( 6 ) ) .gt. hmax )
                if( fail )ifail( 6 ) = ifail( 6 ) + 1
              end if
c best fit values to be returned
              if( .not. fail )then
                do i = 1, mfit
                  pars( i, 1 ) = h34( i )
                end do
                pars( mfit + 2, 1 ) = chisq
              end if
            end if
            if( .not. fail )then
c put current fit values into common block
              cont = h34( 1 )
              amp = h34( 2 )
              mean = h34( 3 )
              sigma = h34( 4 )
              if( mfit .gt. 4 )then
                h3 = h34( 5 )
                if( mfit .eq. 6 )h4 = h34( 6 )
c get velocity of fitted function maximum
                call groots( j, vm )
c test for bi-modality on best fit only
                if( iboot .eq. 0 )bimod = j .gt. 1
              else
                vm = h34( 3 )
              end if
c abandon fit for negative amplitudes or negative sigma (in bootstraps)
              if( ( h34( 2 ) .lt. 0.d0 ) .or. ( h34( 4 ) .lt. 0.d0 )
     +  .or. ( vm .lt. vval( 1 ) ) .or. ( vm .gt. vval( nsizev ) ) )then
                ifail( 7 ) = ifail( 7 ) + 1
                fail = .true.
                do i = 1, mfit
                  pars( i, 1 ) = 0
                end do
                pars( mfit + 1, 1 ) = vdef
                pars( mfit + 2, 1 ) = 0
              end if
c save the velocity of the best fit function maximum
              if( ( .not. fail ) .and.
     +            ( iboot .eq. 0 ) )pars( mfit + 1, 1 ) = vm
            end if
c skip bootstrap iterations if fail flag is up or no bootstraps requested
            if( fail .or. ( nboot .eq. 0 ) )then
              iboot = nboot + 1
            else
              if( iboot .eq. 0 )then
c compute and save residuals
                f = -1
                do i = 1, nv
                  v = vals( 1, i )
                  bstfit( i ) = GHfit( v )
c remember channel of line peak
                  if( bstfit( i ) .gt. f )then
                    m = i
                    f = bstfit( i )
                  end if
                  resids( i ) = vals( 2, i ) - bstfit( i )
                end do
              end if
c save these fitted parameters
              if( iboot .gt. 0 )then
                do i = 1, mfit
                  parms( iboot, i ) = h34( i )
                end do
                parms( iboot, mfit + 1 ) = vm
              end if
c prepare for next fit, unless this is the last
              if( iboot .lt. nboot )then
c restore best guess as first guess
                do i = 1, mfit
                  h34( i ) = pars( i, 1 )
                end do
c rearrange residuals to create new pseudo data
                j = m - nboot / 2 + iboot
                if( iboot .ge. nboot / 2 )j = j + 1
                if( j .le. 0 )j = j + nv
                do i = 1, nv
                  if( j .gt. nv )j = 1
                  vals( 2, i ) = bstfit( i ) + resids( j )
                  j = j + 1
                end do
              end if
              iboot = iboot + 1
c end if ( no bootstraps ) block 
            end if
c end loop over bootstraps
          end do
          if( .not. fail )then
            if( nboot .gt. 0 )then
c compute and save uncertainties
              do j = 1, mfit + 1
c                s2 = 0
c                do i = 1, nboot
c                  s2 = s2 + ( parms( i, j ) - pars( j, 1 ) )**2
c                end do
c                pars( j, 2 ) = sqrt( s2 / dble( nboot ) )
                call biwght( parms( 1, j ), nboot, f, sg )
                pars( j, 2 ) = sg
              end do
c save bi-modality flag
              if( bimod )pars( mfit + 1, 2 ) = -pars( mfit + 1, 2 )
            else
              if( bimod )pars( 2, 1 ) = -pars( 2, 1 )
            end if
          end if
          if( vbs )then
            print *, 'return status for sumsl', iw( 1 ), fail
            if( .not. fail )then
              print '( a, 3i8, f8.2 )',
     +                   'fun calls, grad calls, iterations, and chisq',
     +                                iw( 6 ), iw( 30 ), iw( 31 ), chisq
              print '( a, i8, 4f11.4 )',
     +                                  'guess', nv, noise, smx, amn, sg
              print '( a, i8, 7f11.4 )', '  fit', nv,
     +                                 ( pars( j, 1 ), j = 1, mfit + 1 )
              if( nboot .gt. 0 )print '( a, 8x, 7f12.4 )', ' errs',
     +                                 ( pars( j, 2 ), j = 1, mfit + 1 )
            end if
          end if
c end if ( nv > 10 ) block
        end if
c end if ( max > Imin * noise ) block
      end if
      return
      end

      subroutine GHfun2( m, xc, nf, fc, uiparm, urparm, ufparm )
c  Copyright (C) 2020, Jerry Sellwood
      implicit none
c
c calling arguments
      external ufparm
      integer nf, m, uiparm
      real*8 xc( m ), fc, urparm
c
c common block
c
      include 'comcube.h'
c
c local variables
      integer i
      real*8 alpha, fv, H3, H4, pred, ser, v, w
c
      fc = 0
c sum over line profile
      do i = 1, nv
        v = vals( 1, i )
        fv = vals( 2, i )
        w = ( v - xc( 3 ) ) / xc( 4 )
        alpha = exp( -.5d0 * w**2 ) / sqrt( 2. * pi )
        ser = 1
        if( m .gt. 4 )then
          H3 = ( 2.d0 * w**3 - 3.d0 * w ) / sqrt( 3.d0 )
          ser = ser + xc( 5 ) * H3
          if( m .eq. 6 )then
            H4 = ( 4.d0 * w**4 - 12.d0 * w**2 + 3.d0 ) / sqrt( 24.d0 )
            ser = ser + xc( 6 ) * H4
          end if
        end if
        pred = xc( 1 ) + xc( 2 ) * alpha * ser / xc( 4 )
        fc = fc + ( fv - pred )**2
      end do
c penalty to keep the continuum value near its true value
      fc = fc + 1.d9 * ( xc( 1 ) - cton )**4
      return
      end

      subroutine GHgrd2( m, xc, nf, gc, uiparm, urparm, ufparm )
c  Copyright (C) 2020, Jerry Sellwood
      implicit none
c
c calling arguments
      external ufparm
      integer nf, m, uiparm
      real*8 xc( m ), gc( m ), urparm
c
c common block
c
      include 'comcube.h'
c
c local variables
      integer i, j
      real*8 alpha, fv, H3, H3p, H4, H4p, pred, ser, v, w
c
      do j = 1, m
        gc( j ) = 0
      end do
c sum over line profile
      do i = 1, nv
        v = vals( 1, i )
        fv = vals( 2, i )
        w = ( v - xc( 3 ) ) / xc( 4 )
        alpha = exp( -.5d0 * w**2 ) / sqrt( 2. * pi )
        ser = 1
        if( m .gt. 4 )then
          H3 = ( 2.d0 * w**3 - 3.d0 * w ) / sqrt( 3.d0 )
          H3p = ( 6. * w**2 - 3. ) / sqrt( 3. )
          ser = ser + xc( 5 ) * H3
          if( m .eq. 6 )then
            H4 = ( 4.d0 * w**4 - 12.d0 * w**2 + 3.d0 ) / sqrt( 24.d0 )
            H4p = ( 16. * w**3 - 24. * w ) / sqrt( 24. )
            ser = ser + xc( 6 ) * H4
          end if
        end if
        pred = xc( 1 ) + xc( 2 ) * alpha * ser / xc( 4 )
        gc( 1 ) = gc( 1 ) - 2. * ( fv - pred )
        gc( 2 ) = gc( 2 ) - 2. * ( fv - pred ) * alpha * ser / xc( 4 )
        if( m .eq. 6 )then
          gc( 3 ) = gc( 3 ) - 2. * ( fv - pred ) * xc( 2 ) * alpha *
     +          ( w * ser - xc( 5 ) * H3p - xc( 6 ) * H4p ) / xc( 4 )**2
          gc( 4 ) = gc( 4 ) - 2. * ( fv - pred ) * xc( 2 ) * alpha *
     +          ( ( w**2 - 1 ) * ser - 
     +              w * ( xc( 5 ) * H3p + xc( 6 ) * H4p ) ) / xc( 4 )**2
          gc( 5 ) = gc( 5 ) - 
     +               2. * ( fv - pred ) * xc( 2 ) * alpha * H3 / xc( 4 )
          gc( 6 ) = gc( 6 ) - 
     +               2. * ( fv - pred ) * xc( 2 ) * alpha * H4 / xc( 4 )
        else if( m .eq. 5 )then
          gc( 3 ) = gc( 3 ) - 2. * ( fv - pred ) * xc( 2 ) * alpha *
     +                          ( w * ser - xc( 5 ) * H3p ) / xc( 4 )**2
          gc( 4 ) = gc( 4 ) - 2. * ( fv - pred ) * xc( 2 ) * alpha *
     +           ( ( w**2 - 1 ) * ser - w * xc( 5 ) * H3p ) / xc( 4 )**2
          gc( 5 ) = gc( 5 ) - 
     +               2. * ( fv - pred ) * xc( 2 ) * alpha * H3 / xc( 4 )
        else
          gc( 3 ) = gc( 3 ) - 2. * ( fv - pred ) * xc( 2 ) * alpha *
     +                                          ( w * ser ) / xc( 4 )**2
          gc( 4 ) = gc( 4 ) - 2. * ( fv - pred ) * xc( 2 ) * alpha *
     +                                   ( w**2 - 1 ) * ser / xc( 4 )**2
        end if
      end do
c penalty to keep the continuum value near its true value
      gc( 1 ) = gc( 1 ) + 4.d9 * ( xc( 1 ) - cton )**3
      return
      end

      subroutine GHparm
      print *, 'GHparm called!'
      return
      end
