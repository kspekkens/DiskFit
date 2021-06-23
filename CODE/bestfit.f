      subroutine bestfit( params, eparams )
c copies best fitting parameters and values to more useful areas
c
c   Created by KS
c   Re-created by lifting code from bestfitp and bestfitv - JAS June 2012
c
      include 'commons.h'
c
c calling arguments
      real*8 params( md ), eparams( md )
c
c local variables
      integer i, j, k
      real arg, bn, dIo, dre, dn, expon, f0, f1, rad, sumsq
      real*8 diskel, diskpa, rw, welm, wpm
c
      if( lphot )then
c separate fitval array into disk, bar and bulge
        k = 0
        if( ldisk )then
          do j = 1, nellip
            id( j ) = fitval( k + j )
            eid( j ) = efitval( k + j )
          end do
          k = k + nellip
        end if
        if( lnax )then
          do j = 1, nbar
            if( ( j .lt. nminr ) .or. ( j .gt. nmaxr ) )then
              ibar( j ) = 0.0
              eibar( j ) = 0.0
            else
              i = j - nminr + 1
              ibar( j ) = fitval( k + i )
              eibar( j ) = efitval( k + i )
            end if
          end do
          k = k + nbar
        end if
        if( lbulge )then
          do j = 1, nbulge
            ibulge( j ) = fitval( k + j )
            eibulge( j ) = efitval( k + j )
          end do
          k = k + nbulge
        end if
      end if
      if( lvels )then
c separate fitval array into rotation, radial, non-axisymm and systemic:
        k = 0
        do j = 1, nellip
          id( j ) = fitval( k + j )
          eid( j ) = efitval( k + j )
        end do
        k = k + nellip
        if( lradial )then
          do j = 1, nellip
            if( ( j .lt. nminr ) .or. ( j .gt. nmaxr ) )then
              irad( j ) = 0.
              eirad( j ) = 0.
            else
              i = j - nminr + 1
              irad( j ) = fitval( k + i )
              eirad( j ) = efitval( k + i )
            end if
          end do
          k = k + nradial
        end if
        if( lnax )then
          do j = 1, nellip
            if( ( j .lt. nminr ) .or. ( j .gt. nmaxr ) )then
              ibirad( j ) = 0
              eibirad( j ) = 0
              ibitan( j ) = 0
              eibitan( j ) = 0
            else
              i = j - nminr + 1
              ibitan( j ) = fitval( k + i )
              eibitan( j ) = efitval( k + i )
              ibirad( j ) = fitval( k + nnasymm + i )
              eibirad( j ) = efitval( k + nnasymm + i )
            end if
          end do
          k = k + 2 * nnasymm
        end if
        if( lsystemic )then
          ivsys = fitval( ntot )
          eivsys = efitval( ntot )
          k = k + 1
        end if
      end if
      if( k .ne. ntot )then
        print *, 'Mismatch between k and ntot'
        print *, 'Contact code administrator'
        call crash( 'bestfit' )
      end if
c set other parameters to their final values:
      if( lphot )then
        k = 0
        if( lpa )then
          newpa = params( k + 1 )
          enewpa = eparams( k + 1 )
          k = k + 1
        end if
        if( leps )then
          neweps  = params( k + 1 )
          eneweps = eparams( k + 1 )
          k = k + 1
        end if
        if( lcentre )then
          newxcen = params( k + 1 )
          enewxcen = eparams( k + 1 )
          newycen = params( k + 2 )
          enewycen = eparams( k + 2 )
          k = k + 2
        end if
        if( lnax )then
          if( lphib )then
            newbar_pa = params( k + 1 )
            enewbar_pa = eparams( k + 1 )
            k = k + 1
          end if
          if( lepsb )then
            newbar_eps = params( k + 1 )
            enewbar_eps = eparams( k + 1 )
            k = k + 1
          end if
        end if
        if( lbulge )then
          if( lbleps )then
            newbulge_l = params( k + 1 )
            enewbulge_l = eparams( k + 1 )
            k = k + 1
          else
            newbulge_l = bulge_el
            enewbulge_l = 0
          end if
          if( lsersn )then
            newbulge_n = params( k + 1 )
            enewbulge_n = eparams( k + 1 )
            k = k + 1
          else
            newbulge_n = bulge_n
            enewbulge_n = 0
          end if
          if( lr_e )then
            newr_bulge = params( k + 1 )
            enewr_bulge = eparams( k + 1 )
            k = k + 1
          else
            newr_bulge = r_bulge
            enewr_bulge = 0
          end if
c if a bulge is fitted, compute the bulge profile out to nellip
c also compute uncertainties on bulge profile
c for 1 < n < 10 Graham 01 gives k = 1.9992n - 0.3271
          bn = 1.9992 * newbulge_n - 0.3271
          expon = 1. / newbulge_n
          do i = 1, nellip
            rad = sma( i ) / newr_bulge
            arg = -bn * ( rad**expon - 1. )
            arg = min( arg, 20. )
            blgprof( i ) = ibulge( 1 ) * exp( arg )
            dIo = blgprof( i ) / ibulge( 1 )
            dre = blgprof( i ) * bn /
     +           ( newbulge_n * newr_bulge ) * rad**expon
            if( rad .gt. 0.01 )then
              dn  = blgprof( i ) *
     +             ( ( bn * rad**expon * alog( rad ) ) / newbulge_n**2
     +             + 1.9992 * arg / bn )
            else
              dn = blgprof( i ) * 1.9992 * arg / bn
            end if
            sumsq = dIo**2 * eibulge(1)**2 + dre**2 * enewr_bulge**2
     +           + dn**2 * enewbulge_n**2
            eblgprof( i ) = sqrt( sumsq )
          end do
        end if
      end if
      if( lvels )then
c set other parameters to their final values:
        k = 0
        if( lpa )then
          diskpa = params( k + 1 )
          newpa = diskpa
          enewpa = eparams( k + 1 )
          k = k + 1
        else
          diskpa = pa
        end if
        if( leps )then
          diskel = params( k + 1 )
          neweps = diskel
          eneweps = eparams( k + 1 )
          k = k + 1
        else
          diskel = eps
        end if
        if( lcentre )then
          newxcen = params( k + 1 )
          enewxcen = eparams( k + 1 )
          newycen = params( k + 2 )
          enewycen = eparams( k + 2 )
          k = k + 2
        end if
        if( lnax )then
          newphib = params( k + 1 )
          enewphib = eparams( k + 1 )
          k = k + 1
        end if
        if( lwarp )then
          if( lrwarp )then
            rw = params( k + 1 )
            enewrw = eparams( k + 1 )
            k = k + 1
          else
            rw = rwarp
            enewrw = 0
          end if
          newrw = rw
          if( lwepsm )then
            welm = params( k + 1 )
            enewem = eparams( k + 1 )
            k = k + 1
          else
            welm = wepsm
            enewem = 0
          end if
          newem = welm
          if( lwpm )then
            wpm = params( k + 1 )
            enewpm = eparams( k + 1 )
            k = k + 1
          else
            wpm = wphim
            enewpm = 0
          end if
          newpm = wpm
        end if
      end if
      if( k .ne. nd )then
        print *, 'Mismatch between k and nd'
        print *, 'Contact code administrator'
        call crash( 'bestfit' )
      end if
      return
      end
