      subroutine bestfitp( params, eparams, n, chi2, iterl )
c copies data from fitting arrays to more useful areas, and creates model
c   and residuals arrays for output
c
c   Created by KS
c   Polished by JAS October 2009
c   Adapted for photometry option by JAS March 2011
c   Output bulge intensity profile and uncertainties- KS March 2012
c   Simplified use of mask - JAS March 2012
c
      include 'commons.h'
c
c calling arguments
      integer iterl, n
      real*8 chi2, params( n ), eparams( n )
c
c local variables
      integer i, ij, ix, iy, j, k, km
      real arg, barmod, bn, bulgemod, diskmod, expon, f0, f1, resd, tmod
      real rad, dIo, dre, dn, sumsq
c
c separate fitval array into disk, bar and bulge
      k = 0
      if( ldisk )then
        do j = 1, nellip
          id( j ) = fitval( k + j )
          eid( j ) = efitval( k + j )
        end do

c        print *, 'disk intensity array'
c        print '( 8f10.2 )', ( id( j ), j = 1, nellip )

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

c        print *, 'bar intensity array'
c        print '( 8f10.2 )', ( fitval( j ), j = k + 1, k + nbar )

        k = k + nbar
      end if
      if( lbulge )then
        do j = 1, nbulge
          ibulge( j ) = fitval( k + j )
          eibulge( j ) = efitval( k + j )
        end do

c        print *, 'bulge intensity'
c        print '( 8f10.2 )', ( fitval( j ), j = k + 1, k + nbulge )

        k = k + nbulge
      end if
      if( k .ne. ntot )then
        print*, 'Mismatch between k and ntot in bestfit'
        print*, 'Contact code administrator'
        call crash( 'bestfitp' )
      end if
c set other parameters to their final values:
      final_chi2 = chi2
      final_iter = iterl
      j = 0
      if( lpa )then
        newpa = params( j + 1 )
        enewpa = eparams( j + 1 )
        j = j + 1
      end if
      if( leps )then
        neweps  = params( j + 1 )
        eneweps = eparams( j + 1 )
        j = j + 1
      end if
      if( lcentre )then
        newxcen = params( j + 1 )
        enewxcen = eparams( j + 1 )
        newycen = params( j + 2 )
        enewycen = eparams( j + 2 )
        j = j + 2
      end if
      if( lnax )then
        if( lphib )then 
          newbar_pa = params( j + 1 )
          enewbar_pa = eparams( j + 1 )
          j = j + 1
        end if
        if( lepsb )then
          newbar_eps = params( j + 1 )
          enewbar_eps = eparams( j + 1 )
          j = j + 1
        end if
      end if
      if( lbulge )then
        if( lbleps )then
          newbulge_l = params( j + 1 )
          enewbulge_l = eparams( j + 1 )
          j = j + 1
        else
          newbulge_l = bulge_el
          enewbulge_l = 0
        end if
        if( lsersn )then
          newbulge_n = params( j + 1 )
          enewbulge_n = eparams( j + 1 )
          j = j + 1
        else
          newbulge_n = bulge_n
          enewbulge_n = 0
        end if
        if( lr_e )then
          newr_bulge = params( j + 1 )
          enewr_bulge = eparams( j + 1 )
          j = j + 1
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
          arg = -bn * ( ( rad )**expon - 1. )
          arg = min( arg, 20. ) 
          blgprof( i ) = ibulge( 1 ) * exp( arg )
          dIo = blgprof( i ) / ibulge( 1 )
          dre = blgprof( i ) * bn / 
     +      ( newbulge_n * newr_bulge ) * rad**expon
          if( rad .gt. 0.01 )then
            dn  = blgprof( i ) * ( 
     +       ( bn * rad**expon * alog( rad ) )/newbulge_n**2 
     +       + 1.9992 * arg / bn )
          else
            dn = blgprof( i ) * 1.9992 * arg / bn
          end if
          sumsq = dIo**2 * eibulge(1)**2 + dre**2 * enewr_bulge**2 
     +            + dn**2 * enewbulge_n**2 
          eblgprof( i ) = sqrt( sumsq )
        end do
      end if
c
c make model for all pixels within the mask - the weights were preset in lfracs
c
c redetermine the maximum number of ellipses contributing to any pixel
      km = 0
      do iy = 1, yrange
        do ix = 1, xrange
          do j = 1, ntot
            do i = 1, mk
              if( iwgt( i, ix, iy ) .eq. j )km = max( km, i )
            end do
          end do
        end do
      end do
c initialize in case any are never set
      diskmod = 0
      barmod = 0
      bulgemod = 0
      tmod = 0
c work over pixels
      do iy = 1, yrange
        do ix = 1, xrange
c select only pixels within the mask
          if( inmask( ix, iy ) )then
            diskmod = 0
            k = 0
c sum model over components
            if( ldisk )then
              do j = k + 1, k + nellip
                ij = 0
                do i = 1, km
                  if( iwgt( i, ix, iy ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )diskmod = diskmod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
              end do
              k = k + nellip
            end if
            tmod = diskmod
c bar
            if( lnax )then
              barmod = 0
              do j = k + 1, k + nbar
                ij = 0
                do i = 1, km
                  if( iwgt( i, ix, iy ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )barmod = barmod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
              end do
              k = k + nbar
              tmod = tmod + barmod
            end if
c bulge
            if( lbulge )then
              bulgemod = 0
              do j = k + 1, k + nbulge
                ij = 0
                do i = 1, km
                  if( iwgt( i, ix, iy ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )bulgemod = bulgemod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
              end do
              k = k + nbulge
              tmod = tmod + bulgemod
            end if
c save model values
            model( ix, iy ) = tmod
            diskint( ix, iy ) = diskmod
            barint( ix, iy ) = barmod
            bulgeint( ix, iy ) = bulgemod
c difference from observed value only for good pixels
            if( lgpix( ix, iy ) )then
              resd = sdat( ix, iy ) - tmod
              res( ix, iy ) = resd
              x2res( ix, iy ) = resd / sigma( ix, iy )
            end if
          end if
        end do
      end do
c clear the array
      do j = 1, nellip
        ringpts( j ) = 0
      end do
      if( ldisk )then
c work over all pixels
        do iy = 1, yrange
          do ix = 1, xrange
            if( lgpix( ix, iy ) )then
c calculate the weights
              f0 = wgt( 2, ix, iy )
              f0 = max( 0., f0 )
              f0 = min( 1., f0 )
              f1 = 1. - f0
c ringpts sums the pixel weights that contribute to each ring
              k = iwgt( 1, ix, iy )
              ringpts( k ) = ringpts( k ) + f1
              ringpts( k + 1 ) = ringpts( k + 1 ) + f0
            end if
          end do
        end do
      end if
      return
      end
