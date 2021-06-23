      subroutine bestmod( params )
c make models for all pixels and residuals within the mask for output
c    - the weights were preset in lfracs
c
c   Split off from bestfitp and besttfitv - JAS June 2012
c
      include 'commons.h'
c
c calling argument
      real*8 params( md )
c
c local variables
      integer i, ii, ij, ix, iy, j, k, km
      real amod, barmod, biradmod, bitanmod, bulgemod, diskmod, f0, f1
      real radmod, resd, rotmod, tmod, w
c
      if( lphot )then
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
      end if
      if( lvels )then
c initialize in case these are never set
        radmod = 0
        bitanmod = 0
        biradmod = 0
c
        if( l2D )then
          if( lseeing )then
c all weights have been set, need only determine maximum value of km
            km = 0
            do iy = 1, yrange
              do ix = 1, xrange
                if( inmask( ix, iy ) )then
                  do i = 1, mk
                    if( iwgt( i, ix, iy ) .gt. 0 )km = max( km, i )
                  end do
                end if
              end do
            end do
          else
c set all weights
            fronly = .false.
            call setwgt( params, .true., km )
          end if
c work over pixels within the mask
          do iy = 1, yrange
            do ix = 1, xrange
              if( inmask( ix, iy ) )then
c sum model velocity over components
                amod = 0
                do j = 1, nellip
                  ij = 0
                  do i = 1, km
                    if( iwgt( i, ix, iy ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )amod = amod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
                end do
                k = nellip
                tmod = amod
                rotmod = amod
c radial velocities
                if( lradial )then
                  amod = 0
                  do j = k + 1, k + nradial
                    ij = 0
                    do i = 1, km
                      if( iwgt( i, ix, iy ) .eq. j )ij = i
                    end do
                    if( ij .gt. 0 )amod = amod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
                  end do
                  k = k + nradial
                  tmod = tmod + amod
                  radmod = amod
                end if
c non-symmetric velocities
                if( lnax )then
                  amod = 0
                  do j = k + 1, k + nnasymm
                    ij = 0
                    do i = 1, km
                      if( iwgt( i, ix, iy ) .eq. j )ij = i
                    end do
                    if( ij .gt. 0 )amod = amod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
                  end do
                  k = k + nnasymm
                  tmod = tmod + amod
                  bitanmod = amod
                  amod = 0
                  do j = k + 1, k + nnasymm
                    ij = 0
                    do i = 1, km
                      if( iwgt( i, ix, iy ) .eq. j )ij = i
                    end do
                    if( ij .gt. 0 )amod = amod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
                  end do
                  k = k + nnasymm
                  tmod = tmod + amod
                  biradmod = amod
                end if
c add systemic velocity
                if( lsystemic )then
                  ij = 0
                  do i = 1, km
                    if( iwgt( i, ix, iy ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )tmod = tmod + wgt( ij, ix, iy ) * ivsys
                end if
c difference from observed value for good pixels only
                if( lgpix( ix, iy ) )then
                  if( lsystemic )then
                    resd = sdat( ix, iy ) - tmod
                  else
                    resd = sdat( ix, iy ) - vsys - tmod
                  end if
                  res( ix, iy ) = resd
                  x2res( ix, iy ) = resd / sigma( ix, iy )
                end if
                model( ix, iy ) = tmod
                rotvels( ix, iy ) = rotmod
                radvels( ix, iy ) = radmod
                bivtan( ix, iy ) = bitanmod
                bivrad( ix, iy ) = biradmod
              end if
            end do
          end do
        else
c determine maximum value of km
          km = 0
          do i = 1, inp_pts
            if( inmsk1( i ) )then
              do j = 1, mk
                if( iwgt1( j, i ) .gt. 0 )km = max( km, j )
              end do
            end if
          end do
          do ii = 1, inp_pts
            if( inmsk1( ii ) )then
c sum model velocity over components
              amod = 0
              do j = 1, nellip
                ij = 0
                do i = 1, km
                  if( iwgt1( i, ii ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )amod = amod + wgt1( ij, ii ) *fitval( j )
              end do
              k = nellip
              tmod = amod
              rotmod = amod
c radial velocities
              if( lradial )then
                amod = 0
                do j = k + 1, k + nradial
                  ij = 0
                  do i = 1, km
                    if( iwgt1( i, ii ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )amod = amod +
     +                                      wgt1( ij, ii ) * fitval( j )
                end do
                k = k + nradial
                tmod = tmod + amod
                radmod = amod
              end if
c non-symmetric velocities
              if( lnax )then
                amod = 0
                do j = k + 1, k + nnasymm
                  ij = 0
                  do i = 1, km
                    if( iwgt1( i, ii ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )amod = amod +
     +                                      wgt1( ij, ii ) * fitval( j )
                end do
                k = k + nnasymm
                tmod = tmod + amod
                bitanmod = amod
                amod = 0
                do j = k + 1, k + nnasymm
                  ij = 0
                  do i = 1, km
                    if( iwgt1( i, ii ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )amod = amod +
     +                                      wgt1( ij, ii ) * fitval( j )
                end do
                k = k + nnasymm
                tmod = tmod + amod
                biradmod = amod
              end if
c add systemic velocity
              if( lsystemic )then
                ij = 0
                do i = 1, km
                  if( iwgt1( i, ii ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )tmod = tmod + wgt1( ij, ii ) * ivsys
              end if
c difference from observed value
              if( lsystemic )then
                resd = sdat1( ii ) - tmod
              else
                resd = sdat1( ii ) - vsys - tmod
              end if
              res1( ii ) = resd
              x2res1( ii ) = resd / sigma1( ii )
              model1( ii ) = tmod
              rotvels1( ii ) = rotmod
              radvels1( ii ) = radmod
              bivtan1( ii ) = bitanmod
              bivrad1( ii ) = biradmod
            end if
          end do
        end if
      end if
c clear the array
      do j = 1, nellip
        ringpts( j ) = 0
      end do
      if( lphot .and. ldisk )then
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
      if( lvels )then
        if( l2D )then
          do iy = 1, yrange
            do ix = 1, xrange
c need only good pixels
              if( lgpix( ix, iy ) )then
                do j = 1, km
                  k = iwgt( j, ix, iy )
                  if( ( k .gt. 0 ) .and. ( k .le. nellip ) )then
                    w = max( 0., wgt( j, ix, iy ) )
                    w = min( 1., w )
                    ringpts( k ) = ringpts( k ) + w
                  end if
                end do
              end if
            end do
          end do
        else
          do ii = 1, inp_pts
c need only good pixels
            if( lgpix1( ii ) )then
              do j = 1, km
                k = iwgt1( j, ii )
                if( ( k .gt. 0 ) .and. ( k .le. nellip ) )then
                  w = max( 0., wgt1( j, ii ) )
                  w = min( 1., w )
                  ringpts( k ) = ringpts( k ) + w
                end if
              end do
            end if
          end do
        end if
      end if
      return
      end
