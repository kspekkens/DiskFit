      subroutine bestmod( params )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c make models for all pixels and residuals within the mask for output
c    - the weights were preset in lfracs
c
c   Split off from bestfitp and besttfitv - JAS June 2012
c   Fixed bug in calculating ringpts for vels - JAS Dec 2014
c   Updated to f90 - JAS Jan 2015
c   Minor change to use of vsys - JAS Jul 2015
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
        allocate ( barint( xrange, yrange ) )
        allocate ( bulgeint( xrange, yrange ) )
        allocate ( diskint( xrange, yrange ) )
        allocate ( model( xrange, yrange ) )
        allocate ( res( xrange, yrange ) )
        allocate ( x2res( xrange, yrange ) )
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
          allocate ( bivrad( xrange, yrange ) )
          allocate ( bivtan( xrange, yrange ) )
          allocate ( model( xrange, yrange ) )
          allocate ( res( xrange, yrange ) )
          allocate ( rotvels( xrange, yrange ) )
          allocate ( radvels( xrange, yrange ) )
          allocate ( x2res( xrange, yrange ) )
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
                else
                  tmod = tmod + vsys
                end if
c difference from observed value for good pixels only
                if( lgpix( ix, iy ) )then
                  resd = sdat( ix, iy ) - tmod
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
          allocate ( bivrad( inp_pts, 1 ) )
          allocate ( bivtan( inp_pts, 1 ) )
          allocate ( model( inp_pts, 1 ) )
          allocate ( res( inp_pts, 1 ) )
          allocate ( rotvels( inp_pts, 1 ) )
          allocate ( radvels( inp_pts, 1 ) )
          allocate ( x2res( inp_pts, 1 ) )
c determine maximum value of km
          km = 0
          do i = 1, inp_pts
            if( inmask( i, 1 ) )then
              do j = 1, mk
                if( iwgt( j, i, 1 ) .gt. 0 )km = max( km, j )
              end do
            end if
          end do
          do ii = 1, inp_pts
            if( inmask( ii, 1 ) )then
c sum model velocity over components
              amod = 0
              do j = 1, nellip
                ij = 0
                do i = 1, km
                  if( iwgt( i, ii, 1 ) .eq. j )ij = i
                end do
                if(
     +            ij .gt. 0 )amod = amod + wgt( ij, ii, 1 ) *fitval( j )
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
                    if( iwgt( i, ii, 1 ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )amod = amod +
     +                                    wgt( ij, ii, 1 ) * fitval( j )
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
                    if( iwgt( i, ii, 1 ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )amod = amod +
     +                                    wgt( ij, ii, 1 ) * fitval( j )
                end do
                k = k + nnasymm
                tmod = tmod + amod
                bitanmod = amod
                amod = 0
                do j = k + 1, k + nnasymm
                  ij = 0
                  do i = 1, km
                    if( iwgt( i, ii, 1 ) .eq. j )ij = i
                  end do
                  if( ij .gt. 0 )amod = amod +
     +                                    wgt( ij, ii, 1 ) * fitval( j )
                end do
                k = k + nnasymm
                tmod = tmod + amod
                biradmod = amod
              end if
c add systemic velocity
              if( lsystemic )then
                ij = 0
                do i = 1, km
                  if( iwgt( i, ii, 1 ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )tmod = tmod + wgt( ij, ii, 1 ) * ivsys
              else
                tmod = tmod + vsys
              end if
c difference from observed value
              resd = sdat( ii, 1 ) - tmod
              res( ii, 1 ) = resd
              x2res( ii, 1 ) = resd / sigma( ii, 1 )
              model( ii, 1 ) = tmod
              rotvels( ii, 1 ) = rotmod
              radvels( ii, 1 ) = radmod
              bivtan( ii, 1 ) = bitanmod
              bivrad( ii, 1 ) = biradmod
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
c reset all weights
        fronly = .true.
        call setwgt( params, .true., km )
        if( l2D )then
          do iy = 1, yrange
            do ix = 1, xrange
c need only good pixels
              if( lgpix( ix, iy ) )then
c need only interpolation weights
                do j = 1, 2
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
            if( lgpix( ii, 1 ) )then
              do j = 1, km
                k = iwgt( j, ii, 1 )
                if( ( k .gt. 0 ) .and. ( k .le. nellip ) )then
                  w = max( 0., wgt( j, ii, 1 ) )
                  w = min( 1., w )
                  ringpts( k ) = ringpts( k ) + w
                end if
              end do
            end if
          end do
        end if
c revert to the usual weights
        fronly = .false.
        call setwgt( params, .true., km )
      end if
      return
      end
