      subroutine bestfitv( params, eparams, n, chi2, iterl )
c copies data from fitting arrays to more useful areas, and creates model
c   and residuals arrays for output
c
c   Created by KS
c   Polished by JAS October 2009
c   Revised for warp models by JAS July 2011
c   Compute model for all pixels within the mask - JAS March 2012
c
      include 'commons.h'
c
c calling arguments
      integer iterl, n
      real*8 chi2, params( n ), eparams( n )
c
c local variables
      integer i, ii, ij, ix, iy, j, k, km
      real amod, biradmod, bitanmod, radmod, resd, rotmod, tmod, w
      real*8 diskel, diskpa, rw, welm, wpm
c
c separate fitval array into rotation, radial, non-axisymm and systemic:
      k = 0
      do j = 1, nellip
        id( j ) = fitval( k + j )
        eid( j ) = efitval( k + j )
      end do

c      print *, 'circular speed array'
c      print '( 8f10.2 )', ( id( j ), j = 1, nellip )

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
      if( k .ne. ntot )then
        print*, 'Mismatch between k and ntot in bestfit'
        print*, 'Contact code administrator'
        call crash( 'bestfit' )
      end if
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
              if( ij .gt. 0 )amod = amod + wgt1( ij, ii ) * fitval( j )
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
c set other parameters to their final values:
      final_chi2 = chi2
      final_iter = iterl
      j = 0
      if( lpa )then
        diskpa = params( j + 1 )
        newpa = diskpa
        enewpa = eparams( j + 1 )
        j = j + 1
      else
        diskpa = pa
      end if
      if( leps )then
        diskel = params( j + 1 )
        neweps = diskel
        eneweps = eparams( j + 1 )
        j = j + 1
      else
        diskel = eps
      end if
      if( lcentre )then
        newxcen = params( j + 1 )
        enewxcen = eparams( j + 1 )
        newycen = params( j + 2 )
        enewycen = eparams( j + 2 )
        j = j + 2
      end if
      if( lnax )then
        newphib = params( j + 1 )
        enewphib = eparams( j + 1 )
        j = j + 1
      end if
      if( lwarp )then
        if( lrwarp )then
          rw = params( j + 1 )
          enewrw = eparams( j + 1 )
          j = j + 1
        else
          rw = rwarp
          enewrw = 0
        end if
        newrw = rw
        if( lwepsm )then
          welm = params( j + 1 )
          enewem = eparams( j + 1 )
          j = j + 1
        else
          welm = wepsm
          enewem = 0
        end if
        newem = welm
        if( lwpm )then
          wpm = params( j + 1 )
          enewpm = eparams( j + 1 )
          j = j + 1
        else
          wpm = wphim
          enewpm = 0
        end if
        newpm = wpm
      end if
c sum contributions from each pixel to the orbital velocities
c   formerly done by subroutine weightrings
c recompute weights for each pixel in the map
      fronly = .true.
      call setwgt( params, .false., km )
c clear the array
      do j = 1, nellip
        ringpts( j ) = 0
      end do
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
      return
      end
