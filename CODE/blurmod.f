      subroutine blurmod( km )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c This subroutine blurs the pre-calculated weight array to make the
c   appropriate seeing correction as described in the documentation
c
c   Created by JAS March 2011
c   Updated to f90 - JAS Jan 2015
c   Eliminated bug in setting iwgt array - JAS Jun 2015
c
c common blocks
      include 'commons.h'
c
c calling argument
      integer km
c
c local arrays
      integer, allocatable :: jwgt( :, :, : )
      real, allocatable :: bwgt( :, :, : )
c
c local variables:
      integer i, ix, iy, j, jw, jx, jy, k, kl, l, m, mxx, mxy, mnx, mny
      logical newe
      real wtn
c
      if( .not. l2D )then
        print *, 'Seeing corrections not programmed for pixel list data'
        call crash( 'BLURMOD' )
      end if
c copy old values to working arrays - needed for pixels that are unaffected
      allocate ( jwgt( mk, xrange, yrange ) )
      allocate ( bwgt( mk, xrange, yrange ) )
      do j = 1, yrange
        do i = 1, xrange
          do k = 1, mk
            jwgt( k, i, j ) = iwgt( k, i, j )
            bwgt( k, i, j ) = wgt( k, i, j )
          end do
        end do
      end do
      kl = km
c ignore blurring near the edge
      mnx = iblur / 2 + 1
      mxx = xrange - iblur / 2
      mny = jblur / 2 + 1
      mxy = yrange - jblur / 2
c work over blurred pixels
      do iy = mny, mxy
        do ix = mnx, mxx
c need only good pixels
          if( lgpix( ix, iy ) )then
c initialize blurred weights and note current # of ellipses with non-zero wts
            m = 0
            do l = 1, mk
              if( iwgt( l, ix, iy ) .gt. 0 )m = l
              bwgt( l, ix, iy ) = 0
            end do
c sum contributions from neighbours
            do i = 1, nblur
              jx = ix + ibl( i )
              jy = iy + jbl( i )
c work over all non-zero unblurred weights of this neighbour
              do j = 1, km
                jw = iwgt( j, jx, jy )
                if( jw .gt. 0 )then
c compute fragment needed for current pixel
                  wtn = wblur( i ) * wgt( j, jx, jy )
                  newe = .true.
c add this fragment to the appropriate ellipse
                  do l = 1, m
c add to a pre-existing ellipse, if possible
                    if( jwgt( l, ix, iy ) .eq. jw )then
                      newe = .false.
                      bwgt( l, ix, iy ) = bwgt( l, ix, iy ) + wtn
                    end if
                  end do
c new ellipse for this pixel
                  if( newe )then
                    m = m + 1
                    jwgt( m, ix, iy ) = jw
                    bwgt( m, ix, iy ) = wtn
                    kl = max( kl, m )
                  end if
                end if
              end do
            end do
          end if
        end do
      end do
c copy back from working arrays
      do j = 1, yrange
        do i = 1, xrange
          do l = 1, kl
            iwgt( l, i, j ) = jwgt( l, i, j )
            wgt( l, i, j ) = bwgt( l, i, j )
          end do
        end do
      end do
      deallocate ( bwgt )
      deallocate ( jwgt )
      km = kl
      return
      end
