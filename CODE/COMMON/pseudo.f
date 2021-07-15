      subroutine pseudo( jdum )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c This subroutine uses a bootstrap technique to generate a realization
c    of the velocity field that will be minimized.
c
c Created by KS
c   Included in DiskFit by JS Nov 11
c   Updated to f90 - JAS Jan 2015
c   Made swap always in range 1 - inp_pts - JAS Aug 2015
c
      include 'commons.h'
c
      integer jdum
c
c externals
      integer nearpt
      real*8 ran1_dbl
c
c local arrays
      integer, allocatable :: indlist( : )
      real, allocatable :: newx2res( : )
c
c local variables
      integer ii, newind, swap
      logical ll
      real find, val
c
      if( l2D )then
        print *, 'Should be using pseudp for 2D maps'
        call crash( 'PSEUDO' )
      end if
      allocate( indlist( inp_pts ) )
      allocate( newx2res( inp_pts ) )
c swap relative residuals at random
      do ii = 1, inp_pts
        if( lgpix( ii, 1 ) )then
          ll = .true.
          do while ( ll )
            swap = int( ran1_dbl( jdum ) * inp_pts ) + 1
            ll = .not. lgpix( swap, 1 )
          end do
          newx2res( ii ) = x2res( swap, 1 )
        end if
      end do
c junc>1 requires quasi-coherent residuals
      if( junc .gt. 1. )then
        find = 1. / junc
c set index values to 0 for dependent points
        do ii = 1, inp_pts
          indlist( ii ) = 1
          if( lgpix( ii, 1 ) )then
            val = ran1_dbl( jdum )
            if( val .gt. find )indlist( ii ) = 0
          end if
        end do
c assign dependent points to value of nearest independent point
        do ii = 1, inp_pts
          if( lgpix( ii, 1 ) )then
            if( indlist( ii ) .eq. 0 )then
              newind = nearpt( indlist, inp_pts, ii )
              newx2res( ii ) = newx2res( newind )
            end if
          end if
        end do
      end if
c generate realization of velocity field + errors
      do ii = 1, inp_pts
        if( lgpix( ii, 1 ) )then
          sdat( ii, 1 ) =
     +                  model( ii, 1 ) + newx2res( ii ) * sigma( ii, 1 )
        end if
      end do
      return
      end
