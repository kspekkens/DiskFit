      subroutine nsevls( work, n )
c Copyright (C) 2020, Jerry Sellwood
c
c Routine collect many likely noise values from a data cube
      use aacube
      implicit none
c
c calling argument
      integer n
      real work( n )
c
c common blocks
c
      include 'comcube.h'
c
c local variables
      integer i, j, k, l
c
      l = 0
c 2 lowest velocity planes
      do k = 1, 2
        do j = 1, nsizey
          do i = 1, nsizex
            l = l + 1
            work( l ) = dcube( i, j, k )
          end do
        end do
      end do
c 2 highest velocity planes
      do k = nsizev - 1, nsizev
        do j = 1, nsizey
          do i = 1, nsizex
            l = l + 1
            work( l ) = dcube( i, j, k )
          end do
        end do
      end do
      do k = 3, nsizev - 2
c left y-v edge
        do j = 1, 2
          do i = 1, nsizex
            l = l + 1
            work( l ) = dcube( i, j, k )
          end do
        end do
c right y-v edge
        do j = nsizey - 1, nsizey
          do i = 1, nsizex
            l = l + 1
            work( l ) = dcube( i, j, k )
          end do
        end do
        do j = 3, nsizey - 2
c front x-v edge
          do i = 1, 2
            l = l + 1
            work( l ) = dcube( i, j, k )
          end do
c back x-v edge
          do i = nsizex - 1, nsizex
            l = l + 1
            work( l ) = dcube( i, j, k )
          end do
        end do
      end do
      if( l .ne. n )then
        print *, l, n
        call crash( 'Logic error in NSEVLS' )
      end if
      return
      end
