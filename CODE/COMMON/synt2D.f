      subroutine synt2D( temp, nx, ny )
c  Copyright (C) 2020, Jerry Sellwood
      implicit none
c routine to perform the double Fourier synthesis of the transformed data
c   on 2-D raster.  The results overwrite the input array and are one
c   quarter the size
c
c Uses single precision FFTPAK routines
c
c calling argument
      integer nx, ny
      real temp( 4 * nx * ny )
c
c local allocatable arrays
      real, allocatable :: trig( : )
      real, allocatable :: w( : )
c
c local variables
      integer ix, iy, j, m, mtrig
      real scale
c
c allocate space
      j = 2 * nx * ny
      allocate ( w( j ) )
      j = 4 * max( nx, ny ) + 15
      allocate ( trig( j ) )
c initialize
      mtrig = 2 * nx
      call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
c set pointer
      j = 1
c Fourier synthesis in x first
      do iy = 1, 2 * ny
        call sfftb1( mtrig, temp( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c re-order data for synthesis in y - ignore 2nd half of x-values
      j = 0
      do ix = 1, nx
        m = ix
        do iy = 1, 2 * ny
          j = j + 1
          w( j ) = temp( m )
          m = m + 2 * nx
        end do
      end do
c initialize if the vector length is new
      if( mtrig .ne. 2 * ny )then
        mtrig = 2 * ny
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
      end if
c set pointer
      j = 1
c Fourier synthesis in y
      do ix = 1, nx
        call sfftb1( mtrig, w( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c copy results to final area
      m = 0
      scale = 1. / real( 4 * nx * ny )
      do iy = 1, ny
        j = iy
        do ix = 1, nx
          m = m + 1
          temp( m ) = w( j ) * scale
          j = j + 2 * ny
        end do
        m = m + nx
      end do
      return
      end
