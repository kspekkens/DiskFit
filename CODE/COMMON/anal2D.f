      subroutine anal2D( temp, nx, ny )
c  Copyright (C) 2020, Jerry Sellwood
c
c routine to perform the double Fourier analysis of the data array on a
c   2-D raster.
c
c Uses single precision FFTPAK routines
      implicit none
c
c calling arguments
      integer nx, ny
      real temp( 4 * nx * ny )
c
c local allocatable arrays
      real, allocatable :: trig(:)
      real, allocatable :: w(:)
c
c local variables
      integer i, ix, iy, j, m, mshb, mtrig
c
c allocate space
      mshb = 4 * nx * ny
      allocate ( w( mshb ) )
c arrange input data and add padding for FFT in y
      j = 0
      do ix = 1, nx
        m = ix
        do iy = 1, ny
          j = j + 1
          w( j ) = temp( m )
          w( j + ny ) = 0
          m = m + 2 * nx
        end do
        j = j + ny
      end do
c initialize if the vector length is new
      i = 4 * max( nx, ny ) + 15
      allocate ( trig( i ) )
      mtrig = 2 * ny
      call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
c Fourier analysis in y
      j = 1
      do i = 1, nx
        call sfftf1( mtrig, w( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c re-arrange data and add more padding for FFT in x
      j = 0
      do iy = 1, mtrig
        m = iy
        do ix = 1, nx
          j = j + 1
          temp( j ) = w( m )
          temp( j + nx ) = 0
          m = m + 2 * nx
        end do
        j = j + nx
      end do
c initialize if the vector length is new
      if( 2 * nx .ne. mtrig )then
        mtrig = 2 * nx
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
      end if
c Fourier analysis in x
      j = 1
      do i = 1, 2 * ny
        call sfftf1( mtrig, temp( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
      return
      end
