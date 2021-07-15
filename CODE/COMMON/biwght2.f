      subroutine biwght2( val, n, mean, sigma )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c estimates the mean and dispersion of a vector of real*8 values using the
c   biweight algorithm that ignores outliers (see Beers, Flynn, & Gebhardt
c   1990, AJ, v100, p32) 
c
c calling arguments
      integer n
      real*8 val( n ), mean, sigma
c
c local allocatable array
      real*8, allocatable :: x(:)
c
c local variables
      integer i, j
      real*8 den, med, mad, num, u
c
c copy values
      allocate ( x( n ) )
      do i = 1, n
        x( i ) = val( i )
      end do
c find median value
      call srtmrg2( x, n )
      j = max( ( n + 1 ) / 2, n / 2 )
      med = x( j )
c find median abs deviation (mad)
      do i = 1, n
        x( i ) = abs( val( i ) - med )
      end do
      call srtmrg2( x, n )
      mad = x( j )
      deallocate ( x )
c compute biweight estimates of mean and dispersions
      num = 0
      den = 0
      do i = 1, n
c use biweight if possible
        if( mad .gt. 0.d0 )then
          u = ( val( i ) - med ) / ( 9.d0 * mad )
c ignore outliers
          if( abs( u ) .lt. 1.d0 )then
            num = num + ( val( i ) - med )**2 * ( 1.d0 - u**2 )**4
            den = den + abs( ( 1.d0 - u**2 ) * ( 1.d0 - 5.d0 * u**2 ) )
          end if
        else
c compute normal mean and st dev
          num = num + val( i )**2
          den = den + val( i )
        end if
      end do
      u = n
      if( mad .gt. 0.d0 )then
        sigma = sqrt( u * num ) / den
        mean = med
      else
        mean = den / u
        sigma = ( num - u * mean**2 ) / ( u - 1.d0 )
        sigma = sqrt( max( sigma, 0.d0 ) )
      end if
      return
      end
