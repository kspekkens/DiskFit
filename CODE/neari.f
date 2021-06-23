      integer function neari( a, n, x )
c utility routine to find the nearest value in a 1D array to the given value
c   assumes the values in the input array are sequenced and works by bisection
      implicit none
c
c calling arguments
      integer n
      real a( n ), x
c
c local variables
      integer i, j, k
      real d
c
      neari = 1
      if( n .gt. 1 )then
        k = 1
        j = n
c find bracketing values
        do while ( k + 1 .lt. j )
          i = ( k + j ) / 2
          if( a( i ) .lt. x )then
            k = i
          else
            j = i
          end if
        end do
c find absolute nearest
        j = max( 1, k - 1 )
        k = min( k + 1, n )
        d = abs( a( j ) - x )
        neari = j
        do i = j + 1, k
          if( abs( a( i ) - x ) .lt. d )then
            neari = i
            d = abs( a( i ) - x )
          end if
        end do
      end if
      return
      end
