      real function mask( ix, iy )
c  Copyright (C) 2020, Jerry Sellwood
c
c returns a value for the given pixel that is the maximum I/N in the column
c   of velocities, except that zero is returned if no three consecutive
c   velocity channels exceed the threshold specified by Imin.
c
c This strategy is adapted from section 3.6 of Walter et al (2008)
      use aacube
      implicit none
c
c calling arguments
      integer ix, iy
c
c common block
c
      include 'comcube.h'
c
c local variables
      integer i
      logical l
      real a( 3 ), b
c
c check relative intensities over the full velocity range
      mask = 0
      a( 2 ) = dcube( ix, iy, 1 ) / noise
      a( 3 ) = dcube( ix, iy, 2 ) / noise
      b = max( a( 2 ), a( 3 ) )
      do i = 3, nsizev
        a( 1 ) = a( 2 )
        a( 2 ) = a( 3 )
        a( 3 ) = dcube( ix, iy, i ) / noise
        b = max( b, a( 3 ) )
c set mask = 1 when any 3 consecutive channels have S/N > 2
        l = ( a( 1 ) .gt. Imin ) .and.
     +      ( a( 2 ) .gt. Imin ) .and. ( a( 3 ) .gt. Imin )
        if( l )mask = 1
      end do
      mask = mask * b
      return
      end

