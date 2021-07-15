      integer function nearpt( indlist, m, inpind )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c returns index of pt in indlist whose corresponding xval,yval
c    are closest to that in inpind, and whose indlist value is nonzero
c
c created by KS
c   included in diskfit by JAS Nov 11
c   Updated to f90 - JAS Jan 2015
c
      include 'commons.h'
c
c calling arguments
      integer indlist( m ), inpind, m
c
c local variables
      integer ii, minind
      real min, diff, xref, yref
c
      min = 1000
      minind = 0
      xref = xval( inpind )
      yref = yval( inpind )
c
      do ii = 1, inp_pts
        if( lgpix( ii, 1 ) )then
          diff = sqrt( ( xval( ii ) - xref )**2 +
     +                 ( yval( ii ) - yref )**2 )
          if( ( diff .lt. min ) .and. ( indlist( ii ) .ne. 0 ) )then
            min = diff
            minind = ii
          end if
        end if
      end do
      nearpt = minind
      return
      end
