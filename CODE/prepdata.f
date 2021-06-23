      subroutine prepdata
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c a driver routine to read and prepare data for the fit
c
c   Created by JAS Sep 11
c   Updated to f90 - JAS Jan 2015
c
      include 'commons.h'
c
c local variables
      character*1 a
      integer i, j, k
      real gtdst
c
c read in the data
      if( lvels )call readvels
      if( lphot )call readimg
c check input values and space
      if( lFITS )then
        if( ( lox .lt. 1 ) .or. ( lox - 1 + xrange .gt. nsizex ) .or.
     +      ( loy .lt. 1 ) .or. ( loy - 1 + yrange .gt. nsizey ) )then
          print *, lox, loy, lox + xrange, loy + yrange
          print *, 'Selected region partly outside input data'
          print *, 'Image/map size is', nsizex, nsizey
          call crash( 'prepdata' )
        end if
      end if
c select the data to fit
      if( lvels )call prepvels( gtdst )
      if( lphot )call prepimg( gtdst )
c return space that is no longer needed
      if( allocated( ldat ) )deallocate ( ldat )
      if( allocated( ldate ) )deallocate ( ldate )
      if( allocated( x2res ) )deallocate ( x2res )
c check pixel count
      if( pixct .eq. 0 )then
        print *, 'prepdata returned 0 pixels!!!'
        print *, 'Check your boxes, initials guesses, etc'
        call crash( 'prepdata' )
      end if
c check for sensible size of outer ellipse
      do i = 1, nellip
        if( gtdst + 1. .lt. sma( i ) )then
          print *, 'Your largest ellipse has a sma of', sma( nellip )
          if( i .lt. nellip )print *,
     +                       'and ellipse', i, ' has a sma of', sma( i )
          print *, 'while selected pixels extend only to', gtdst
          k = 0
    2     k = k + 1
          if( k .gt. 5 )then
            print *, 'Too many tries'
            call crash( 'prepdata' )
          end if
          print *, 'Enter q to quit or c to continue'
          read '( a )', a
c convert to uppercase if needed
          j = ichar( a )
          if( ( j .gt. 96 ) .and. ( j .lt. 123 ) )then
            a = char( j - 32 )
          end if
          if( a .eq. 'Q' )call crash( 'prepdata' )
          if( a .eq. 'C' )go to 1
          go to 2
        end if
      end do
    1 continue
      return
      end
