      subroutine warpset( diskel, diskpa, rw, wpm, welm )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c sets the values of position angle and inclination of each ellipse
c   for the selected warp model
c
c created by JAS July 2011
c
      include 'commons.h'
c
c calling arguments
      real*8 diskel, diskpa, rw, welm, wpm
c
c local variables
      integer i
      real f
c
c work over all ellipses
      do i = 1, nellip
        if( sma( i ) .le. rw )then
          wel( i ) = diskel
          wsi( i ) = sini
          wba( i ) = b_over_a
          wphi( i ) = diskpa
          wcp( i ) = cphi
          wsp( i ) = sphi
        else
c quadratic function from rw to disk edge
          f = ( sma( i ) - rw ) / ( sma( nellip ) - rw )
          wel( i ) = diskel + f**2 * welm
          wphi( i ) = diskpa + f**2 * wpm
          if( wel( i ) .lt. 0. )then
            wel( i ) = 1. - 1. / ( 1. - wel( i ) )
            wphi( i ) = wphi( i ) + .5 * pi
          end if
          wba( i ) = 1.d0 - wel( i )
          if( wba( i ) .gt. q )then
            wsi( i ) =
     +           sqrt( 1.d0 - ( wba( i )**2 - q**2 ) / ( 1.d0 - q**2 ) )
          else
            wsi( i ) = 1
          end if
          wcp( i ) = real( dcos( dble( wphi( i ) ) ) )
          wsp( i ) = real( dsin( dble( wphi( i ) ) ) )
        end if
      end do
      return
      end
