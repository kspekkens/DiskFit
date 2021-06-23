      subroutine setpars( params )
c set initial guesses for non-linear parameters
c   called from diskfit and bootstrap
c
c  Created by JAS - March 2012
c
      include 'commons.h'
c
c calling argument
      real*8 params( md )
c
c local variable
      integer j
c
      j = 0
      if( lpa )then
        params( j + 1 ) = pa
        j = j + 1
      end if
      if( leps )then
        params( j + 1 ) = eps
        j = j + 1
      end if
      if( lcentre )then
        params( j + 1 ) = xcen
        params( j + 2 ) = ycen
        j = j + 2
      end if
      if( lvels )then
c set initial pa of non-axisymmetric distortions
c Two changes from v1:
c   - phib is rel. to disk major axis phibprime is rel. to sky coords
c     to be consistent with notation in SS07
c   - initial guess is now an input parameter
        if( lnax )then
          if( lphib )then
            params( j + 1 ) = phib
            j = j + 1
          end if
        end if
c set initial warp parameters
        if( lwarp )then
          if( lrwarp )then
            params( j + 1 ) = rwarp
            j = j + 1
          end if
          if( lwepsm )then
            params( j + 1 ) = wepsm
            j = j + 1
          end if
          if( lwpm )then
            params( j + 1 ) = wphim
            j = j + 1
          end if
        end if
      end if
      if( lphot )then
c set initial guesses for bar pa and projected eccentriciy
        if( lnax )then
          if( lphib )then
            params( j + 1 ) = bar_pa
            j = j + 1
          end if
          if( lepsb )then
            params( j + 1 ) = bar_eps
            j = j + 1
          end if
        end if
c set initial guesses for bulge parameters
        if( lbulge )then
          if( lbleps )then
            params( j + 1 ) = bulge_el
            j = j + 1
          end if
          if( lsersn )then
            params( j + 1 ) = bulge_n
            j = j + 1
          end if
          if( lr_e )then
            params( j + 1 ) = r_bulge
            j = j + 1
          end if
        end if
      end if
      if( j .ne. nd )then
        print *, 'mismatched number of parameters'
        call crash( 'setpars' )
      end if
      return
      end
