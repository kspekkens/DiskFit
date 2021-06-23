      subroutine setwgt( params, everyp, km )
c driving routine to interpret parameter list and set weights
c   Extracted from func - JAS March 2012
c
c common block
      include 'commons.h'
c
c calling arguments
      integer km
      logical everyp
      real*8 params( md )
c
c local variables
      integer i, ix, iy, j
      logical warpon
      real xcn, xpos, ycn, ypos
      real*8 barel, barpa, diskel, diskpa, naxphi, rw, welm, wpm
c
c interpret calling parameter list
      j = 0
      if( lpa )then
        diskpa = params( j + 1 )
        j = j + 1
      else
        diskpa = pa
      end if
      if( leps )then
        diskel = params( j + 1 )
        j = j + 1
      else
        diskel = eps
      end if
      cphi = cos( diskpa )
      sphi = sin( diskpa )
      b_over_a = 1.d0 - diskel
      if( lcentre )then
        xcn = params( j + 1 )
        ycn = params( j + 2 )
        j = j + 2
      else
        xcn = xcen
        ycn = ycen
      end if
      if( lvels )then
        if( lnax )then
          if( lphib )then
            naxphi = params( j + 1 )
            j = j + 1
          else
            naxphi = phib
          end if
          cosphib = cos( naxphi )
          sinphib = sin( naxphi )
          cos2phib = 1.d0 - 2.d0 * ( sinphib**2 )
          sin2phib = 2.d0 * sinphib * cosphib
        end if
c set sini - q = 0 is hardwired for now
        if( b_over_a .gt. q )then
          sini = sqrt( 1.d0 - ( b_over_a**2 - q**2 ) / ( 1.d0 - q**2 ) )
        else
          sini = 1
        end if
      end if
c options for photometry
      if( lphot )then
        if( lnax )then
          if( lphib )then
            barpa = params( j + 1 )
            j = j + 1
          else
            barpa = bar_pa
          end if
          if( lepsb )then
            barel = params( j + 1 )
            j = j + 1
          else
            barel = bar_eps
          end if
          cosphib = cos( barpa )
          sinphib = sin( barpa )
          b_over_a2 = 1.d0 - barel
        end if
        if( lbulge )then
          if( lbleps )then
            bulgel = params( j + 1 )
            j = j + 1
          else
            bulgel = bulge_el
          end if
          if( lsersn )then
            bulgen = params( j + 1 )
            j = j + 1
          else
            bulgen = bulge_n
          end if
          if( lr_e )then
            rbulge = params( j + 1 )
            j = j + 1
          else
            rbulge = r_bulge
          end if
        end if
      end if
c simple warp model
      if( lwarp )then
c inner radius of warp
        if( lrwarp )then
          rw = params( j + 1 )
          j = j + 1
        else
          rw = rwarp
        end if
c maximum eccentricity of warp
        if( lwepsm )then
          welm = params( j + 1 )
          j = j + 1
        else
          welm = wepsm
        end if
c max position angle of warp
        if( lwpm )then
          wpm = params( j + 1 )
          j = j + 1
        else
          wpm = wphim
        end if
      end if
c check that the number of parameters set agrees with the number expected
      if( j .ne. nd )then
        print *, 'logical error in number of parameters'
        call crash( 'setwgt' )
      end if
c set warp data arrays if needed
      warpon = lwarp .and. ( rw .lt. sma( nellip ) )
      if( warpon )call warpset( diskel, diskpa, rw, wpm, welm )
c compute weights for pixels in the map
      km = 0
      if( l2D )then
        do iy = 1, yrange
          ypos = yval( iy ) - ycn
          do ix = 1, xrange
            xpos = xval( ix ) - xcn
c need all pixels if requested, only good pixels otherwise
            if( everyp .or. lgpix( ix, iy ) )then
              if( lphot )call setwgtp(
     +             xpos, ypos, km, iwgt( 1, ix, iy ), wgt( 1, ix, iy ) )
              if( lvels )then
                if( warpon )then
                  call setwgtw(
     +             xpos, ypos, km, iwgt( 1, ix, iy ), wgt( 1, ix, iy ) )
                else
                  call setwgtv(
     +             xpos, ypos, km, iwgt( 1, ix, iy ), wgt( 1, ix, iy ) )
                end if
              end if
            end if
          end do
        end do
      else
        do i = 1, inp_pts
          ypos = yval( i ) - ycn
          xpos = xval( i ) - xcn
c no seeing corrections possible
          if( everyp .or. lgpix1( i ) )then
            if( lphot )call setwgtp(
     +                     xpos, ypos, km, iwgt1( 1, i ), wgt1( 1, i ) )
            if( lvels )then
              if( warpon )then
                call setwgtw(
     +                     xpos, ypos, km, iwgt1( 1, i ), wgt1( 1, i ) )
              else
                call setwgtv(
     +                     xpos, ypos, km, iwgt1( 1, i ), wgt1( 1, i ) )
              end if
            end if
          end if
        end do
      end if
      if( km .le. 0 )then
        print *, 'Something went wrong in setwgt'
        call crash( 'func' )
      end if
c blur the weights to allow for seeing
      if( lseeing )call blurmod( km )
c check that array was large enough
      if( km .gt. mk )then
        print *, 'value of km', km, ' returned from setwgt or blurmod'
        print *, 'current max allowed value =', mk
        print *, 'increase parameter mk in commons.h and recompile'
        call crash( 'setwgt' )
      end if
      return
      end
