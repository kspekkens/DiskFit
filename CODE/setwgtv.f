      subroutine setwgtv( xpos, ypos, km, iwt, wt )
c Sets the weight factors for a point at xpos, ypos relative to the
c   current disk center.  Called both from func & from writemod
c
c   Created by JAS by extracting code from func - October 2009
c   Modified by JAS March 2011 to reduce memory needed
c
c common block
      include 'commons.h'
c
c calling arguments
      integer km, iwt( mk )
      real xpos, ypos, wt( mk )
c
c local variables
      integer i, j, jj, k, kk
      real ai, costheta, cos2theta, cosmthetab, frac0, frac1
      real sintheta, sin2theta, sinmthetab, vwgt, xel, yel
c
c clear the wgt array
      do j = 1, mk
        iwt( j ) = 0
        wt( j ) = 0
      end do
      jj = 0
c compute semi-major axis of ellipse through this pixel
      xel = xpos * cphi + ypos * sphi
      yel = -xpos * sphi + ypos * cphi
      ai = sqrt( xel**2 + ( yel / b_over_a )**2 )
      if( ai .gt. 1.e-3 )then
        costheta = xel / ai
        sintheta = yel / ( b_over_a * ai )
      else
        costheta = 1
        sintheta = 0
      endif
      cos2theta = 1. - 2. * sintheta**2
      sin2theta = 2. * sintheta * costheta
c find nearest interior ellipse by bi-section - assumed in ascending order
c   this is substantially more efficient than the previous "call closest"
c
c setting j = nellip here requires k to be less than nellip, which is what
c   is needed for the following code.
c it would be OK to set j = nellip + 1 if desired, and then pixels outside
c   the largest ellipse would be flagged by k = nellip
c
c assume velocity fits always fit for mean orbital speed
      k = 1
      if( linter0 )k = 0
      j = nellip
      do while ( k + 1 .lt. j )
        i = ( k + j ) / 2
        if( sma( i ) .lt. ai )then
          k = i
        else
          j = i
        end if
      end do
c calculate the weights for the mean orbital speed
      j = 0
      vwgt = sini * costheta
c
c The following strategy differs from that in the velfit version 1 release
c   because that version created large residuals from pixels outside the
c   area between the smallest and largest ellipses.  In these areas the
c   model predictions were tapered to zero, whereas the observed
c   velocities are generally large.
c
c This code now reverts to the old strategy adopted by EIB, ASR & JAS,
c   which is to predict velocities by extrapolation from the last (or
c   first) two ellipses.  The model will therefore be constrained to
c   favor a constant velocity gradient from the last two ellipses that
c   extrapolates most consistently with the data.
c
c One could also opt to exclude altogether pixels that are not bracketed
c   by ellipses, but then care must be taken to ensure that those pixels
c   do not contribute to the chi^2 returned by this routine.  (See the
c   comment in the calculation of chi^2 below.)  I do not favor this option
c   because chi^2 will be reduced simply by inclining the model more so
c   that fewer pixels contribute.  (JAS, May-1-09)
c
c The following code implements linear interpolation between rings and
c   linear extrapolation to values inside the smallest or outside the
c   largest rings
c
      kk = k
      if( k .gt. 0 )then
        frac0 = ( ai - sma( k ) ) / ( sma( k + 1 ) - sma( k ) )
        frac1 = 1. - frac0
      else
c interpolate to zero at center inside innermost ring
        frac0 = ai / sma( 1 )
        frac1 = 0
      end if
      iwt( jj + 1 ) = j + kk
      iwt( jj + 2 ) = j + kk + 1
      if( fronly )then
        wt( jj + 1 ) = frac1
        wt( jj + 2 ) = frac0
      else
        wt( jj + 1 ) = vwgt * frac1
        wt( jj + 2 ) = vwgt * frac0
      end if
      jj = jj + 2
c
c The above discussion of what to do about pixels that are not bracketed
c   by ellipses is important only for the circular speed.  Contributions
c   from pixels exterior to this region can safely be excluded for radial
c   or non-axsymmetric flow fits since the dominant motion is presumed to
c   be purely orbital everywhere else. (JAS)
c
      j = j + nellip
c redefine k if the point is outside the largest ellipse
      if( ai .gt. sma( nellip ) )k = nellip
c adjust weights if the number of ellipses for these models is restricted
      if( lradial .or. lnax )then
c taper to zero one steps inside innermost ring
        if( k .lt. nminr )then
          frac1 = 0
          if( k .lt. nminr - 1 )frac0 = 0
        end if
c taper to zero one steps outside outermost ring
        if( k .ge. nmaxr )then
          frac0 = 0
          if( k .gt. nmaxr )frac1 = 0
        end if
      end if
c compute radial weights:
      if( lradial )then
        vwgt = sini * sintheta
c linear interpolation
        kk = k - nminr + 1
        if( ( kk .gt. -1 ) .and. ( kk .le. nradial ) )then
          if( kk .gt. 0 )then
            jj = jj + 1
            iwt( jj ) = j + kk
            wt( jj ) = vwgt * frac1
          end if
          if( kk .le. nradial )then
            jj = jj + 1
            iwt( jj ) = j + kk + 1
            wt( jj ) = vwgt * frac0
          end if
        end if
        j = j + nradial
      end if
c compute non-axisymm weights:
      if( lnax )then
        if( order .eq. 1 )then
          cosmthetab = costheta * cosphib + sintheta * sinphib
          sinmthetab = sintheta * cosphib - costheta * sinphib
        else
          cosmthetab = cos2theta * cos2phib + sin2theta * sin2phib
          sinmthetab = sin2theta * cos2phib - cos2theta * sin2phib
        end if
c tangential component first:
        vwgt = -sini * costheta * cosmthetab
        kk = k - nminr + 1
c linear interpolation
        if( ( kk .gt. -1 ) .and. ( kk .le. nnasymm ) )then
          if( kk .gt. 0 )then
            jj = jj + 1
            iwt( jj ) = j + kk
            wt( jj ) = vwgt * frac1
          end if
          if( kk .le. nnasymm )then
            jj = jj + 1
            iwt( jj ) = j + kk + 1
            wt( jj ) = vwgt * frac0
          end if
        end if
        j = j + nnasymm
c radial component next:
        vwgt = -sini * sintheta * sinmthetab
c linear interpolation
        if( ( kk .gt. -1 ) .and. ( kk .le. nnasymm ) )then
          if( kk .gt. 0 )then
            jj = jj + 1
            iwt( jj ) = j + kk
            wt( jj ) = vwgt * frac1
          end if
          if( kk .le. nnasymm )then
            jj = jj + 1
            iwt( jj ) = j + kk + 1
            wt( jj ) = vwgt * frac0
          end if
        end if
        j = j + nnasymm
      end if
c systemic velocity weights - they are all 1
      if( lsystemic )then
        jj = jj + 1
        iwt( jj ) = ntot
        wt( jj ) = 1
      end if
      km = max( jj, km )
      return
      end
