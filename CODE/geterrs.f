      subroutine geterrs( aparams, afitval, lfrac, eparams, m )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c make uncertainty estimates from bootstrap iterations
c   uses parameter sdev returned from Numerical Recipes routine _moment
c   called from bootstrap and bootlace
c
c  Created by cutting code from bootstrap.f by JAS - June 2012
c  Made indata an allocatable array - JAS Aug 2015
c
      include 'commons.h'
c
c calling arguments
      integer m
      real*8 afitval( m, mtot ), aparams( m, md ), eparams( md )
      real*8 lfrac( 3, m )
c
c local arrays
      logical langle( md )
      real*8, allocatable :: indata(:)
c
c local variables
      integer i, j, k
      logical lc
      real*8 a, adev, ave, b, curt, sdev, skew, var
c
c set flags for non-angle parameters
      do j = 1, nd
        langle( j ) = .false.
      end do
c flag angle parameters that need more spohisticated treatment
      j = 0
      if( lpa )then
        j = j + 1
        langle( j ) = .true.
      end if
      if( leps )j = j + 1
      if( lcentre )j = j + 2
      if( lphot )then
        if( lnax )then
          if( lphib )then
            j = j + 1
            langle( j ) = .true.
          end if
          if( lepsb )j = j + 1
        end if
        if( lbulge )then
          if( lbleps )j = j + 1
          if( lsersn )j = j + 1
          if( lr_e )j = j + 1
        end if
      else
        if( lnax )then
         if( lphib )then
           j = j + 1
           langle( j ) = .true.
         end if
        end if
        if( lwarp )then
          if( lrwarp )j = j + 1
          if( lwepsm )j = j + 1
          if( lwpm )j = j + 1
        end if
      end if
      if( j .ne. nd )then
        print *, 'mismatched number of parameters'
        call crash( 'geterrs' )
      end if
      allocate ( indata( nunc ) )
c compute mean and standard deviation of outputs
      do j = 1, nd
        do i = 1, nunc
          indata( i ) = aparams( i, j )
        end do
c allow for restricted range of angular parameters
        if( langle( j ) )then
          a = 0
          b = 0
          do i = 1, nunc
            indata( i ) = mod( indata( i ), pi )
            if( indata( i ) .lt. 0. )indata( i ) = indata( i ) + pi
            a = a + cos( 2. * indata( i ) )
            b = b + sin( 2. * indata( i ) )
          end do
c pilot estimate of mean
          ave = .5 * atan2( b, a )
          if( ave .lt. 0. )ave = ave + pi 
c ensure distribution is symmetric about this mean
          a = .5 * pi
          do i = 1, nunc
            if( indata( i ) - ave .gt. a )indata( i ) = indata( i ) - pi
            if( indata( i ) - ave .lt.-a )indata( i ) = indata( i ) + pi
          end do
        end if
        call moment_dbl(
     +                  indata, nunc, ave, adev, sdev, var, skew, curt )
        eparams( j ) = sdev
      end do
c uncertainties in fitted values
      do j = 1, ntot
        do i = 1, nunc
          indata( i ) = afitval( i, j )
        end do
        call moment_dbl(
     +                  indata, nunc, ave, adev, sdev, var, skew, curt )
        efitval( j ) = sdev
      end do
c separate efitval array into disk, bar and bulge
      k = 0
      if( ldisk )then
        do j = 1, nellip
          eid( j ) = efitval( k + j )
        end do
        k = k + nellip
      end if
      if( lphot )then
        if( lnax )then
          do j = 1, nbar
            if( ( j .lt. nminr ) .or. ( j .gt. nmaxr ) )then
              eibar( j ) = 0.0
            else
              i = j - nminr + 1
              eibar( j ) = efitval( k + i )
            end if
          end do
          k = k + nbar
        end if
        if( lbulge )then
          do j = 1, nbulge
            eibulge( j ) = efitval( k + j )
          end do
          k = k + nbulge
        end if
      end if
      if( lvels )then
        if( lnax )then
          do j = 1, nellip
            if( ( j .lt. nminr ) .or. ( j .gt. nmaxr ) )then
              eibirad( j ) = 0
              eibitan( j ) = 0
            else
              i = j - nminr + 1
              eibitan( j ) = efitval( k + i )
              eibirad( j ) = efitval( k + nnasymm + i )
            end if
          end do
          k = k + 2 * nnasymm
        end if
        if( lradial )then
          do j = 1, nellip
            if( ( j .lt. nminr ) .or. ( j .gt. nmaxr ) )then
              eirad( j ) = 0.
            else
              i = j - nminr + 1
              eirad( j ) = efitval( k + i )
            end if
          end do
          k = k + nradial
        end if
        if( lsystemic )then
          eivsys = efitval( ntot )
          k = k + 1
        end if
      end if
      if( k .ne. ntot )then
        print *, 'Mismatch between k and ntot'
        print *, 'Contact code administrator'
        call crash( 'geterrs' )
      end if
c uncertainties in light fractions
      if( lphot )then
        edskfrac = 0
        ebarfrac = 0
        eblgfrac = 0
        do k = 1, 3
          if( k .eq. 1 )lc = ldisk .and. ( lnax .or. lbulge )
          if( k .eq. 2 )lc = lnax
          if( k .eq. 3 )lc = lbulge
          if( lc )then
            do i = 1, nunc
              indata( i ) = lfrac( k, i )
            end do
            call moment_dbl(
     +                  indata, nunc, ave, adev, sdev, var, skew, curt )
            if( k .eq. 1 )edskfrac = sdev
            if( k .eq. 2 )ebarfrac = sdev
            if( k .eq. 3 )eblgfrac = sdev
          end if
        end do
      end if
      return
      end
