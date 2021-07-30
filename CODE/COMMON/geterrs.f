      subroutine geterrs( aparams, afitval, n, lfrac, eparams, m )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c make uncertainty estimates from bootstrap iterations
c   uses parameter sdev returned from Numerical Recipes routine _moment
c   called from bootstrap and bootlace
c
c  Created by cutting code from bootstrap.f by JAS - June 2012
c  Made indata an allocatable array - JAS Aug 2015
c  Fixed bug in error on angles     - JAS Mar 2017
c  Fix-up for angles -> fix-angles  - JAS Jul 2017
c  Use the bi-weight to compute uncertainties- JAS Aug 2017
c  Reordered subscripts of calling args - JAS Dec 2017
c
      include 'commons.h'
c
c calling arguments
      integer m, n
      real*8 afitval( n, m ), aparams( n, m ), eparams( md )
      real*8 lfrac( 3, m )
c
c local array
      real*8, allocatable :: indata(:)
c
c local variables
      integer i, j, k
      logical lc
      real*8 ave, sdev
c
      allocate ( indata( nunc ) )
c compute mean and standard deviation of outputs
      do j = 1, nd
        do i = 1, nunc
          indata( i ) = aparams( j, i )
        end do
        call biwght2( indata, nunc, ave, sdev )
        eparams( j ) = sdev
      end do
c uncertainties in fitted values
      do j = 1, ntot
        do i = 1, nunc
          indata( i ) = afitval( j, i )
        end do
        call biwght2( indata, nunc, ave, sdev )
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
            call biwght2( indata, nunc, ave, sdev )
            if( k .eq. 1 )edskfrac = sdev
            if( k .eq. 2 )ebarfrac = sdev
            if( k .eq. 3 )eblgfrac = sdev
          end if
        end do
      end if
      return
      end
