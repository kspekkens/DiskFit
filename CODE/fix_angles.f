      subroutine fix_angles( aparams, m )
c Copyright (C) 2017, Jerry Sellwood and Kristine Spekkens
c
c eliminates possible artificial bimodality in the distribution of angle
c   parameters
c   called from bootstrap and bootlace
c
c  Created by cutting code from geterrs.f by JAS - July 2017
c
      include 'commons.h'
c
c calling arguments
      integer m
      real*8 aparams( m, md )
c
c local arrays
      integer msect( md )
      logical langle( md )
      real*8, allocatable :: indata(:)
c
c local variables
      integer i, j
      real*8 a, am, ave, b
c
c set default values
      do j = 1, nd
        langle( j ) = .false.
        msect( j ) = 1
      end do
c set flags for angle parameters
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
            msect( j ) = order
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
            msect( j ) = order
          end if
        end if
        if( lwarp )then
          if( lrwarp )j = j + 1
          if( lwepsm )then
            j = j + 1
            langle( j ) = .true.
          end if
          if( lwpm )j = j + 1
        end if
      end if
      if( j .ne. nd )then
        print *, 'mismatched number of parameters'
        call crash( 'fix_angles' )
      end if
      allocate ( indata( nunc ) )
c deal with possible artificial bimodality of angle parameters 
      do j = 1, nd
        if( langle( j ) )then
          do i = 1, nunc
            indata( i ) = aparams( i, j )
          end do
          am = msect( j )
          a = 0
          b = 0
          do i = 1, nunc
            indata( i ) = mod( indata( i ), pi )
            a = a + cos( am * indata( i ) )
            b = b + sin( am * indata( i ) )
          end do
c pilot estimate of mean
          ave = atan2( b, a ) / am
          if( ave .lt. 0. )ave = ave + 2. * pi / am
c ensure distribution is symmetric about this mean
          a = .5 * pi / am
          do i = 1, nunc
            if( indata( i ) - ave .gt.  a )indata( i ) =
     +                                     indata( i ) - pi / am
            if( indata( i ) - ave .lt. -a )indata( i ) =
     +                                     indata( i ) + pi / am
          end do
c replace values
          do i = 1, nunc
            aparams( i, j ) = indata( i )
          end do
        end if
      end do
      return
      end
