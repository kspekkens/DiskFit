      subroutine warpinit( params )
c Copyright (C) 2020, Jerry Sellwood and Kristine Spekkens
c
c re-creates the working arrays of the best fitting warp
c
c   Created by JAS Aug 2020
c
      include 'commons.h'
c
c calling argument
      real*8 params( md )
c
c local variables
      integer k
      real*8 diskel, diskpa, rw, welm, wpm
c
      if( .not. lwarp )then
        print *, 'warp not on'
        call crash( 'WARPINIT' )
      end if
c get fitted or default parameters
      k = 0
      if( lpa )then
        diskpa = params( k + 1 )
        k = k + 1
      else
        diskpa = pa
      end if
      if( leps )then
        diskel = params( k + 1 )
        k = k + 1
      else
        diskel = eps
      end if
c these params are not needed for this routine
      if( lcentre )k = k + 2
      if( lnax .and. lphib )k = k + 1
      if( lrwarp )then
        rw = params( k + 1 )
        k = k + 1
      else
        rw = rwarp
      end if
      if( lwepsm )then
        welm = params( k + 1 )
        k = k + 1
      else
        welm = wepsm
      end if
      if( lwpm )then
        wpm = params( k + 1 )
        k = k + 1
      else
        wpm = wphim
      end if
      if( k .ne. nd )then
        print *, 'Mismatch between k and nd'
        print *, 'Contact code administrator'
        call crash( 'WARPINIT' )
      end if
c set up the warp arrays
      call warpset( diskel, diskpa, rw, wpm, welm )
      return
      end
