      subroutine vGHfit( ix, iy, pars, ifail )
c  Copyright (C) 2020, Jerry Sellwood
      use aacube
      implicit none
c
c calling arguments
      integer ix, iy, ifail( 7 )
      real pars( 8, 2 )
c
c common block
c
      include 'comcube.h'
c
c local array
      real, allocatable, save :: iv( : )
c
c local variables
      integer i
      logical firstc
      save firstc
      data firstc / .true. /
c
      if( firstc )then
        allocate ( iv( nsizev ) )
        firstc = .false.
      end if
c get intensity at each velocity
      do i = 1, nsizev
        iv( i ) = dcube( ix, iy, i )
      end do
c fit line profile
      call fitvGH( iv, nsizev, pars, ifail )
      return
      end

