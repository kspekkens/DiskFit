      module alldata
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c allocatable arrays that are accessible globally
      integer, allocatable, save :: iwgt( :, :, : )
      logical, allocatable, save :: inmask( :, : )
      logical, allocatable, save :: lgpix( :, : )
      real, allocatable, save :: barint( :, : )
      real, allocatable, save :: bivrad( :, :) 
      real, allocatable, save :: bivtan( :, : )
      real, allocatable, save :: bulgeint( :, : )
      real, allocatable, save :: diskint( :, : )
      real, allocatable, save :: dmask( :, : )
      real, allocatable, save :: ldat( :, : )
      real, allocatable, save :: ldate( :, : )
      real, allocatable, save :: model( :, : )
      real, allocatable, save :: res( :, : )
      real, allocatable, save :: radvels( :, : )
      real, allocatable, save :: rotvels( :, : )
      real, allocatable, save :: sdat( :, : )
      real, allocatable, save :: sdate( :, : )
      real, allocatable, save :: sigma( :, : )
      real, allocatable, save :: xval( : )
      real, allocatable, save :: yval( : )
      real, allocatable, save :: wgt( :, :, : )
      real, allocatable, save :: x2res( :, : )
      end module alldata
