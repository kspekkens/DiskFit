      subroutine crash( string )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
      implicit none
c provides an exit when a problem has been identified
c
      character*(*) string
c
      print '( ''crash called from routine '', a )', string
      stop
      end
