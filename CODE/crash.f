      subroutine crash( string )
      implicit none
c
      character*(*) string
c
      print '( ''crash called from routine '', a )', string
      stop
      end
