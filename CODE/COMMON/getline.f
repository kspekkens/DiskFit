      subroutine getline( iunit, line, istat, ic )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c utility routine to read in the next significant line of ASCII data from a
c   specified logical unit.   The # character is interpreted as the beginning
c   of a comment and it is eliminated together with all subsequent characters
c   blank lines and lines containing only comments are skipped
c
c   Created by JAS March 2010
c
      implicit none
c
c calling arguments
      character*(*) line
      integer ic, istat, iunit
c
c external
      integer lnblnk
c
c local variables
      integer i, j, k
c
c read, and keep count of, the next line from this file
    1 ic = ic + 1
      read( iunit, '( a )', iostat = istat )line
c successful read
      if( istat .eq. 0 )then
c skip completely blank lines
        k = lnblnk( line )
        if( k .eq. 0 )go to 1
c find comment delimeter
        j = index( line, '#' )
c skip comment lines
        if( j .eq. 1 )go to 1
c clear away comment
        if( j .gt. 0 )then
          do i = j, k
            line( i:i ) = ' '
          end do
        end if
      else
c read failed for some reason
        print *, 'istat =', istat, ' from read in getline'
      end if
      return
      end
