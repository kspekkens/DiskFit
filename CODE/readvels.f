      subroutine readvels
c reads in a dataset to xval, yval, vels2, evel2
c
c   Orignal version requires a text file - KS
c   Input fits files - RZS March 2009
c   Created option to read either - JAS October 2009
c
      include 'commons.h'
c
c input appropriate type of file
      if( lFITS )then
        call rvels_FITS
      else if( ltext )then
        call rvels_text
      else
        print *, 'Are vels data in neither a FITS nor a text file?'
        call crash( 'readvels' )
      end if
      return
      end
