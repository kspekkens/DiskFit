      subroutine timestamp( str )
c Copyright (C) 2020, Jerry Sellwood
      implicit none
c
c calling argument
      character*(*) str
c
c local variables
      character date*8, time*10
c
c f90 external
      call date_and_time( date, time )
c format the text string
      str = date( 1:4 ) // '/' // date( 5:6 ) // '/' // date( 7:8 ) //
     + ' ' // time( 1:2 ) // ':' // time( 3:4 ) // ':' // time( 5:6 )
      return
      end
