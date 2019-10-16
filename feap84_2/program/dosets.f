c$Id:$
      subroutine dosets(inp,otp,res,sav,plt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Declare characters as 128                        12/01/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Add initial character to default files using input
c               filename

c      Inputs:
c         inp  - Input   filename without inital character
c         otp  - Output  filename without inital character
c         res  - Restart filename without inital character
c         sav  - Save    filename without inital character
c         plt  - Plot    filename without inital character

c      Outputs:
c         inp  - Input   filename with inital character
c         otp  - Output  filename with inital character
c         res  - Restart filename with inital character
c         sav  - Save    filename with inital character
c         plt  - Plot    filename with inital character
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character inp*128,otp*128,res*128,sav*128,plt*128

      save

      if(otp.eq.' ') then
        otp      = inp
        otp(1:1) = 'O'
      endif
      if(res.eq.' ') then
        res      = inp
        res(1:1) = 'R'
      endif
      if(sav.eq.' ') then
        sav      = inp
        sav(1:1) = 'R'
      endif
      if(plt.eq.' ') then
        plt = inp
        plt(1:1) = 'P'
      endif

      end
