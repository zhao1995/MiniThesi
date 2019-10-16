c$Id:$
      subroutine plftyp(pstyp,nel,iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'nel.eq.20' option                           15/02/2007
c       2. Add 'nel.eq.14/15' option                        06/09/2007
c       3. Add user options for pstyp < 0                   31/08/2008
c       4. Add 'nel.eq.64' cubic brick option               07/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set the plot for this element

c      Inputs:
c         pstyp   - Element topology
c         nel     - Number of element nodes
c         iel     - Element number

c      Output:
c         inord   - Number of plot
c-----[--.----+----.----+----.-----------------------------------------]
       implicit   none

       integer    pstyp,nel,iel

c      One dimensional plot sets

       if(    pstyp.eq.1) then
         if(nel.eq.1) then
           call pltpt1(iel)
         elseif(nel.eq.3) then
           call pltln3(iel)
         else
           call pltln2(iel)
         endif

c      Two dimensional plot sets

       elseif(pstyp.eq.2) then
         if(nel.eq.3) then
           call pltri3(iel)
         elseif(nel.eq.6 .or. nel.eq.7) then
           call pltri6(iel)
         elseif(nel.eq.8 .or. nel.eq.9) then
           call plqud8(iel)
         elseif(nel.eq.16) then
           call pltq16(iel)
         else
           call plqud4(iel)
         endif

c      Three dimensional plot sets

       elseif(pstyp.eq.3) then
         if(nel.eq.4) then
           call pltet4(iel)
         elseif(nel.eq.10 .or. nel.eq.11) then
           call pltet10(iel)
         elseif(nel.eq.14 .or. nel.eq.15) then
           call pltet10(iel)
         elseif(nel.eq.20 .or. nel.eq.27) then
           call plbk27(iel)
         elseif(nel.eq.64) then
           call plbkpqr(3,iel)
         else
           call plbrk8(iel)
         endif

c      User plots

       elseif(pstyp.lt.0) then

         call uftyplib(pstyp,nel,iel)

       endif

       end
