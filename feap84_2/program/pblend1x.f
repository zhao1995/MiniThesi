c$Id:$
      subroutine pblend1x(nn,nr,ni,ndm, fxim,nty,x,nflag,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct format spacings                          21/10/2007
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Construct one dimensional interpolation using blending

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      logical    nflag,prt,prth
      integer    i,k, nr,ni, nn,ndm, nty(*)
      real*8     n1i,n2i, rnr, fxim(ndm,0:nr), x(ndm,*)

      save

      nn  = ni - 1
      rnr = 1.d0/dble(nr)

      if(prt) then
        call prtitl(prth)
        write(iow,2000) (i,i=1,ndm)
      endif
      do i = 0,nr
        nn = nn + 1
        if(nflag .or. nty(nn).ge.0) then
          n2i     = dble(i)*rnr
          n1i     = 1.d0 - n2i
          nty(nn) = 0
          do k = 1,ndm
c           x(k,nn) = fxim(k,i) - n1i*fxim(k,0) - n2i*fxim(k,nr)
            x(k,nn) = fxim(k,i)
          end do ! k
          if(prt) then
            write(iow,2001) nn,(x(k,nn),k=1,ndm)
          endif
        end if
      end do ! i

c     Formats

2000  format('   B l e n d e d   C o o r d i n a t e s'//
     &       '     Node',4(i4,'-Coordinate':))

2001  format(i8,1x,1p,4e15.5)

      end
