c$Id:$
      subroutine faceset(d,ie,ix,ic,ir,intel,nn,hn1, hisflg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct dimension of ic to 'numnp'               02/09/2007
c       2. Change nih1,nih2,nih3 to h1ni,h2ni,h3ni int4     17/04/2011
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Set data list for element interfaces
c     Inputs:
c       ix(*)  - Nodes connected to each element
c       ic(*)  - List of pointers for element connections to nodes
c       ir(*)  - List of elements connected to nodes
c       hisflg - Determine history storage if true

c     Outputs:
c       intel  - Interface data storage
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'facset.h'
      include   'hdatam.h'
      include   'ieldat.h'
      include   'iofile.h'
      include   'sdata.h'

      logical    matchf, hisflg
      integer    i, n,nn, mi,mj,il,ixi, nec1,nec2, ftyp1,ftyp2,factyp
      integer    hn1, h1ni,h2ni,h3ni
      integer    ie(nie,*),ix(nen1,numel),ic(0:numnp),ir(*),intel(8,*)
      integer    intnod(5)
      real*8     d(ndd,*),ul1(1),ul2(1),xl1(1),xl2(1),tl1(1),tl2(1)
      real*8     s(1),p(1)

      save

      nn      = 0
      hn1     = 0
      hnimax  = 0
      hni3max = 0
      write(iow,2000)
      do n = 1,numel
        ftyp1 = factyp(ix(1,n), nec1)
        do mi = 1,nec1
          call setfac(ftyp1,nec1, mi,ixi,ix(1,n))
          do il = ic(ixi-1)+1,ic(ixi)
            ftyp2 = factyp(ix(1,ir(il)), nec2)
            if(ftyp1.eq.ftyp2) then
              do mj = 1,nec2
                if(matchf(ftyp1,nec2, mj,ix(1,ir(il)))) then
                  if(n.lt.ir(il)) then
                    write(iow,2001) n,ir(il),ichk

                    nn = nn + 1
                    intel(1,nn) = n
                    intel(2,nn) = ir(il)
                    intel(3,nn) = ichk(1)
                    intel(4,nn) = ichk(2)
                    intel(5,nn) = ichk(3)

                    intnod(1)   = ichk(1)
                    intnod(2)   = ichk(2)
                    intnod(3)   = ichk(3)
                    intnod(4)   = n
                    intnod(5)   = ir(il)

c                   Check for history terms on interfaces

                    if(hisflg) then
                      ma1  = ix(nen1,n)
                      ma2  = ix(nen1,ir(il))
                      iel1 = ie(nie-1,ma1)
                      iel2 = ie(nie-1,ma2)
                      h1ni = 0
                      h2ni = 0
                      h3ni = 0
                      nsts = 1
                      call ielmlib(d(1,ma1),ul1,xl1,ix(1,n),tl1,
     &                             d(1,ma2),ul2,xl2,ix(1,ir(il)),tl2,
     &                             intnod,s,p,1)

                      hnimax  = max(hnimax ,h1ni)
                      hni3max = max(hni3max,h3ni)
                      intel(6,nn) = hn1
                      intel(7,nn) = hn1         + h1ni
                      intel(8,nn) = intel(7,nn) + h1ni
                      hn1         = hn1 + h1ni*2 + h3ni
                    endif ! hisflg
                    go to 100
                  endif
                endif
              end do ! mj
            endif
          end do ! il

c         Boundary segment:

          do il = nn,1,-1
            if(intel(2,il).eq.n .and. intel(3,il).eq.ichk(1)
     &                          .and. intel(4,il).eq.ichk(2)
     &                          .and. intel(5,il).eq.ichk(3)) then
              go to 100
            endif
          end do ! il
          nn          = nn + 1
          intel(1,nn) = n
          intel(2,nn) = 0
          intel(3,nn) = ichk(1)
          intel(4,nn) = ichk(2)
          intel(5,nn) = ichk(3)
          write(iow,2001) (intel(i,nn),i=1,5)

          intnod(1)   = ichk(1)
          intnod(2)   = ichk(2)
          intnod(3)   = ichk(3)
          intnod(4)   = n
          intnod(5)   = 0

c         Check for history terms on interfaces

          if(hisflg) then
            ma1  = ix(nen1,n)
            ma2  = ma1
            iel1 = ie(nie-1,ma1)
            iel2 = 0
            h1ni = 0
            h2ni = 0
            h3ni = 0
            nsts = 1
            call ielmlib(d(1,ma1),ul1,xl1,ix(1,n),tl1,
     &                   d(1,ma2),ul2,xl2,ix(1,ir(il)),tl2,
     &                   intnod,s,p,1)

            hnimax      = max(hnimax ,h1ni)
            hni3max     = max(hni3max,h3ni)
            intel(6,nn) = hn1
            intel(7,nn) = hn1         + h1ni
            intel(8,nn) = intel(7,nn) + h1ni
            hn1         = hn1         + h1ni*2 + h3ni
          endif ! hisflg
100       continue
        end do ! mi
      end do ! n

c     Place marker end

      if(hisflg) then
        write(iow,*) ' Total matching faces =',nn
        nn = nn + 1
        do n = 1,8
          intel(n,nn) = 0
        end do ! n
        intel(6,nn) = hn1
      endif ! hisflg

c     Formats

2000  format(/7x,'M a t c h e d   F a c e s   f o r   E l e m e n t s'//
     &       10x,'  Matching Elements         Face Node Numbers'/
     &       10x,'1-Element   2-Element   1-Node    2-Node    3-Node')

2001  format(5x,2i12,3i10)

      end
