c$Id:$
      subroutine faclset(d,ie,ix,intel,nn,hn1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change nih1 to h1ni, nih3 to h3ni for int4 use   17/04/2011
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Set data list for element interfaces
c     Inputs:
c       ix(*)  - Nodes connected to each element

c     Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'facset.h'
      include   'hdatam.h'
      include   'ieldat.h'
      include   'iofile.h'
      include   'sdata.h'

      logical    errck, pinput, readin
      integer    i, n,nn, mi,mj,nec1,nec2, ftyp1,ftyp2,factyp
      integer    hn1, h1ni,h3ni, el1,el2
      integer    ie(nie,*),ix(nen1,numel),intel(8,*),intnod(5)
      real*8     d(ndd,*),ul1(1),ul2(1),xl1(1),xl2(1),tl1(1),tl2(1)
      real*8     s(1),p(1), td(4)

      save

      nn      = 0
      hn1     = 0
      hnimax  = 0
      hni3max = 0
      write(iow,2000)
      readin = .true.
      do while (readin)
        errck = pinput(td,4)
        el1 = nint(td(1))
        if(el1.gt.0) then
          ftyp1 = factyp(ix(1,el1), nec1)
          el2   = nint(td(3))
          if(el2.gt.0) then
            ftyp2 = factyp(ix(1,el2), nec2)
          else
            ftyp2 = 0
          endif

c         Same face types: Match

          if(ftyp1.eq.ftyp2) then

c           One-d elements

            if(ftyp1.eq.1) then
              ichk(1) = nint(td(2))
              ichk(2) = nint(td(4))
              ichk(3) = 0
              do mi = 1,nec1
                if(ix(mi,el1).eq.ichk(1)) then
                  do mj = 1,nec2
                    if(ix(mj,el2).eq.ichk(2)) then
                      nn = nn + 1
                      intel(1,nn) = el1
                      intel(2,nn) = el2
                      intel(3,nn) = ichk(1)
                      intel(4,nn) = ichk(2)
                      intel(5,nn) = ichk(3)

                      intnod(1)   = ichk(1)
                      intnod(2)   = ichk(2)
                      intnod(3)   = ichk(3)
                      intnod(4)   = el1
                      intnod(5)   = el2

                      if(ior.lt.0) then
                        write(*,2001) (intel(i,nn),i=1,5)
                      endif
                      write(iow,2001) (intel(i,nn),i=1,5)

c                     Check for history terms on interfaces

                      ma1  = ix(nen1,el1)
                      ma2  = ix(nen1,el2)
                      iel1 = ie(nie-1,ma1)
                      iel2 = ie(nie-1,ma2)
                      nsts = 1
                      call ielmlib(d(1,ma1),ul1,xl1,ix(1,el1),tl1,
     &                             d(1,ma2),ul2,xl2,ix(1,el2),tl2,
     &                             intnod,s,p,1)

                      hnimax  = max(hnimax ,h1ni)
                      hni3max = max(hni3max,h3ni)
                      intel(6,nn) = hn1
                      intel(7,nn) = hn1         + h1ni
                      intel(8,nn) = intel(7,nn) + h1ni
                      hn1         = hn1 + h1ni*2 + h3ni

                      go to 100
                    endif
                  end do ! mj
                endif
              end do ! mi
            endif

c         Boundary segment

          elseif(el2.eq.0) then
            if(ftyp1.eq.1) then
              ichk(1) = nint(td(2))
              ichk(2) = 0
              ichk(3) = 0

              nn = nn + 1
              intel(1,nn) = el1
              intel(2,nn) = el2
              intel(3,nn) = ichk(1)
              intel(4,nn) = ichk(2)
              intel(5,nn) = ichk(3)

              intnod(1)   = ichk(1)
              intnod(2)   = ichk(2)
              intnod(3)   = ichk(3)
              intnod(4)   = el1
              intnod(5)   = el2

              if(ior.lt.0) then
                write(*,2001) (intel(i,nn),i=1,5)
              endif
              write(iow,2001) (intel(i,nn),i=1,5)

c             Check for history terms on interfaces

              ma1  = ix(nen1,el1)
              ma2  = ma1
              iel1 = ie(nie-1,ma1)
              iel2 = 0
              nsts = 1
              call ielmlib(d(1,ma1),ul1,xl1,ix(1,el1),tl1,
     &                     d(1,ma2),ul2,xl2,ix(1,el1),tl2,
     &                     intnod,s,p,1)

              hnimax      = max(hnimax ,h1ni)
              hni3max     = max(hni3max,h3ni)
              intel(6,nn) = hn1
              intel(7,nn) = hn1         + h1ni
              intel(8,nn) = intel(7,nn) + h1ni
              hn1         = hn1         + h1ni*2 + h3ni
            endif

          endif

100       continue

        else
          readin = .false.
        endif
      end do ! while

      if(ior.lt.0) then
        write(*,*) ' Total matching faces =',nn
      endif
      write(iow,*) ' Total matching faces =',nn

c     Place marker end

      nn = nn + 1
      do n = 1,8
        intel(n,nn) = 0
      end do ! n
      intel(6,nn) = hn1

c     Formats

2000  format(/7x,'M a t c h e d   F a c e s   f o r   E l e m e n t s'//
     &       10x,'  Matching Elements         Face Node Numbers'/
     &       10x,'1-Element   2-Element   1-Node    2-Node    3-Node')

2001  format(5x,2i12,3i10)

      end
