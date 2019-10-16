c$Id:$
      subroutine pspinin(spinv,spnum,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    09/03/2009
c       1. Add 'prt' parameter to control outputs           08/04/2009
c       2. Add 'edge' inputs to tables                      04/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Input data for spin about an axis

c     Input:
c        spnum           - Number for spin data
c        prt             - Print flag

c     Output:
c        spinv(12,spnum) - Spin data for number 'spnum'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'edgdat.h'
      include   'iofile.h'
      include   'pconstant.h'
      include   'sdata.h'

      character  text*15
      logical    prt, pcomp,errck,tinput, nflag
      integer    spnum, j
      real*8     spinv(12,*), td(6), nn

c     Input data for spin about axes

      espfl = .false.
      nflag = .false.
      text = 'xxxx'
      do while (.not.pcomp(text,'    ',4))
        errck = tinput(text,1,td,6)
        if(pcomp(text,'cent',4)) then
          do j = 1,3
            spinv(j,spnum) = td(j)
          end do ! j
        elseif(pcomp(text,'norm',4)) then
          nn = sqrt(td(1)**2 + td(2)**2 + td(3)**2)
          if(nn.gt.0.0d0) then
            do j = 1,3
              spinv(j+3,spnum) = td(j)/nn
            end do ! j
          else
            write(ilg,3000)
            write(iow,3000)
            call plstop()
          endif
          nflag = .true.
        elseif(pcomp(text,'rate',4)) then
          spinv(7,spnum) = td(1)
        elseif(pcomp(text,'prop',4)) then
          spinv(8,spnum) = td(1)
        elseif(pcomp(text,'velo',4)) then
          do j = 1,min(ndf,3)
            spinv(j+8,spnum) = td(j)
          end do ! i
        elseif(pcomp(text,'edge',4)) then
          do j = 1,3
            spinv(j+3,spnum) = 0.0d0
          end do ! i
          j                = nint(td(1))
          spinv(j+3,spnum) = 1.0d0
          spinv(12 ,spnum) = td(2)
          espfl = .true.
        endif
      end do !

c     Output result

      if(prt) then
        write(iow,2000) spnum,(spinv(j,spnum),j=1,7),
     &                  nint(spinv(8,spnum)),
     &                 (spinv(j,spnum),j=9,8+min(ndf,3))
        if(ior.lt.0) then
          write(*,2000) spnum,(spinv(j,spnum),j=1,7),
     &                  nint(spinv(8,spnum)),
     &                 (spinv(j,spnum),j=9,8+min(ndf,3))
        endif
        if(spinv(12,spnum).ne.0.0d0) then
          write(iow,2001) spinv(12,spnum)
          if(ior.lt.0) then
            write(*,2001) spinv(12,spnum)
          endif
        endif
      endif

c     Convert spin rate from degrees to radians

      spinv(7,spnum) = spinv(7,spnum)*pi/180.0d0

c     Check for input error

      if(espfl .and. nflag) then
        write(iow,3001)
        call plstop()
      endif

c     Formats

2000  format(/5x,'S p i n   D a t a   V a l u e s'//
     &  8x,'Number',i4,' Geometric data and spin rate'/
     & 10x,'Center  : x =',1p,1e11.3,' y =',1p,1e11.3,' z =',1p,1e11.3/
     & 10x,'Normal  : x =',1p,1e11.3,' y =',1p,1e11.3,' z =',1p,1e11.3/
     & 10x,'Rate    : w =',1p,1e11.3,' degrees/time'/
     &  8x,'Translation data'/10x,'Prop Num:    ',i4/
     & 10x,'Velocity: x =',1p,1e11.3, ' y =',1p,1e11.3:,
     &              ' z =',1p,1e11.3)
2001  format(5x,'Edge Value  =',1p,1e11.3)

3000  format(5x,'**ERROR** Zero length normal vector for spin')
3001  format(5x,'**ERROR** Can not specify NORM and EDGE')

      end
