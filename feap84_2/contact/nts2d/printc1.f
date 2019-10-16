c$Id:$
      subroutine printc1 (npair,nel1,nod1,ix1,ix2,ch2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: PRint Node To Segment status

c      Purpose: Print information for the basic NTS conatc

c      Inputs :
c         npair   - # of current pair
c         nel1    - Current contact element of surf. 1
c         nod1    - Current node of contact element nel1
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c                 - On the listing file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_keyh.h'
      include  'c_geom.h'
      include  'c_pair.h'
      include  'iofile.h'

      integer   npair,nel1,nod1,ix1(dnope1,*),ix2(dnope2,*)
      integer   masts,istgt,istgn,istfr,ns,n1,n2
      real*8    ch2(*), s21,c21,d21,csi,gn,area,dgt,fn,dgnfn,presn
      real*8    ft,dgtft,dgnft,dgte,dgtp,prestg,totf,totp,ftx,fty
      real*8    fnx,fny,s21c,c21c,totfx,totfy,costf,sintf

      save

      call cdebug0 ('      printc1',-1)

c     Get data for geometry

      masts = nint(ch2(p1(1)))
      istgt = nint(ch2(p1(3)))
      istgn = nint(ch2(p1(4)))
      s21   = ch2(p1(5))
      c21   = ch2(p1(6))
      d21   = ch2(p1(7))
      csi   = ch2(p1(8))
      gn    = ch2(p1(9))
      area  = ch2(p1(11))
      dgt   = ch2(p1(12))

      ns    = ix1(nod1,nel1)
      n1    = ix2(1,masts)
      n2    = ix2(2,masts)

c     Get data for normal contact

      fn    = -ch2(p1(51))
      dgnfn =  ch2(p1(52))

c     Compute data for normal contact

      presn = fn/area

c     Get data for friction

      if (iffric.eq.1) then
        ft    = -ch2(p1(53))
        dgtft =  ch2(p1(54))
        dgnft =  ch2(p1(55))
        dgte  =  ch2(p1(56))
        dgtp  =  ch2(p1(57))
        istfr =  nint(ch2(p1(58)))

c       Compute data for friction

        prestg =  ft/area
        totf   =  sqrt(fn**2 + ft**2)
        totp   =  sqrt(presn**2+prestg**2)
        ftx    =  ft * c21
        fty    =  ft * s21

        if (istgn.eq.1) then
          fnx  =  fn * s21
          fny  = -fn * c21
        else
          s21c =  ch2(p1(13))
          c21c =  ch2(p1(14))
          fnx  =  fn * s21c
          fny  = -fn * c21c
        endif

        totfx  =  ftx + fnx
        totfy  =  fty + fny
        if (totf.ne.0.d0) then
          costf =  totfx/totf
          sintf =  totfy/totf
        else
          costf = 0.d0
          sintf = 0.d0
         endif
      endif

c     Printout

      write (iow,2000) nel1,rnpair,npair
      write (iow,2001) nod1,masts,ns,n1,n2,istgt,istgn
      write (iow,2002)
      write (iow,2003) s21,c21,d21,csi,gn,area,dgt
      write (iow,2004)
      write (iow,2003) fn,presn,dgnfn
      if (ior.lt.0) then
        write (*,2000) nel1,rnpair,npair
        write (*,2001) nod1,masts,ns,n1,n2,istgt,istgn
        write (*,2002)
        write (*,2003) s21,c21,d21,csi,gn,area,dgt
        write (*,2004)
        write (*,2003) fn,presn,dgnfn
      endif

      if (iffric.eq.1) then
        write (iow,2005)
        write (iow,2006) ft,prestg,dgtft,dgnft,dgte,dgtp,istfr
        write (iow,2007)
        write (iow,2003) totf,totp,costf,sintf
        if (ior.lt.0) then
          write (*,2005)
          write (*,2006) ft,prestg,dgtft,dgnft,dgte,dgtp,istfr
          write (*,2007)
          write (*,2003) totf,totp,costf,sintf
        endif
      endif

2000  format (//
     &         5x,'contact element # ',i6,16x,'pair # ',i6,
     &            ' internal pair # ',i6/
     &         4x,' local node',' mast  elem',' node slave',
     &            ' nod mast 1',' nod mast 2','  status gt',
     &            '  status gn')
2001  format (4x,7i11)
2002  format (/4x,'      sin t','      cos t',' mast lengt',
     &            '        csi','         gn','  cont area',
     &            '    sliding')
2003  format (4x,7e11.3)
2004  format (/4x,'         fn',' norm press','      dgnfn')
2005  format (/4x,'         ft',' tang press','      dgtft',
     &            '      dgnft','       dgte','       dgtp',
     &            '  stat fric')
2006  format (4x,6e11.3,i11)
2007  format (/4x,'  tot force','  tot press',' orient cos',
     &            ' orient sin' )

      end
