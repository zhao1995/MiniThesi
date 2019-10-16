c$Id:$
      subroutine aug2 (alpha,ch2,rnorm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996           1.00
c               Robert L. Taylor         October 12, 1996           1.01

c      Acronym: AUGMENTation

c      Purpose: Augment contact forces

c      Inputs :
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_keyh.h'
      include  'augdat.h'
      include  'iofile.h'

      integer   istgn
      real*8    alpha,ch2(*),rnorm, fn,augfna,area,gn,augfn,kn

      save

      call cdebug0 ('      aug2',-1)

c     check for opening

      istgn = nint(ch2(p1(4)))
      if (istgn.ge.0) then
        fn     = ch2(p1(51))
        augfna = ch2(p1(151))
        area   = ch2(p1(11))
        gn     = ch2(p1(9))
        augg   = max(abs(gn),augg)

        augfn = augfna + (fn-augfna)*alpha
        rnorm = rnorm + gn**2

        if (ifdb) then
          if (gn.ne.0.d0) then
             kn = fn/(gn*area)
          else
             kn = 0.d0
          endif
          write (*,*) 'AUG2 augfn kn gn dfn',augfn,kn,gn,augfn-augfna
          write (iow,*) 'augfn kn gn dfn',augfn,kn,gn,augfn-augfna
        endif
      else
        augfn = 0.d0
        if (ifdb) then
          write (*,*) 'open gap'
        endif
      endif
      ch2(p1(151)) = augfn

      end
