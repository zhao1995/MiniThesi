c$Id:$
      subroutine coutcis(ncom)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_dict.h'
      include   'iofile.h'

      integer    ncom, typ,fep,opp,scp,sop, ii,jj

c     Commands

      write(iow,4000) (cis(typ(ncom,ii)),ii=1,nty(ncom))

c     Features

      do ii = 1,nfe(ncom)
        write(iow,4001) cis(fep(ncom,ii)),
     &                 (cis(opp(ncom,ii,jj)),jj=1,nop(ncom))
      end do ! ii

c     Subcommands

      do ii = 1,nsc(ncom)
        write(iow,4002) cis(scp(ncom,ii)),
     &                 (cis(sop(ncom,ii,jj)),jj=1,nso(ncom))
      end do ! ii

c     Formats

4000  format(/5x,'Commands: '/(7x,6('  ',a)))
4001  format(/7x,'Feature: ',a/
     &       /9x,'Feature Options:'/(9x,6('  ',a)))
4002  format(/11x,'Sub-Command: ',a/
     &       /13x,'Sub-Command Options:'/(13x,6('  ',a)))

      end
