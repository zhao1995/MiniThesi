c$Id:$
      subroutine crmate (nmate,cm0,cm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read SURfaces data

c      Purpose: Input of contact surface data

c      Inputs :
c         nmate   - Number of the MATErial

c      Outputs:
c         cm0(*)  - Contact material control data
c         cm(*)   - Contact materials data storage
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_comnd.h'
      include  'c_dict.h'
      include  'iofile.h'
      include  'print.h'

      logical   whfl
      character getlab*70,label*70,type*4,cc*4,c2*4
      integer   nmate, ncom,ntype,nfeat,nopti,nsubc,nsopt,kr,k,labl
      real*8    cm0(nr0,n0c2:*),cm(*), tydat(15),td(14)

      save

      call cdebug0 ('    crmate',-1)

c     Print command title

      if (prt) then
        write(iow,2200) nmate
        label = getlab(2,labl)
        if (labl.ne.0) write (iow,2201) label
      endif

c     Store material #, offset for material data

      ncom      = 2
      cm0(1,-1) = nmate
      cm0(2,-1) = ofsmate

c     TYPE declaration

      call crtype (ncom,nmate,type,ntype,tydat)
      cm0(1,0) = ntype
      do kr = 2,c_nr0
        cm0(kr,0) = tydat(kr-1)
      end do
      if (prt) then
        write (iow,2210) type
        if (ntype.eq.1) then
          write (iow,2211)
        elseif (ntype.eq.2) then
          write (iow,2212)
        elseif (ntype.eq.3) then
          write (iow,2213)
        endif
      endif

c     Read data

      whfl = .true.
      do while (whfl)
        call crdata (ncom,cc,c2,nfeat,nopti,nsubc,nsopt,td)

c       FEATURE - deposit values in command table

        if (nfeat.gt.0) then
          cm0(1,nfeat) = nfeat
          cm0(2,nfeat) = nopti
          do k = 1,14
            cm0(k+2,nfeat) = td(k)
          end do
          if (prt) then
            write (iow,2220) cc
            if ((nfeat.eq.1)) then
c              write (iow,2221) c2
            endif
          endif

c       SUB-COMMAND - perform subcommand

        elseif (nsubc.gt.0) then
          if (prt) then
            write (iow,2230) cc
            if (nopti.gt.0) then
              write (iow,2231) c2
            endif
          endif
c         call cr... ()

c       TYPE DATA --> read type data

        elseif (nfeat.eq.-3) then
          if (prt) then
            write (iow,2240)
          endif
          backspace (ior)
          if (ntype.eq.1) then
            call crmat01 (cm(ofsmate))
          elseif (ntype.eq.2) then
            call crmat02 (cm(ofsmate))
          elseif (ntype.eq.3) then
            call cumater (cm(ofsmate))
          endif

c       Blank line found - search again

        elseif (nfeat.eq.-1) then
          continue

c       New command found - go back to PCONT

        elseif (nfeat.eq.-2) then
          backspace (ior)
          whfl = .false.
        endif
      end do

c     Compute new offset for next material set

      ofsmate   = ofsmate + c_lmv

2200  format (/5x,'C o n t a c t   M a t e r i a l   D a t a'//
     &         5x,'Data Set for Material Number      ',i4/)

2201  format (/5x,'Material Label: ',a)

2210  format ( 5x,'Contact Material Type             ',a)
2211  format (/5x,'Non-deformable Material with ',
     &            'Coulomb Friction')
2212  format (/5x,'Non-Linear Friction Material')
2213  format (/5x,'User Material Model')

2220  format (/5x,'  Material Feature                ',a)

2230  format (/5x,'  Material sub-command            ',a)
2231  format ( 5x,'    Sub-command option            ',a)

2240  format ( 5x,'  Type Declaration Data:          ')

      end
