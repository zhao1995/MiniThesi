c$Id:$
      subroutine crtype (ncom,nn,type,ntype,td)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read TYPE declaration

c      Purpose: Input of contact type

c      Inputs :
c         ncom    - internal number of command
c         nn      - user number for command

c      Outputs:
c         type    - type string
c         ntype   - # of type
c         td      - # input variables list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dict.h'
      include  'iofile.h'

      logical   errck,gettxtd,pcomp,whfl
      character type*4, tx(2)*15
      integer   ncom,nn,ntype, typ,ke, ii
      real*8    td(15)

      save

      call cdebug0 ('    crtype',-1)

      errck = gettxtd(tx(1),1,td,15,'skip')
      type  = tx(1)(1:4)

c     Set type comparing with dictionary

      whfl = .true.
      ke   = 0
      do while (whfl)
        ke = ke + 1
        if (ke.gt.nty(ncom)) then
          write (ilg,3001) type,cwd(ncom),nn,
     &          (cis(typ(ncom,ii)),ii=1,nty(ncom))
          write (iow,3001) type,cwd(ncom),nn,
     &          (cis(typ(ncom,ii)),ii=1,nty(ncom))
          write (  *,3001) type,cwd(ncom),nn,
     &          (cis(typ(ncom,ii)),ii=1,nty(ncom))
          call plstop()
        elseif (pcomp(type,cis(typ(ncom,ke)),4))then
          ntype = ke
          whfl  = .false.
        endif
      end do

c     Formats

3001  format(/' *ERROR* CRTYPE: Illegal type found : ',a4/
     &        '         for the contact command    : ',a4,i5/
     &        '         OPTIONS: ',8(a:,', ')/(18x,8(a:,', ')))

      end
