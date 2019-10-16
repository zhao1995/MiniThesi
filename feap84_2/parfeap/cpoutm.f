c$Id:$
      subroutine cpoutm (cs0,cm0,cp0,ics,hic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    26/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Control for output of parallel contact data to file

c     Inputs:
c       cs0(*)    - Surface  table
c       cm0(*)    - Material table
c       cp0(*)    - Pair     table
c       ics(*)    - Contact facet data
c       hic(*)    - History correspondence table

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'      ! ios
      include   'setups.h'      ! ndomn
      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_contac.h'
      include   'c_geom.h'      ! nsurf1, nsurf2

      include   'pointer.h'
      include   'comblk.h'

      logical    dflags, dflagm, allfls, allflm, first
      character  type(9)*8
      integer    n
      integer    ofss, nepss,dnopes,nopes, ntypes
      integer    ofsm, nepsm,dnopem,nopem, ntypem
      integer    ics(*), hic((c_lp1+c_lp3),*)
      real*8     cs0(nr0,n0c1:nc01,*), cm0(nr0,n0c2:nc02,*)
      real*8     cp0(nr0,n0c3:nc03,*)

      save

      data       type / 'LINE    ','TRIANGLE','QUAD    ','BEAM    ',
     &                  'POINT   ','RIGID   ','NURBS   ','TSPLINE ',
     &                  'PART    ' /

c     Pair outputs

      first = .true.
      do n = 1, numcp

        call setcomp(n,cs0,cm0,cp0,hic)

c       Slave surface pointers

        ofss   = nint(abs(cs0(2,-1,nsurf1)))
        nepss  = nint(abs(cs0(3,-1,nsurf1)))
        dnopes = nint(abs(cs0(4,-1,nsurf1)))
        nopes  = nint(abs(cs0(2, 0,nsurf1)))
        ntypes = nint(abs(cs0(1, 0,nsurf1)))

c       Master surface pointers

        ofsm   = nint(abs(cs0(2,-1,nsurf2)))
        nepsm  = nint(abs(cs0(3,-1,nsurf2)))
        dnopem = nint(abs(cs0(4,-1,nsurf2)))
        nopem  = nint(abs(cs0(2, 0,nsurf2)))
        ntypem = nint(abs(cs0(1, 0,nsurf2)))

c       Check Slave surface for output limits

        call cpoutmc(ics(ofss),dnopes,nopes,nepss,
     &               ndomn, mr(np(254)),mr(np(111)),mr(np(112)),
     &               dflags, allfls )

c       Check Master surface for output limits

        call cpoutmc(ics(ofsm),dnopem,nopem,nepsm,
     &               ndomn, mr(np(254)),mr(np(111)),mr(np(112)),
     &               dflagm, allflm )  ! XX Debug changed to m

c       Set all flags for current partition

        allfls = allfls .or. dflagm
        allflm = allflm .or. dflags

c       Slave surface output

        if(dflags .or. allfls) then

c         Output beginning of contact data

          if(first) then
            write(ios,'(/a)') 'CONTACT'
            first = .false.
          endif

          call cpoutms(nsurf1,type(ntypes),ics(ofss),dnopes,nopes,nepss,
     &                 ndomn, mr(np(254)),mr(np(111)),mr(np(112)),
     &                 dflags, allfls )
        endif

c       Master surface output

        if(dflagm .or. allflm) then

          call cpoutms(nsurf2,type(ntypem),ics(ofsm),dnopem,nopem,nepsm,
     &                 ndomn, mr(np(254)),mr(np(111)),mr(np(112)),
     &                 dflagm, allflm )

c         Pair output

          call coutmp(n,cp0)

        endif
      end do ! n

c     Output end of contact data

      if(.not.first) then
        write(ios,'(/a)') 'END CONTACT'
      endif

      end
