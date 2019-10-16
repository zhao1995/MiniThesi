c$Id:$
      subroutine pedgin()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'curv'e inputs for boundary conditions       13/11/2008
c       2. Add load table option for eforc and edisp        10/01/2009
c       3. Change 'np(27)' to 'point' on peforc call loads  04/06/2010
c       4. Use setext to assign file extenders & set number 20/12/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control routine of data inputs based on edge coordinate

c      Inputs:
c         none      - Data retrieved through common blocks

c      Outputs:
c         none      - Data stored in pointers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'cdata.h'
      include  'conval.h'
      include  'edgdat.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'pload1.h'
      include  'pointer.h'
      include  'print.h'
      include  'p_point.h'
      include  'sdata.h'
      include  'trdata.h'
      include  'comblk.h'

      character fext*8, type*4
      logical   oprt,oprth
      integer   l1,isd

      save

      data      isd  / 16 /

c     Set edge angle values

      oprt  = prt
      oprth = prth
      if(eanfl) then                  ! Set for 'eang'
        do l1 = 0,neang-1
          call setext('eang', l1, fext, .false.)
          call pinpfl('PEDGIN',fext, type, 1)
          call peforc(hr(np(43)),hr(np(45)),mr(np(190)),
     &                ndm,1,numnp,type,prt,prth,'Angle')
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge boundary conditions

      if(ebcfl) then                  ! Set for 'ebou'
        do l1 = 0,nebcs-1
          call setext('ebou', l1, fext, .false.)
          call pinpfl('PEDGIN',fext, type, 1)
          call pedges(hr(np(43)),mr(np(31)+ndf*numnp),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'B.C.')
          call pinpfl('PEDGIN',fext, type, 2)

        end do ! l1
      endif

c     Set edge displacement values

      if(edifl) then                  ! Set for 'edis'
        do l1 = 0,nedis-1
          call setext('edis', l1, fext, .false.)
          call pinpfl('PEDGIN',fext, type, 1)
          if(ldflg) then
            if(np(266).ne.0) then
              call pldexp(2,mr(np(265)),mr(np(266)),hr(np(267)),
     &                    hr(np(26)))
            else
              call pzero(hr(np(26)),nneq)
            endif
            point = np(26)
          else
            point = np(27)+nneq
          endif
          call peforc(hr(np(43)),hr(point),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Displ')
          if(ldflg) then
            call pldseta(mr(np(265)),hr(point), 2, ndf,numnp)
          endif
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge force values

      if(efcfl) then                  ! Set for 'efor'
        do l1 = 0,nefrc-1
          call setext('efor', l1, fext, .false.)
          call pinpfl('PEDGIN',fext, type, 1)
          if(ldflg) then
            if(np(266).ne.0) then
              call pldexp(1,mr(np(265)),mr(np(266)),hr(np(267)),
     &                    hr(np(26)))
            else
              call pzero(hr(np(26)),nneq)
            endif
            point = np(26)
          else
            point = np(27)
          endif
          call peforc(hr(np(43)),hr(point),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Force')
          if(ldflg) then
              call pldseta(mr(np(265)),hr(point), 1, ndf,numnp)
          endif
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge proportional load numbers

      if(eprfl) then                  ! Set for 'epro'
        do l1 = 0,nepro-1
          call setext('epro', l1, fext, .false.)
          call pinpfl('PEDGIN',fext, type, 1)
          type = 'set'
          call pedges(hr(np(43)),mr(np(29)),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Prop')
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge base conditions

      if(ebsfl) then                  ! Set for 'ebas'
        do l1 = 0,nebas-1
          call setext('ebas', l1, fext, .false.)
          call pinpfl('PEDGIN',fext, type, 1)
          call pedges(hr(np(43)),mr(np(125)),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Base')
          call pinpfl('PEDGIN',fext, type, 2)

        end do ! l1
      endif

c     Set curve boundary inputs

      if(curfl) then                  ! Set for 'curv'
        do l1 = 0,ncurv-1
          call setext('curv', l1, fext, .false.)
          call pinpfl('PEDGIN',fext, type, 1)
c                        x      ,   id          ,   xs
          call pcurve(hr(np(43)),mr(np(31)+nneq),hr(np(161)),
     &                mr(np(162)),numsd,isd,ndm,ndf,numnp)
c                        is
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1

      endif

      prt  = oprt
      prth = oprth

      end
