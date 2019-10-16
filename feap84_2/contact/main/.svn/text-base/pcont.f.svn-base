c$Id:$
       subroutine pcont (cs0,cm0,cp0,c90,c100,c110,cm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Assign length to cc = 4 characters, c2 = 8       17/04/2007
c       2. Check if np(140) exists before deleting CTEM5    15/01/2008
c       3. Move call to ptranf after the check for [tran]   04/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Equivalent to PMESH

c      Purpose: Data input routine for contact description

c      Inputs:
c         cs0(*)  - Contact surfaces control data
c         cm0(*)  - Contact material control data
c         cm(*)   - Contact materials data storage
c         cp0(*)  - Contact pair control data

c      Outputs:
c         ICS     - Stored within np(133) pointer definition, length
c                   depends on commands specified
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_comnd.h'
      include  'c_dict.h'
      include  'cdata.h'
      include  'chdata.h'
      include  'iofile.h'
      include  'iosave.h'
      include  'print.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp,gettxtd,errck,readfl,savefl, setvar,palloc
      character tx(2)*15,cc*4,c2*8
      integer   lics0,lics1, ncom,nn,nd,k, offs1
      real*8    cs0(nr0,n0c1:nc01,*) ,cm0(nr0,n0c2:nc02,*)
      real*8    cp0(nr0,n0c3:nc03,*) ,c90(nr0,n0c9:nc09,*), td(14)
      real*8    c100(nr0,n0c10:nc010,*),c110(nr0,n0c11:nc011,*),cm(*)

      save

      call cdebug0 ('  pcont',-1)

c     clear command array

      do k = 1, c_ncc
        cck(k) = 0
      end do ! k

c     Print main title

      if (prt) then
        call prtitl(prth)
        write (iow,2000)
      endif
      write (ilg,2001)

c     Scan data
c     WARNING - data check for command number performed in PNUMC
c     WARNING - backspace performed to use always gettxtd to allow
c               parameters use everywhere

      errck = gettxtd(tx,2,td,0,'skip')
      cc = tx(1)(1:4)
      c2 = tx(2)(1:8)

      do while (.not.pcomp(cc,'end',3))

c       [surf]ace - contact surface geometry data input

        if (pcomp(cc,cwd(1),4)) then
          write(ilg,4001) cc,c2
          ncom = 1
          cck(ncom) = cck(ncom) + 1
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

c         Input surface data

          lics1  = 0
          offs1  = ofssurf - 1
          call crsurf (nn,cs0(1,n0c1,cck(ncom)),lics1)
          if(lics1.gt.0) then
            lics0 = lics0 + lics1
            setvar = palloc(133,'ICS  ',max(1,lics0),1)
            do k = 0,lics1-1
              mr(np(133)+k+offs1) = mr(np(140)+k)
            end do ! k
          endif
          if(np(140).ne.0) then
            setvar = palloc(140,'CTEM5',0,1)  ! Allocated in crsurf
          endif

c       [mate]rial,ma: Data input for contact material set ma

        elseif (pcomp(cc,cwd(2),4)) then
          write(ilg,4001) cc,c2
          ncom = 2
          cck(ncom) = cck(ncom) + 1
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

c         input material data

          call crmate (nn,cm0(1,n0c2,cck(ncom)),cm)

c       [pair]

        elseif (pcomp(cc,cwd(3),4)) then
          write(ilg,4001) cc,c2
          ncom = 3
          cck(ncom) = cck(ncom) + 1
          backspace (ior)
          errck = gettxtd(tx,1,td,2,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))
          nd = nint(td(2))
          nd = min(nd,ndm)

c         Input record for surface type

          call crpair (nn,nd,cp0(1,n0c3,cck(ncom)))

c       [auto] auto surface data generation

        elseif (pcomp(cc,cwd(4),4)) then
          write(ilg,4001) cc,c2
          ncom = 4
          cck(ncom) = cck(ncom) + 1
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

c         Generate the boundary slidelines

          setvar = palloc(140,'CTEM5',6*numnp,1)
          lics1  = 0
          nd     = max(nn,1)
          offs1  = ofssurf - 1
          call crauto (nd,cs0(1,n0c1,cck(1)+1),mr(np(140)),lics1)
          lics0 = lics0 + lics1
          setvar = palloc(133,'ICS  ',max(1,lics0),1)
          do k = 0,lics1-1
            mr(np(133)+k+offs1) = mr(np(140)+k)
          end do ! k
          setvar = palloc(140,'CTEM5',0,1)

c         Increment the number of slide lines by 'nd'

          cck(1) = cck(1) + nd

c       [read] data from a new file

        elseif (pcomp(cc,cwd(5),4)) then
          write(ilg,4001) cc,c2
          errck = readfl(c2)

c       [save] data in a file

        elseif (pcomp(cc,cwd(6),4)) then
          write(ilg,4001) cc,c2
          errck = savefl(c2)

c       [tran]sform

        elseif (pcomp(cc,cwd(7),4)) then
          write(ilg,4001) cc,c2
          call ptranf(c2,prt)

c       [help]

        elseif (pcomp(cc,cwd(8),4)) then
cc        call phelp (c2,wd,ed,10,'PCONT')
          continue

c     cc      - character*4 name of command for help
c     wd      - character*4 name for list of commands
c     ed      - integer for manual levels (0-3, 0=basic -> 3 = all)
c     nwd     - number of words in list
c     subname - character*(*) calling subprogram name (in quotes)

c       [unu1] user commands

        elseif (pcomp(cc,cwd(9),4)) then
          write(ilg,4001) cc,c2
          ncom = 9
          cck(ncom) = cck(ncom) + 1
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))
          call crunu1 (nn,c90(1,n0c9,cck(ncom)))

c       [unu2] user commands

        elseif (pcomp(cc,cwd(10),4)) then
          write(ilg,4001) cc,c2
          ncom = 10
          cck(ncom) = cck(ncom) + 1
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))
          call crunu2 (nn,c100(1,n0c10,cck(ncom)))

c       [unu3] user commands

        elseif (pcomp(cc,cwd(11),4)) then
          write(ilg,4001) cc,c2
          ncom = 11
          cck(ncom) = cck(ncom) + 1
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))
          call crunu3 (nn,c110(1,n0c11,cck(ncom)))
        endif

        errck = gettxtd(tx,2,td,0,'skip')
        cc = tx(1)(1:4)
        c2 = tx(2)(1:8)
      end do ! while

c     [end] of contact data inputs
c     Lsave check performed in PNUMC

c     Formats

2000  format (/'   C O N T A C T   I N P U T   D A T A'/)
2001  format (/'  CONTACT INPUT DATA'/'  ------------------')

4001  format (5x,a,' ',a)

      end
