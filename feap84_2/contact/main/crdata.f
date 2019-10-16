c$Id:$
      subroutine crdata (ncom,cc,c2,nfeat,nopti,nsubc,nsopt,td)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Assign length of cc and c2 to 4 characters       17/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read TYPE declaration

c      Purpose: Input of contact type

c      Inputs :
c         ncom    - # of command

c      Outputs:
c         cc      - first input string  (feature or sub-command)
c         c2      - second input string (feat. or sub-comm. option)
c         string  - text input string
c         nfeat   - # of FEATure
c         nopti   - # of OPTIon
c         nsubc   - # of SUB-Command
c         nsopt   - # of Sub-command OPTion
c         td      - # input variables list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dict.h'
      include  'chdata.h'
      include  'iofile.h'

      logical   errck,gettxtd,pcomp,whfl,whf2,whf3
      character cc*4,c2*4, tx(2)*15
      integer   ncom,nfeat,nopti,nsubc,nsopt
      integer   fep,opp,scp,sop,ke,ks,ko,kc, ii
      real*8    td(14)

      save

      call cdebug0 ('    crdata',-1)

c     Set defaults

      nfeat = 0
      nopti = 0
      nsubc = 0
      nsopt = 0
      do ke = 1,14
        td(ke) = 0.d0
      end do

c     Get data string (absolutely in this way because of command labels)

      errck = gettxtd(tx,2,td,0,'skip')
      cc = tx(1)(1:4)
      c2 = tx(2)(1:4)

c     Blank line found

      if (pcomp(cc,'    ',4)) then
        nfeat = -1

c     Decode string looking for feature and option; sub-command and
c     sub-command option; or new command comparing with dictionary.

      else
        whfl = .true.
        ke = 0
        do while (whfl)
          ke = ke+1
          if (ke.gt.nfe(ncom)) then

c           Feature not found, check if it is sub-command

            ks = 0
            whf2 = .true.
            do while (whf2)
              ks = ks+1
              if (ks.gt.nsc(ncom)) then

c               Sub-command not found, check if it is new command

                kc = 0
                whf3 = .true.
                do while (whf3)
                  kc = kc+1

c                 New command not found: should be type declaration

                  if (kc.gt.c_ncc) then
                    nfeat = -3
                    whf3  = .false.
                    whf2  = .false.
                    whfl  = .false.

c                 New command found

                  elseif (pcomp(cc,cwd(kc),4))then
                    nfeat = -2
                    whf3  = .false.
                    whf2  = .false.
                    whfl  = .false.
                  endif
                end do

c             Sub-command found

              elseif (pcomp(cc,cis(scp(ncom,ks)),4))then
                nsubc = ks
                whf2  = .false.
                whfl  = .false.

c               Check sub-command option

                if (.not.pcomp(c2,'    ',4)) then
                  whf3 = .true.
                  ko   = 0
                  do while (whf3)
                    ko = ko+1

c                   Sub-command option not found  --> input error

                    if (ko.gt.nso(ncom)) then
                      write (ilg,3002) c2,cc,cwd(ncom),
     &                      (cis(sop(ncom,nsubc,ii)),ii=1,nso(ncom))
                      write (iow,3002) c2,cc,cwd(ncom),
     &                      (cis(sop(ncom,nsubc,ii)),ii=1,nso(ncom))
                      call plstop()

c                   Sub-command option found

                    elseif (pcomp(c2,cis(sop(ncom,nsubc,ko)),4))then
                      nsopt = ko
                      whf3  = .false.
                    endif
                  end do
                endif
              endif
            end do

c         Feature found

          elseif (pcomp(cc,cis(fep(ncom,ke)),4))then
            nfeat = ke
            whfl  = .false.

c           Check feature option

            if (.not.pcomp(c2,'    ',4)) then
              ko = 0
              whf3 = .true.
              do while (whf3)
                ko = ko+1

c               Option not found  --> input error

                if (ko.gt.nop(ncom)) then
                  write (ilg,3003) c2,cc,cwd(ncom),
     &                      (cis(opp(ncom,nfeat,ii)),ii=1,nop(ncom))
                  write (iow,3003) c2,cc,cwd(ncom),
     &                      (cis(opp(ncom,nfeat,ii)),ii=1,nop(ncom))
                  call plstop()

c               Feature option found

                elseif (pcomp(c2,cis(opp(ncom,nfeat,ko)),4))then
                  nopti = ko
                  whf3 = .false.
                endif
              end do
            endif
          endif
        end do
      endif

c     If feature or sub-command found then get numerical data

      if ((nfeat.gt.0) .or. (nsubc.gt.0)) then
        backspace (ior)
        errck = gettxtd(tx,2,td,14,'skip')
      endif

c     Formats

3002  format(' *ERROR* CRDATA:'/
     &       '          Illegal Sub-command Option     : ',a/
     &       '          for the sub-command            : ',a/
     &       '          of the contact command         : ',a/
     &       '         OPTIONS: ',8(a:,', ')/(18x,8(a:,', ')))

3003  format(' *ERROR* CRDATA:'/
     &       '          Illegal Feature Option         : ',a/
     &       '          for the Feature                : ',a/
     &       '          of the contact command         : ',a/
     &       '         OPTIONS: ',8(a:,', ')/(18x,8(a:,', ')))

      end
