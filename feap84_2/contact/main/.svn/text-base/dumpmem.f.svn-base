c$Id:$
      subroutine dumpmem (debf,nn,vect,cs0,ics,cm0,cm,cp0,
     &                    ch1,ch2,ch3,lp1,lp3,hic,screen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Modify print formats to fit screen               09/07/2010
c       2. Initialize 'c' on for all calls                  05/10/2011
c       3. Check for rank before inputs from *              14/03/2013
c       4. Align formats 2018 & 2019 better                 15/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: DUMP of MEMory

c      Purpose: Write contents of contact memory arrays

c      Inputs :
c         debf    - DEBug File
c         nn      - # data set to dump
c         vect(*) - string command for dump
c         cs0(*)  - Contact surfaces control data
c         ics(*)  - Contact element nodal connection array
c         cm0(*)  - Contact material control data
c         cm(*)   - Contact materials data storage
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         lp1     - Length of p1 (static)
c         lp3     - Length of p3 (static)
c         hic(*)  - HIstory Correspondence table
c         screen  - flag for screen output

c      Outputs:
c                 - on file fort.90
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'iofile.h'
      include  'print.h'
      include  'setups.h'

      logical   screen, whfl,pcomp,prnam
      character vect*(*)
      character ctemp*10,csys*10,ctyp*10,cfea(10)*10,cusr*10, c*1
      integer   debf,nn,ics(*),hic(*),lp1,lp3, fep,typ
      integer   nf,nl,ks,kr,kc,km,kp,ke,kn,i0,kd,ofs,neps,nope
      integer   ofh1,ofh3,lh1,lh3,ns1,ns2,ndv,i1,i3,ih1,ih3,kv
      integer   nset,kp1,kp3,rnsurf,rnmat,rnpair,kk,dnope
      real*8    cs0(nr0,n0c1:nc01,*),cm0(nr0,n0c2:nc02,*),cm(*)
      real*8    cp0(nr0,n0c3:nc03,*),ch1(*),ch2(*),ch3(*)

      save

      data      csys,cusr /'    system','      user'/

c     Set c to blank for all inputs

      c  = ' '
c     Dump of SURF command table

      if ((pcomp(vect,'surf',4)) .or. (pcomp(vect,'C0',2))) then

c       Check if dump is limited to a specific surf

        if (nn.eq.0) then
          nf = 1
          nl = numcs

c       Search set

        else
          do kc = 1, cck(1)
            if (cs0(1,-1,kc).eq.nn) then
              nf = kc
              nl = kc
            endif
          end do ! kc
        endif

c       Perform dump

        c  = ' '
        ks = nf
        do while (ks.le.nl)
          rnsurf = nint(cs0(1,-1,ks))
          ctyp       = ' '
          ctyp(7:10) = cis(typ(1,nint(cs0(1,0,ks))))
          do kc = 1,nfe(1)
            ctemp       = ' '
            ctemp(7:10) = cis(fep(1,kc))
            cfea(kc)    = ctemp
          end do ! kc
          if (ifdb) then
            write (debf,2000) rnsurf,ks
            write (debf,2018) (csys,kc=1,abs(n0c(1))),ctyp,
     &                        (cfea(kc),kc=1,nfe(1)),(cusr,kc=1,nuc(1))
            write (debf,2019) (kc,kc= n0c(1),nc0(1))
          endif
          if (screen) then
            write (iow,2000) rnsurf,ks
            write (iow,2018) (csys,kc=1,abs(n0c(1))),ctyp,
     &                       (cfea(kc),kc=1,nfe(1)),(cusr,kc=1,nuc(1))
            write (iow,2019) (kc,kc= n0c(1),nc0(1))
          endif
          if ((ior.lt.0) .and. (screen)) then
            write (*,2000) rnsurf,ks
            write (*,2018) (csys,kc=1,abs(n0c(1))),ctyp,
     &                     (cfea(kc),kc=1,nfe(1)),(cusr,kc=1,nuc(1))
            write (*,2019) (kc,kc= n0c(1),nc0(1))
          endif
          do kr = 1,nr0
            if (ifdb) then
              write (debf,2001) (cs0(kr,kc,ks),kc = n0c(1),nc0(1))
            endif
            if (screen) then
              write (iow,2001) (cs0(kr,kc,ks),kc = n0c(1),nc0(1))
            endif
            if ((ior.lt.0) .and. (screen)) then
              write (*,2001) (cs0(kr,kc,ks),kc = n0c(1),nc0(1))
            endif
          end do ! kr
          if ((ior.lt.0) .and. (screen) .and.
     &        (.not.pcomp(c,'r',1))          ) then
4001        call pprint ('  Enter = next, r = run to last,')
            call pprint(' g # = go to #, s = stop ===> ')
            if(rank.eq.0) read (*,2017,err=4001) c,kk
            if (pcomp(c,'s',1)) then
               screen = .false.
               ifdb   = .false.
               ks = nl+1
            elseif (pcomp(c,'g',1)) then
              ks = max(1,kk)
            else
              ks = ks+1
            endif
          else
            ks = ks+1
          endif
        end do ! while
      endif


c     Dump of MATE command table

      if ((pcomp(vect,'mate',4)) .or. (pcomp(vect,'C0',2))) then

c       Check if dump is limited to a specific mate

        if (nn.eq.0) then
          nf = 1
          nl = numcm

c         Search set

        else
          do kc = 1, cck(2)
            if (cm0(1,-1,kc).eq.nn) then
              nf = kc
              nl = kc
            endif
          end do ! kc
        endif

c       Perform dump

        km = nf
        do while (km.eq.nl)
          rnmat = nint(cm0(1,-1,km))
          ctyp       = ' '
          ctyp(7:10) = cis(typ(2,nint(cm0(1,0,km))))
          do kc = 1,nfe(2)
            ctemp       = ' '
            ctemp(7:10) = cis(fep(2,kc))
            cfea(kc)    = ctemp
          end do ! kc
          if (ifdb) then
            write (debf,2002) rnmat,km
            write (debf,2018) (csys,kc=1,abs(n0c(2))),ctyp,
     &                        (cfea(kc),kc=1,nfe(2)),(cusr,kc=1,nuc(2))
            write (debf,2019) (kc,kc= n0c(2),nc0(2))
          endif
          if (screen) then
            write (iow,2002) rnmat,km
            write (iow,2018) (csys,kc=1,abs(n0c(2))),ctyp,
     &                       (cfea(kc),kc=1,nfe(2)),(cusr,kc=1,nuc(2))
            write (iow,2019) (kc,kc= n0c(2),nc0(2))
          endif
          if ((ior.lt.0) .and. (screen)) then
            write (*,2002) rnmat,km
            write (*,2018) (csys,kc=1,abs(n0c(2))),ctyp,
     &                     (cfea(kc),kc=1,nfe(2)),(cusr,kc=1,nuc(2))
            write (*,2019) (kc,kc= n0c(2),nc0(2))
          endif
          do kr = 1,nr0
            if (ifdb) then
              write (debf,2001) (cm0(kr,kc,km),kc = n0c(2),nc0(2))
            endif
            if (screen) then
              write (iow,2001) (cm0(kr,kc,km),kc = n0c(2),nc0(2))
            endif
            if ((ior.lt.0) .and. (screen)) then
              write (*,2001) (cm0(kr,kc,km),kc = n0c(2),nc0(2))
            endif
          end do ! kr
          if ((ior.lt.0) .and. (screen) .and.
     &        (.not.pcomp(c,'r',1))          ) then
4002        call pprint ('  Enter = next, r = run to last,')
            call pprint(' g # = go to #, s = stop ===> ')
            if(rank.eq.0) read (*,2017,err=4002) c,kk
            if (pcomp(c,'s',1)) then
               km = nl+1
               screen = .false.
               ifdb   = .false.
            elseif (pcomp(c,'g',1)) then
              km = max(1,kk)
            else
              km = km+1
            endif
          else
            km = km+1
          endif
        end do ! while
      endif


c     Dump of PAIR command table

      if ((pcomp(vect,'pair',4)) .or. (pcomp(vect,'C0',2))) then

c       Check if dump is limited to a specific pair

        if (nn.eq.0) then
          nf = 1
          nl = numcp

c         Search set

        else
          do kc = 1, cck(3)
            if (cp0(1,-1,kc).eq.nn) then
              nf = kc
              nl = kc
            endif
          end do ! kc
        endif

c       Perform dump

        kp = nf
        do while (kp.le.nl)
          rnpair = nint(cp0(1,-1,kp))
          ctyp       = ' '
          ctyp(7:10) = cis(typ(3,nint(cp0(1,0,kp))))
          do kc = 1,nfe(3)
            ctemp       = ' '
            ctemp(7:10) = cis(fep(3,kc))
            cfea(kc)    = ctemp
          end do ! kc
          if (ifdb) then
            write (debf,2003) rnpair,kp
            write (debf,2018) (csys,kc=1,abs(n0c(3))),ctyp,
     &                        (cfea(kc),kc=1,nfe(3)),(cusr,kc=1,nuc(3))
            write (debf,2019) (kc,kc= n0c(3),nc0(3))
          endif
          if (screen) then
            write (iow,2003) rnpair,kp
            write (iow,2018) (csys,kc=1,abs(n0c(3))),ctyp,
     &                       (cfea(kc),kc=1,nfe(3)),(cusr,kc=1,nuc(3))
            write (iow,2019) (kc,kc= n0c(3),nc0(3))
          endif
          if ((ior.lt.0) .and. (screen)) then
            write (*,2003) rnpair,kp
            write (*,2018) (csys,kc=1,abs(n0c(3))),ctyp,
     &                     (cfea(kc),kc=1,nfe(3)),(cusr,kc=1,nuc(3))
            write (*,2019) (kc,kc= n0c(3),nc0(3))
          endif
          do kr = 1,nr0
            if (ifdb) then
              write (debf,2004) kr,(cp0(kr,kc,kp),kc = n0c(3),nc0(3))
            endif
            if (screen) then
              write (iow,2004) kr,(cp0(kr,kc,kp),kc = n0c(3),nc0(3))
            endif
            if ((ior.lt.0) .and. (screen)) then
              write (*,2004) kr,(cp0(kr,kc,kp),kc = n0c(3),nc0(3))
            endif
          end do ! kr
          if ((ior.lt.0) .and. (screen) .and.
     &        (.not.pcomp(c,'r',1))          ) then
4003        call pprint ('  Enter = next, r = run to last,')
            call pprint(' g # = go to #, s = stop ===> ')
            if(rank.eq.0) read (*,2017,err=4003) c,kk
            if (pcomp(c,'s',1)) then
               kp = nl+1
               screen = .false.
               ifdb   = .false.
            elseif (pcomp(c,'g',1)) then
              kp = max(1,kk)
            else
              kp = kp+1
            endif
          else
            kp = kp+1
          endif
        end do ! while
      endif

c     Dump of material data (CM)

      if (pcomp(vect,'CM',2)) then

c       Check if dump is limited to a specific mate

        if (nn.eq.0) then
          nf = 1
          nl = numcm

c         Search set

        else
          do kc = 1, cck(2)
            if (cm0(1,-1,kc).eq.nn) then
              nf = kc
              nl = kc
            endif
          end do ! kc
        endif

c       Perform dump

        km = 1
        do while (km.le.numcm)
          rnmat = nint(cm0(1,-1,km))

c         Find information

          ofs  = nint(cm0(2,-1,km))

c         Perform minimal printout

          if (ifdb) then
            write (debf,2007) rnmat,km
          endif
          if (screen) then
            write (iow,2007) rnmat,km
          endif
          if ((ior.lt.0) .and. (screen)) then
            write (*,2007) rnmat,km
          endif

c         Perform dump of material data if requested

          if ((km.ge.nf) .and. (km.le.nl)) then
            i0 = ofs - 1
            if (ifdb) then
              write (debf,2008) (cm(i0+kd),kd = 1,c_lmv)
            endif
            if (screen) then
              write (iow,2008) (cm(i0+kd),kd = 1,c_lmv)
            endif
            if ((ior.lt.0) .and. (screen)) then
              write (*,2008) (cm(i0+kd),kd = 1,c_lmv)
            endif
          endif
          if ((ior.lt.0) .and. (screen) .and.
     &        (.not.pcomp(c,'r',1))          ) then
4004        call pprint ('  Enter = next, r = run to last,')
            call pprint(' g # = go to #, s = stop ===> ')
            if(rank.eq.0) read (*,2017,err=4004) c,kk
            if (pcomp(c,'s',1)) then
               km = numcm+1
            elseif (pcomp(c,'g',1)) then
              km = max(1,kk)
            else
              km = km+1
            endif
          else
            km = km+1
          endif
        end do ! while
      endif

c     Dump of element node (ICS)

      if (pcomp(vect,'ICS',3)) then

c       Check if dump is limited to a specific surf

        if (nn.eq.0) then
          nf = 1
          nl = numcs
        else
          do kc = 1, cck(1)
            if (cs0(1,-1,kc).eq.nn) then
              nf = kc
              nl = kc
            endif
          end do ! kc
        endif

c       Perform dump

        ks = 1
        do while (ks.le.numcs)
          rnsurf = nint(cs0(1,-1,ks))

c         Find information

          ofs  = nint(cs0(2,-1,ks))
          neps = nint(cs0(3,-1,ks))
          dnope= nint(cs0(4,-1,ks))
          nope = nint(cs0(2,0,ks))

c         Perform minimal printout

          if (ifdb) then
            write (debf,2005) rnsurf,ks,neps,nope,dnope-nope
          endif
          if (screen) then
            write (iow,2005) rnsurf,ks,neps,nope,dnope-nope
          endif
          if ((ior.lt.0) .and. (screen)) then
            write (*,2005) rnsurf,ks,neps,nope,dnope-nope
          endif

c         Perform dump of nodes if requested

          if ((ks.ge.nf) .and. (ks.le.nl)) then
            do ke = 1,neps
              i0 = ofs + (ke-1)*dnope - 1
              if (ifdb) then
                write (debf,2006) ke,(ics(i0+kn),kn = 1,dnope)
              endif
              if (screen) then
                write (iow,2006) ke,(ics(i0+kn),kn = 1,dnope)
              endif
              if ((ior.lt.0) .and. (screen)) then
                write (*,2006) ke,(ics(i0+kn),kn = 1,dnope)
              endif
            end do ! ke
          endif
          if ((ior.lt.0) .and. (screen) .and.
     &        (.not.pcomp(c,'r',1))          ) then
4005        call pprint ('  Enter = next, r = run to last,')
            call pprint(' g # = go to #, s = stop ===> ')
            if(rank.eq.0) read (*,2017,err=4005) c,kk
            if (pcomp(c,'s',1)) then
               ks = numcs+1
            elseif (pcomp(c,'g',1)) then
              ks = max(1,kk)
            else
              ks = ks+1
            endif
          else
            ks = ks+1
          endif
        end do ! while
      endif

c     Dump of history corresponence vector HIC

      if (pcomp(vect,'HIC',3)) then

c       Check if dump is limited to a specific pair

        if (nn.eq.0) then
          nf = 1
          nl = numcp

c         Search set

        else
          do kc = 1, cck(3)
            if (cp0(1,-1,kc).eq.nn) then
              nf = kc
              nl = kc
            endif
          end do ! kc
        endif

c       Perform dump

        kp = 1
        do while (kp.le.numcp)
          rnpair = nint(cp0(1,-1,kp))
          i0 =  (kp-1)*(c_lp1+c_lp3)
          if ((kp.ge.nf) .and. (kp.le.nl)) then
            if (ifdb) then
              write (debf,2014) rnpair,kp
              write (debf,2015) (hic(i0+kv),kv=1,c_lp1)
              write (debf,2016)
            endif
            if (screen) then
              write (iow,2014) rnpair,kp
              write (iow,2015) (hic(i0+kv),kv=1,c_lp1)
              write (iow,2016)
            endif
            if ((ior.lt.0) .and. (screen)) then
              write (*,2014) rnpair,kp
              write (*,2015) (hic(i0+kv),kv=1,c_lp1)
              write (*,2016)
            endif
            i0 = i0+c_lp1
            if (ifdb) then
              write (debf,2015) (hic(i0+kv),kv=1,c_lp3)
            endif
            if (screen) then
              write (iow,2015) (hic(i0+kv),kv=1,c_lp3)
            endif
            if ((ior.lt.0) .and. (screen)) then
              write (*,2015) (hic(i0+kv),kv=1,100)
            endif
          endif
          if ((ior.lt.0) .and. (screen) .and.
     &        (.not.pcomp(c,'r',1))          ) then
4006        call pprint ('  Enter = next, r = run to last,')
            call pprint(' g # = go to #, s = stop ===> ')
            if(rank.eq.0) read (*,2017,err=4006) c,kk
            if (pcomp(c,'s',1)) then
               kp = numcp+1
            elseif (pcomp(c,'g',1)) then
              kp = max(1,kk)
            else
              kp = kp+1
            endif
          else
            kp = kp+1
          endif
        end do ! while
      endif

c     Dump of history vector CH

      if (pcomp(vect,'CH',2)) then

c       Check if dump is limited to a specific pair

        if (nn.eq.0) then
          nf = 1
          nl = numcp
        else

c         Search set

          do kc = 1, cck(3)
            if (cp0(1,-1,kc).eq.nn) then
              nf = kc
              nl = kc
            endif
          end do ! kc
        endif

c       Perform dump

        do kp = 1,numcp

          rnpair = nint(cp0(1,-1,kp))

c         Find information

          ofh1 = nint(cp0(2,-1,kp))
          ofh3 = nint(cp0(3,-1,kp))
          lh1  = nint(cp0(4,-1,kp))
          lh3  = nint(cp0(5,-1,kp))
          nset = nint(cp0(6,-1,kp))
          ndv  = nint(cp0(1,0,kp))
          ns1  = nint(cp0(2,0,kp))
          ns2  = nint(cp0(3,0,kp))

c         Perform minimal printout

          if (ifdb) then
            write (debf,2009) rnpair,kp,ns1,ns2,lh1,lh3,ofh1,ofh3,
     &                        nset,ndv
          endif
          if (screen) then
            write (iow,2009) rnpair,kp,ns1,ns2,lh1,lh3,ofh1,ofh3,
     &                       nset,ndv
          endif
          if ((ior.lt.0) .and. (screen)) then
            write (*,2009) rnpair,kp,ns1,ns2,lh1,lh3,ofh1,ofh3,
     &                     nset,ndv
          endif

c         Perform dump of history data if requested

          if ((kp.ge.nf) .and. (kp.le.nl)) then
            ih1 = (c_lp1+c_lp3) * (kp-1)
            ih3 = ih1 + c_lp1

            ks = 1
            do while (ks.le.nset)

c             Set pointer for history data set

              i1 = ofh1 - 1 + lh1*(ks-1)
              i3 = ofh3 - 1 + lh3*(ks-1)

              if (ifdb) then
                write (debf,2010) rnpair,ks
              endif
              if (screen) then
                write (iow,2010) rnpair,ks
              endif
              if ((ior.lt.0) .and. (screen)) then
                write (*,2010) rnpair,ks
              endif

c             Loop on all ch1 variables

              do kv = 1, lh1

c               Search corresponding name of variable

                whfl = .true.
                kp1 = 0
                do while (whfl)
                  kp1 = kp1+1
                  if (kp1.gt.lp1) then
                    prnam = .false.
                    whfl  = .false.
                  elseif (hic(ih1+kp1).eq.kv) then
                    prnam = .true.
                    whfl  = .false.
                  endif
                end do ! while

                if (ifdb) then
                  if (prnam) then
                    write (debf,2011) ch1(i1+kv),ch2(i1+kv),kv,kp1,
     &                                w1(kp1,ndv)
                  else
                    write (debf,2011) ch1(i1+kv),ch2(i1+kv),kv
                  endif
                endif
                if (screen) then
                  if (prnam) then
                    write (iow,2011) ch1(i1+kv),ch2(i1+kv),kv,kp1,
     &                                w1(kp1,ndv)
                  else
                    write (iow,2011) ch1(i1+kv),ch2(i1+kv),kv
                  endif
                endif
                if ((ior.lt.0) .and. (screen)) then
                  if (prnam) then
                    write (*,2011) ch1(i1+kv),ch2(i1+kv),kv,kp1,
     &                                w1(kp1,ndv)
                  else
                    write (*,2011) ch1(i1+kv),ch2(i1+kv),kv
                  endif
                endif
              end do ! kv

              if (ifdb) then
                write (debf,2012)
              endif
              if (screen) then
                write (iow,2012)
              endif
              if ((ior.lt.0) .and. (screen)) then
                write (*,2012)
              endif

c             Loop on all ch3 variables

              do kv = 1, lh3

c               Search corresponding name of variable

                whfl = .true.
                kp3 = 0
                do while (whfl)
                  kp3 = kp3+1
                  if (kp3.gt.lp3) then
                    prnam = .false.
                    whfl  = .false.
                  elseif (hic(ih3+kp3).eq.kv) then
                    prnam = .true.
                    whfl  = .false.
                  endif
                end do ! while
                if (ifdb) then
                  if (prnam) then
                    write (debf,2013) ch3(i3+kv),kv,kp3,w3(kp3,ndv)
                  else
                    write (debf,2013) ch3(i3+kv),kv
                  endif
                endif
                if (screen) then
                  if (prnam) then
                    write (iow,2013) ch3(i3+kv),kv,kp3,w3(kp3,ndv)
                  else
                    write (iow,2013) ch3(i3+kv),kv
                  endif
                endif
                if ((ior.lt.0) .and. (screen)) then
                  if (prnam) then
                    write (*,2013) ch3(i3+kv),kv,kp3,w3(kp3,ndv)
                  else
                    write (*,2013) ch3(i3+kv),kv
                  endif
                endif
              end do ! kv
              if ((ior.lt.0) .and. (screen) .and.
     &            (.not.pcomp(c,'r',1))          ) then
4007            call pprint ('  Enter = next, r = run to last,')
                call pprint(' g # = go to #, s = stop ===> ')
                if(rank.eq.0) read (*,2017,err=4007) c,kk
                if (pcomp(c,'s',1)) then
                  ks = nset+1
                elseif (pcomp(c,'g',1)) then
                  ks = max(1,kk)
                else
                  ks = ks+1
                endif
              else
                ks = ks+1
              endif
            end do ! while
          endif
        end do ! kp
      endif

c     Formats

2000  format (/1x,'  Array CS0 - surface #    ',i5/
     &         1x,'     Internal surface #    ',i5/)
2001  format ( 1x,1p,7e10.2)
2002  format (/1x,'  Array CM0 - material #   ',i5/
     &         1x,'     Internal material #   ',i5/)
2003  format (/1x,'  Array CP0 - pair #       ',i5/
     &         1x,'     Internal pair #       ',i5/)
2004  format ( i5,1p,7e10.2/ 5x,1p,7e10.2)
2005  format (/1x,'  Array ICS - surface #    ',i5/
     &         1x,'     Internal surface #    ',i5/
     &         1x,'  # of elements            ',i5/
     &         1x,'  # of nodes per element   ',i5/
     &         1x,'  # of dims  per element   ',i5/)
2006  format ( 1x,i7,',,',20i7)

2007  format (/1x,'  Array CM - material #    ',i5/
     &         1x,'    Internal material #    ',i5/)
2008  format ( 1x,1p,5e15.7)

2009  format (/1x,'  Array CH1 - pair #       ',i5/
     &         1x,'     Internal pair #       ',i5/
     &         1x,'  First  surface #         ',i5/
     &         1x,'  Second surface #         ',i5/
     &         1x,'  # of terms in ch1 & ch2  ',i5/
     &         1x,'  # of terms in ch3        ',i5/
     &         1x,'  Offset for    ch1 & ch2  ',i5/
     &         1x,'  Offset for    ch3        ',i5/
     &         1x,'  # of history set         ',i5/
     &         1x,'  # of contact element     ',i5)
2010  format (/1x,' Dump of set  #            ',i5,
     &            ' Internal set # ',i5/
     &         1x,'         History set CH1',
     &            '         History set CH2',
     &            '  # in CH1','   Absol #',
     &            '  Name')
2011  format ( 1x,1p,2e24.15,2i10,2x,a8)
2012  format (/1x,24x,
     &            '         History set CH3',
     &            '  # in CH3','   Absol #',
     &            '  Name')
2013  format ( 1x,24x,1p,e24.15,2i10,2x,a8)
2014  format (/1x,'  Array HIC - pair #       ',i5/
     &         1x,'     Internal pair #       ',i5//
     &         1x,'  Corresp. for CH1 & CH2   ')
2015  format ( 1x,20i4)
2016  format (/1x,'  Corresp. for CH3         ')

2017  format ( a,i10)
2018  format ( 1x,7a10)
2019  format ( 1x,7i10)

      end
