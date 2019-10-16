c$Id:$
      subroutine aslid2b(con, slid, numnp, ns)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d
c               Establish slideline array

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'iofile.h'

      integer             slidn
      common     /aslids/ slidn

      integer    numnp,nslid
      integer    i,ii, jj, n,ns

      integer    con(2,numnp), slid(2,*)

      call cdebug0 ('    aslid2b',-1)

c     Establish slide lines

      jj    = 0

      do  n = 1,numnp
       if(con(1,n).gt.0) then

c        Start list
         jj         = jj + 1
         ns         = jj
         slid(1,jj) = n
         i          = con(1,n)
         ii         = n

c        Fill list

         do while (i.ne.n)
           if(ii.eq.0 .or. ii.gt.numnp) then
             write(  *,*) 'ERROR in ASLID2b',n,ii,jj
             write(iow,*) 'ERROR in ASLID2b',n,ii,jj
             call iprint(con ,2,numnp,2,'Error CON ')
             call iprint(slid,2,   jj,2,'Error SLID')
             call plstop()
           endif
           jj         = jj + 1
           slid(1,jj) = con(1,ii)
           slid(2,jj) = con(2,ii)
           i          = con(1,ii)
           con(1,ii)  = 0
           ii         = i
         end do ! while

c        Terminate list

         slid(2,ns) = jj

        endif
      end do ! n

      if(ifdb) then
        call iprint(slid,2,jj,2,'SLID')
      endif

      nslid = jj
      slidn = jj

      ns = 0
      jj = 1
100   ns = ns + 1
      jj = slid(2,jj) + 1
      if(jj.lt.nslid) go to 100

      end
