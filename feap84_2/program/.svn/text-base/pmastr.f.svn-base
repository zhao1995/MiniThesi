c$Id:$
      subroutine pmastr(ixt,rlink,ntyp,x,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Construct table of master to slave nodes (small deform)

c     Inputs:
c        ntyp(*)      - Active node type
c        x(ndm,*)     - Nodal coordinate array
c        prt          - Print flag: Output results if 'true'

c     Outputs:
c        ixt(*)       - Array to hold list of slave/master nodes
c        rlink(ndf,*) - List of master-slave node links
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'sdata.h'

      character txt*15
      logical   prt, errck, tinput, pcomp, readfl, clflg
      integer   i, idir, n, mas,slv, ixt(*),rlink(ndf,*),ntyp(*)
      real*8    dotx, tmn, xmn, gap, x(ndm,*),td(15)

      save

      readfl = .true.
      gap    =  1.d-08

c     Perform set for master-slave values

      if(prt) then
        write(iow,2000) head, (i,i=1,ndf)
      endif
      do while (readfl)

c       Input data item

        errck = tinput(txt,1,td,15)

c       Set gap to different value

        if(pcomp(txt,'gap ',4)) then

          gap = td(1)
          write(iow,2001) gap

c       Find number for master node

        elseif(.not.pcomp(txt,'    ',4)) then
          clflg = .false.
          do n = 1,numnp
            if(ntyp(n).ge.0) then
              tmn = dotx(td(1),x(1,n),ndm)
              if(clflg) then
                if(tmn.lt.xmn) then
                  xmn = tmn
                  mas = n
                endif
              else
                xmn   =  tmn
                mas   =  n
                clflg = .true.
              endif
            endif
          end do ! n
        endif

c       Slave node specification

        if(pcomp(txt,'node',4) .or. pcomp(txt,'slav',4)) then

c         Find slave node

          clflg = .false.
          do n = 1,numnp
            if(ntyp(n).ge.0) then
              tmn = dotx(td(ndm+1),x(1,n),ndm)
              if(clflg) then
                if(tmn.lt.xmn) then
                  xmn = tmn
                  slv = n
                endif
              else
                xmn   =  tmn
                slv   =  n
                clflg = .true.
              endif
            endif
          end do ! n

          if(clflg) then
            if(ixt(slv).eq.0 .and. ixt(mas).eq.0) then
              ixt(slv) = mas
              do i = 1,ndf
                rlink(i,slv) = nint(td(ndm+ndm+i))
              end do ! i
              if(prt) then
                write(iow,2002) mas,slv,(rlink(i,slv),i=1,ndf)
              endif
            else
              write(ilg,3000) slv,mas,ixt(mas)
              write(iow,3000) slv,mas,ixt(mas)
              call plstop()
            endif
          endif

c       Slave surface specification

        elseif(pcomp(txt,'surf',4)) then

          idir = nint(td(ndm+1))
          do n = 1,numnp
            if(abs(x(idir,n) - x(idir,mas)).lt.gap) then
              if(n.ne.mas) then
                if(ixt(n).eq.0 .and. ixt(mas).eq.0) then
                  ixt(n) = mas
                  do i = 1,ndf
                    rlink(i,n) = nint(td(ndm+1+i))
                  end do ! i
                  if(prt) then
                    write(iow,2002) mas,n,(rlink(i,n),i=1,ndf)
                  endif
                else
                  write(iow,3000) n,mas,ixt(mas)
                  call plstop()
                endif
              endif
            endif
          end do ! n

c       Exit loop

        elseif(pcomp(txt,'    ',4)) then

          readfl = .false.

        endif

      end do ! while

c     Formats

2000  format(1x,19a4,a3//5x,'M a s t e r   t o   S l a v e   N o d e',
     &       3x,'D e f i n i t i o n s'//
     &       5x,'Master   Slave',7(i3,'-link')/(21x,7(i3,'-link')))

2001  format(/5x,'Search Gap Tolerance =',1p,1e13.5/)

2002  format(2x,9i8/(18x,7i8))

3000  format(/' *ERROR* PMASTR: Slave node',i8,' can not be assigned',
     &        ' to master node',i8/'    since it was previously',
     &        ' assigned as slave to master node',i8/)
      end
