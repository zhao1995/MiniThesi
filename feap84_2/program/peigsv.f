c$Id:$
      subroutine peigsv( lct, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Save eigen-pairs on disk for use in another problem.

c     Inputs:
c              lct  - Filename for pairs
c              isw  - 1 - write eigenpairs to   "filename"
c                     2 - read  eigenpairs from "filename"

c     Outputs:
c              Eigen-pairs on file: "filename"
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include  'cdata.h'
      include  'comblk.h'
      include  'evdata.h'
      include  'iofile.h'
      include  'pointer.h'

      character lct*15
      logical   exst,palloc
      integer   neqold,isw,i

      save

c     Write a file

      if(isw.eq.1) then
        open (unit = 35, file = lct , form = 'unformatted')

        write(35) mf,mq,neq,neqold
        write(35) (hr(np(76)+ i),i=0,mq-1),(hr(np(77)+i),i=0,mq*neq-1)
        close(35)

        write(iow,2000) lct
        if(ior.lt.0) then
          write(*,2000) lct
        endif

c     Read a set of eigen pairs

      elseif(isw.eq.2) then
        inquire(file = lct , exist = exst )
        if(exst) then
          open (unit = 35, file = lct , form = 'unformatted')

          read (35) mf,mq,neqold
          if(neq.ne.neqold) then
            if(ior.lt.0) then
              write(*,3000) neq,neqold
            else
              write(ilg,3000) neq,neqold
              write(iow,3000) neq,neqold
              call plstop()
            endif
          else
c           Allocate space for eigenpairs if necessary
            exst = palloc( 76,'EVAL',mq    , 2 )
            exst = palloc( 77,'EVEC',mq*neq, 2 )
            npev = np(77)
            read(35) (hr(np(76)+ i),i=0,mq-1),
     &               (hr(np(77)+i),i=0,mq*neq-1)
            close(35)
            write(iow,2001) lct
            if(ior.lt.0) then
              write(*,2001) lct
            endif
          endif
        else
          if(ior.lt.0) then
            write(*,3001) lct
          else
            write(ilg,3001) lct
            write(iow,3001) lct
            call plstop()
          endif
        endif
      endif

c     Formats

2000  format(/5x,'Eigenpairs saved on file:',a/)

2001  format(/5x,'Eigenpairs read from file:',a/)

3000  format(' *ERROR* PEIGSV: Number equations differs from current',
     &       ' problem'/
     &       '         Current neq =',i9,': Old neq =',i9)
3001  format(' *ERROR* PEIGSV: File:',a,' does not exist, respecify')

      end
