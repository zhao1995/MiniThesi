c$Id:$
      subroutine uimport(lct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: User function to import displacements into FEAP

c      Inputs:
c         lct    - File name for imports

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   pcomp,cfr,test
      character lct*(*),fnamr*21
      integer   i,kk
      real*8    ul(16)

      save

c     Compute number dof's and open import file

      fnamr = lct
      if(.not.pcomp(fnamr,'    ',4)) then
        call opnfil(fnamr,fnamr,-2,ios,cfr)
        write(iow,2000) fnamr
        if(ior.lt.0) then
          write(*,2000) fnamr
        endif

c       Read all records from file & save in F0 array displacement part.

        if(cfr) then
          test = .true.
          do while (test)
            read(ios,*,end=1501,err=1501) kk,(ul(i),i=1,ndf)
            if(kk.le.0. or. kk.gt.numnp) then
              test = .false.
            else
              fp(1) = np(31) + (kk-1)*ndf - 1           ! ID - array
              fp(2) = np(28) + (kk-1)*ndf - 1 + nneq    ! F0 - array
              do i = 1,ndf
                if(mr(fp(1)+i).le.0) then
                   hr(fp(2)+i) = ul(i)
                endif
              end do
            endif
          end do
1501      close(ios)
        endif
      else
        write(*,*) ' *ERROR* specify file name for import'
      endif

2000  format('   Import DISPLACEMENTS from: ',a/)

      end
