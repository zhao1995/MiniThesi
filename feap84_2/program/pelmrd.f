c$Id:$
      subroutine pelmrd(ix,rben,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add set of last element number to last_elm       29/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input all element connections

c      Inputs:
c         prt        - Print results if ture

c      Outputs:
c         ix(nen1,*) - Element connection data
c         rben(*)    - Rigid body/flexible element indicator
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'chdata.h'
      include   'comfil.h'
      include   'sdata.h'
      include   'iofile.h'
      include   'ioincl.h'
      include   'p_ptname.h'
      include   'region.h'
      include   'rigid2.h'

      logical    prt, pcomp, incfl
      character  inam*4,fnam*15,gnam*4
      integer    i,j,n,noge, ibc,ix(nen1,numel),rben(*)

c     Check for an include record
c     [include filename <noge> - Input 'filename': noge = omit field
c     [nogeneration            - Omit generation field

      read(ior,1000) record
      call pstrip(xxx,record,1)
      call acheck(xxx,yyy,15,80,80)
      read(yyy,1002)  inam,fnam,gnam
      if(pcomp(inam,'incl',4)) then
        call pincld(fnam)
        noge  =  1
        incfl = .true.
        if(pcomp(gnam,'noge',4)) then
          noge  =  0
        endif
      elseif(pcomp(inam,'noge',4)) then
        noge  =  0
        incfl = .false.
      else
        backspace(ior)
        noge  =  1
        incfl = .false.
      endif

c     List directed input of all the data

      do n = 1,numel
        if(noge.eq.0) then
          read(ior,*) j,    ix(nen1,j),(ix(i,j),i=1,nen)
        else
          read(ior,*) j,ibc,ix(nen1,j),(ix(i,j),i=1,nen)
        endif
      end do ! n
      last_elm = numel

c     Set region and rigid body indicators

      do n = 1,numel
        rben(n)      = nrigid
        ix(nen1-1,n) = nreg
      end do ! n

c     Output data
      if(prt) then
        call prtitl(prt)
        write(iow,2000)  (i,i=1,nen)
        do n = 1,numel
          write(iow,2001) n,ix(nen1,n),ix(nen1-1,n),(ix(i,n),i=1,nen)
        end do ! n
      endif

      irecrd(isf) = irecrd(isf) + numel

c     Clear the include file if necessary

      if(incfl) then
        call pincld('end')
      endif

c     Formats

1000  format(a)
1002  format(a4,11x,a15,a4)

2000  format(5x,'E l e m e n t s'//3x,'Elmt Mat Reg',8(i3,' Node':),
     &      (15x,8(i3,' Node':)))
2001  format(i7,2i4,8i8:/(15x,8i8))

      end
