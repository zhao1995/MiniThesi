c$Id:$
      subroutine prtlis(x,b,ttim,prop,ndm,ndf,nlist,list,ii,
     &                  ndf1,ndf2,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal values for real solutions based on list

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         b(*)      - Current value of solution
c         ttim      - Value of solution time
c         prop      - Value of total proportional load
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         nlist     - Number items in list
c         list(*)   - List of node numbers to output
c         ii        - Type of output: 1 = displacement; 2 = velocity;
c                                     3 = acceleration; 4 = eigenvector
c         ndf1      - First component to output
c         ndf2      - Last  component to output
c         prth      - Output title/header data if true

c      Outputs:
c         None      - Outputs to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   prth,tot
      character nd*4,cd*6,di(4)*6,fmt1*30,fmt2*30
      integer   ndm,ndf,nlist, n,nn, i,ii, ndf1,ndf2, kount, list(*)
      real*8    ttim,prop, x(ndm,*),b(ndf,*)

      save

      data      nd   /'Node'/,cd /' Coord'/
      data      di   /' Displ',' Veloc',' Accel',' EigV.'/
      data      fmt1 /'(3x,a4,3(i6,a6)/(7x,6(i6,a6)))'/
      data      fmt2 /'(i7,1p,3e12.4/(7x,1p,6e12.4))'/

      tot   = (ndf+ndm).le.6
      if(.not.tot) then
        write(fmt1(8:8),'(i1)') ndm
        write(fmt2(8:8),'(i1)') ndm
      endif
      kount = 0
      do nn = 1,nlist
        n = list(nn)
        kount = kount - 1
        if(kount.le.0) then
          call prtitl(prth)
          if(ii.le.3) then
            write(iow,2000) ttim,prop
            if(ior.lt.0) then
              write(*,2000) ttim,prop
            endif
          else
            write(iow,2001) ttim,prop
            if(ior.lt.0) then
              write(*,2001) ttim,prop
            endif
          endif
          if(tot) then
            write(iow,2002) (i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
            if(ior.lt.0) then
              write(*,2002) (i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
            endif
          else
            write(iow,fmt1) nd,(i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
            if(ior.lt.0) then
              write(*,fmt1) nd,(i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
            endif
          endif
          kount = 48000000
        endif
        if(tot) then
          write(iow,2003) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
          if(ior.lt.0) then
            write(*,2003) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
          endif
        else
          write(iow,fmt2) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
          if(ior.lt.0) then
            write(*,fmt2) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
          endif
        endif
      end do ! nn

c     Formats

2000  format('  N o d a l   O u t p u t s',17x,
     &  'Time',e18.5/44x,'Prop. Ld.',1pe13.5/1x)

2001  format('  N o d a l   O u t p u t s',17x,
     &  'Time',e18.5/43x,'Eigenvalue',1pe13.5/1x)

2002  format('   Node',6(i6,a6):/(7x,6(i6,a6):))

2003  format(i7,1p,6e12.4:/(7x,1p,6e12.4:))

      end
