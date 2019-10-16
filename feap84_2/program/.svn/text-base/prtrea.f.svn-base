c$Id:$
      subroutine prtrea(r,x,u,ndm,ndf,n1,n2,n3,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dabs' to 'abs'                           17/11/2006
c       2. Add exchange of sums for 'glob'al option         13/04/2007
c       3. Add computation of moments for 2 and 3-d probs   12/04/2013
c          in reference and current configurations.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal reactions for current solution

c      Inputs:
c         r(*)      - Current value of reactions
c         x(ndm,*)  - Nodal coordinates of mesh
c         u(ndf,*)  - Nodal displacements of mesh
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last node to output
c         n3        - Increment to n1
c         prth      - Output title/header data if true

c      Outputs:
c         None      - Outputs to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'xtout.h'
      include  'pfeapb.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   prth,nzprt
      integer   ndm,ndf,n1,n2,n3, i,k,n,  count, nxt1
      real*8    x(ndm,*),u(ndf,*),r(ndf,*),rsum(30),asum(30),psum(30)
      real*8    tbuf(60),msum(3),usum(3)

      save

      do k = 1,ndf
        psum(k) = 0.d0
        rsum(k) = 0.d0
        asum(k) = 0.d0
      end do ! k
      msum = 0.d0
      usum = 0.d0
      do i = 1,numnp
        do k = 1,ndf
          rsum(k) = rsum(k) - r(k,i)
          asum(k) = asum(k) + abs(r(k,i))
        end do ! k
      end do ! i
      count = 0
      nxt1  = max(1,nxt)
      fp(1) = np(190) - 1
      do n = n1,n2,n3
        if( (mr(fp(1)+n).ge.0) .and. (nxt.eq.0 .or.
     &      (abs(x(nxt1,n)-xt).le.xtol) )  ) then
          nzprt = .false.
          do k = 1,ndf
            psum(k) = psum(k) - r(k,n)
            if(r(k,n).ne.0.0d0) then
              nzprt = .true.
            endif
          end do ! k
          if(nzprt) then
            if(ndm.eq.2 .and. ndf.ge.2) then
              msum(3) = msum(3) + x(2,n)*r(1,n) - x(1,n)*r(2,n)
              usum(3) = usum(3) + u(2,n)*r(1,n) - u(1,n)*r(2,n)
            elseif(ndm.eq.3 .and. ndf.ge.3) then
              msum(1) = msum(1) + x(3,n)*r(2,n) - x(2,n)*r(3,n)
              msum(2) = msum(2) + x(1,n)*r(3,n) - x(3,n)*r(1,n)
              msum(3) = msum(3) + x(2,n)*r(1,n) - x(1,n)*r(2,n)
              usum(1) = usum(1) + u(3,n)*r(2,n) - u(2,n)*r(3,n)
              usum(2) = usum(2) + u(1,n)*r(3,n) - u(3,n)*r(1,n)
              usum(3) = usum(3) + u(2,n)*r(1,n) - u(1,n)*r(2,n)
            endif
            count = count - 1
            if(count.le.0) then
              call prtitl(prth)
              write(iow,2000) (k,k=1,ndf)
              if(ior.lt.0.and.pfr) then
                write(*,2000) (k,k=1,ndf)
              endif
              count = 50000000
            endif
            if(ior.lt.0.and.pfr) then
              write(*,2001) n,(-r(k,n),k=1,ndf)
            endif
            write(iow,2001) n,(-r(k,n),k=1,ndf)
          endif
        endif
      end do ! n

c     Print sum checks

      write(iow,2002) ' Pr.Sum',(psum(k),k=1,ndf)
      write(iow,2003) '   Sum ',(rsum(k),k=1,ndf)
      write(iow,2003) '  |Sum|',(asum(k),k=1,ndf)
      write(iow,2003) '  M_ref',(msum(k),k=1,3)
      write(iow,2003) '  M_cur',(msum(k)+usum(k),k=1,3)
      if(ior.lt.0.and.pfr) then
        write(*,2002) ' Pr.Sum',(psum(k),k=1,ndf)
        write(*,2003) '   Sum ',(rsum(k),k=1,ndf)
        write(*,2003) '  |Sum|',(asum(k),k=1,ndf)
        write(*,2003) '  M_ref',(msum(k),k=1,3)
        write(*,2003) '  M_cur',(msum(k)+usum(k),k=1,3)
      endif

c     Exchange of sums for parallel run

      if(pfeap_on) then
        write(iow,2004) (k,k=1,ndf)
        if(ior.lt.0.and.pfr) then
          write(*,2004) (k,k=1,ndf)
        endif
        call pfeapsr(psum,tbuf,ndf,.true.)
        call pfeapsr(rsum,tbuf,ndf,.true.)
        call pfeapsr(asum,tbuf,ndf,.true.)
        call pfeapsr(msum,tbuf,  3,.true.)
        call pfeapsr(usum,tbuf,  3,.true.)
        write(iow,2003) ' Pr.Sum',(psum(k),k=1,ndf)
        write(iow,2003) '   Sum ',(rsum(k),k=1,ndf)
        write(iow,2003) '  |Sum|',(asum(k),k=1,ndf)
        write(iow,2003) '  M_ref',(msum(k),k=1,3)
        write(iow,2003) '  M_cur',(msum(k)+usum(k),k=1,3)
        if(ior.lt.0.and.pfr) then
          write(*,2003) ' Pr.Sum',(psum(k),k=1,ndf)
          write(*,2003) '   Sum ',(rsum(k),k=1,ndf)
          write(*,2003) '  |Sum|',(asum(k),k=1,ndf)
          write(*,2003) '  M_ref',(msum(k),k=1,3)
          write(*,2003) '  M_cur',(msum(k)+usum(k),k=1,3)
        endif
      endif

c     Formats

2000  format('  N o d a l    R e a c t i o n s'//
     &  '   Node',6(i8,' dof'):/(7x,6(i8,' dof'):))

2001  format(i7,1p,6e12.4:/(7x,1p,6e12.4:))

2002  format(/a,1p,6e12.4:/(7x,1p,6e12.4:))

2003  format( a,1p,6e12.4:/(7x,1p,6e12.4:))

2004  format(/'  G l o b a l    N o d a l    S u m s'//
     &  '   Node',6(i8,' dof'):/(7x,6(i8,' dof'):))

      end
