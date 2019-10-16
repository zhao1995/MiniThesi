c$Id:$
      subroutine prtstr(x,dp,ds,ndm,numnp,n1,n2,n3,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set point for hist plots                         20/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal projected stresses and principal values

c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         dp(numnp,*) - Principal values at nodes
c         ds(numnp,*) - Stress values at nodes
c         ndm         - Spatial dimension of mesh
c         numnp       - Number of nodes in mesh
c         n1          - Number of first node to output
c         n2          - Number of last  node to output
c         n3          - Increment to node from n1
c         prth        - Output title/header data if true

c      Outputs:
c         None        - Outputs to file/screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'eldatp.h'
      include  'fdata.h'
      include  'pfeapb.h'
      include  'p_point.h'
      include  'strnum.h'
      include  'xtout.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   cknon0,vnon0,prth
      integer   ndm,numnp,n1,n2,n3, i, n, ista, count, nxt1
      integer   lnod,gnod, bserchi
      real*8    x(ndm,*),dp(numnp,*),ds(numnp,*)

      save

c     Determine how many non-zero nodal values there are

      ista = 0
      do n = 1,abs(istv)
        if(cknon0(ds(1,n),numnp) ) ista = n
      end do ! n

      if(ista.gt.0 .or. histpltfl) then
        count = 0
        nxt1  = max(1,nxt)
        do n = n1,n2,n3
          if( nxt.eq.0 .or. abs(x(nxt1,n)-xt).le.xtol ) then
            vnon0 = .false.
            do i = 1,ista
              if(ds(n,i).ne.0.0d0) vnon0 = .true.
            end do ! i
            if(vnon0 .or. histpltfl) then
              count = count - 1
              if(count.le.0) then
                call prtitl(prth)
                write(iow,2000)
                write(iow,2001) (i,i=1,3),(i,i=1,ista)
                if(histpltfl) then
                  write(iow,2002) (i,i=1,hplmax)
                endif
                if(ior.lt.0.and.pfr) then
                  write(*,2000)
                  write(*,2001) (i,i=1,3),(i,i=1,ista)
                  if(histpltfl) then
                    write(*,2002) (i,i=1,hplmax)
                  endif
                endif
                count = 50000000
              endif
              if(.not.pfeap_on) then
                write(iow,2003) n,(dp(n,i),i=1,7),(ds(n,i),i=1,ista)
                if(ior.lt.0.and.pfr) then
                  write(*,2003) n,(dp(n,i),i=1,7),(ds(n,i),i=1,ista)
                endif
                if(histpltfl) then
                  point = np(305) + n - 2
                  gnod  = (hplmax - 1)*numnp + 1
                  write(iow,2005) (hr(point+i),i=1,gnod,numnp)
                  if(ior.lt.0.and.pfr) then
                    write(*,2005) (hr(point+i),i=1,gnod,numnp)
                  endif
                endif
              else
                if(pfeap_gnod) then
                  lnod = bserchi(mr(np(244)),numpn, n)
                  gnod = mr(np(244)+lnod-1)
                else
                  lnod = n
                  gnod = mr(np(244)+n-1)
                endif
                if(gnod.gt.0) then
                  write(iow,2004) lnod,(dp(lnod,i),i=1,4),gnod,
     &                                 (dp(lnod,i),i=5,7),
     &                                 (ds(lnod,i),i=1,ista)
                  if(ior.lt.0.and.pfr) then
                    write(*,2004) lnod,(dp(lnod,i),i=1,4),gnod,
     &                                 (dp(lnod,i),i=5,7),
     &                                 (ds(lnod,i),i=1,ista)
                  endif
                  if(histpltfl) then
                    point = np(305) + numnp*(lnod-1) - 1
                    gnod  = hplmax*numnp - 1
                    write(iow,2005) (hr(point+i),i=1,gnod,numnp)
                    if(ior.lt.0.and.pfr) then
                      write(*,2005) (hr(point+i),i=1,gnod,numnp)
                    endif
                endif
                endif
              endif
            endif
          endif
        end do ! n
      else
        if(ior.lt.0.and.pfr) write(*,*) 'All values zero'
        write(iow,*) 'All values zero'
      endif

c     Formats

2000  format('   N o d a l   P r o j e c t i o n s')
2001  format(/'   Node', 3(i3,'-Pr.Value'),'  1-Pr.Angle'/
     & 8x,'   I_1 Value   J_2 Value   J_3 Value':/(8x,6(i6,' Value'):))
2002  format((8x,6(i4,' History'):))

c2003  format(/i8,1p,4e16.8/8x,1p,3e12.4/(8x,1p,6e16.8:))
2003  format(/i8,1p,4e12.4/8x,1p,3e12.4/(8x,1p,6e12.4:))
2004  format(/i8,1p,4e12.4/i8,1p,3e12.4/(8x,1p,6e12.4:))
2005  format(8x,1p,6e12.4)

      end
