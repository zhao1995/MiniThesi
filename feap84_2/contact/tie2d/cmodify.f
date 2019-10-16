c$Id:$
      subroutine cmodify(ism,ixs,id, nel,mel,ndf,nst,tanm,resn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor             6 March 2003            1.0

c      Acronym: Contact MODIFY for 2D Tied
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nel,mel,ndf,nst
      integer    i,j, n,ns, m,nm, il, im
      integer    ixs(nel),id(ndf,*), ism(2)
      real*8     tanm(nst,nst),resn(nst)

      il = (nel+mel)*ndf
      im = il
      do m = 1,2
        ns = ixs(m)
        do n = 1,2
          nm = ism(n)
          if(nm.gt.0) then
            do i = 1,ndf
              if(id(i,ns).ne.0 .and. id(i,nm).ne.0) then  ! fixed dof
                resn(il+i) = 0.0d0
                do j = 1,im
                  tanm(j,il+i) = 0.0d0
                  tanm(il+i,j) = 0.0d0
                end do ! j
                tanm(il+i,il+i) = 1.0d0
              endif
            end do ! i
          endif
        end do ! n
        il = il + ndf
      end do ! m

      end
