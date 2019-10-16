c$Id:$
      subroutine pleige(ct,ie,ix,x,dr,eval,evec,cs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot element eigenvectors

c      Inputs:
c         ct(3)     - Parameters
c         ie(nie,*) - Element directives
c         ix(nen1,*)- Element connection
c         x(ndm,*)  - Nodal coordinates
c         eval(*)   - Element Eigenvalues
c         evec(*)   - Element Eigenvectors
c         cs        - Scale factor

c      Temporary storage
c         dr(3,*)   - Deformed position

c      Outputs:
c         none      - Graphics to screen and PostScript file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'iofile.h'
      include  'pdata3.h'
      include  'sdata.h'

      integer   i,ii,jj,k1,k2,k3
      real*8    cs,dumv

      integer   ie(*),ix(nen1,*)
      real*8    ct(*),x(ndm,*),dr(3,*),eval(*),evec(ndf,nen,*)

      save

c     Add mode shape to last element coordinates

      k1 = nint(abs(ct(1)))
      if(k1.eq.0) k1 = 1
      if(ior.lt.0) then
        write(*,2005) k1,eval(k1)
      endif
      if(ct(3).eq.0.0d0) then
        call pleigt(eval(k1))
      endif

      dumv = 0.0d0
      do i = 1,nen
        k3 = ix(i,numel)
        if(k3.gt.0) then
          do k2 = 1,ndm
            dumv = max(abs(evec(k2,i,k1)),dumv)
          end do ! k2
        end if
      end do ! i

      if(dumv.gt.0.0d0) then
        call plopen
        if(ct(2).gt.0.0d0) then
          k3 = nint(ct(2))
          call pppcol(k3,0)
        else
          call pppcol(k1+1,0)
        end if
        dumv = cs/dumv
        do ii = 1,2
          jj = 0
          do i = 1,nen
            k3 = ix(i,numel)
            if(k3.gt.0) then
              jj = i
              do k2 = 1,ndm
                dr(k2,i) = x(k2,k3) + dumv*evec(k2,i,k1)
              end do ! k2
              do k2 = ndm+1,3
                dr(k2,i) = 0.0d0
              end do ! k2
            end if
          end do ! i

c         Plot element

          k1 = nie*ix(nen1,numel) - 1
          if(hide) then
            k3 = 1
          else
            k3 = 3
          endif
          call plot9(ie(k1),ix(1,numel),dr,3,jj,k3)
          call pppcol(1,1)
          dumv = 0.0d0
        end do ! ii
      else
        if(ior.lt.0) then
          write(*,*) '   *WARNING* Eigenvector zero at nodes'
        endif
      endif

c     Format

2005  format(5x,'Eigenvalue',i4,' = ',1p,e13.5)

      end
