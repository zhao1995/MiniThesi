c$Id:$
      subroutine updatl(eq,id,ixt,du,u,ul,f)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'ndf' and 'numnp' from argument list,     27/12/2007
c          add include 'cdata' and 'sdata' to retrieve.
c       2. Add triad rotation treatment                     11/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Control program for update to solution values
c                Considers: (1) Restrained (specified) boundaries;
c                           (2) Sloping 2-d boundaries;
c                           (3) Euler angle 3-d boudaries

c      Inputs:
c        eq(ndf,*)   - Active equation numbers
c        id(ndf,*)   - Boundary condition codes
c        ixt(*)      -
c        du(*)       - Solution increments
c        u(ndf,*)    - Previous solution values
c        f(ndf,*)    - Boundary displacement values
c        ndf         - Number degree of freedoms at each node
c        numnp       - Number of nodes

c      Working
c        du(ndf)     - Solution update increments

c      Outputs:
c        u(ndf,*)    - Updated solution values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'corset.h'
      include   'counts.h'
      include   'iofile.h'
      include   'mdata.h'
      include   'part0.h'
      include   'ptest.h'
      include   'p_point.h'
      include   'sdata.h'
      include   'tdatb.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    i,j,n
      real*8     cc
      integer    eq(ndf,*),id(ndf,*),ixt(*)
      real*8     du(*),u(ndf,numnp,*),ul(ndf,3),f(ndf,numnp)
      real*8     ubl(30),ang

c     If test solution step set 'cc' to small value 'testva'

      if( testfl .and. niter.le.1 ) then
        cc = testva
      else
        cc = 1.d0 ! full step
      endif

      do n = 1,numnp
        do i = 1,ndf
          ubl(i) = u(i,n,1)
        end do ! i

c       Check for sloping angle boundaries

        if(anglefl) then
          if(hr(np(45)+n-1).ne.0.0d0) then
            ang = hr(np(45)+n-1)
            call upang(dal,ang,ubl,ndf,1)
            if(ndf.ge.6) then
              call upang(ral,ang,ubl,ndf,1)
            endif
          endif
        endif

c       Check for Euler angle boundaries

        if(eulerfl) then
          if(hr(np(242)+3*n-3).ne.0.0d0 .or.
     &       hr(np(242)+3*n-2).ne.0.0d0 .or.
     &       hr(np(242)+3*n-1).ne.0.0d0) then
            call upeul(dal,hr(np(242)+3*n-3),ubl,ndf,1)
            if(ndf.ge.6) then
              call upeul(ral,hr(np(242)+3*n-3),ubl,ndf,1)
            endif
          endif
        endif

c       Check for triad angle boundaries

        if(triadfl) then
          point = np(274) + 9*n - 9
          if(hr(point).ne.-100.0d0) then
            call uptriad(dal,hr(point),ubl,ndf,1)
            if(ndf.ge.6 .and. ral(1).ne.0) then
              call uptriad(ral,hr(point),ubl,ndf,1)
            endif
          endif
        endif

c       Extract local displacements

        do i = 1,ndf
          if(ndfp(i).eq.npart) then
            j = eq(i,n)
            if (j.gt.0) then

c             For active dof compute values from solution
c             where 'du(j)' is increment of 'u' for active dof 'j'.

              ul(i,1) = du(j)*cc1
              ul(i,2) = du(j)*cc2
              ul(i,3) = du(j)

c           Set for restrained dof (not a rigid body)

            elseif(ixt(n).le.0 .and. id(i,n).gt.0) then

c             Set value from input for restrained dof

              ul(i,1) = (f(i,n) - ubl(i))*cc

c             Incremental boundary solution

              if(cc3.ne.0.0d0) then
                ul(i,3) = ul(i,1)*cc3
              else
                ul(i,3) = 0.0d0
                if(ul(i,1).ne.0.0d0) then
                  write(iow,*)' WARNING - infinite acceleration'
                endif
              endif
              ul(i,2) = ul(i,3)*cc2
            else
              ul(i,1) = 0.0d0
              ul(i,2) = 0.0d0
              ul(i,3) = 0.0d0
            endif
          endif
        end do ! i

c       Check for sloping angle boundaries

        if(anglefl) then
          if(hr(np(45)+n-1).ne.0.0d0) then
            call upang(dal,hr(np(45)+n-1),ul,ndf,2)
            if(ndf.ge.6) then
              call upang(ral,hr(np(45)+n-1),ul,ndf,2)
            endif
          endif
        endif

c       Check for Euler angle boundaries

        if(eulerfl) then
          if(hr(np(242)+3*n-3).ne.0.0d0 .or.
     &       hr(np(242)+3*n-2).ne.0.0d0 .or.
     &       hr(np(242)+3*n-1).ne.0.0d0) then
            call upeul(dal,hr(np(242)+3*n-3),ul,ndf,2)
            if(ndf.ge.6) then
              call upeul(ral,hr(np(242)+3*n-3),ul,ndf,2)
            endif
          endif
        endif

c       Check for triad angle boundaries

        if(triadfl) then
          point = np(274) + 9*n - 9
          if(hr(point).ne.-100.0d0) then
            call uptriad(dal,hr(point),ul,ndf,2)
            if(ndf.ge.6) then
              call uptriad(ral,hr(point),ul,ndf,2)
            endif
          endif
        endif

c       Perform update
        do i = 1,ndf
          if(ndfp(i).eq.npart) then
            u(i,n,1) = u(i,n,1) + ul(i,1)
            u(i,n,2) = u(i,n,2) + ul(i,2)
            u(i,n,3) =            ul(i,3)
          endif
        end do ! i
      end do ! n

      end
