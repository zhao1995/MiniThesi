c$Id:$
      subroutine gnqtq (kset,x,u,ix1,ix2,xs,ch1,ch2,ch3,ids,im,
     &                  surpoin,mue,test)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: Compute gap etc. with quadrilateral segments (3d)

c      Inputs :
c         kset    - Node number
c         x(*)    - Nodal coordinates
c         u(*)    - Current nodal solution vectors
c         ix1(*)  - Element nodal connection list for surface 2
c         ix2(*)  - Element nodal connection list for surface 2
c         xs(*)   - Coordinates of slave node S
c         ch1(*)  - Contact history variables (t_n)
c         ids(3)  - Equation/BC indicators for slave node
c         surpoin - Surface points
c         mue     - Friction coefficient
c         test    - Flag

c      Outputs:
c         im(4)   - Master setment nodes
c         ch2(*)  - Contact history variables (t_n+1)
c         ch3(*)  - Contact history variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_geom.h'
      include   'c_keyh.h'
      include   'c_pair.h'
      include   'c_tole.h'
      include   'sdata.h'
      include   'tdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    test
      integer    kset,nod2,masts,istgn,istgt,im(4),mastsold
      integer    iter,imax,nrold,nrneben,istgnold, i,nr,nodes
      integer    ix1(dnope1,*),ix2(dnope2,*),surpoin(*),ids(3)
      real*8     gapn,gn,nbar,xi1,xi2,det,a11,a12,a22,aslave
      real*8     temp1,temp2, a1(3),a2(3),aij(2,2),xm(3,4)
      real*8     ni(4),niab(4,2),niabab(4,2,2),n(3),nvec(3),g(3)
      real*8     t1(3),t2(3),xp(3),xi(2),xold(3),xnew(3)
      real*8     x(ndm,*),u(ndf,*),xs(3),ch1(*),ch2(*),ch3(*),mue

      save

      data       imax / 25 /

      call cdebug0 ('    gnqtq',-1)

      call as(kset,x,u,ix1,surpoin,aslave, nodes)

      ch3(p3(16)) = aslave
      istgn       = nint(ch2(p1(4)))

      if(istgn.eq.0) then

        if(test) then
          ch3(p3(2)) = 1.0d0
        endif
        test = .false.
        return

      elseif (istgn.eq.1) then

        masts = nint(ch2(p1 (1)  ))
        xi1   = ch2(p1(24)  )
        xi2   = ch2(p1(24)+1)

        do nod2 = 1,nope2
          im(nod2) = ix2(nod2,masts)
          do i = 1,3
            xm(i,nod2) = x(i,im(nod2)) + u(i,im(nod2))
          end do ! i
        end do ! nod2

c       Perform constrained closest point projection

        call cnproj( ids, xs,xm, xi,xp,t1,t2,g, iter)

        if(ifdb .and. iter.ge.imax) then
          write(*,*) ' No convergence in GNQTQ, slave node =', nodes
        endif

c       Compute normal vector to facet (times twice area)

        n(1) = t1(2)*t2(3) - t1(3)*t2(2)
        n(2) = t1(3)*t2(1) - t1(1)*t2(3)
        n(3) = t1(1)*t2(2) - t1(2)*t2(1)
        gapn = g(1)*n(1) + g(2)*n(2) + g(3)*n(3)
        nbar = 1.0d0/( n(1)*n(1) + n(2)*n(2) + n(3)*n(3) )

c       If gap negative check postion

c       if( gapn*sqrt(nbar) .le. tlopen ) then

c         CONTACT!!!

          aij(1,1) = t1(1)*t1(1) + t1(2)*t1(2) + t1(3)*t1(3)
          aij(1,2) = t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)
          aij(2,2) = t2(1)*t2(1) + t2(2)*t2(2) + t2(3)*t2(3)

          if (ifdb .and. max(abs(xi(1)),abs(xi(2))).gt.1.d0+tlouts) then
            write( *,3000) kset,masts,ttim
            write(99,3000) kset,masts,ttim
          endif

c         Contact point in segment

          nbar =  sqrt(nbar)
          gapn =  gapn*nbar
c         gn   = -abs(gapn)
          gn   =  gapn
          a11  = aij(1,1)
          a12  = aij(1,2)
          a22  = aij(2,2)
          xi1  = xi(1)
          xi2  = xi(2)
          do i = 1,3
            nvec(i) = nbar*n(i)
            a1(i)   = t1(i)
            a2(i)   = t2(i)
            xnew(i) = xp(i)
          end do ! i
c       else
c         gn = abs(gapn)
c       end if

c       Set history variables

        ch3(p3( 9)  ) = gn
        ch3(p3(18)  ) = nvec(1)
        ch3(p3(18)+1) = nvec(2)
        ch3(p3(18)+2) = nvec(3)
        ch2(p1(19)  ) = a1(1)
        ch2(p1(19)+1) = a1(2)
        ch2(p1(19)+2) = a1(3)
        ch3(p3(20)  ) = a2(1)
        ch3(p3(20)+1) = a2(2)
        ch3(p3(20)+2) = a2(3)
        ch2(p1(24)  ) = xi1
        ch2(p1(24)+1) = xi2

        if(ifdb .and. indb.ge.1) then
          write(99,3001) kset,masts,gn,xi1,xi2
        endif

      elseif((istgn.eq.2).or.(istgn.eq.3))then

        nr      = nint(ch2(p1( 1)  ))
        nrneben = nint(ch2(p1(25)))

        do i = 1,3
          xm(i,1) = x(i,nr)      + u(i,nr)
          xm(i,2) = x(i,nrneben) + u(i,nrneben)
          a1(i)   = xm(i,2)      - xm(i,1)
          a2(i)   = xs(i)        - xm(i,1)
        end do ! i

c       Projection of slave node

        a11 =  a1(1)*a1(1) + a1(2)*a1(2) + a1(3)*a1(3)
        xi1 = (a1(1)*a2(1) + a1(2)*a2(2) + a1(3)*a2(3))/a11
        do i = 1,3
          nvec(i) = a2(i)   - xi1*a1(i)
          xnew(i) = xm(i,1) + xi1*a1(i)
        end do ! i
        gn =  sqrt(nvec(1)*nvec(1) + nvec(2)*nvec(2) + nvec(3)*nvec(3))
        do i = 1,3
          nvec(i) =  nvec(i)/gn
        end do ! i

c       Set history variables

        ch3(p3( 9)  ) = gn
        ch3(p3(18)  ) = nvec(1)
        ch3(p3(18)+1) = nvec(2)
        ch3(p3(18)+2) = nvec(3)
        ch2(p1(19)  ) = a1(1)
        ch2(p1(19)+1) = a1(2)
        ch2(p1(19)+2) = a1(3)
        ch2(p1(24)  ) = xi1
        if(ifdb .and. indb.ge.1) then
          write(99,3002) kset,nr,nrneben,gn,xi1
        endif

      elseif(istgn.eq.4)then

        nr = nint(ch2(p1(1)))
        do i = 1,3
          xnew(i) = x(i,nr) + u(i,nr)
          nvec(i) = xs(i)  - xnew(i)
        end do ! i
        gn = -sqrt(nvec(1)*nvec(1) + nvec(2)*nvec(2) + nvec(3)*nvec(3))
        do i = 1,3
          nvec(i) = -nvec(i)/gn
        end do ! i

c       Set history variables

        ch3(p3( 9)  ) = gn
        ch3(p3(18)  ) = nvec(1)
        ch3(p3(18)+1) = nvec(2)
        ch3(p3(18)+2) = nvec(3)

        if(ifdb .and. indb.ge.1) then
          write(99,3004) kset,nr,gn
        endif

      endif

c     if(gn.gt.tlopen) then
c       gn   =  0.d0
c       if(test) ch3(p3(2)) = 1.0d0
c       test = .false.
c       if(ifdb) write(*,*)' No longer have contact with segment',masts
c     else
        if(test) ch3(p3(2)) = 0.0d0
c       test = .true.
c     endif

c     Friction case

      if (mue.ne.0.d0) then

        istgt    = 0
        istgnold = nint(ch1(p1(4)))

        if(istgnold.eq.1) then

          mastsold = nint(ch1(p1(1)))
          do nod2 = 1,4
            im(nod2) = ix2(nod2,mastsold)
            do i = 1,3
              xm(i,nod2) = x(i,im(nod2)) + u(i,im(nod2))
            end do ! i
          end do ! nod2

          xi(1) = ch1(p1(24)  )
          xi(2) = ch1(p1(24)+1)

          call formfkt(xi,ni,niab,niabab)

          do i = 1,3
            xold(i) = ni(1)*xm(i,1) + ni(2)*xm(i,2)
     &              + ni(3)*xm(i,3) + ni(4)*xm(i,4)
          end do ! i

        elseif((istgnold.eq.2).or.(istgnold.eq.3)) then

          nrold      = nint(ch1(p1(1)))
          do i = 1,3
            xm(i,1) = x(i,nrold) + u(i,nrold)
          end do ! i
          do i = 1,3
            xold(i) = xm(i,1) + ch1(p1(24))*ch1(p1(19)+i-1)
          end do ! i

        elseif(istgnold.eq.4) then

          nrold = nint(ch1(p1(1)))
          do i = 1,3
            xold(i) = x(i,nrold) + u(i,nrold)
          end do ! i
        endif

c       TANGENTIAL GAP => TANGENTIAL CONTACT

        if (istgn.eq.1) then

          if (istgnold.eq.1) then

            if((xi1.ne.ch1(p1(24))).and.(xi2.ne.ch1(p1(24)+1)))then

              istgt = 1
              if (masts.ne.ch1(p1(1))) then

                det         =  a11*a22 - a12*a12
                temp1       =  xnew(1)*a1(1)
     &                      +  xnew(2)*a1(2)
     &                      +  xnew(3)*a1(3)
     &                      -  xold(1)*a1(1)
     &                      -  xold(2)*a1(2)
     &                      -  xold(3)*a1(3)
                temp2       =  xnew(1)*a2(1)
     &                      +  xnew(2)*a2(2)
     &                      +  xnew(3)*a2(3)
     &                      -  xold(1)*a2(1)
     &                      -  xold(2)*a2(2)
     &                      -  xold(3)*a2(3)
                ch2(p1(10)) = (a22*temp1 - a12*temp2)/det
                ch2(p1(11)) = (a11*temp2 + a12*temp1)/det

c             Same master segment as in last time step

              else

                ch2(p1(10)) = xi1 - ch1(p1(24)  )
                ch2(p1(11)) = xi2 - ch1(p1(24)+1)

              endif
            endif

          elseif((istgnold.eq.2).or.(istgnold.eq.3)
     &                          .or.(istgnold.eq.4)) then

            det         =  a11*a22 - a12*a12
            temp1       =  xnew(1)*a1(1)
     &                  +  xnew(2)*a1(2)
     &                  +  xnew(3)*a1(3)
     &                  - xold(1)*a1(1)
     &                  -  xold(2)*a1(2)
     &                  -  xold(3)*a1(3)
            temp2       =  xnew(1)*a2(1)
     &                  +  xnew(2)*a2(2)
     &                  +  xnew(3)*a2(3)
     &                  -  xold(1)*a2(1)
     &                  -  xold(2)*a2(2)
     &                  -  xold(3)*a2(3)
            ch2(p1(10)) = (a22*temp1 - a12*temp2)/det
            ch2(p1(11)) = (a11*temp2 + a12*temp1)/det
            istgt       = 1

          elseif(istgnold.eq.0) then

            istgt = 0

          endif

        elseif((istgn.eq.2).or.(istgn.eq.3))then

          if(istgnold.eq.1) then

c           gt(i) = xnew(i)-xold(i)

            do i = 1,3
              ch2(p1(26)+i-1) = xnew(i) - xold(i)
            end do ! i
            istgt = 1

          elseif((istgnold.eq.2).or.(istgnold.eq.3)) then

            if((nr.eq.ch1(p1(1))).and.(nrneben.eq.ch1(p1(25)))) then

              ch2(p1(10)) = xi1 - ch1(p1(24))
              istgt       = 1
            else

c             gt(i) = xnew(i)-xold(i)

              do i = 1,3
                ch2(p1(26)+i-1) = xnew(i) - xold(i)
              end do ! i
              istgt = 1
            endif

          elseif(istgnold.eq.4) then

c           gt(i) = xnew(i)-xold(i)

            do i = 1,3
              ch2(p1(26)+i-1) = xnew(i) - xold(i)
            end do ! i
            istgt = 1

          elseif(istgnold.eq.0) then

            istgt = 0

          endif

        elseif(istgn.eq.4) then

          if(istgnold.eq.1) then

            do i = 1,3
              ch2(p1(26)+i-1) = xnew(i) - xold(i)
            end do ! i
            istgt = 1

          elseif((istgnold.eq.2).or.(istgnold.eq.3)) then

            do i = 1,3
              ch2(p1(26)+i-1) = xnew(i) - xold(i)
            end do ! i
            istgt = 1

          elseif(istgnold.eq.4) then

            if (nr.eq.ch1(p1(1))) then
              istgt = 0
            else
              do i = 1,3
                ch2(p1(26)+i-1) = xnew(i) - xold(i)
              end do ! i
              istgt = 1
            endif

          elseif(istgnold.eq.0) then

            istgt = 0

          endif

        elseif(istgn.eq.0) then

          istgt = 0

        endif

c       Save sliding type

        ch2(p1(3)) = istgt

      endif

c     Format

3000  format(' Slave node',i5,' out of master segment',i5,
     &       ': Time =',1p,1e12.4)

3001  format(' ISTGN =1 SLAVE =',i5,' MASTS =',i5/
     &       ' GN,XI ',1p,3e12.4)

3002  format(' ISTGN =2 SLAVE =',i5,' NR/NEBEN =',2i5/
     &       ' GN,XI ',1p,2e12.4)

3004  format(' ISTGN =4 SLAVE =',i5,' NR =',i5/
     &       ' GN ',1p,1e12.4)

      end
