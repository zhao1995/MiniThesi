c$Id:$
      subroutine stfqtq1 (cp0,cm,ch2,ch1,ch3,tanm,resv,xm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Reset loop order for sequential indexing         29/10/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: Stiffness for node to surface

c      Inputs:
c         cp0(*)  - Table parameters
c         cm(*)   - Material
c         ch2(*)  - History parameters for t_n+1
c         ch1(*)  - History parameters for t_n
c         ch3(*)  - History parameters element
c         xm(*)   - Deformed position of facet nodes

c      Outputs:
c         tanm    - Element tangent
c         resv    - Element residual
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_contac.h'
      include   'c_keyh.h'
      include   'c_mate.h'
      include   'c_pair.h'
      include   'counts.h'
      include   'iofile.h'
      include   'ndata.h'
      include   'comblk.h'

      integer    i,j,k,l,istgn,inew,jnew,istgt,istgnold
      real*8     cp0(nr0,n0c3:*),cm(*),ch1(*),ch2(*),ch3(*)
      real*8     tanm(18,18),resv(3,6),xm(3,4)
      real*8     det,epsn,epst,gn,fn,mue,temp1,temp2,temp3,temp4,temp5
      real*8     a1(3),a2(3),naabl(2,2),mx(2,2),mxinv(2,2)
      real*8     aijinv(2,2),aij(2,2),nvec(3),a1oben(3),a2oben(3)
      real*8     ninew(5),niablnew(5,2),xabab(3,2,2)
      real*8     ni(4),niabl(4,2),niabab(4,2,2),xi(2),new(2,5,3)
      real*8     delxi(2,3,5),delgn(3,5),ddgn(3,5,3,5),xabl(3,2)
      real*8     ttvec(3),ttnorm,dela12,dela11,dela22
      real*8     deltttr(2,3,5),delt1(2,3,5),deldelxi(2,3,5,3,5)
      real*8     ftrial,diffxi(2),tttrial(2),gtold(3)
      real*8     dela1(3,3,5),dela2(3,3,5),delta(3,3)

      save

      data       delta / 1.d0, 3*0.0d0, 1.d0, 3*0.0d0, 1.d0 /

      call cdebug0 ('     stfqtq1',-1)

c     Extract state

      istgn    = nint(ch2(p1(4)))
      istgnold = nint(ch1(p1(4)))

c     Get geometrical variables

      gn       = ch3(p3(9))

      nvec(1)  = ch3(p3(18))
      nvec(2)  = ch3(p3(18)+1)
      nvec(3)  = ch3(p3(18)+2)

      a1(1)    = ch2(p1(19)  )
      a1(2)    = ch2(p1(19)+1)
      a1(3)    = ch2(p1(19)+2)

      a2(1)    = ch3(p3(20)  )
      a2(2)    = ch3(p3(20)+1)
      a2(3)    = ch3(p3(20)+2)

      aij(1,1) = a1(1)*a1(1) + a1(2)*a1(2) + a1(3)*a1(3)
      aij(1,2) = a1(1)*a2(1) + a1(2)*a2(2) + a1(3)*a2(3)
      aij(2,1) = aij(1,2)
      aij(2,2) = a2(1)*a2(1) + a2(2)*a2(2) + a2(3)*a2(3)

      xi(1)    = ch2(p1(24)  )
      xi(2)    = ch2(p1(24)+1)

      call formfkt(xi,ni,niabl,niabab)

      do j=1,2
        do i=1,3
          xabl(i,j) = niabl(1,j)*xm(i,1)
     &              + niabl(2,j)*xm(i,2)
     &              + niabl(3,j)*xm(i,3)
     &              + niabl(4,j)*xm(i,4)
          do k=1,2
            xabab(i,j,k) = niabab(1,j,k)*xm(i,1)
     &                   + niabab(2,j,k)*xm(i,2)
     &                   + niabab(3,j,k)*xm(i,3)
     &                   + niabab(4,j,k)*xm(i,4)
          end do ! k
        end do ! i
      end do ! j

c     Compute var. of xi

      do j=1,2
        do i=1,2
          naabl(i,j) = nvec(1)*xabab(1,i,j)
     &               + nvec(2)*xabab(2,i,j)
     &               + nvec(3)*xabab(3,i,j)
        end do ! i
      end do ! j

      do j = 1,2
        do i = 1,2
          mx(i,j) = aij(i,j) - gn*naabl(i,j)
        end do ! i
      end do ! j

      det           =  1.d0/(mx(1,1)*mx(2,2)-mx(1,2)*mx(2,1))
      mxinv(1,1)    =  mx(2,2)*det
      mxinv(2,2)    =  mx(1,1)*det
      mxinv(1,2)    = -mx(1,2)*det
      mxinv(2,1)    = -mx(2,1)*det

      det           =  1.d0/(aij(1,1)*aij(2,2)-aij(1,2)*aij(2,1))
      aijinv(1,1)   =  aij(2,2)*det
      aijinv(2,2)   =  aij(1,1)*det
      aijinv(1,2)   = -aij(1,2)*det
      aijinv(2,1)   = -aij(2,1)*det

      ninew(1)      =  1.d0
      niablnew(1,1) =  0.d0
      niablnew(1,2) =  0.d0

      do i=1,4
        ninew(i+1)      = -ni(i)
        niablnew(i+1,1) =  niabl(i,1)
        niablnew(i+1,2) =  niabl(i,2)
      end do

      do j=1,5
        do i=1,3
          delgn(i,j)   = nvec(i)*ninew(j)
          delxi(1,i,j) =
     &             mxinv(1,1)*(a1(i)*ninew(j)+gn*nvec(i)*niablnew(j,1))
     &           + mxinv(1,2)*(a2(i)*ninew(j)+gn*nvec(i)*niablnew(j,2))
          delxi(2,i,j) =
     &             mxinv(2,1)*(a1(i)*ninew(j)+gn*nvec(i)*niablnew(j,1))
     &           + mxinv(2,2)*(a2(i)*ninew(j)+gn*nvec(i)*niablnew(j,2))
        end do ! i
      end do ! j

      do k=1,3
        do j=1,5
          do i=1,2
            new(i,j,k) = niablnew(j,i)*nvec(k)
     &                 + naabl(i,1)*delxi(1,k,j)
     &                 + naabl(i,2)*delxi(2,k,j)
          end do ! i
        end do ! j
      end do ! k

      do i=1,3
        a1oben(i) = aijinv(1,1)*a1(i) + aijinv(1,2)*a2(i)
        a2oben(i) = aijinv(2,1)*a1(i) + aijinv(2,2)*a2(i)
      end do ! i

      do l=1,5
        do k=1,5
          do j=1,3
            do i=1,3
              ddgn(i,k,j,l)=
     &           - ninew(l)*(a1oben(j)*new(1,k,i)
     &                     + a2oben(j)*new(2,k,i))
     &                     -(delxi(1,i,k)*niablnew(l,1)
     &                     + delxi(2,i,k)*niablnew(l,2))*nvec(j)
            end do ! i
          end do ! j
        end do ! k
      end do ! l

c     Get penalty for normal stiffness

      epsn = cp0(3,2)
      fn   = epsn*gn
      if(ifsolm.eq.2) then                   ! Lagrange multiplier case
        fn        =  fn + ch2(p1(31))
        resv(1,6) = -gn*ch3(p3(16))
        do i = 1,3
          do j = 1,5
            inew          = 3*(j-1)+i
            tanm(inew,16) = delgn(i,j)*ch3(p3(16))
            tanm(16,inew) = tanm(inew,16)
          end do ! j
        end do ! i
      endif
      ch3(p3(14)) = fn

c     Multiply by area factor to form stiffness

      temp1 = ch3(p3(16))*epsn
      temp2 = ch3(p3(16))*fn
      do i = 1,3
        do j = 1,5
          inew = 3*(j-1)+i
          do k = 1,3
            do l = 1,5
              jnew = 3*(l-1)+k
              tanm(inew,jnew) = temp1*delgn(i,j)*delgn(k,l)
     &                        + temp2*ddgn(i,j,k,l)
            end do ! l
          end do ! k
        end do ! j
      end do ! i

c     Form residual

      do j=1,5
        do i = 1,3
          resv(i,j) = -temp2*delgn(i,j)
        end do ! i
      end do ! j

      if (iffric.eq.1) then

        mue         = cm(1)
        istgt       = nint(ch2(p1(3)))
        ch3(p3( 5)) = 0.0d0
        ch3(p3(15)) = 0.0d0

        if (istgt.eq.1) then

          if(istgnold.ne.0) then

            do i = 1,3
              gtold(i) = ch1(p1(26)+i-1)
            end do

            diffxi(1)  =  ch2(p1(10))
     &                 +  a1oben(1)*gtold(1)
     &                 +  a1oben(2)*gtold(2)
     &                 +  a1oben(2)*gtold(3)
            diffxi(2)  =  ch2(p1(11))
     &                 +  a2oben(1)*gtold(1)
     &                 +  a2oben(2)*gtold(2)
     &                 +  a2oben(2)*gtold(3)

            epsn       =  cp0(3,2)
            epst       =  cp0(4,2)

            tttrial(1) = -epst*(diffxi(1)*aij(1,1) + diffxi(2)*aij(2,1))
            tttrial(2) = -epst*(diffxi(1)*aij(1,2) + diffxi(2)*aij(2,2))

            do i = 1,3
              do j = 1,5
                do k = 1,3
                  dela1(k,i,j) = delta(k,i)*niablnew(j,1)
     &                         + xabab(k,1,1)*delxi(1,i,j)
     &                         + xabab(k,1,2)*delxi(2,i,j)
                  dela2(k,i,j) = delta(k,i)*niablnew(j,2)
     &                         + xabab(k,2,1)*delxi(1,i,j)
     &                         + xabab(k,2,2)*delxi(2,i,j)
                end do ! k

                dela11         = (dela1(1,i,j)*a1(1)
     &                         +  dela1(2,i,j)*a1(2)
     &                         +  dela1(3,i,j)*a1(3))*2.d0
                dela22         = (dela2(1,i,j)*a2(1)
     &                         +  dela2(2,i,j)*a2(2)
     &                         +  dela2(3,i,j)*a2(3))*2.d0

                dela12         =  dela1(1,i,j)*a2(1)
     &                         +  dela1(2,i,j)*a2(2)
     &                         +  dela1(3,i,j)*a2(3)
     &                         +  dela2(1,i,j)*a1(1)
     &                         +  dela2(2,i,j)*a1(2)
     &                         +  dela2(3,i,j)*a1(3)

                deltttr(1,i,j) = -epst*(delxi(1,i,j)*aij(1,1)
     &                                + delxi(2,i,j)*aij(2,1)
     &                                + diffxi(1)*dela11
     &                                + diffxi(2)*dela12)

                deltttr(2,i,j) = -epst*(delxi(1,i,j)*aij(1,2)
     &                                + delxi(2,i,j)*aij(2,2)
     &                                + diffxi(1)*dela12
     &                                + diffxi(2)*dela22)
              end do ! j
            end do ! i

            do i=1,3
              ttvec(i) = tttrial(1)*a1oben(i) + tttrial(2)*a2oben(i)
            end do ! i

            ttnorm = sqrt(ttvec(1)**2 + ttvec(2)**2 + ttvec(3)**2)

            ftrial = ttnorm - mue*abs(epsn*gn)

            do i=1,3
              do j=1,5
                do k=1,3
                  do l=1,5

                    temp1 =- a1(1)*(dela1(1,i,j)*delxi(1,k,l)
     &                            + dela2(1,i,j)*delxi(2,k,l))
     &                     - a1(2)*(dela1(2,i,j)*delxi(1,k,l)
     &                            + dela2(2,i,j)*delxi(2,k,l))
     &                     - a1(3)*(dela1(3,i,j)*delxi(1,k,l)
     &                            + dela2(3,i,j)*delxi(2,k,l))
     &                     - a1(k)*(niablnew(l,1)*delxi(1,i,j)
     &                            + niablnew(l,2)*delxi(2,i,j))
     &                     + dela1(1,k,l)*(ninew(j)*delta(i,1)
     &                                   - a1(1)*delxi(1,i,j)
     &                                   - a2(1)*delxi(2,i,j))
     &                     + dela1(2,k,l)*(ninew(j)*delta(i,2)
     &                                   - a1(2)*delxi(1,i,j)
     &                                   - a2(2)*delxi(2,i,j))
     &                     + dela1(3,k,l)*(ninew(j)*delta(i,3)
     &                                   - a1(3)*delxi(1,i,j)
     &                                   - a2(3)*delxi(2,i,j))
     &                     + dela1(1,i,j)*(ninew(l)*delta(k,1)
     &                                   - a1(1)*delxi(1,k,l)
     &                                   - a2(1)*delxi(2,k,l))
     &                     + dela1(2,i,j)*(ninew(l)*delta(k,2)
     &                                   - a1(2)*delxi(1,k,l)
     &                                   - a2(2)*delxi(2,k,l))
     &                     + dela1(3,i,j)*(ninew(l)*delta(k,3)
     &                                   - a1(3)*delxi(1,k,l)
     &                                   - a2(3)*delxi(2,k,l))
     &                     + gn*(nvec(k)*(niabab(l-1,1,1)*delxi(1,i,j)
     &                                  + niabab(l-1,1,2)*delxi(2,i,j))
     &                         + nvec(i)*(niabab(j-1,1,1)*delxi(1,k,l)
     &                                  + niabab(j-1,1,2)*delxi(2,k,l)))

                    temp2 =- a2(1)*(dela1(1,i,j)*delxi(1,k,l)
     &                            + dela2(1,i,j)*delxi(2,k,l))
     &                     - a2(2)*(dela1(2,i,j)*delxi(1,k,l)
     &                            + dela2(2,i,j)*delxi(2,k,l))
     &                     - a2(3)*(dela1(3,i,j)*delxi(1,k,l)
     &                            + dela2(3,i,j)*delxi(2,k,l))
     &                     - a2(k)*(niablnew(l,1)*delxi(1,i,j)
     &                            + niablnew(l,2)*delxi(2,i,j))
     &                     + dela2(1,k,l)*(ninew(j)*delta(i,1)
     &                                   - a1(1)*delxi(1,i,j)
     &                                   - a2(1)*delxi(2,i,j))
     &                     + dela2(2,k,l)*(ninew(j)*delta(i,2)
     &                                   - a1(2)*delxi(1,i,j)
     &                                   - a2(2)*delxi(2,i,j))
     &                     + dela2(3,k,l)*(ninew(j)*delta(i,3)
     &                                   - a1(3)*delxi(1,i,j)
     &                                   - a2(3)*delxi(2,i,j))
     &                     + dela2(1,i,j)*(ninew(l)*delta(k,1)
     &                                   - a1(1)*delxi(1,k,l)
     &                                   - a2(1)*delxi(2,k,l))
     &                     + dela2(2,i,j)*(ninew(l)*delta(k,2)
     &                                   - a1(2)*delxi(1,k,l)
     &                                   - a2(2)*delxi(2,k,l))
     &                     + dela2(3,i,j)*(ninew(l)*delta(k,3)
     &                                   - a1(3)*delxi(1,k,l)
     &                                   - a2(3)*delxi(2,k,l))
     &                     + gn*(nvec(k)*(niabab(l-1,2,1)*delxi(1,i,j)
     &                                  + niabab(l-1,2,2)*delxi(2,i,j))
     &                         + nvec(i)*(niabab(j-1,2,1)*delxi(1,k,l)
     &                                  + niabab(j-1,2,2)*delxi(2,k,l)))

                    deldelxi(1,i,j,k,l) = mxinv(1,1)*temp1
     &                                  + mxinv(1,2)*temp2

                    deldelxi(2,i,j,k,l) = mxinv(2,1)*temp1
     &                                  + mxinv(2,2)*temp2

                  end do ! l
                end do ! k
              end do ! j
            end do ! i

            if (ftrial.le.0.d0) then
            else
              temp1 = mue*epsn/ttnorm
              temp2 = mue*abs(epsn*gn)/ttnorm
              temp3 = tttrial(1)*aijinv(1,1) + tttrial(2)*aijinv(1,2)
              temp4 = tttrial(1)*aijinv(2,1) + tttrial(2)*aijinv(2,2)
              do k=1,2
                temp5 = tttrial(k)/(ttnorm*ttnorm)
                do i=1,3
                  do j=1,5
                    delt1(k,i,j)=-temp1*delgn(i,j)*tttrial(k)
     &                          + temp2*(deltttr(k,i,j)
     &                          + temp5*(( ttvec(1)*dela1(1,i,j)
     &                                   + ttvec(2)*dela1(2,i,j)
     &                                   + ttvec(3)*dela1(3,i,j)
     &                                   - deltttr(1,i,j))*temp3
     &                                 + ( ttvec(1)*dela2(1,i,j)
     &                                   + ttvec(2)*dela2(2,i,j)
     &                                   + ttvec(3)*dela2(3,i,j)
     &                                   - deltttr(2,i,j))*temp4))
                  end do ! j
                end do ! i
              end do ! k

              do j = 1,5
                do i = 1,3
                  do k = 1,2
                    deltttr(k,i,j) = delt1(k,i,j)
                  end do ! k
                end do ! i
              end do ! j

              do i = 1,2
                tttrial(i) = mue*abs(epsn*gn)*tttrial(i)/ttnorm
              end do ! i
            endif

            do i = 1,3
              do j = 1,5
                inew = 3*(j-1) + i
                do k = 1,3
                  do l = 1,5
                    jnew = 3*(l-1) + k
                    tanm(inew,jnew) = tanm(inew,jnew) - ch3(p3(16))
     &                              * (deltttr(1,k,l)*delxi(1,i,j)
     &                              +  deltttr(2,k,l)*delxi(2,i,j)
     &                              +  tttrial(1)*deldelxi(1,i,j,k,l)
     &                              +  tttrial(2)*deldelxi(2,i,j,k,l))

                  end do ! l
                end do ! k
              end do ! j
            end do ! i

            do j = 1,5
              do i = 1,3
                resv(i,j) = resv(i,j) + ch3(p3(16))
     &                                * (tttrial(1)*delxi(1,i,j)
     &                                +  tttrial(2)*delxi(2,i,j))
              end do ! i
            end do ! j

            do i=1,3
              ttvec(i) = tttrial(1)*a1oben(i) + tttrial(2)*a2oben(i)
            end do ! i

            ch3(p3( 5)) = ftrial
            ch3(p3(15)) = sqrt(ttvec(1)**2 + ttvec(2)**2 + ttvec(3)**2)

          elseif(istgnold.eq.0) then

            ch3(p3( 5)) = 0.0d0
            ch3(p3(15)) = 0.0d0

          endif

        endif                  !istgt.eq.1
      endif                    !mue.ge.0.d0

      end
