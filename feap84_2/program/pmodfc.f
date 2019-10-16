c$Id:$
      subroutine pmodfc(eval,phi,ixt,imf,rlam,f,ftn,mf,ndf,
     &                  ndm,numnp,afd,afl,afu,wf,nmod,mnorm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute residual and tangent for modal force

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'eltran.h'
      include  'modreg.h'

      integer   mf,ndf,ndm,numnp,nmod, i,j,k,n,nn
      integer   ixt(*), imf(ndf,*)
      real*8    eval(*),phi(neqmf,mf)
      real*8    rlam(9,6,*),f(ndf,numnp),ftn(ndf,numnp,*)
      real*8    afd(*),afl(mf,*),afu(mf,*),wf(mf,*),mnorm
      real*8    fmax,y(3),f1(3,3),fh(3,3),lam1(3,3),lamh(3,3),lamn(3,3)

      save

c     Extract rotation matrices

      call quamat(rlam(1,3,modbod(nmod)),lam1)
      call quamat(rlam(1,1,modbod(nmod)),lamn)
      do i = 1,ndm
        do j = 1,ndm
          lamh(i,j) = 0.5d0*(lamn(i,j) + lam1(i,j))
        end do ! j
      end do ! i

c     Set residual from modal terms and diagonal tangent

      do i = 1,mf
        wf(i,3) = -eval(i)*wf(i,6) - wf(i,8)
        afd(i)  =  eval(i)*ctan(1) + ctan(3)
        do j = 1,3
          afu(i,j) = 0.d0
          afl(i,j) = 0.d0
        end do ! j
      end do ! i

      do i = 1,ndf
        do j = 1,numnp
          f(i,j) = 0.5d0*(ftn(i,j,3) + ftn(i,j,4))
        end do ! j
      end do ! i

c     Loop over nodes and construct modal residual loads

      do n = 1,numnp
        if(ixt(n).eq.modbod(nmod)) then

          fmax = 0.d0
          do i = 1,ndm
            fmax = max(fmax,abs(f(i,n)))
          end do ! i

          if(fmax.gt.0.d0) then

c           Modal residual: precompute Lambda * F terms

            call pzero(y,ndm)
            do i = 1,ndm
              do j = 1,ndm
                y(i) = y(i) + lamh(j,i)*f(j,n)
              end do ! j
            end do ! i

c           Modal tangent: precompute F-hat * Lambda terms

            do j = 1,3
              fh(1,j)=(lamh(3,j)*f(2,n)-lamh(2,j)*f(3,n))*ctan(1)
              fh(2,j)=(lamh(1,j)*f(3,n)-lamh(3,j)*f(1,n))*ctan(1)
              fh(3,j)=(lamh(2,j)*f(1,n)-lamh(1,j)*f(2,n))*ctan(1)

              f1(1,j)=(lam1(3,j)*f(2,n)-lam1(2,j)*f(3,n))*ctan(1)
              f1(2,j)=(lam1(1,j)*f(3,n)-lam1(3,j)*f(1,n))*ctan(1)
              f1(3,j)=(lam1(2,j)*f(1,n)-lam1(1,j)*f(2,n))*ctan(1)
            end do ! j

            do i = 1,mf

              do j = 1,3

                nn = imf(j,n)

                if(nn.gt.0) then

c                 Modal residual for loads

                  wf(i,3) = wf(i,3) + phi(nn,i)*y(j)

c                 Modal tangent for loads

                  do k = 1,3
                    afu(i,k) = afu(i,k) + phi(nn,i)*f1(k,j)
                    afl(i,k) = afl(i,k) + phi(nn,i)*fh(k,j)
                  end do ! k
                endif

              end do ! j

            end do ! i

            mnorm = 0.d0
            do i = 1,mf
              mnorm = mnorm + wf(i,3)*wf(i,3)
            end do ! i

          endif

        endif

      end do ! n

      end
