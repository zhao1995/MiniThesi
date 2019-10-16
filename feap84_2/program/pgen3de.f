c$Id:$
      subroutine pgen3de(x, ix,ix2 ,nen2, nume2,numn2,nseg, prt)

      implicit   none

      include   'bdata.h'
      include   'iofile.h'
      include   'sdata.h'

      logical    prt,genfl
      integer    nen2,nume2,numn2,nseg, ne,nn,nd, n,ns, i,ii
      integer    ix2(nen2,*),ix(nen1,*), iq(8)
      real*8     x(ndm,*),dx(3,3)

c     Output header

      genfl = .true.
      if(prt) then
        write(iow,2000) head,(i,i=1,8)
      endif

c     Set default generation ordering for solid elements

      do i = 1,8
        iq(i) = i
      end do ! i

c     Loop over segments

      ne = 0
      nn = 0
      do ns = 1,nseg
        do n = 1,nume2

c         Count elements and set material number

          ne = ne + 1
          ix(nen1,ne) = ix2(nen2,n)

c         1-node to 2-node element case

          if(ix2(2,n).eq.0) then

            ix(1,   ne) = ix2(1,n) + nn
            ix(2,   ne) = ix2(2,n) + nn + numn2
            nd          = 2

c         2-node to 4-node element case

          elseif(ix2(3,n).eq.0) then

            ix(1,   ne) = ix2(1,n) + nn
            ix(2,   ne) = ix2(2,n) + nn
            ix(3,   ne) = ix2(2,n) + nn + numn2
            ix(4,   ne) = ix2(1,n) + nn + numn2
            nd          = 4

c         4-node to 8-node element case

          else

            if(genfl) then

              genfl = .false.
              do i = 1,3
                dx(i,1) = x(i,ix2(2,n)+nn      ) - x(i,ix2(1,n)+nn)
                dx(i,2) = x(i,ix2(4,n)+nn      ) - x(i,ix2(1,n)+nn)
                dx(i,3) = x(i,ix2(1,n)+nn+numn2) - x(i,ix2(1,n)+nn)
              end do ! i

c             Check jacobian of element

              if( dx(1,1)*(dx(2,2)*dx(3,3)-dx(2,3)*dx(3,2))
     &          + dx(1,2)*(dx(2,3)*dx(3,1)-dx(2,1)*dx(3,3))
     &          + dx(1,3)*(dx(2,1)*dx(3,2)-dx(2,2)*dx(3,1))
     &          .lt. 0.0d0) then
                do i = 1,4
                  iq(i  ) = i+4
                  iq(i+4) = i
                end do ! i
              endif
            endif

            ii = nen2 - 1
            do i = 1,ii
              ix(iq(i   ),ne) = ix2(i,n) + nn
              ix(iq(i+ii),ne) = ix2(i,n) + nn + numn2
            end do ! i
            nd = 8

          endif

c         Output connection list

          if(prt) then
            write(iow,2001) ne,ix(nen1,ne),(ix(i,ne),i=1,nd)
          endif
        end do ! n
        nn = nn + numn2
      end do ! ns

c     Formats

2000  format(/1x,20a4//'   E l e m e n t   C o n n e c t i o n s'//
     &      '    Elem  Mate',8(i3,'-Node'))

2001  format(i8,i6,8i8)

      end
