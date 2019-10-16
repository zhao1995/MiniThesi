c$Id:$
      subroutine pload(id,eq,u,f1,dr,prop,flg,afl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add input of load table values                   10/01/2009
c       2. Add argument for ns and nf on call to pldadd     05/03/2009
c       3. Add id(1,2) on call to pldadd                    07/04/2009
c       4. Separate 'id' and 'eq' on argument list          29/04/2009
c       5. Check that ldtab exists before call to pldadd    05/04/2011
c       6. Set periodic values for rank zero only           21/05/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form nodal load vector for current time

c      Inputs:
c         id(*)    - Boundary conditions
c         eq(*)    - Equation numbers for degree of freedom
c         u(*)     - Current solution state
c         prop     - Total proportional load level
c         flg      - Flag: Form residual if true; else reactions
c         afl      - Flag: Assemble tangent for surface loads

c      Outputs:
c         f1(*)    - Total nodal load for t_n+1
c         dr(*)    - Total reaction/residual
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'compas.h'
      include  'complx.h'
      include  'ddata.h'
      include  'fdata.h'
      include  'idptr.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'prld1.h'
      include  'rdat1.h'
      include  'sdata.h'
      include  'setups.h'
      include  'tdato.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   afl,flg
      integer   i,j,n, ipro, id(*), eq(*)
      real*8    dfn, prop,thn, f1(nneq,*),dr(*),u(*),step0,step1

      save

      fl(11) = .false.

c     Set force vectors for t_n+1

      step1 = steps(npart)
      step0 = 1.d0 - step1

      do i = 1,ndf
        if(ndfp(i).eq.npart) then
          do n = i,nneq,ndf

c                   F
            fp(1) = np(27) + n - 1
            fp(2) = fp(1) + nneq
c                   FU
            fp(3) = np(28) + n - 1
            fp(4) = fp(3) + nneq
c                   F0
            fp(5) = fp(4) + nneq
            fp(6) = fp(5) + nneq

c                       FPRO
            ipro = mr(np(29)+n-1)
            if(ipro.eq.0) then
              if(id(n).eq.0) then
                f1(n,1) = hr(fp(1))*prop  + hr(fp(5)) + hr(fp(3))
              elseif(id(n).gt.0) then
                f1(n,1) = hr(fp(2))*prop  + hr(fp(6)) + hr(fp(4))
              endif
              f1(n,3)   = hr(fp(1))*prop + hr(fp(5)) + hr(fp(3))
            else
              if(id(n).eq.0) then
                f1(n,1) = hr(fp(1))*prldv(ipro) + hr(fp(5)) + hr(fp(3))
              elseif(id(n).gt.0) then
                f1(n,1) = hr(fp(2))*prldv(ipro) + hr(fp(6)) + hr(fp(4))
              endif
              f1(n,3) = hr(fp(1))*prldv(ipro) + hr(fp(5)) + hr(fp(3))
            endif
            f1(n,1) = f1(n,2)*step0 + f1(n,1)*step1
            f1(n,3) = f1(n,4)*step0 + f1(n,3)*step1

          end do ! n
        endif
      end do ! i

c     Add loads from ldtab

      if(np(265).ne.0) then
        call pldadd(mr(np(265)),mr(np(266)),hr(np(267)), f1,id,1,2,
     &              prop)
      endif

c     Master-Slave force transformations for couples

      if(np(167).ne.0) then
        call rfclnk(f1,hr(np(43)),mr(np(100)),mr(np(167)))
      endif

c     Set boundary displacements for periodic boundary cases

      if(perflg.and.rank.eq.0) then
        call pperdis(mr(id31),hr(np(43)),f1)
        if(np(257).ne.0) then
          call psetper(mr(np(257)),hr(np(43)),hr(np(40)))
        endif
      endif

c     Initialize residual/reaction

      if(flg) then
        do i = 1,max(neq,nneq)*ipc
          dr(i) = 0.0d0
        end  do ! i
        if(compre) then
          rnorm1 = 0.0d0
          rnormn = 0.0d0
        endif
      endif

c     Form surface type load and tangent load arrays

      call ploadl(eq,mr(np(20+npart)),mr(np(34)),f1,
     &            hr(na),hr(nal),hr(nau),hr(np(35)),hr(np(36)),
     &            hr(np(43)),hr(np(44)),u,hr(np(41)),
     &            prop,ndf,ndm,flg,afl)

c     Compute interpolated load vector

      thn = 1.0d0 - theta(3)

      do i = 1,ndf
        if(ndfp(i).eq.npart) then
          do n = i,nneq,ndf

            j = eq(n)
            if(j.gt.0) then
              if(id(n).eq.0) then
                if(flg) then
                  dfn   = theta(3)*f1(n,1) + thn*f1(n,2)
                  dr(j) = dr(j) + dfn
                  if(compre) then
                    rnorm1 = rnorm1 + abs(dfn)
                    rnormn = rnormn + 1.d0
                  endif
                else
                  dr(n) = theta(3)*f1(n,1) + thn*f1(n,2)
                endif
              endif
            endif

          end do ! n
        endif
      end do ! i

      end
