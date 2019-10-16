c$Id:$
      subroutine psolve(stype,b,fp,factor,solve,cfr, prnt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Solver driver
c      Inputs:
c         stype  - Solver matrix type: 0 = diagonal
c                                     <0 = direct
c                                     >0 = iterative
c         b(*)   - Vector for solutions
c         fp(4)  - Pointers for arrays
c         factor - Factor option
c         solve  - Solve option
c         cfr    - Unsymmtric flag
c         prnt   - Print data if true

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'cdata.h'
      include   'compas.h'
      include   'endata.h'
      include   'eqsym.h'
      include   'fdata.h'
      include   'iofile.h'
      include   'ndata.h'
      include   'part0.h'
      include   'part7.h'
      include   'pathn.h'
      include   'pscal.h'
      include   'rdata.h'
      include   'rdat0.h'
      include   'rdat1.h'
      include   'rjoint.h'
      include   'setups.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'p_int.h'

      logical    factor,solve,cfr,prnt, setvar,palloc,flags(5)
      integer    stype, i
      real*4     etime, tt,tary(2)
      real*8     rn, dot, b(*), rsdi

      save

      if(solver) then

c       Factor tangent matrix

        if(factor) then
          if(stype.le.-1) then
            tt    = etime(tary)
            tdiff = tary(1)

c           Call out-of-core or in-core direct solver

            if(stype.eq.-1) then

c             Store AU and AL on disk in block form

              call xnumb( mr(fp(4)), neq, nnr-neq, maxbl )
              setvar = palloc( 92,'OINB',3*maxbl,1)

              if(scale(npart)) then
                call pscalc(mr(np(94)),mr(np(93)),hr(fp(1)),hr(fp(3)),
     &                      hr(fp(2)),hr(np(234+npart)),neq, cfr)
              endif

              call xsave(mr(np(94)),mr(np(93)),hr(fp(2)),hr(fp(3)),
     &                   fau,fal,mr(fp(4)), iuau, iual, neq, nnr-neq,
     &                   mr(np(92)), maxbl, cfr)

c             Factor tangent matrix

              call xdatri(fau,fal,iuau,iual,mr(fp(4)),neq,hr(fp(1)),
     &                    hr(fp(2)),hr(fp(3)),mr(np(92)),maxbl,cfr)

            else
              if(scale(npart)) then
c               Sparse solver
                if(stype.eq.-2) then
                  call pscals(mr(np(93)+neq),mr(np(94)),hr(fp(1)),
     &                        hr(np(234+npart)),neq)
c               Profile solver
                elseif(stype.eq.-3) then
                  call pscala(hr(fp(1)),hr(fp(3)),hr(fp(2)),mr(fp(4)),
     &                        hr(np(234+npart)),neqs,neq)
                endif
              endif
              if(neqr.lt.neq) then
                call datrim(hr(fp(3)),hr(fp(2)),hr(fp(1)),
     &                      mr(fp(4)),neqr,neqs,neq)
              else
                call datri (hr(fp(3)),hr(fp(2)),hr(fp(1)),
     &                      mr(fp(4)),neqs,neq)
              endif
            endif

            tt    = etime(tary)
            tdiff = tary(1) - tdiff  ! save for timing solutions
            if(prnt) then
              write(iow,2001) tary
              if(ior.lt.0) then
                write(*,2001) tary
              endif
            endif
          endif
          if(stype.ne.0) then
            fl(4) = .false.  ! Tangent formed
            if(monofl) then
              flp(4,5) = .false.
            endif
          endif
        endif

c       Solve equations

        if(solve) then

c         Diagonal solution

          if(stype.eq.0) then

            aengy = 0.0d0
            do i = 0,neq-1
              if(hr(fp(1)+i) .ne. 0.0d0) then
                rsdi   = b(i+1)
                b(i+1) = b(i+1)/hr(fp(1)+i)
                aengy  = aengy + rsdi*b(i+1)
              endif
            end do ! i

c         Direct solution

          elseif(stype.le.-1) then

            if(neqr.lt.neq) then
              call dasolm(hr(fp(3)),hr(fp(2)),hr(fp(1)),b,
     &                    mr(fp(4)),neqr,neqs,neq,aengy,scale(npart))
            else
              call dasol (hr(fp(3)),hr(fp(2)),hr(fp(1)),b,
     &                    mr(fp(4)),neqs,neq,aengy,scale(npart))
            endif

c         Conjugate gradient solver

          else

            if(stype.ge.2) then
              fp(3) = np(68)
              fp(4) = np(80)
            else
              fp(3) = 1
              fp(4) = 1
            endif

            setvar = palloc(111,'TEMP1',5*neq,2) ! temporary storage
            fp(5) = np(111)
            do i = 6,9
              fp(i) = fp(i-1) + neq
            end do ! i

            if(prnt .and. ior.lt.0) write(*,*) 'START CG: SOLVER'

            call conjgd(hr(fp(1)),hr(fp(2)),hr(fp(5)),hr(fp(3)),
     &                  mr(fp(4)),b,mr(np(93)),mr(np(94)),hr(fp(6)),
     &                  hr(fp(7)),hr(fp(8)),hr(fp(9)),
     &                  neq,icgits,itol,rn,rn0,stype)

            aengy  = dot(b,hr(fp(9)),neq)
            setvar = palloc(111,'TEMP1', 0, 2) ! destroy storage

            if(prnt) then
              tt    = etime(tary)
              write(iow,2002) tary
              if(ior.lt.0) then
                write(*,2002) tary
              endif
            endif

          endif
        endif

c     User supplied solver routine

      else
        flags(1) = .false.
        flags(2) =  factor
        flags(3) =  cfr
        flags(4) =  solve
        flags(5) = .false.
        call usolve(flags,b)
        if(factor) then
          fl(4) = .false.  ! Tangent formed
        endif
      endif

c     Formats

2001  format('   End Triangular Decomposition',28x,'t=',0p,2f9.2)
2002  format('   End Conjugate Gradient Solution',25x,'t=',0p,2f9.2)

      end
