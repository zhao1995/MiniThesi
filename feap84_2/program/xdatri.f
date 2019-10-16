c$Id:$
      subroutine xdatri(fau, fal, iunau, iunal, jp, neq,
     &                  ad, au, al, iblk, maxbl, flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute triangular factors of matrix stored in
c               profile form in blocks on disk.

c      Inputs:
c         fau       - Name of file containing upper blocks
c         fal       - Name of file containing lower blocks
c         iunau     - Logical unit number for upper blocks
c         iunal     - Logical unit number for lower blocks
c         jp(*)     - Pointer to end of row/columns in profile
c         neq       - Number of equations
c         ad(*)     - Diagonals of matrix
c         au(*)     - Upper profile of matrix
c         al(*)     - Lower profile of matrix
c         iblk(3,*) - Block information
c                     iblk(1,I) = first equation in block I
c                     iblk(2,I) = last  equation in block I
c                     iblk(3,I) = Number of block needed to start
c         maxblk    - Number of blocks to store profile
c         flg       - if .true. equations are unsymmetric

c      Outputs:
c         ad(*)     - Reciprocal diagonals of factored matrix
c         au(*)     - Upper profile of factored matrix
c         al(*)     - Lower profile of factored matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'counts.h'
      include  'fdata.h'
      include  'iofile.h'

      logical   flg
      character fau*(*), fal*(*)
      integer   iunau, iunal, neq, isize, isizn, maxbl
      integer   i,j,ij,in, jr,jh, ih,is,is1, jrh,jrl,jru,idh
      integer   ista, istau,istal, ifig,jp(neq), iblk(3, maxbl)
      real*4    etime, tt,tary(2)
      real*8    dot, dsred, dured, tol, dimx, dimn, dfig, dd, daval
      real*8    al(*), au(*), ad(neq)

      save

c     N.B.  Tol should be set to approximate half-word precision.

      data tol/0.5d-07/

c     Set initial values for conditioning check and solve first eq.

      dimx = abs(ad(1))
      dimn = 0.0d0
      do j = 1,neq
        dimn = max(dimn,abs(ad(j)))
      end do ! j
      dfig = 0.0d0
      if(ad(1).ne.0.0d0) ad(1) = 1.d0/ad(1)

c     Loop through columns to perform triangular decomposition

      do ij = 1, maxbl
        isize = jp(iblk(2,ij)) - jp(iblk(1,ij)-1)
        call rfile(iunau, fau, ij, isize, au)
        istau  = iblk(1, ij) - 1
        do in = iblk(3, ij), ij-1
          isizn = jp(iblk(2,in)) - jp(iblk(1,in)-1)
          call rfile(iunal, fal, in, isizn, al)
          istal = iblk(1, in) - 1
          do j = iblk(1, ij), iblk(2, ij)
            if(jp(j).gt.jp(j-1)) then
              jr = jp(j-1) + 1
              is = j - jp(j) + jr
              if(is.gt.istal) then
                is1 = is
              else
                is1 = istal + 1
                jr  = jr + is1 - is
              endif
              do i = is1, iblk(2, in)
                jr = jr + 1
                ih = min(jp(i)-jp(i-1),i-is+1)
                if(ih.gt.0) then
                  jru = jr    - jp(istau)
                  jrh = jru   - ih
                  idh = jp(i) - jp(istal) - ih + 1
                  au(jru) = au(jru) - dot(au(jrh),al(idh),ih)
                endif
              end do ! i
            endif
          end do ! j
        end do ! in

        if(flg) then
          call rfile(iunal, fal, ij, isize, al)
          istal  = iblk(1, ij) - 1
          if(iblk(1, ij).ne.iblk(1, iblk(3, ij))) then
            call wfile(iunau, fau, ij, isize, au)
          endif
          do in = iblk(3,ij), ij-1
            isizn = jp(iblk(2,in)) - jp(iblk(1,in)-1)
            call rfile(iunau, fau, in, isizn, au)
            istau = iblk(1, in) - 1
            do j = iblk(1, ij), iblk(2,ij)
              if(jp(j).gt.jp(j-1)) then
                jr = jp(j-1) + 1
                is = j - jp(j) + jr
                if(is.gt.istau) then
                  is1 = is
                else
                  is1 = istau + 1
                  jr  = jr + is1 - is
                endif
                do i = is1, iblk(2, in)
                  jr = jr + 1
                  ih = min(jp(i)-jp(i-1),i-is+1)
                  if(ih.gt.0) then
                    jrl = jr    - jp(istal)
                    jrh = jrl   - ih
                    idh = jp(i) - jp(istau) - ih + 1
                    al(jrl) = al(jrl) - dot(al(jrh),au(idh),ih)
                  endif
                end do ! i
              endif
            end do ! j
          end do ! in

          if(istau.ne.istal) then
             call rfile(iunau, fau, ij, isize, au)
          endif

c       For symmetric equation copy au into al

        else
          call pmove( au, al, isize )
        endif

        ista = iblk(1, ij)
        istau = jp(ista-1)
        do j = ista, iblk(2, ij)
          if(jp(j).ge.jp(j-1)) then
            jr = jp(j-1) + 1
            jh = jp(j)   - jr
            if(j - jh .lt. ista) jr = jp(j) - j + ista
            is = j - jh
            is1 = max(is, ista)
            if(ad(j).eq.0.0d0) call datest(au(jr),jh,daval)
            do i = is1, j-1
              jr = jr + 1
              ih = min(jp(i)-jp(i-1),i-is+1)
              if(ih.gt.0) then
                jru = jr    - istau
                jrh = jru   - ih
                idh = jp(i) - istau - ih + 1
                au(jru) = au(jru) - dot(au(jrh),al(idh),ih)
                if(flg) then
                  al(jru) = al(jru) - dot(al(jrh),au(idh),ih)
                else
                  al(jru) = au(jru)
                endif
              endif
            end do ! i

c           Reduce diagonal

            dd = ad(j)
            if(jh.ge.0) then
              jr  = jp(j) - jh
              jrh = j     - jh - 1
              jrl = jr    - istau
              if(flg) then
                ad(j) = ad(j) - dured(al(jrl),au(jrl),ad(jrh),jh+1)
              else
                ad(j) = ad(j) - dsred(au(jrl),ad(jrh),jh+1)
                call pmove(au(jrl),al(jrl),jh+1)
              endif

c             Check for possible errors and print warnings

              if(abs(ad(j)).lt.tol*abs(dd))      write(iow,2000) j
     &                                          ,nstep,niter
              if(dd.lt.0.0d0.and.ad(j).gt.0.0d0) write(iow,2001) j
     &                                          ,nstep,niter
              if(dd.gt.0.0d0.and.ad(j).lt.0.0d0) write(iow,2001) j
     &                                          ,nstep,niter
              if(ad(j) .eq.  0.0d0)              write(iow,2002) j
     &                                          ,nstep,niter

c             Complete rank test for a zero diagonal case

              if(dd.eq.0.0d0.and.jh.gt.0) then
                if(abs(ad(j)).lt.tol*daval)      write(iow,2003) j
              endif
              if(ior.lt.0) then
                if(abs(ad(j)).lt.tol*abs(dd))      write(*,2000) j
     &                                          ,nstep,niter
                if(dd.lt.0.0d0.and.ad(j).gt.0.0d0) write(*,2001) j
     &                                          ,nstep,niter
                if(dd.gt.0.0d0.and.ad(j).lt.0.0d0) write(*,2001) j
     &                                          ,nstep,niter
                if(ad(j) .eq.  0.0d0)              write(*,2002) j
     &                                          ,nstep,niter

c               Complete rank test for a zero diagonal case

                if(dd.eq.0.0d0.and.jh.gt.0) then
                  if(abs(ad(j)).lt.tol*daval)      write(*,2003) j
                endif
              endif
            endif

c           Store reciprocal of diagonal

            if(ad(j).ne.0.0d0) then
              dimx  = max(dimx,abs(ad(j)))
              dimn  = min(dimn,abs(ad(j)))
              dfig  = max(dfig,abs(dd)/abs(ad(j)))
              ad(j) = 1.d0/ad(j)
            endif
          endif
        end do ! j
        if(flg) call wfile(iunal, fal, ij, isize, al)
        call wfile(iunau, fau, ij, isize, au)
        if(pfr .and. ior.lt.0) then
          tt = etime(tary)
          write(*,2005) ij, iblk(2,ij), tary
        endif
      end do ! ij

c     Print conditioning information

      dd = 0.0d0
      if(dimn.ne.0.0d0) dd = dimx/dimn
      if(dfig.gt.0.0d0) then
        ifig = nint(log10(dfig) + 0.6d0)
      else
        ifig = -1
      endif
      if(pfr) write(iow,2004) dimx,dimn,dd,ifig
      if(pfr.and.ior.lt.0) write(*,2004) dimx,dimn,dd,ifig

c     Formats

2000  format(' *WARNING* Lost at least 7 digits reducing equation',i8/
     &       ' Step:',i7,' Iteration:',i6)

2001  format(' *WARNING* Sign of diagonal changed for equation',i8/
     &       ' Step:',i7,' Iteration:',i6)

2002  format(' *WARNING* Reduced diagonal is zero for equation',i8/
     &       ' Step:',i7,' Iteration:',i6)

2003  format(' *WARNING* Rank failure, zero diagonal in equation',i8/
     &       ' Step:',i7,' Iteration:',i6)

2004  format(' Condition check: D-max',1p,1e11.4,'; D-min',1p,1e11.4,
     & '; Ratio',1p,1e11.4/' Maximum no. diagonal digits lost:',i3)

2005  format('    --> Block',i4,' completed: No. Equations =',i8,
     &       '    t=',2f9.2)

      end
