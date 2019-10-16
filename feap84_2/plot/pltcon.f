c$Id:$
      subroutine pltcon(x,ie,ix,ip,u,
     &                  nie,ndm,ndf,nen1,nen0,ic,mc,lc,mmc,label)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add user contour plots for meshes                31/08/2008
c       2. Add 'nel' on call to ufacelib                    08/12/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot of mesh contours: With inter-element smoothing

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ie(nie,*) - Assembly data for material sets
c         ix(nen1,*)- Element nodal connections
c         ip(8,*)   - Sorted element order for each quadrant
c         u(*)      - Solution state
c         nie       - Dimension of ie array
c         ndm       - Dimension of x array
c         ndf       - Number dof/node
c         nen1      - Dimension of ix array
c         nen0      - Number nodes on plot face
c         ic        - Component number to plot
c         mc        - Number of contour lines: < 0 fill; > 0 line
c         lc        - Dimensioning information on component to plot
c         mmc       - Type of plot
c         label     - Flag, put labels on plots if true

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'chdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'pbody.h'
      include  'pdata1.h'
      include  'pdata2.h'
      include  'pdata6.h'
      include  'pdatay.h'
      include  'pdatri.h'
      include  'plcon1.h'
      include  'plflag.h'
      include  'ppers.h'
      include  'prange.h'
      include  'prmptd.h'
      include  'psdat1.h'
      include  'rpdata.h'

      character y*1
      logical   errck,cont,pinput,label,labl
      integer   nie, ndm, ndf, nen0,nen1, ic, mc, lc, mmc
      integer   i, j, n, ma, nc, nnc, nerr, nume, numf,numfac
      integer   iuf, iutot, ns, nsy, ii, ne,nel, pstyp
      real*8    vmx, vmn

      integer   ie(nie,*),ix(nen1,*),ip(8,*),iplt(50),icl(30)
      integer   ilq(4),iq9(4,4),iq16(4,9),it6(3,4)
      real*8    xl(3,29),xq(3,4),x(ndm,*),u(*),v(29),vc(12),vq(4)

      save

      data      it6 / 1, 4, 6,   4, 2, 5,   4, 5, 6,   6, 5, 3/
      data      iq9 / 1, 5, 9, 8,   5, 2, 6, 9,   8, 9, 7, 4,
     &                9, 6, 3, 7/
      data      iq16/ 1, 5,13,12,   5, 6,14,13,   6, 2, 7,14,
     &               12,13,16,11,  13,14,15,16,  14, 7, 8,15,
     &               11,16,10, 4,  16,15, 9,10,  15, 8, 3, 9/

c     Contour plot routine for elements: lines if mc > 0;
c                                        fills if mc < 0
      cont = .true.
      labl = label
      call pzerol ( tvc , .true. , 81 )
1     if(mc.gt.0) then
        nc    = max(1,min(mc,12))
        nlabi = 0
        dx1   = .024d0/scale
        vflg  = ipb.eq.0
        nerr  = 0
11      if(ior.lt.0) write(*,2001) nc
        nnc = min(8,nc)
        if (prompt .and. .not.defalt) then
          if(ior.lt.0) call pprint('   >')
          errck = pinput(vc,nnc)
          nerr  = nerr+1
          if(nerr.gt.5) then
            if(ior.lt.0) return
            call plstop
          endif
          if(errck) go to 11
        else
          vc(1) = 0.0d0
          vc(2) = 0.0d0
        endif
        if(nc.gt.1 .and. vc(1).eq.0.0d0 .and. vc(2).eq.0.0d0) then
          vc(1)   = rmn +    (rmx - rmn)/(nc+1)
          vc(nc)  = rmn + nc*(rmx - rmn)/(nc+1)
          do i = 2,nc-1
            vc(i) = vc(1) + (vc(nc)-vc(1))*(i-1)/(nc-1)
          end do ! i
        else
          if(nc.gt.8) then
            nnc = min(12,nc) - 8
            if(ior.lt.0) call pprint('   >')
            errck = pinput(vc(9),nnc)
            if(errck) go to 11
          endif
        endif
        if (prompt) then
          if(pfr) write(iow,2000) (vc(i),i=1,nc)
          if(ior.lt.0 .and. .not.defalt) then
            write(*,2000) (vc(i),i=1,nc)
          endif
        endif

c       Input label and color for first contour

13      if (prompt .and. .not.defalt) then
          if (ior.lt.0) then
            call pprint(' Input number for first contour label > ')
          endif
          errck = pinput(vlu(1),1)
          if(errck) go to 13
          nlabi = max(1,min(int(vlu(1)),12)) - 1
          if(nlabi+nc.gt.12) then
            if(ior.lt.0) write(*,2007)
            nlabi = 12 - nc
          endif
        else
          nlabi = 0
        endif

c     Inputs for filled plots

      else
        cont  = .false.
        nc    = 11
        if(ipb.ge.0) then
15        if(rangfl) then
            vc(1) = rangmn
            vc(2) = rangmx
          elseif(prompt .and. .not.defalt) then
            if(ior.lt.0) then
              write(xxx,2008) rmn,rmx
              call pprint(xxx)
            endif
            errck = pinput(vc,2)
            if(errck) go to 15
          else
            vc(1) = 0.0d0
            vc(2) = 0.0d0
          endif
          if(vc(1).eq.vc(2)) then
            vc( 1) = rmn +       (rmx - rmn)/12.0d0
            vc(11) = rmn + 11.d0*(rmx - rmn)/12.0d0
          else
            vc(11) = vc(2)
          endif
          do i = 2,10
            vc(i) = vc(1) + (vc(11)-vc(1))*(i-1)*0.1d0
          end do ! i
          if(prompt) then
            if(pfr) write(iow,2000) (vc(i),i=1,nc)
            if(ior.lt.0 .and. .not.defalt) then
              write(*,2000) (vc(i),i=1,nc)
            endif
          endif
        endif
      endif

c     If interactive, offer chance to change inputs

      if(ior.lt.0 .and. prompt .and. .not.defalt) then
        call pprint(' Input values correct? (y or n, c = cancel) > ')
20      read (*,1000,err=21,end=22) y
        goto  23
21      call  errclr ('PLTCON')
        goto  20
22      call  endclr ('PLTCON',y)
23      if(y.eq.'c' .or.y.eq.'C') return
        if(y.ne.'y'.and.y.ne.'Y') go to 1
      endif

c     Open plot and find max/min of plot variable

      call plopen
      do nsy = 1,nsym
        lsym = isym(nsy)
        nume = nfac(lsym)
        call pltsym(x,ndm,numnp,lsym)
        j   = ic
        xmx = x(1,1)
        ymx = x(2,1)
        xmn = x(1,1)
        ymn = x(2,1)
        vmn = u(j)
        vmx = u(j)
        do i = 1,numnp
          xmx = max(x(1,i),xmx)
          ymx = max(x(2,i),ymx)
          xmn = min(x(1,i),xmn)
          ymn = min(x(2,i),ymn)
          vmn = min(vmn,u(j))
          vmx = max(vmx,u(j))
          j   = j + ndf
        end do ! i
        if(xmx.ne.xmn) xmx = 8.2d0/(xmx-xmn)
        if(ymx.ne.ymn) ymx = 8.2d0/(ymx-ymn)
        if(vmn.eq.0.0d0 .and. vmx.eq.0.0d0) then
          write(iow,2005)
          if(ior.lt.0) write(*,2005)
          return
        endif

c       Open plot and loop through elements

        call pzero(xl,3*max(4,nen0))
        ic = max(1,min(ic,ndf))
        if (nsy .eq.1 ) then
          psmx  = vmn - 1.
          psmn  = vmx + 1.
        endif
        do ne = 1,nume

c         Get plot order for each element

          n  = ip(lsym,ne)

c         Plot active regions: material number: maplt; all if maplt = 0;

          ma    = ix(nen1,n)
          pstyp = ie(1,ma)
          if((ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2) .and.
     &       (pstyp.ne.0) .and.
     &       ((ma.gt.0 .and. maplt.eq.0) .or. ma.eq.maplt) ) then
            ma = ie(nie-1,ma) ! iel

c           Determine maximum number of nodes on element

            do i = nen0,1,-1
              if(ix(i,n).ne.0) then
                nel = i
                exit
              endif
            end do ! i
            nel = min(nel,nen)
            if(nel.gt.2 .and. ix(1,n).eq.ix(nel,n)) then
              nel       = nel - 1
            endif

c           Get plot order for each element

            call plftyp(pstyp,nel,ma)

c           Perspective or standard feap elements

            if(kpers.eq.1 .or. pstyp.gt.0) then
              call pltord(ix(1,n),ma, iuf,iplt)

c             Set values of vlu(1) and vlu(2)

              iutot  = 0
              vlu(1) = vmx
              vlu(2) = vmn
              do i = 1,iuf
                ns = iplt(i)
c               if(ns.le.nen0) then
                if(ns.le.nel) then
                  ii = ix(iplt(i),n)
                  if(ii.gt.0) then
                    iutot    = iutot + 1
                    xl(1,ns) = x(1,ii)
                    xl(2,ns) = x(2,ii)
                    if(ndm.ge.3) xl(3,ns) = x(3,ii)
                    j      = ndf*(ii-1) + ic
                    v(ns)  = u(j)
                    vlu(1) = min(vlu(1),v(ns))
                    vlu(2) = max(vlu(2),v(ns))

c                   Plot min/max for graphics

                    if(psmn.gt.v(ns)) then
                      psmn    = v(ns)
                      xpsn(1) = xl(1,ns)
                      xpsn(2) = xl(2,ns)
                      xpsn(3) = xl(3,ns)
                    endif
                    if(psmx.lt.v(ns)) then
                      psmx    = v(ns)
                      xpsx(1) = xl(1,ns)
                      xpsx(2) = xl(2,ns)
                      xpsx(3) = xl(3,ns)
                    endif
                  endif
                endif
              end do ! i
              if(iutot.gt.3) then

c               Linear Triangle

                if(nel.eq.3) then
                  if(cont) then
                    call pltecn(xl,v,vc,nc)
                  else
                    call pltcor(3,icl,v,vc,nc)
                    call pltefl(3,icl,xl,v,vc,nc)
                  endif

c               Linear Quadrilateral

                elseif(nel.eq.4) then
                  call pltcor(nel,icl,v,vc,nc)
                  call pltqfl(icl,xl,v,vc,nc,cont)

c               Quadratic Triangle

                elseif(nel.eq.6 .or. nel.eq.7) then
                  do j = 1,4
                    do i = 1,3
                      ii      = it6(i,j)
                      xq(1,i) = xl(1,ii)
                      xq(2,i) = xl(2,ii)
                      xq(3,i) = xl(3,ii)
                      vq(i)   = v(ii)
                    end do ! i
                    xq(1,4) = xq(1,1)
                    xq(2,4) = xq(2,1)
                    xq(3,4) = xq(3,1)
                    vq(4)   = vq(1)
                    if(cont) then
                      call pltecn(xq,vq,vc,nc)
                    else
                      call pltcor(3,ilq,vq,vc,nc)
                      call pltefl(3,ilq,xq,vq,vc,nc)
                    endif
                  end do ! j

c               Quadratic quad

                elseif((nel.eq.8 .and. iuf.ne.16) .or. nel.eq.9) then
                  if(nel.eq.8) then
                    do i = 1,3
                      xl(i,9) = - 0.25d0*(xl(i,1) + xl(i,2)
     &                                  + xl(i,3) + xl(i,4))
     &                          + 0.50d0*(xl(i,5) + xl(i,6)
     &                                  + xl(i,7) + xl(i,8))
                    end do ! i
                    v(9) = - 0.25d0*(v(1) + v(2) + v(3) + v(4))
     &                     + 0.50d0*(v(5) + v(6) + v(7) + v(8))
                    nel = 9
                  else
                    ii = ix(9,n)
                    do i = 1,ndm
                      xl(i,9) = x(i,ii)
                    end do ! i
                    v(9)  = u(ndf*(ii-1) + ic)
                  endif
                  nel = 9
                  call pltcor(nel,icl,v,vc,nc)
                  do j = 1,4
                    do i = 1,4
                      ii     = iq9(i,j)
                      ilq(i) = icl(ii)
                      xq(1,i) = xl(1,ii)
                      xq(2,i) = xl(2,ii)
                      xq(3,i) = xl(3,ii)
                      vq(i)   = v(ii)
                    end do ! i
                    call pltqfl(ilq,xq,vq,vc,nc,cont)
                  end do ! j

c               Cubic quad

                elseif(nel.eq.16) then
                  do j = 13,16
                    ii = ix(j,n)
                    do i = 1,ndm
                      xl(i,j) = x(i,ii)
                    end do ! i
                    v(j) = u(ndf*(ii-1) + ic)
                  end do ! j
                  call pltcor(nel,icl,v,vc,nc)
                  do j = 1,9
                    do i = 1,4
                      ii      = iq16(i,j)
                      ilq(i)  = icl(ii)
                      xq(1,i) = xl(1,ii)
                      xq(2,i) = xl(2,ii)
                      xq(3,i) = xl(3,ii)
                      vq(i)   = v(ii)
                    end do ! i
                    call pltqfl(ilq,xq,vq,vc,nc,cont)
                  end do ! j

c               If all else fails plot a quad

                else
                  call pltcor(nel,icl,v,vc,nc)
                  call pltqfl(icl,xl,v,vc,nc,cont)
                endif

c               Draw border around element

                call pppcol(0,1)
                if(.not.cont .and. ipb.eq.0 ) then
                  call plotl(xl(1,iplt(1)),
     &                       xl(2,iplt(1)),
     &                       xl(3,iplt(1)),3)
                  do i = 2,iuf-1
                    j = iplt(i)
                    if(j.le.nel) then
                      call plotl(xl(1,j),xl(2,j),xl(3,j),2)
                    endif
                  end do ! i
                  call plotl(xl(1,iplt(1)),
     &                       xl(2,iplt(1)),
     &                       xl(3,iplt(1)),2)
                else if(edgfl .and. iuf.eq.5) then
                  if(cont) call pppcol(5,1)
                  call pedge4(nel,ix,xl,nen1,n)
                endif

c             Line elements

              else
                if(ix(nel,ne).eq.ix(nel-1,ne)) then
                  nel = nel - 1
                endif
                call pltlfl(nel,xl,v,vc,nc)
              endif

c           User plot types

            elseif(pstyp.lt.0) then

              call ufacelib(pstyp,nel,ipu,numfac)

c             Set values of vlu(1) and vlu(2)

              do numf = 1,numfac
                iutot  = 0
                vlu(1) = vmx
                vlu(2) = vmn
                do i = 1,4
                  ns = ipu(i,numf)
                  if(ns.le.nel) then
                    ii = ix(ns,n)
                    if(ii.gt.0) then
                      iutot    = iutot + 1
                      xl(1,iutot) = x(1,ii)
                      xl(2,iutot) = x(2,ii)
                      if(ndm.ge.3) xl(3,iutot) = x(3,ii)
                      j      = ndf*(ii-1) + ic
                      v(iutot)  = u(j)
                      vlu(1) = min(vlu(1),v(iutot))
                      vlu(2) = max(vlu(2),v(iutot))

c                     Plot min/max for graphics

                      if(psmn.gt.v(iutot)) then
                        psmn    = v(iutot)
                        xpsn(1) = xl(1,iutot)
                        xpsn(2) = xl(2,iutot)
                        xpsn(3) = xl(3,iutot)
                      endif
                      if(psmx.lt.v(iutot)) then
                        psmx    = v(iutot)
                        xpsx(1) = xl(1,iutot)
                        xpsx(2) = xl(2,iutot)
                        xpsx(3) = xl(3,iutot)
                      endif
                    endif
                  endif
                end do ! i
                do i = 1,3
                  xl(i,iutot+1) = xl(i,1)
                end do ! i
                call pltcor(  4,icl,v,vc,nc)
                call pltqfl(icl,xl,v,vc,nc,cont)
              end do ! numf

c             Draw border around element

              call pppcol(0,1)
c             ma  = ie(nie-1,n) ! iel
              iuf = inord(ma)
              if(.not.cont .and. ipb.eq.0 ) then
                call plotl(x(1,ix(ipord(1,ma),n)),
     &                     x(2,ix(ipord(1,ma),n)),
     &                     x(3,ix(ipord(1,ma),n)),3)
                do i = 2,iuf
                  call plotl(x(1,ix(ipord(i,ma),n)),
     &                       x(2,ix(ipord(i,ma),n)),
     &                       x(3,ix(ipord(i,ma),n)),2)
                end do ! i
                call plotl(x(1,ix(ipord(1,ma),n)),
     &                     x(2,ix(ipord(1,ma),n)),
     &                     x(3,ix(ipord(1,ma),n)),2)
              endif
            endif
          endif

        end do ! ne

c       Put on labels

        if(labl) then
          if(cont) then
            call pltctx(vc,lc,nlabi,nc,mmc)
          else
            call pltftx(vc,-mc,mmc)
          endif
          labl = .false.
        end if

        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

c     Formats

1000  format(a)
2000  format('   ------ Contour Values for Plot ------'/(3x,1p,5e15.6))
2001  format(' Input',i3,' Contour Values for Plot - 8 Values/Line')
2005  format(' *WARNING* No plot - all zero values')
2007  format(' *WARNING* Initial label reset to fit screen')
2008  format(' Input Min/Max (Default:',1p,e9.2,'/',1p,e9.2,'): >')

      end
