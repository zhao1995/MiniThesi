c$Id:$
      subroutine pltlns(x,ix,u,numn,ndm,ndf,nen1,nume,ic,mc,mmc,label)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot of mesh contours: With inter-element smoothing

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ix(nen1,*)- Element nodal connections
c         u(*)      - Solution state
c         numn      - Number of solution states
c         ndm       - Dimension of x array
c         ndf       - Number dof/node
c         nen1      - Dimension of ix array
c         nume      - Number of segments on surface
c         ic        - Component number to plot
c         mc        - Number of contour lines: < 0 fill; > 0 line
c         mmc       - Type of plot
c         label     - Flag, put labels on plots if true

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'pbody.h'
      include  'pdata1.h'
      include  'pdata2.h'
      include  'pdatri.h'
      include  'prange.h'
      include  'prmptd.h'
      include  'psdat1.h'
      include  'rpdata.h'
      include  'pdatay.h'

      character y*1
      logical   vflg,errck,cont,pinput,label,labl
      integer   numn,ndm,ndf,nen1,nume,ic,mc,mmc
      integer   i, j, k, n, nc, nnc, nerr, nlabi
      integer   ns, nsy, ii, ix(nen1,*),iln(2)
      real*8    dx1, vmx, vmn
      real*8    xl(3,29),x(ndm,*), u(*),v(9),vc(12),vlu(2)

      save

c     Contour plot routine for elements: lines if mc > 0;
c                                        fills if mc < 0

      cont = .true.
      labl = label
1     if(mc.gt.0) then
        nc    = max(1,min(mc,12))
        nlabi = 0
        dx1   = .024d0/scale
        vflg  = (idev.eq.2.and.ipb.eq.0).or.(idev.eq.1.and.ipb.ne.0)
        nerr=0
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
          if(ior.lt.0 .and. idev.eq.1 .and. .not.defalt ) then
            write(*,2000) (vc(i),i=1,nc)
          endif
        endif

c       Input label and color for first contour

13      if(prompt .and. .not.defalt) then
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
              write(*,2008) rmn,rmx
              call pprint('   >')
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
            if(ior.lt.0 .and. .not.defalt ) then
              write(*,2000) (vc(i),i=1,nc)
            endif
          endif
        endif
      endif

c     If interactive, offer chance to change inputs

      if(ior.lt.0 .and. prompt .and. .not.defalt ) then
        call pprint(' Input values correct? (y or n, c = cancel) > ')
20      read (*,1000,err=21,end=22) y
        goto  23
21      call  errclr ('PLTLNS')
        goto  20
22      call  endclr ('PLTLNS',y)
23      if(y.eq.'c' .or.y.eq.'C') return
        if(y.ne.'y'.and.y.ne.'Y') go to 1
      endif

c     Set a thick line

      call plopen

c     Find max/min of plot variable

      do nsy = 1,nsym

        lsym = isym(nsy)
        call pltsym(x,ndm,numnp,lsym)
        iln(1) = 1
        iln(2) = 5
        call plline(iln)
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
        end do ! i
        do i = 1,nume+1
          vmn = min(vmn,u(j))
          vmx = max(vmx,u(j))
          j   = j + ndf
        end do ! i
        if(xmx.ne.xmn) xmx = 8.2d0/(xmx-xmn)
        if(ymx.ne.ymn) ymx = 8.2d0/(ymx-ymn)

c       Open plot and loop through elements

        call pzero(xl,3*max(4,nen))

        psmx  = vmn*(1.d0 - 1.d-8)
        psmn  = vmx*(1.d0 + 1.d-8)

        do n = 1,nume

c         Set values of vl and vu

          vlu(1) = vmx
          vlu(2) = vmn
          ns = 0
          do i = 1,2
            ii = ix(i,n)
            if(ii.gt.0) then
              ns = ns + 1
              xl(1,ns) = x(1,ii)
              xl(2,ns) = x(2,ii)
              xl(3,ns) = x(3,ii)
              if(n+i-1 .le. numn) then  ! Check for closed surfaces
                k = n + i - 1
              else
                k = 1
              endif
              v(ns)    = u(k)
              vlu(1)   = min(vlu(1),v(ns))
              vlu(2)   = max(vlu(2),v(ns))

c             Plot min/max for graphics

              if(psmn.gt.v(ns)) then
                psmn = v(ns)
                xpsn(1) = xl(1,ns)
                xpsn(2) = xl(2,ns)
                xpsn(3) = xl(3,ns)
              endif
              if(psmx.lt.v(ns)) then
                psmx = v(ns)
                xpsx(1) = xl(1,ns)
                xpsx(2) = xl(2,ns)
                xpsx(3) = xl(3,ns)
              endif
            endif
          end do ! i
          xl(1,ns+1) = xl(1,1)
          xl(2,ns+1) = xl(2,1)
          xl(3,ns+1) = xl(3,1)
          v(ns+1)    = v(1)

          call pltlfl(ns,xl,v,vc,nc)

        end do ! n

c       Set  default for line thickness

        iln(1) = 1
        iln(2) = 1
        call plline(iln)

c       Put on labels

        if(labl) then
          call pltftx(vc,-mc,mmc)
          labl = .false.
        end if

        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

c     Formats

1000  format(a)
2000  format('   ------ Contour Values for Plot ------'/(3x,5e15.6))
2001  format(' Input',i3,' Contour Values for Plot - 8 Values/Line')
2007  format(' ** WARNING ** Initial label reset to fit screen')
2008  format(' Input Minimum/Maximum Values for Fill Plot')

      end
