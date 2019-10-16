c$Id:$
      subroutine bpltcon(x,ix,ip,u,unumnp,ndm,ndf,nen,nen1,ic,label)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot of beam surface mesh contours: With inter-element
c               smoothing

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ix(nen1,*)- Element nodal connections
c         ip(8,*)   - Sorted element order for each quadrant
c         u(*)      - Solution state
c         ndm       - Dimension of x array
c         ndf       - Number dof/node
c         nen       - Size of connections in ix array
c         nen1      - Dimension of ix array
c         ic        - Component number to plot
c         label     - Flag, put labels on plots if true

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include   'fdata.h'
      include   'iofile.h'
      include   'pbody.h'
      include   'pdata2.h'
      include   'pdatay.h'
      include   'pdatri.h'
      include   'prmptd.h'
      include   'psdat1.h'
      include   'rpdata.h'

      character y*1
      logical   zoom, tvc(9,9),vflg,cont,errck,offer,pinput,label,labl
      integer   ndm, ndf, nen, nen1, ic, mc, mmc, icolor
      integer   i, j, n, nc, nlabi, nume, unumnp
      integer   iu, iutot, ma, ns, nsy, ii, ne
      real*8    dx1, vl, vu, vmx, vmn

      integer   ix(nen1,*),ip(8,*),iplt(5)
      real*8    xl(3,6),x(ndm,*),u(*),v(6),vc(12)

      save

      data       iplt /1,2,3,4,1/

c     Contour plot routine for elements: lines if mc > 0;
c                                        fills if mc < 0
      labl = label
      cont = .false.
      vflg = .false.
      mc   = 1 ! Contour # 1
      mmc  = 1 ! Label S T R E S S
      call pzerol ( tvc , .true. , 81 )
1     offer = .false.

c     Inputs for filled plots

      nc    = 11
15    offer = idev.ne.2
      if(prompt .and. .not.defalt) then
        if(ior.lt.0) then
          if(idev.eq.1 ) then
            write(*,2008)
          elseif(idev.eq.2 ) then
            write(*,2009) rmn,rmx
          endif
          call pprint('   ')
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
        if(ior.lt.0 .and. idev.eq.1 .and. .not.defalt ) then
          write(*,2000) (vc(i),i=1,nc)
        endif
      endif

c     If interactive, offer chance to change inputs

      if(ior.lt.0 .and. offer .and. prompt .and. .not.defalt ) then
        call pprint(' Input values correct? (y or n, c = cancel) > ')
20      read (*,1000,err=21,end=22) y
        goto  23
21      call  errclr ('PLTCON')
        goto  20
22      call  endclr ('PLTCON',y)
23      if(y.eq.'c' .or.y.eq.'C') return
        if(y.ne.'y'.and.y.ne.'Y') go to 1
      endif

c     Find max/min of plot variable

      call plopen

      do nsy = 1,nsym

        lsym = isym(nsy)
        nume = nfac(lsym)
        call pltsym(x,ndm,unumnp,lsym)
        j   = ic
        xmx = x(1,1)
        ymx = x(2,1)
        xmn = x(1,1)
        ymn = x(2,1)
        vmn = u(j)
        vmx = u(j)
        do i = 1,unumnp
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
          write(iow,3000)
          if(ior.lt.0) write(*,3000)
          return
        endif

c       Open plot and loop through elements

        call pzero(xl,3*max(4,nen))
        ic = max(1,min(ic,ndf))
        if (nsy .eq.1 ) then
          psmx  = vmn - 1.
          psmn  = vmx + 1.
        endif

        do ne = 1,nume

c         Check for active regions

          ma = ix(nen1,n)
          if((ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2) .and.
     &      ((ma.gt.0 .and. maplt.eq.0) .or. maplt.eq.ma)  ) then

c           Get plot order for each element

            n  = ip(lsym,ne)

c           Plot active regions: matl number: maplt; all if maplt = 0;

            iu    = 5
            iutot = iu

c           Check if element is in window and set values of vl and vu

            vl = vmx
            vu = vmn
            ns = 0
            do i = 1,iu
              ii = ix(iplt(i),n)
              if(ii.gt.0) then
                ns = ns + 1
                xl(1,ns) = x(1,ii)
                xl(2,ns) = x(2,ii)
                if(ndm.ge.3) xl(3,ns) = x(3,ii)
                j = ndf*(ii-1) + ic
                v(ns) = u(j)
                vl = min(vl,v(ns))
                vu = max(vu,v(ns))

c               Plot min/max for graphics

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
            if(ns.gt.3 .and. zoom(xl,3) ) then
              call pltris(ic,nc,n,ns,iutot,ndm,ndf,nen,nen1,nlabi,
     &                    icolor,ix,x,xl,v,vc,dx1,vl,vu,tvc,cont,vflg)
            endif

          endif ! region check
        end do ! ne

c       Put on labels

        if(labl) then
          call pltftx(vc,mc,mmc)
          labl = .false.
        end if

        call pltsym(x,ndm,unumnp,lsym)

      end do ! nsy

c     Formats

1000  format(a)
2000  format('   ------ Contour Values for Plot ------'/(3x,5e15.6))
2008  format(' Input Minimum/Maximum Values for Fill Plot')
2009  format(' Input Min/Max (Default:',1p,e9.2,'/',1p,e9.2,')')
3000  format(' *WARNING* No plot - all zero values')

      end
