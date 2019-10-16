c$Id:$
      subroutine pltwir(x,ie,ix,ip,u,
     &                  nie,ndm,ndf,nen1,nen0,ic,mc,mmc,label)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change pstyp.gt.0 to pstyp.ne.0                  31/08/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot of mesh contours: Wire frame plots

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
      include  'pdatay.h'
      include  'pdatri.h'
      include  'plcon1.h'
      include  'plflag.h'
      include  'prange.h'
      include  'prmptd.h'
      include  'psdat1.h'
      include  'rpdata.h'

      character y*1
      logical   errck,cont,pinput,label,labl
      integer   nie, ndm, ndf, nen0,nen1, ic, mc, mmc
      integer   i, j, n, ma, nc, nume
      integer   iu, iutot, ns, nsy, ii, ne,nel, pstyp
      real*8    vmx, vmn

      integer   ie(nie,*),ix(nen1,*),ip(8,*),iplt(50)
      real*8    xl(3,29),xq(3,4),x(ndm,*),u(*),v(29),vc(12),vq(4)

      save

c     Contour plot routine for elements: lines

      cont = .true.
      labl = label

c     Inputs for filled plots

      cont  = .false.
      nc    = 11
1     if(ipb.ge.0) then

15      if(rangfl) then
          vc(1) = rangmn
          vc(2) = rangmx
        elseif(prompt .and. .not.defalt) then
          if(ior.lt.0) then
            write(xxx,2003) rmn,rmx
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

c     If interactive, offer chance to change inputs

      if(ior.lt.0 .and. prompt .and. .not.defalt) then
        call pprint(' Input values correct? (y or n, c = cancel) > ')
20      read (*,'(a)',err=21,end=22) y
        goto  23
21      call  errclr ('PLTWIR')
        goto  20
22      call  endclr ('PLTWIR',y)
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
          write(iow,2001)
          if(ior.lt.0) write(*,2001)
          return
        endif

c       Loop through elements

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

          ma = ix(nen1,n)
          pstyp = ie(1,ma)
          if((ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2) .and.
     &       (pstyp.ne.0) .and.
     &      ((ma.gt.0 .and. maplt.eq.0) .or. ma.eq.maplt) ) then
            ma = ie(nie-1,ma)

c           Get maximum number of nodes on element

            do i = nen0,1,-1
              if(ix(i,n).ne.0) then
                nel = i
                exit
              endif
            end do ! i
            nel = min(nel,nen)

c           Get plot order for each element

            call plftyp(pstyp,nel,ma)
            call pltord(ix(1,n),ma, iu,iplt)

c           Set values of vlu(1) and vlu(2)

            iutot  = 0
            vlu(1) = vmx
            vlu(2) = vmn
            do i = 1,iu
              ns = iplt(i)
              if(ns.le.nel) then
                ii = ix(iplt(i),n)
                if(ii.gt.0) then
                  iutot    = iutot + 1
                  xl(1,iutot) = x(1,ii)
                  xl(2,iutot) = x(2,ii)
                  if(ndm.ge.3) xl(3,iutot) = x(3,ii)
                  j      = ndf*(ii-1) + ic
                  v(iutot)  = u(j)
                  vlu(1) = min(vlu(1),v(iutot))
                  vlu(2) = max(vlu(2),v(iutot))

c                 Plot min/max for graphics

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

c           Line elements for edges

            do i = 1,iutot-1

              vq(1)   = v(i)
              vq(2)   = v(i+1)

              xq(1,1) = xl(1,i)
              xq(2,1) = xl(2,i)
              xq(3,1) = xl(3,i)

              xq(1,2) = xl(1,i+1)
              xq(2,2) = xl(2,i+1)
              xq(3,2) = xl(3,i+1)
              call pltlfl(2,xq,vq,vc,nc)
            end do ! i
          endif

        end do ! ne

c       Put on labels

        if(labl) then
          call pltftx(vc,-mc,mmc)
          labl = .false.
        end if

        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

c     Formats

2000  format('   ------ Contour Values for Plot ------'/(3x,1p,5e15.6))
2001  format(' *WARNING* No plot - all zero values')
2003  format(' Input Min/Max (Default:',1p,e9.2,'/',1p,e9.2,'): >')

      end
