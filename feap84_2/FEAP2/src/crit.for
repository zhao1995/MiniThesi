      subroutine pccrit(fail,crit1,crit2,numnp)
c----------------------------------------------------------------------
c     Failure criterion at nodes and layer boundaries
c----------------------------------------------------------------------
      USE pcrit
      implicit double precision (a-h,o-z)
      dimension fail(numnp,ncs),crit1(numnp,2),crit2(numnp,2)
c....  test criterion:  overshooting at layer boundary ics
      do i = 1,numnp
        n = crit1(i,1)
        i0 = 0
        f0 = 0
        if(icc.eq.0)then
          do ii = 1,ncs
            f  = fail(i,ii)
c....       delamination (Hashin) see book Reddy p. 128
            if(f.gt.f0) then
              i0 = ii
              f0 = f
            end if
          end do
        else if(icc.ne.0)then
          if(fail(i,icc).gt.0) then
            i0 = icc
            f0 = fail(i,icc)
          end if
        end if
c
c.... setup crit2
c
        if(n.eq.0) then
          crit2(i,1) = i0
          crit2(i,2) = f0
        else if(n.ne.0) then
          crit2(i,1) = crit1(i,1)
          crit2(i,2) = crit1(i,2)
        end if
      end do
      return
      end
c
      subroutine pltstr1(dt,st,numnp,ncs)
c----------------------------------------------------------------------
c     calculate average of st field at nodes
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*)
      do ii = 1,numnp
        dh = dt(ii)
        if(dh.ne.0.0d0) then
          do kk = 1,ncs
            st(ii,kk) = st(ii,kk)/dh
          end do
        end if
      end do
      return
      end
c
      subroutine prtcrit(fail,crit2,numnp,ncs,ipc)
c----------------------------------------------------------------------
c     Print Failure criterion at nodes for all layer boundaries
c----------------------------------------------------------------------
      USE bdata
      USE fdata
      USE iofile
      implicit double precision (a-h,o-z)
      dimension fail(numnp,ncs),crit2(numnp,2)
      if(ipc.eq.1) then
        kount = 0
        do n = 1,numnp
          if(crit2(n,1).ne.0) then
            kount = kount - 1
            if(kount.le.0) then
              write(iow,2000) o,head
c             if(ior.lt.0.and.pfr) then
              if(ior.lt.0) then
                write(*,2000) o,head
              end if
              kount = 50
            end if
            write(iow,2001) n,crit2(n,1),crit2(n,2)
c           if(ior.lt.0.and.pfr) then
            if(ior.lt.0) then
              write(*,2001) n,crit2(n,1),crit2(n,2)
            end if
          end if
        end do
      else
        kount = 0
        do n = 1,numnp
          kount = kount - 1
          if(kount.le.0) then
            write(iow,2002) o,head
c           if(ior.lt.0.and.pfr) then
            if(ior.lt.0) then
              write(*,2002) o,head
            end if
            kount = 50
          end if
          write(iow,2003) n,(fail(n,i),i=1,ncs)
c         if(ior.lt.0.and.pfr) then
          if(ior.lt.0) then
            write(*,2003) n,(fail(n,i),i=1,ncs)
          end if
        end do
      end if
c
2000  format(a1,19a4,a3/'   n o d a l   c r i t e r i o n'/
     1 ' node','  crack at boundary  value of criterion')
2001  format(i5,g12.5,2x,g12.5)
2002  format(a1,19a4,a3/'   n o d a l   c r i t e r i o n'/
     1 ' node','  value of f at boundaries ics = 1,ncs')
2003  format(i5,1p6e12.5/(5x,1p6e12.5))
      return
      end
c
      subroutine pinitv(crit1,crit2,v,ndf,numnp)
c----------------------------------------------------------------------
c     initialize of v(ndf-2) to v(ndf)
c
c     ilay must be equal jlay  !!!!!!!!!!!!!!!!!!
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension crit1(numnp,2),crit2(numnp,2),v(ndf,numnp)
      do i = 1,numnp
        ics1 = crit1(i,1)
        ics2 = crit2(i,1)
        dnorm = dot(v(ndf-2,i),v(ndf-2,i),3)
        if ((ics1.eq.0).and.(ics2.ne.0).and.(dnorm.lt.1.d-20)) then
          idf = 3 * (ics2+1)
          do j = 1,3
            v(ndf-3+j,i) = v(idf-3+j,i)
          end do
        end if
      end do

      return
      end
