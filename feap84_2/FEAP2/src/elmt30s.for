      subroutine elmt30(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c----------------------------------------------------------------------+
c
c.... loading element
c     purpose: adding surface loads on 4 and 9 node elements via qload
c    
c   
c     Input
c      via qload
c
c----------------------------------------------------------------------+
      USE cdata
      USE eldata
      USE iofile
      USE pdata6
      USE qload
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),tl(*),d(*),ul(ndf,*),s(nst,*),p(*),ix(*),
     1          shp(3,9),sg(16),tg(16),wg(16)
      dimension yl(3,9),tr(3,3) 
      dimension ipord4(5),ipord9(9)
      dimension h1(*),h2(*),h3(*)    
      
      data ipord4  /1,2,3,4,1/
      data ipord9  /1,5,2,6,3,7,4,8,1/
c.... transfer to correct processor
      go to (1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,22),isw
      return
c.... input material properties
1     continue
                   write(iow,4001) 
      if(ior.lt.0) write(*  ,4001) 
c
2     continue
      if(isw.eq.2) then
        ielno = 30
        if(nel.eq.4) then 
          inord(ielno) = 5
          do ii = 1,5
            ipord(ii,ielno) = ipord4(ii)
          end do
        else if(nel.eq.9) then
          inord(ielno) = 9        
          do ii = 1,9
            ipord(ii,ielno) = ipord9(ii)
          end do
        end if
      end if
      return
c.... Surface load via qload
22    continue
c.... Setting up local coordinate system
      call trans30(xl,yl,tr,ndm,nel)
      if(nel.eq.4.d0) then
        ll=2
      else if(nel.eq.9.d0) then
        ll=3
      end if
      call pgauss(ll,lint,sg,tg,wg)
      do l=1,lint
        call shape(sg(l),tg(l),yl,shp,xsj,ndm,nel,ix,.true.)
        da=xsj*wg(l)
c....   Calculating nodal forces
        call qload30(shp,xl,ul,tr,aqloa,numel,ndf,nst,p,s,da)
      end do
      
c.... formats
4001  format(a1,20a4//'  Loading Element'//
     1 ' no additional input necessary   '/)
      end
c
      subroutine trans30(xl,yl,t,ndm,nel)
c---------------------------------------------------------------------------
c.... compute the transformation array and surface coords.
c     y=Tx 
c     comment ww: restricted to plain elements
c
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x0(3),xl(ndm,nel),yl(3,nel),t(3,3)
c
c.... compute the center (0,0)
      do i = 1,3
        x0(i) = 0.25d0*(xl(i,1) + xl(i,2) + xl(i,3) + xl(i,4))
      end do

c     basic via diagonals        
c.... compute the inplane direction cosines (bisect diagonals)
      t=0
      do i = 1,3
        t(1,i) = xl(i,3) - xl(i,1)
        t(2,i) = xl(i,2) - xl(i,4)
      end do    
      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)
      do i = 1,3
        v1 = t(1,i)/dl1
        v2 = t(2,i)/dl2
        t(1,i) = v1 + v2
        t(2,i) = v1 - v2
      end do    
      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)

c.... vectors t_1 and t_2
      do i = 1,3
        t(1,i) = t(1,i)/dl1
        t(2,i) = t(2,i)/dl2
      end do       

c.... compute the normal to the surface
      t(3,1) = t(1,2)*t(2,3) - t(2,2)*t(1,3)
      t(3,2) = t(1,3)*t(2,1) - t(2,3)*t(1,1)
      t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
      
c.... compute the projected middle surface coordinates
      yl = 0
      do i = 1,nel
        do j = 1,3
          do k = 1,3
            yl(j,i) = yl(j,i) + t(j,k)*(xl(k,i) - x0(k))
          end do ! k
        end do ! j
      end do ! i

c.... set offset coordinates to zero if small compared to plan size
      htol =  0.0d0
      do i = 1,nel
        htol = max(htol,abs(yl(1,i)),abs(yl(2,i)))
      end do 
      htol = htol*1.e-7
      do i = 1,nel
        if(abs(yl(3,i)) .le. htol) yl(3,i) = 0.0d0
      end do  

      return
      end
c
      subroutine qload30(shp,xl,ul,tr,q,numel,ndf,nst,p,s,da)
c-----------------------------------------------------------------------
c.... Input properties
c     q(n,1) = element n, load q_x
c     q(n,2) = element n, load q_y
c     q(n,3) = element n, load q_z
c     q(n,4) = element n, iltyp, 0 = global, 1 = local, 2=follower load 
c
c     comment: follower load only with (local) q_z
c
c-----------------------------------------------------------------------
c     
      USE eldata
      USE prlod
      USE qload
c
      implicit double precision (a-h,o-z)
      dimension q(numel,10),q2(3)
      dimension shp(3,*),xl(3,*),ul(ndf,*),p(*),s(nst,*),tr(3,3)
      dimension xu1(3,3)

      if (mqloa.eq.1) return
      iltyp = q(n,4)

      q2 = 0.d0

      if (iltyp.lt.2) then
c....   Conservative loads 
        if (iltyp.eq.0) then ! using q in global coordinates
          do i=1,3
            q2(i) = q(n,i)
          end do

        else if (iltyp.eq.1) then ! using q in local coordinates
          do i=1,3
            do j=1,3
              q2(i) = q2(i)+tr(j,i)*q(n,j)
            end do
          end do
        end if

c....   Setting up deformation-independent load vector
        ir0 = 0
        do i = 1,nel
          fact = propq*shp(3,i)*da
          do j = 1,ndf
            p(ir0+j) = p(ir0+j) + fact*q2(j)
          end do
          ir0 = ir0 + ndf
        end do

      else if (iltyp.eq.2) then
c....   Follower loads

        q2(3) = q(n,3)

c....   Setting up deformation-dependent load vector
        xu1 = 0.d0
        do i = 1,3
          do j = 1,2
            do k = 1,nel
              xu1(i,j) = xu1(i,j) + (xl(i,k)+ul(i,k))*shp(j,k)
            end do
          end do
        end do

        call vecp(xu1(1,1),xu1(1,2),xu1(1,3))

        ir0 = 0
        do inode = 1,nel
          fact = propq*q2(3)*shp(3,inode)
           do i = 1,3
            p(ir0+i) = p(ir0+i) + xu1(i,3)*fact
          end do

c....   Setting up sub-tangent stiffness matrix
          jc0 = 0
          do jnode = 1,nel
            s1 = (xu1(1,1)*shp(2,jnode) - xu1(1,2)*shp(1,jnode))*fact
            s2 = (xu1(2,1)*shp(2,jnode) - xu1(2,2)*shp(1,jnode))*fact
            s3 = (xu1(3,1)*shp(2,jnode) - xu1(3,2)*shp(1,jnode))*fact
            s(ir0+2,jc0+3) = s(ir0+2,jc0+3) + s1
            s(ir0+3,jc0+1) = s(ir0+3,jc0+1) + s2
            s(ir0+1,jc0+2) = s(ir0+1,jc0+2) + s3
            s(ir0+3,jc0+2) = s(ir0+3,jc0+2) - s1
            s(ir0+1,jc0+3) = s(ir0+1,jc0+3) - s2
            s(ir0+2,jc0+1) = s(ir0+2,jc0+1) - s3
            jc0 = jc0 + ndf
          end do
          ir0 = ir0 + ndf
        end do
      end if

      return
      end
