      subroutine elmt14(d,ul,xl,ix,t,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c----------------------------------------------------------------------
c.... contact node-to-node element                    ELNLE2_1
c.... contact in one direction possible
c.... gap: (x_2+u_2) - (x_1+u_1) > 0
c.... penalty + augmented lagrange formulation
c
c     Input Record 1.
c       idf   - d.o.f. number for applying contact constraint
c               (idf is the same as the coordinate direction)
c       pen   - penalty value for imposing contact constraint
c       tol14 - tolerance on gap (default 0.0)
c.... w. wagner 11/92  extract from PCFEAP+augmented 9/03
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE hdata
      USE iofile
      USE pdata10
      implicit double precision (a-h,o-z)
      dimension ix(*),ixl(5),d(*),ul(ndf,*),s(nst,*),p(*),
     *          xl(ndm,*),t(*),xll(3,4)
      dimension h1(*),h2(*),h3(*)
      data ixl/1,2,3,4,1/,eps/1.0e-7/
      go to (1,2,3,3,3,3,2,2,2,3,2,2,3), isw
      return
c.... input material data
1     if(ior.lt.0) write(*,3000)
      call dinput(d,3)
      idf = d(1)
                   write(iow,2000) idf,d(2),d(3)
      if(ior.lt.0) write(*,2000)   idf,d(2),d(3)
c...  h-array  1 term  Force
      nh1  =  1
2     return
c.... compute the contact matrix and properties
3     idf = d(1)
      tol14 = d(3)
      ttl = 0.0d0
      dx  = xl(idf,2)-xl(idf,1) 
      du  = ul(idf,2)-ul(idf,1)
c.... gap with respect to coordinates
      if(dx.ge.0.d0) then    
        dg  = dx + du
      else
        dg  = -(dx + du)
      end if
c.... penetration
      if(dg .lt. tol14) then
        tt = d(2)*dg
        i  = ndf + idf
c....   stiffness matrix 
        if(isw .eq. 3) then
          s(idf,idf) = d(2)
          s(idf,i  ) =-d(2)
          s(i  ,idf) =-d(2)
          s(i  ,i  ) = d(2) 
        end if
c....   residual (contact force = penalty fac.* penetration + lambda)
        call contforce(ttl,tt,h2)
        if(dx.ge.0.d0) then  !  with respect to gap definition
          p(idf) = ttl
          p(i  ) =-ttl
        else 
          p(idf) =-ttl
          p(i  ) = ttl
        end if
      end if
c.... output contact forces (all = + -)
      if(isw.eq.4) then
                     write(iow,2001) n,dg,ttl
        if(ior.lt.0) write(*,2001)   n,dg,ttl
      endif
c.... update lagrange multiplier
      if(isw.eq.10) then
        epsdg = d(2)*dg
        call contupdate(epsdg,h2)
      endif
c.... plot  contact forces
      if(isw.eq.13) then
      if(nfp.ne.1) return
        klayf = 1
        if(flfp) then
          xmaxf = max(xmaxf,ttl)
          xminf = min(xminf,ttl)
          ccfp  = abs(ttl)
          cfp  = max(cfp,ccfp)
        else
          if(abs(ttl).lt.eps) return
c....     plot on mesh (on deformed mesh impossible  ds -> 0!)
          sn = xl(2,2)-xl(2,1)
          cs = xl(1,2)-xl(1,1)
          sl = dsqrt(cs*cs + sn*sn)
          sn = sn/sl
          cs = cs/sl
          xll(1,1) = xl(1,1)
          xll(2,1) = xl(2,1)
          xll(1,2) = xl(1,nel)
          xll(2,2) = xl(2,nel)
          xll(1,3) = xl(1,nel) - sn*ttl*cfp
          xll(2,3) = xl(2,nel) + cs*ttl*cfp
          xll(1,4) = xl(1,1)   - sn*ttl*cfp
          xll(2,4) = xl(2,1)   + cs*ttl*cfp
          if(ndm.eq.3) then
            xll(3,1) = xl(3,1)
            xll(3,2) = xl(3,2)
            xll(3,3) = xl(3,2)
            xll(3,4) = xl(3,1)
          endif
          call plot9s(ixl,xll,ndm,4)
        endif
      endif
      return
c
c.... formats
2000  format(5x,'node-on-node penalty contact element'/
     1   10x,'contact d.o.f.   ',i5/10x,'penalty parameter',e12.5/
     2   10x,'gap tolerance    ',e12.5/)
2001  format(5x,'Contact element ',i5,' gap = ',e13.5', Force = ',e13.5)
3000  format(3x,'Input: direction, penalty, gap tol'/5x,' >',\)
      end
c
      subroutine contforce(ttl,tt,hlamb)
      implicit double precision (a-h,o-z)
      ttl = tt + hlamb
      return
      end
c
      subroutine contupdate(epsdg,hlamb)
      implicit double precision (a-h,o-z)
      hlamb = hlamb + epsdg
      return
      end      