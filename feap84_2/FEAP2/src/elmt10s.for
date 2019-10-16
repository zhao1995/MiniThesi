      subroutine elmt10(d,ul,xl,ix,t,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c---------------------------------------------------------------------+
c.... point stiffness and/or mass/damping element                     |
c...  stiffness also for 2 nodes                                      | 
c                                                                     |           
c---------------------------------------------------------------------+
c     d( 1- 6): c_i                                                   |
c     d( 7-12): m_i                                                   |
c     d(13-18): d_i                                                   |
c---------------------------------------------------------------------+
c                                                                     |  
c.... extract from PCFEAP              WW 5/92                        |
c.... add terms for cmas               WW 5/94                        |
c.... add stiffness terms for 2nd node WW 2/03                        |
c.... add damping terms                WW 6/04                        |
c                                                                     |
c---------------------------------------------------------------------+
c     open                                                            |
c     forces from mass and damping terms??                            | 
c                                                                     |               
c---------------------------------------------------------------------+
      USE bdata
      USE eldata
      USE iofile
      implicit double precision(a-h,o-z)
      dimension ix(*),t(*),xl(ndm,*),d(*),ul(ndf,*),s(nst,*),p(*),
     1          sig(6),xx(3)
      dimension h1(*),h2(*),h3(*)
      go to (1,2,3,4,5,3,2,2,2,2,2,12,2), isw
      return
1     if(ior.lt.0) write(*,3000)
      call dinput(d,6)
      if(ior.lt.0) write(*,3001)
      call dinput(d( 7),6)
      if(ior.lt.0) write(*,3002)
      call dinput(d(13),6)
      if(ior.lt.0) write(  *,2000) (i,d(i)   ,i=1,ndf)
                   write(iow,2000) (i,d(i)   ,i=1,ndf)
      if(ior.lt.0) write(  *,2001) (i,d(i+ 6),i=1,ndf)
                   write(iow,2001) (i,d(i+ 6),i=1,ndf)
      if(ior.lt.0) write(  *,2002) (i,d(i+12),i=1,ndf)
                   write(iow,2002) (i,d(i+12),i=1,ndf)
2     return
3     if(ix(2).eq.0) then
c....   compute the point stiffness and residual
        do i = 1,ndf
          s(i,i) = d(i)
          p(i)   = p(i) - d(i)*ul(i,1)
        end do
      else
c....   calculate terms for two nodes
        do i=1,ndf
          i1=i
          i2=i+ndf
          s(i1,i1) =  d(i)
          s(i1,i2) = -d(i)
          s(i2,i1) = -d(i)
          s(i2,i2) =  d(i)
          p(i1)    =  p(i1) - d(i)*( ul(i,1)- ul(i,2))
          p(i2)    =  p(i2) - d(i)*(-ul(i,1)+ ul(i,2))
        end do
      end if 
      return
c.... compute the internal forces, forces from mass and damping terms??
c     sign with repect to normal vector component 2-1
4     if(ix(2).eq.0) then
        do i = 1,ndf
          sig(i) = d(i)*ul(i,1)
        end do   
      else 
        do i = 1,ndf
          sig(i) = d(i)*(ul(i,2)-ul(i,1))
        end do   
      end if
      mct = mct - 1
      if(mct.le.0) then
                     write(iow,2003) o,head,ndf
        if(ior.lt.0) write(*  ,2003) o,head,ndf
          mct = 50
        endif
        call pzero(xx,3)
        xx(1) =xl(1,1)
        xx(2) =xl(2,1)
        if(ndm.eq.3) xx(3) =xl(3,1)
                     write(iow,2004) n,ma,xx,(sig(i),i=1,ndf)
        if(ior.lt.0) write(*  ,2004) n,ma,xx,(sig(i),i=1,ndf)
      return
c.... compute the point mass
5     do i = 1,ndf
        s(i,i) = d(i+6)
        p(i)   = d(i+6)
      end do   
      return
c.... compute the point damping 
12    do i = 1,ndf
        s(i,i) = d(i+12)
        p(i)   = d(i+12)
      end do   
      return
c.... formats
2000  format(5x,'point stiffness/mass/damping element'/
     1      (10x,i5,' d.o.f. stiffness =',e13.5))
2001  format(10x,i5,' d.o.f. mass      =',e13.5)
2002  format(10x,i5,' d.o.f. damping   =',e13.5)
2003  format(a1,20a4,//5x,'S p r i n g   E l e m e n t',//,' el  mat',
     1   1x,'1-coord',1x,'2-coord',1x,'3-coord',5x,'forces  1 -',i2)
2004  format(2i4,3f8.3,6e13.5)
3000  format(3x,'Input: Nodal d.o.f. stiffness values'/5x,' >',\)
3001  format(3x,'Input: Nodal d.o.f. mass      values'/5x,' >',\)
3002  format(3x,'Input: Nodal d.o.f. damping   values'/5x,' >',\)
      end
