      subroutine unika
c-----------------------------------------------------------------------
c
c     Purpose: Draw Uni Karlsruhe Symbol
c
c     W.Wagner	2/01
c-----------------------------------------------------------------------
      USE pdata2
      implicit double precision(a-h,o-z)
      dimension x(25),y(25),xw(4),yw(4),xr(4),yr(4),xr2(4),yr2(4)
      integer ih1(5),ih2(5),ih3(5),ih4(5),ih5(5)
      integer iv1(5),iv2(5),iv3(5),iv4(5),iv5(5)
c     data for whitebox
      data xw/ 3.0,15.3,15.3, 3.0/
      data yw/59.3,59.3,65.9,65.9/
c     data for redbox
      data xr/15.3,49.6,49.6,15.3/
      data yr/59.3,59.3,65.9,65.9/
c     data for coordinates
      data x /
     +15.8,24.0,32.0,38.5,45.1,22.4,28.1,34.5,40.5,46.2,27.2,31.8,
     +37.4,43.3,48.7,30.4,34.6,40.1,46.1,52.3,31.9,36.5,42.4,49.1,56.7/
      data y /
     +69.5,66.6,62.6,57.9,51.6,71.7,69.5,66.7,63.7,60.2,74.3,72.7,
     +70.8,68.8,67.0,76.9,75.7,74.4,73.5,73.1,79.2,78.1,77.4,77.4,78.3/
c     data for nodes
      data ih1 / 1, 2, 3, 4, 5/
      data ih2 / 6, 7, 8, 9,10/
      data ih3 /11,12,13,14,15/
      data ih4 /16,17,18,19,20/
      data ih5 /21,22,23,24,25/
      data iv1 / 1, 6,11,16,21/
      data iv2 / 2, 7,12,17,22/
      data iv3 / 3, 8,13,18,23/
      data iv4 / 4, 9,14,19,24/
      data iv5 / 5,10,15,20,25/
c     data for 2nd redbox
      data xr2/25.6,36.9,35.4,25.6/
      data yr2/59.3,59.3,65.9,65.9/
c
      size =  100.0d0
c.... center
      xpos = 0.28
      ypos = 0.16
c.... white box
      call pppcol(16)
      call dplot(xw(1)/size+xpos,yw(1)/size+ypos,1,0)
      do  i = 2,4
      call dplot(xw(i)/size+xpos,yw(i)/size+ypos,2,0)
      enddo
      call dplot(xw(1)/size+xpos,yw(1)/size+ypos,2,0)
      call clpan

c.... red   box
      call pppcol(2)
      call dplot(xr(1)/size+xpos,yr(1)/size+ypos,1,0)
      do  i = 2,4
      call dplot(xr(i)/size+xpos,yr(i)/size+ypos,2,0)
      enddo
      call dplot(xr(1)/size+xpos,yr(1)/size+ypos,2,0)
      call clpan
      
c.... mesh      
      call pppcol(16)
c.... horizontal lines
      ip = ih1(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = ih1(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = ih2(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = ih2(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = ih3(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = ih3(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = ih4(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = ih4(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = ih5(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = ih5(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

c.... vertical lines
      ip = iv1(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = iv1(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = iv2(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = iv2(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = iv3(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = iv3(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = iv4(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = iv4(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

      ip = iv5(1)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,3,0)
      do  i = 2,5
      ip = iv5(i)
      call dplot(x(ip)/size+xpos,y(ip)/size+ypos,2,0)
      enddo

c.... 2nd red   box
      call pppcol(2)
      call dplot(xr2(1)/size+xpos,yr2(1)/size+ypos,1,0)
      do  i = 2,4
      call dplot(xr2(i)/size+xpos,yr2(i)/size+ypos,2,0)
      enddo
      call dplot(xr2(1)/size+xpos,yr2(1)/size+ypos,2,0)
      call clpan
      
c     text
      xp =  3.4/size + xpos
      yp = 61.2/size + ypos
      ip = 1
      call drawtxt(ip,xp,yp,32,3,5,'B A U')
      xp = 15.8/size + xpos
      call drawtxt(ip,xp,yp,16,3,11,'S T A T I K')
      call pltsiz(1)

      return
      end
