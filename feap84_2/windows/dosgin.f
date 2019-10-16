c$Id:$
      subroutine dosgin(ix,iy,butn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace MSFLIB by DFLIB                          08/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Graphical input of screen coordinates with mouse

c      Inputs:
c        none

c      Outputs:
c        ix       - x screen coordinate
c        iy       - y screen coordinate
c        butn     - Button pressed
c-----[--.----+----.----+----.-----------------------------------------]
      use DFLIB

      implicit  none

      include  'wdata.h'

      character(len=1) :: butn
      integer          :: mevent,shftl,shftr,key,result
      integer          :: ix,iy,x,y

      integer         idxl,idyl,jfill
      common /vgraph/ idxl,idyl,jfill

      save

c     Check for left or right button click

      shftl = mouse$ks_shift .or. mouse$ks_lbutton
      shftr = mouse$ks_shift .or. mouse$ks_rbutton

      mevent = mouse$lbuttondown .or. mouse$rbuttondown

      result = waitonmouseevent(mevent,key,x,y)

      if(key.eq.shftl .or. key.eq.shftr) then
        butn = 'm'
      elseif(key.eq.mouse$ks_lbutton) then
        butn = 'l'
      elseif(key.eq.mouse$ks_rbutton) then
        butn = 'r'
      else
        butn = 'e'
      endif

c     Return coordinates

      ix  = x*idxl
      iy  = 22000 - y*idyl

      end
