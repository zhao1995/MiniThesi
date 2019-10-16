c$Id:$
      subroutine umacr3(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'prt' from argument list                  09/07/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  User interface for adding solution command language
c                instructions.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'umac1.h'
      include  'comblk.h'
      include  'pointer.h'
      include  'cdata.h'
      include  'sdata.h'

      
      logical   pcomp
      character lct*15
      real*8    ctl(3)
      
      integer   i, dof, searched_node
      double precision max_val, max_tmp, dam_node, D_tol, deltaD
      double precision deltaD_curr, add_dam

      save

c     Set command word

      if(pcomp(uct,'mac3',4)) then      ! Usual    form
       uct = 'inda'                    ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation


c     Find node with less than ctl(1) % damage

c      max_val = ctl(1)
c      max_tmp = 0.0d0
c      do i = 1,numnp
c        dof = ndf-1
c        dam_node = hr(np(40)+(i-1)*ndf+dof)
c        if(max_tmp .lt. dam_node .and. dam_node .lt. max_val) then
c	  max_tmp = dam_node
c	  searched_node = i
c	endif
c      end do
      
c     Find node with more highest delta_D..	delta_D =  change of damage
      D_tol = ctl(1)/100.0
      deltaD = 0.0
      dof = ndf-1
      do i = 1,numnp
        dam_node = hr(np(40)+(i-1)*ndf+dof)
        deltaD_curr = hr(np(40)+ndf*numnp+ndf*i+dof)
        if(dam_node.gt.D_tol .and. deltaD_curr.gt.deltaD) then
          deltaD = deltaD_curr
          searched_node = i
        endif
       end do
     
c      write(*,*)"Some values:",
c     1 hr(np(40)+ndf*numnp+ndf*searched_node+dof)
      add_dam = 0.1*(1-hr(np(40)+(searched_node-1)*ndf+dof))
      open(10,file="batch2", status="unknown")
      write(10,*)"BATCh"
      write(10,*)" arcl,,",ndf
      write(10,*)" loop time 2 !40"
      write(10,*)"    time"
      write(10,*)"    nopr"
      write(10,*)"    loop iter 40"
      write(10,*)"      tang,,1"
      write(10,*)"    next iter"
      write(10,*)"    plot cont,5"
      write(10,*)"    prin"
!      write(10,*)"    pvie test"
      write(10,*)"  next time"
      write(10,*)"  inda,,70"
      write(10,*)"END"
      write(10,*)"1"
      write(10,"(I,A,I,A,F)")
     1   searched_node,",",ndf,",",add_dam
!      write(10,"(I,I,F7.5)")searched_node,",",ndf,",",ctl(2)
!      write(10,*)"INTEractive"
      if(hr(np(40)+(searched_node-1)*ndf+dof).lt.0.99) then
         write(10,*)"include,batch2"
      else
         write(10,*)"INTEractive"
      endif
      close(10)
     
      write(*,*)"Searched node:",searched_node
      write(*,*)"Delta Damage:",deltaD
!      pause

      endif
      
      
      end
