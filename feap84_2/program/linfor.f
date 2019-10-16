c$Id:$
      subroutine linfor(id,x,f,ndm,ndf,numnp,prth,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Output nodal values

c     Inputs:
c        id(ndf,*) - Boundary condition indicators
c        x(ndm,*)  - Nodal coordinats
c        ndm       - Mesh coordinate dimension
c        ndf       - Degree of freedoms/node
c        numnp     - Number mesh nodesw
c        prth      - Print header flag
c        isw       - Switch flag: 1 = force; 2 = displacement

c     Outputs:
c        f(ndf,*)  - Nodal force/displacement array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'fdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      character cd*6,fd(2)*6
      logical   prth, errck, pinput, flag
      integer   ndm,ndf,numnp,isw, n1,n2,n3,n4, i,n, count
      integer   id(ndf,numnp)
      real*8    x(ndm,numnp),f(ndf,numnp,2), td(8)

      save

      data      cd   /'-coord'/, fd /'-force','-displ'/

c     Input a record: Terminates on a blank

      count = 0
      if(ior.lt.0) then
        call pprint(' Input: node1,node2,inc,dir,a0,ax,ay,az')
        call pprint('   >')
      endif
      errck = pinput(td,8)

c     Compute the force/displacement using a linear relation on coords.

      n1 = nint(td(1))
      n2 = nint(td(2))
      if(n1.gt.0 .and. max(n1,n2).le.numnp ) then
        if(n2.eq.0) then
          n2 = n1
          n3 = 1
        else
          n3 = max(1,int(abs(td(3))))
          n3 = sign(n3,n2-n1)
        endif
        n4 = max(1,min(ndf,int(abs(td(4)))))

c       Output functional form for value

        call prtitl(prth)
        if(ior.lt.0) then
          write(*,2002) td(5),(td(i),i-5,i=6,5+ndm)
        endif
        write(iow,2002) td(5),(td(i),i-5,i=6,5+ndm)

        do n = n1,n2,n3
          if( mr(np(190)+n-1).ge.0 ) then

            if( (id(n4,n).eq.0. and. isw.eq.1) .or.
     &          (id(n4,n).ne.0. and. isw.eq.2)) then
              f(n4,n,isw) = td(5)
              do i = 1,ndm
                f(n4,n,isw) = f(n4,n,isw) + td(5+i)*x(i,n)
              end do ! i
              flag = .true.
            else
              flag = .false.
            endif
            if( flag ) then
              count = count - 1
              if(count.le.0) then
                call prtitl(prth)
                write(iow,2000) (i,cd,i=1,ndm),n4,fd(isw)
                if(ior.lt.0.and.pfr) then
                  write(*,2000) (i,cd,i=1,ndm),n4,fd(isw)
                endif
                count = 48
              endif
              write(iow,2001) n,(x(i,n),i=1,ndm),f(n4,n,isw)
              if(ior.lt.0.and.pfr) then
                write(*,2001) n,(x(i,n),i=1,ndm),f(n4,n,isw)
              endif
            endif
          endif
        end do ! n
      endif

c     Formats

2000  format('  N o d a l   F o r c e / D i s p l.',//'  Node',6(i6,a6))

2001  format(i6,1p,6e12.5)

2002  format('    Value = ',1p,1e11.4,3(' + ',1p,1e11.4,'*X_',i1:))

      end
