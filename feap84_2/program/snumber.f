c$Id:$
      subroutine snumber(n,nc,snode,lstnum,e2,adj,xadj,s,q,p)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Number nodes in component of graph for small profile &
c              rms wavefront

c     Input:
c     ------

c     n      - number of nodes in graph
c     nc     - number of nodes in component of graph
c     snode  - node at which numbering starts
c     lstnum - count of nodes which have already been numbered
c     e2     - twice tne number of edges in graph = xadj(n+1)-1
c     adj    - adjacency list for all nodes in graph
c            - list of length 2e where e is number of edges in
c              graph and 2e = xadj(n+1)-1
c     xadj   - index vector for adj
c            - nodes adjacent to node i are found in adj(j), where
c              j = xadj(i), xadj(i)+1, ..... , xadj(i+1)-1
c     s      - list giving distance of each node in this
c              component from end node
c     q      - list of nodes which are in this component
c            - also used to store queue of active or preactive nodes
c     p      - undefined

c     Output:
c     -------

c     n      - unchanged
c     nc     - unchanged
c     snode  - unchanged
c     lstnum - count of numbered nodes (input value incremented by nc)
c     e2     - unchanged
c     adj    - unchanged
c     xadj   - unchanged
c     s      - list of new node numbers
c            - new number for node i is s(i)
c     q      - not used
c     p      - not used

c     Notes:
c     ------

c     s also serves as list giving status of nodes during numbering:

c     s(i) gt 0 indicates node i is postactive
c     s(i) =  0 indicates node i is active
c     s(i) = -1 indicates node i is preactive
c     s(i) = -2 indicates node i is inactive

c     p is used to hold priorities for each node

c     Programmer:             Scott Sloan
c     -----------

c     Last modified:          1 March 1991    Scott Sloan
c     --------------

c     Copyright 1989:         Scott Sloan
c     ---------------         Department of Civil Engineering
c                             University of Newcastle
c                             NSW 2308
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    i,j,n
      integer    e2,nc,nn,w1,w2
      integer    nbr
      integer    next,node,prty
      integer    nabor,snode
      integer    addres,lstnum,maxprt
      integer    p(n),q(nc),s(n)
      integer    adj(e2)
      integer    xadj(n+1)

      parameter (w1=1, w2=2)

      save

c     Initialise priorities and status for each node in this component
c     Initial priority = w1*dist - w2*degree     where:
c     w1     = a positive weight
c     w2     = a positive weight
c     degree = initial current degree for node
c     dist   = distance of node from end node

      do i = 1,nc
        node    = q(i)
        p(node) = w1*s(node) - w2*(xadj(node+1) - xadj(node)+1)
        s(node) = -2
      end do ! i

c     Insert starting node in queue and assign it preactive status
c     nn is size of queue

      nn       = 1
      q(nn)    = snode
      s(snode) = -1

c     Loop while queue is not empty

      do while (nn.gt.0)

c       Scan queue for node with max priority

        addres = 1
        maxprt = p(q(1))
        do i = 2,nn
          prty = p(q(i))
          if(prty.gt.maxprt) then
            addres = i
            maxprt = prty
          end if
        end do ! i

c       Next is node to be numbered next

        next = q(addres)

c       Delete node next from queue

        q(addres) = q(nn)
        nn        = nn - 1
        if(s(next).eq.-1) then

c         Node next is preactive, examine its neighbours

          do i = xadj(next),xadj(next+1) - 1

c           Decrease current degree of neighbour by -1

            nbr    = adj(i)
            p(nbr) = p(nbr) + w2

c           Add neighbour to queue if it is inactive
c           Assign it preactive status

            if(s(nbr).eq.-2) then
              nn     = nn + 1
              q(nn)  = nbr
              s(nbr) = -1
            end if
          end do ! i
        end if

c       Store new node number for node next
c       Status for node next is now postactive

        lstnum  = lstnum + 1
        s(next) = lstnum

c       Search for preactive neighbours of node next

        do i = xadj(next),xadj(next+1)-1
          nbr = adj(i)
          if(s(nbr).eq.-1) then

c           Decrease current degree of preactive neighbour by -1
c           Assign neighbour an active status

            p(nbr) = p(nbr) + w2
            s(nbr) = 0

c           Loop over nodes adjacent to preactive neighbour

            do j = xadj(nbr),xadj(nbr+1)-1

c             Decrease current degree of adjacent node by -1

              nabor    = adj(j)
              p(nabor) = p(nabor) + w2
              if(s(nabor).eq.-2) then

c               Insert inactive node in queue with preactive status

                nn       = nn + 1
                q(nn)    = nabor
                s(nabor) = -1
               end if
            end do ! j
          end if
        end do ! i
      end do ! while

      end
