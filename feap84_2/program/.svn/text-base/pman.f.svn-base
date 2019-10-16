c$Id:$
      subroutine pman(name,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Echo request type                                21/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output user information as help

c      Inputs:
c         name      - Name of command to display
c         nn        - Type of command: 1=mesh, 2=solution, 3=plot

c      Outputs:
c         To screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      character name*4, htype(3)*12
      integer   nn

      save

      data      htype /'Mesh','Solution','Plot'/

c     Echo request

      if(ior.lt.0) then
        write(*,2001) htype(nn),name
      endif

c     Set for Mesh

      if(    nn.eq.1) then

        call pman_mesh(name)

c     Set for Solution Language

      elseif(nn.eq.2) then

        call pman_macr(name)

c     Set for Plot

      elseif(nn.eq.3) then

        call pman_plot(name)

c     Error

      else
        write(*,2002) name
      endif

c     Formats

2001  format(/4x,'Help requested for ',a,': ',a/)

2002  format(/4x,'No help available for ',a)

      end

      subroutine pman_mesh(name)

      implicit   none

      character  name*4
      logical    pcomp

      if(pcomp(name,'coor',4)) then

      elseif(pcomp(name,'elem',4)) then

      elseif(pcomp(name,'mate',4)) then

      elseif(pcomp(name,'boun',4)) then

      elseif(pcomp(name,'forc',4)) then

      elseif(pcomp(name,'temp',4)) then

      elseif(pcomp(name,'ene',3)) then

      elseif(pcomp(name,'prin',4)) then

      elseif(pcomp(name,'nopr',4)) then

      elseif(pcomp(name,'titl',4)) then

      elseif(pcomp(name,'bloc',4)) then

      elseif(pcomp(name,'polal',4)) then

      elseif(pcomp(name,'ebou',4)) then

      elseif(pcomp(name,'angl',4)) then

      elseif(pcomp(name,'sloa',4)) then

      elseif(pcomp(name,'cons',4)) then

      elseif(pcomp(name,'sphe',4)) then

      elseif(pcomp(name,'btem',4)) then

      elseif(pcomp(name,'icon',4)) then

      elseif(pcomp(name,'pars',4)) then

      elseif(pcomp(name,'nopa',4)) then

      elseif(pcomp(name,'trib',4)) then

      elseif(pcomp(name,'para',4)) then

      elseif(pcomp(name,'efor',4)) then

      elseif(pcomp(name,'eang',4)) then

      elseif(pcomp(name,'cbou',4)) then

      elseif(pcomp(name,'cfor',4)) then

      elseif(pcomp(name,'cang',4)) then

      elseif(pcomp(name,'foll',4)) then

      elseif(pcomp(name,'slav',4)) then

      elseif(pcomp(name,'rese',4)) then

      elseif(pcomp(name,'sblo',4)) then

      elseif(pcomp(name,'curv',4)) then

      elseif(pcomp(name,'rota',4)) then

      elseif(pcomp(name,'setn',4)) then

      elseif(pcomp(name,'setr',4)) then

      elseif(pcomp(name,'btra',4)) then

      elseif(pcomp(name,'fpro',4)) then

      elseif(pcomp(name,'cpro',4)) then

      elseif(pcomp(name,'regi',4)) then

      elseif(pcomp(name,'tran',4)) then

      elseif(pcomp(name,'damp',4)) then

      elseif(pcomp(name,'mass',4)) then

      elseif(pcomp(name,'stif',4)) then

      elseif(pcomp(name,'csur',4)) then

      elseif(pcomp(name,'ereg',4)) then

      elseif(pcomp(name,'reac',4)) then

      elseif(pcomp(name,'manu',4)) then

      elseif(pcomp(name,'body',4)) then

      elseif(pcomp(name,'glob',4)) then

      elseif(pcomp(name,'shif',4)) then

      elseif(pcomp(name,'disp',4)) then

      elseif(pcomp(name,'edis',4)) then

      elseif(pcomp(name,'cdis',4)) then

      elseif(pcomp(name,'debu',4)) then

      elseif(pcomp(name,'side',4)) then

      elseif(pcomp(name,'face',4)) then

      elseif(pcomp(name,'snod',4)) then

      elseif(pcomp(name,'blen',4)) then

      elseif(pcomp(name,'move',4)) then

      elseif(pcomp(name,'rigi',4)) then

      elseif(pcomp(name,'moda',4)) then

      elseif(pcomp(name,'flex',4)) then

      elseif(pcomp(name,'base',4)) then

      elseif(pcomp(name,'epro',4)) then

      elseif(pcomp(name,'mpro',4)) then

      elseif(pcomp(name,'loop',4)) then

      elseif(pcomp(name,'next',4)) then

      elseif(pcomp(name,'file',4)) then

      elseif(pcomp(name,'cdam',4)) then

      elseif(pcomp(name,'cmas',4)) then

      elseif(pcomp(name,'csti',4)) then

      elseif(pcomp(name,'ebas',4)) then

      elseif(pcomp(name,'cbas',4)) then

      elseif(pcomp(name,'eule',4)) then

      elseif(pcomp(name,'ceul',4)) then

      elseif(pcomp(name,'rfor',4)) then

      elseif(pcomp(name,'lfori',4)) then

      elseif(pcomp(name,'load',4)) then

      elseif(pcomp(name,'swee',4)) then

      elseif(pcomp(name,'spin',4)) then

      elseif(pcomp(name,'peri',4)) then

      elseif(pcomp(name,'expl',4)) then

      endif

      write(*,*) ' '

      end

      subroutine pman_macr(name)

      implicit   none

      character  name*4
      logical    pcomp

      if(pcomp(name,'stre',4)) then

        write(*,*) '   stre,,k1,k2,k3'
        write(*,*) '   stre all'
        write(*,*) '   stre node k1 k2 k3'
        write(*,*) '   stre gnod k1 k2 k3'
        write(*,*) '   stre cnod k1 k2 k3'
        write(*,*) '   stre coor nxt xt xtol'
        write(*,*) '   stre cerr k1,k2 k3'
        write(*,*) '   stre erro k1 k2 k3'
        write(*,*) '   stre cont k1 k2 k3'

      elseif(pcomp(name,'tang',4)) then

        write(*,*) '   tang: form symmetric tangent'
        write(*,*) '   tang,,1: form rhs and solve.'
        write(*,*) '   tang line 1 shift,value: line search on value'
        write(*,*) '   tang eigv 0 shift: with no mass/damping added'
        write(*,*) '   tang nume 0 shift: numerical tangent'

      elseif(pcomp(name,'utan',4)) then

        write(*,*) '   utan: form unsymmetric tangent'
        write(*,*) '   utan,,1: form rhs and solve.'
        write(*,*) '   utan line,1,shift,value: line search on value'

      elseif(pcomp(name,'form',4)) then

        write(*,*) '   form: form rhs residual'
        write(*,*) '   form acce: initial acceleration'
        write(*,*) '   form expl: explicit solution with lumped mass'
        write(*,*) '   form conv: check residual for convergence'
        write(*,*) '   form stre: update history only'

      elseif(pcomp(name,'resi',4)) then

        write(*,*) '   resid: form residual only'

      elseif(pcomp(name,'mass',4)) then

        write(*,*) '   mass lump: Lumped mass matrix used'
        write(*,*) '   mass cons: Consistent  mass matrix used'
        write(*,*) '   mass unsy: Unsymmetric mass matrix used'
        write(*,*) '   mass     : Same as consistent'

      elseif(pcomp(name,'reac',4)) then

        write(*,*) '   reac,,k1 k2 k3: output reactions'
        write(*,*) '   reac all: output all reactions'
        write(*,*) '   reac coor nxt xt xtol: at coordinate xt'
        write(*,*) '   reac node x1 x2 x3: at node (x1,x2,x3)'
        write(*,*) '   reac imag k1 k2 k3: complex imaginary'
        write(*,*) '   reac list k1: reactions for list k1'
        write(*,*) '   reac regi k1 k2: reactions for region k1'
        write(*,*) '   reac file: reactions to file fsav.rea'
        write(*,*) '   reac nopr: prevents any prints'

      elseif(pcomp(name,'chec',4)) then

        write(*,*) '   chec: check mesh for errors'
        write(*,*) '   chec init: check history data base changes'

      elseif(pcomp(name,'erro',4)) then

        write(*,*) '   erro:'
        write(*,*) '   erro stre: stress error for adaptivity'
        write(*,*) '   erro ener: energy error for adaptivity'

      elseif(pcomp(name,'damp',4)) then

        write(*,*) '   damp: consistent damping matrix'

      elseif(pcomp(name,'augm',4)) then

        write(*,*) '   augm,,value: value is factor on penalty'
        write(*,*) '   augm pena value: reset penalty parameter'

      elseif(pcomp(name,'geom',4)) then

        write(*,*) '   geom: geometric stiffness for eigenvalues'
        write(*,*) '   geom on/off: nonlinear geometric stiffness'

      elseif(pcomp(name,'dire',4)) then

        write(*,*) '   direct: profile solver (in-core)'
        write(*,*) '   direct block: profile solver (out of core)'
        write(*,*) '   direct sparse: sparse (in-core) symmetric'

      elseif(pcomp(name,'iter',4)) then

        write(*,*) '   iterative: solver (in-core)'

      elseif(pcomp(name,'expo',4)) then

        write(*,*) '   export: tangent and residual'

      elseif(pcomp(name,'impo',4)) then

        write(*,*) '   import: solution value'

      elseif(pcomp(name,'ntan',4)) then

        write(*,*) '   ntangent mate num: tangents for mate num'
        write(*,*) '   ntangent elem num: tangent for element num'
        write(*,*) '   ntangent off: turn off numerical computes'
        write(*,*) '   ntangent cont pair: numerical tangent for pair'

      elseif(pcomp(name,'base',4)) then

        write(*,*) '   base: base modes for multiple support excitation'

      elseif(pcomp(name,'jint',4)) then

        write(*,*) '   jint grad node: J-integral calculation'

      elseif(pcomp(name,'zzhu',4)) then

        write(*,*) '   zzhu,<off>,<ma>: Zienkiewicz-Zhu projections'
        write(*,*) '        - ma .eq. 0: project over all materials'
        write(*,*) '             .ne. 0: project over material ma'
        write(*,*) '        - off  - return to lumped projections'

      elseif(pcomp(name,'solv',4)) then

        write(*,*) '   solv: solve equations'
        write(*,*) '   solv,line,value: line search for ratios > value'

      elseif(pcomp(name,'dsol',4)) then

        write(*,*) '   dsolve: solve with diagonal matrix'

      elseif(pcomp(name,'scal',4)) then

        write(*,*) '   scal,<on,off>: diagonal scaling of tangent'

      elseif(pcomp(name,'hill',4)) then

        write(*,*) '   hill-mandel computations'
        write(*,*) '   hill tang:  - Compute tangent and stress'
        write(*,*) '   hill stre:  - Compute stress only'
        write(*,*) '   hill read:  - Read record from file'
        write(*,*) '   hill clos:  - Close file'

      elseif(pcomp(name,'tol',3)) then

        write(*,*) '   tol,,eval rval: energy (eval) residual (rval)'
        write(*,*) '   tol ener eval: energy (eval)'
        write(*,*) '   tol emax eval: max energy convergence to emax'
        write(*,*) '   tol iter itol atol: residual & absolute tols'

      elseif(pcomp(name,'dt',2)) then

        write(*,*) '   dt,,value: size of time increment'

      elseif(pcomp(name,'loop',4)) then

        write(*,*) '   loop,,number: perform number loops to next'
        write(*,*) '   loop,infinite: exit controlled by solution'

      elseif(pcomp(name,'next',4)) then

        write(*,*) '   next: termination of loop pair'

      elseif(pcomp(name,'prop',4)) then

        write(*,*) '   prop,,num1: input proportional function num1'
        write(*,*) '   prop,,num1,num2: proportional loads num1 to num2'
        write(*,*) '   prop,off: proportional loading, remove tables'

      elseif(pcomp(name,'time',4)) then

        write(*,*) '   time,,<tmax>: Advance time by dt'
        write(*,*) '                 quit after time > tmax'
        write(*,*) '   time,set,ttim_new: Time set to ttim_new'
        write(*,*) '   time,expl,<tmax>,c: Explicit critical time'
        write(*,*) '                        dt = c * dt_cr'

      elseif(pcomp(name,'prin',4)) then

        write(*,*) '   prin data: data from inputs output'
        write(*,*) '   prin command: solution commands on screen'
        write(*,*) '   prin,less: Prints shorter prompts to screen'
        write(*,*) '   prin <mass,cmas,geom>: output mass diagonal'
        write(*,*) '   prin <iden,lmas>: output identity/lump mass'
        write(*,*) '   prin,<tang,utan>: output tangent diagonals'
        write(*,*) '   prin,resi: output residual array'

      elseif(pcomp(name,'nopr',4)) then

        write(*,*) '   nopr data: data from inputs not output'
        write(*,*) '   nopr command: no solution commands on screen'
        write(*,*) '   nopr,less: Prints longer prompts to screen'

      elseif(pcomp(name,'tran',4)) then

        write(*,*) '   trans xxxx beta gamma alpha: transient solve'
        write(*,*) '       - xxxx = off:  Static solution'
        write(*,*) '       - xxxx = back: backward Euler'
        write(*,*) '       - xxxx = eule: Euler implicit'
        write(*,*) '       - xxxx = newm: Newmark method'
        write(*,*) '       - xxxx = hht:  Hilber-Hughes-Taylor'
        write(*,*) '       - xxxx = expl: Explicit Newmark'
        write(*,*) '       - xxxx = cent: Central difference explicit'

      elseif(pcomp(name,'init',4)) then

        write(*,*) '   init disp: set initial displacements'
        write(*,*) '   init rate: set initial rates'
        write(*,*) '   init velo: set initial velocity (same as rate)'
        write(*,*) '   init acce: set initial accelerations'
        write(*,*) '   init spin omg1 omg2 omg3: - set initial spins'
        write(*,*) '   init mate v1 v2 v3: velocity for material'
        write(*,*) '   init regi v1 v2 v3: velocity for region'

      elseif(pcomp(name,'iden',4)) then

        write(*,*) '   iden,,n1,n2: set dof n1 to n2 to unity'
        write(*,*) '                used as mass eigen matrix'

      elseif(pcomp(name,'newf',4)) then

        write(*,*) '   newf; set F0 to current force vector'
        write(*,*) '   newf zero: Set F0 to zero'

      elseif(pcomp(name,'iden',4)) then

        write(*,*) '   iden,,n1,n2: set dof n1 to n2 to unity'

      elseif(pcomp(name,'back',4)) then

        write(*,*) '   back,,dt: back-up to beginning of time'
        write(*,*) '             step, set new dt'

      elseif(pcomp(name,'debu',4)) then

        write(*,*) '   debug <on,off>: turn on/off debug prints'

      elseif(pcomp(name,'auto',4)) then

        write(*,*) '   auto,time,iopt,fact,repeat: set auto time step'
        write(*,*) '   auto,dt,dtmin,dtmax; set dt range'
        write(*,*) '   auto,mate: set by material behavior'
        write(*,*) '   auto,off: auto time control off '

      elseif(pcomp(name,'if',2)) then

        write(*,*) '   if expression: Start of if/else/endif'
        write(*,*) '     -expression controls use'

      elseif(pcomp(name,'else',4)) then

        write(*,*) '   else expression: Start of if/else/endif'
        write(*,*) '       -expression controls use'

      elseif(pcomp(name,'endi',4)) then

        write(*,*) '   endif: end of if/else/endif'

      elseif(pcomp(name,'jump',4)) then

        write(*,*) '   jump: marker for jump point'

      elseif(pcomp(name,'echo',4)) then

        write(*,*) '   echo <on/off>: batch solution commands appear'
        write(*,*) '                  to screen if on'

      elseif(pcomp(name,'conv',4)) then

        write(*,*) '   conv <on/off>: convergence to log file if on'

      elseif(pcomp(name,'omeg',4)) then

        write(*,*) '   omega   : frequency in radians/second'
        write(*,*) '   omega hz: frequency in Hertz'

      elseif(pcomp(name,'disp',4)) then

        write(*,*) '   disp,all: output all nodal displacements'
        write(*,*) '   disp,,k1,k2,k3: output displ. k1 to k2 step k3'
        write(*,*) '   disp,gnod,k1,k2,k3: output displ. for global'
        write(*,*) '                       nodes k1 to k2 step k3'
        write(*,*) '   disp,coor,k1,xt,xtol; output displ. all nodes'
        write(*,*) '                         where x-k1=xt-xtol'
        write(*,*) '   disp,node,x1,x2,x3: output displ. for node'
        write(*,*) '                       closest to x1,x2,x3'
        write(*,*) '   disp,imag,k1,k2,k3: output all imaginary'
        write(*,*) '   disp,cmpl,k1,k2,k3: output all real/imag'
        write(*,*) '   disp,list,k1: output displacements in list k1'
        write(*,*) '   disp,glob: output all global displacements'
        write(*,*) '   disp,glob,k1,k2,k3: output global k1 to k2 in k3'

      elseif(pcomp(name,'velo',4)) then

        write(*,*) '   velo,all: output all nodal velocities'
        write(*,*) '   velo,,k1,k2,k3: output veloc. k1 to k2 step k3'
        write(*,*) '   velo,gnod,k1,k2,k3: output veloc. for global'
        write(*,*) '                       nodes k1 to k2 step k3'
        write(*,*) '   velo,coor,k1,xt,xtol; output veloc. all nodes'
        write(*,*) '                         where x-k1=xt-xtol'
        write(*,*) '   velo,node,x1,x2,x3: output veloc. for node'
        write(*,*) '                       closest to x1,x2,x3'
        write(*,*) '   velo,imag,k1,k2,k3: output all imaginary'
        write(*,*) '   velo,cmpl,k1,k2,k3: output all real/imag'
        write(*,*) '   velo,list,k1: output velocities in list k1'
        write(*,*) '   velo,glob: output all global velocities'
        write(*,*) '   velo,glob,k1,k2,k3: output global k1 to k2 in k3'

      elseif(pcomp(name,'acce',4)) then

        write(*,*) '   acce,all: output all nodal accelerations'
        write(*,*) '   acce,,k1,k2,k3: output accel. k1 to k2 step k3'
        write(*,*) '   acce,gnod,k1,k2,k3: output accel. for global'
        write(*,*) '                       nodes k1 to k2 step k3'
        write(*,*) '   acce,coor,k1,xt,xtol; output accel. all nodes'
        write(*,*) '                         where x-k1=xt-xtol'
        write(*,*) '   acce,node,x1,x2,x3: output accel. for node'
        write(*,*) '                       closest to x1,x2,x3'
        write(*,*) '   acce,imag,k1,k2,k3: output all imaginary'
        write(*,*) '   acce,cmpl,k1,k2,k3: output all real/imag'
        write(*,*) '   acce,list,k1: output accelerations in list k1'
        write(*,*) '   acce,glob: output all global accelerations'
        write(*,*) '   acce,glob,k1,k2,k3: output global k1 to k2 in k3'

      elseif(pcomp(name,'mesh',4)) then

        write(*,*) '   mesh: Reenter mesh generation phase'
        write(*,*) '   mesh,filename: Read the data from  filename'

      elseif(pcomp(name,'opti',4)) then

        write(*,*) '   optimize:  node numbering then reset profile'
        write(*,*) '   opti off:  turn off profile optimization'
        write(*,*) '   opti cont: optimization during contact solution'
        write(*,*) '   opti hoit: Set to use optimizer by Wilson/Hoit'
        write(*,*) '   opti sloan: Set to use Sloan algorithm (default)'

      elseif(pcomp(name,'plot',4)) then

        write(*,*) '   plot: enter interactive plot mode'
        write(*,*) '   plot <optn k1 k2 k3>: see plot manual'

      elseif(pcomp(name,'subs',4)) then

        write(*,*) '   subs <prin> k1 <k2>: subspace for k1 eigenpairs'
        write(*,*) '                  - k2 number guard vectors'
        write(*,*) '                  - prin ouputs projecte matrices'

      elseif(pcomp(name,'writ',4)) then

        write(*,*) '   write fnam: open write file named fnam'
        write(*,*) '   write disp: write displacements to fnam'
        write(*,*) '   write stre: write nodal streses to fnam'
        write(*,*) '   write eige: write eigenpairs to fnam'
        write(*,*) '   write wind: rewind  fnam'
        write(*,*) '   write clos: close  fnam'

      elseif(pcomp(name,'read',4)) then

        write(*,*) '   read fnam: open read file named fnam'
        write(*,*) '   read disp: read displacements to fnam'
        write(*,*) '   read stre: read nodal streses to fnam'
        write(*,*) '   read eige: read eigenpairs to fnam'
        write(*,*) '   read wind: rewind  fnam'
        write(*,*) '   read clos: close  fnam'

      elseif(pcomp(name,'cont',4)) then

        write(*,*) '   cont check:   Geometry check now'
        write(*,*) '   cont nocheck: Geometry check at each iteration'
        write(*,*) '   cont on:  Enable  contact'
        write(*,*) '   cont off: Disable contact'
        write(*,*) '   cont pena n,pen: n = pair number; pen = penalty'
        write(*,*) '   cont friction  : Contact friction ON'
        write(*,*) '   cont nofriction: Contact friction OFF'

      elseif(pcomp(name,'rest',4)) then

        write(*,*) '   restart ext_name kk: Restart file extension'
        write(*,*) '                        ext_name + kk'

      elseif(pcomp(name,'bfgs',4)) then

        write(*,*) '   bfgs,,k1,k2,k3: BFGS soln k1 = no. steps;'
        write(*,*) '                               k2 = line search tol'
        write(*,*) '                               k3 = bfgs energy tol'

      elseif(pcomp(name,'arcl',4)) then

        write(*,*) '   arcl,,kflag,lflag: set arc length parameters'
        write(*,*) '   arcl,add,k1,tau: add eigvenvector k1, amount tau'
        write(*,*) '   arcl,chec,k1: check bifurcation using eigv. k1'
        write(*,*) '   arcl,off: set arclength to off'

      elseif(pcomp(name,'save',4)) then

        write(*,*) '   save ext_name: save Restart file with file'
        write(*,*) '                  extension ext_name'

      elseif(pcomp(name,'paus',4)) then

        write(*,*) '   pause: Pause on no convergence'

      elseif(pcomp(name,'eige',4)) then

        write(*,*) '   eige <vect k1>: eigenpairs for element k1'
        write(*,*) '   eige <mass k1>: mass eigenvalues for element k1'
        write(*,*) '              k1 <= 0 for last element'

      elseif(pcomp(name,'expl',4)) then

        write(*,*) '   explicit: do explicit solution'

      elseif(pcomp(name,'acti',4)) then

        write(*,*) '   acti all:       activate all regions'
        write(*,*) '   acti ,k1,k2,k3: activate regions: n = k1,k2,k3'
        write(*,*) '   acti init,k1,k2,k3: (also initialize strains)'

      elseif(pcomp(name,'deac',4)) then

        write(*,*) '   deac all:       deactivate all regions'
        write(*,*) '   deac: show current parameters on screen'
        write(*,*) '   deac,,k1,k2,k3: deactivate regions: n = k1,k2,k3'

      elseif(pcomp(name,'zero',4)) then

        write(*,*) '   zero: Zero entire problem to start over'
        write(*,*) '   zero rate: Zero transient vectors (VEL)'
        write(*,*) '   zero regi n1: Zero displacements in region n1'

      elseif(pcomp(name,'epri',4)) then

        write(*,*) '   eprint: output last element tangent and residual'

      elseif(pcomp(name,'moda',4)) then

        write(*,*) '   modal: Solve linear problem by modal method'

      elseif(pcomp(name,'eigv',4)) then

        write(*,*) '   eigv dofs/dof-list: DOFS for comp (1=on; 0=off)'
        write(*,*) '   eigv all nnn: Output eigenvector nnn (all)'
        write(*,*) '   eigv coor k1 xt nnn: output forr x_k1 = xt'
        write(*,*) '   eigv list k1,nnn: Output nnn using list k1'
        write(*,*) '   eigv nnn k1,k2,k3: Output nnn nd k1-k2 @ inc k3'

      elseif(pcomp(name,'rayl',4)) then

        write(*,*) '   rayleigh,,a0,a1: Rayleigh damping specification'
        write(*,*) '   rayl,freq,zeta,w1,w2: frequencies/damping ratio'

      elseif(pcomp(name,'csso',4)) then

        write(*,*) '   cxsoleve,,omega: Complex solution for omega'
        write(*,*) '                    in rad/sec'

      elseif(pcomp(name,'broy',4)) then

        write(*,*) '   broy,,iter: Broyden update for unsymmric prob'

      elseif(pcomp(name,'rect',4)) then

        write(*,*) '   rect: cartesian output of displ and stress'

      elseif(pcomp(name,'cyli',4)) then

        write(*,*) '   cyli: Cylindrical output of displ and stress'

      elseif(pcomp(name,'sphe',4)) then

        write(*,*) '   sphe: Spherical output of displ and stress'

      elseif(pcomp(name,'forc',4)) then

        write(*,*) '   force: force prints'
        write(*,*) '   force all: output all'
        write(*,*) '   force coor dir,x_dir: output for edge coordinate'
        write(*,*) '   force node x(i) i=1,ndm: output at node at x(i)'
        write(*,*) '   force k1 k2 k3: output for k1 to k2 in incr k3'

      elseif(pcomp(name,'tie',3)) then

        write(*,*) '   tie,,off: Untie all mesh connections'

      elseif(pcomp(name,'real',4)) then

        write(*,*) '   real: Set to output real part of complex'

      elseif(pcomp(name,'imag',4)) then

        write(*,*) '   imag: Set to output imaginary part of complex'

      elseif(pcomp(name,'outm',4)) then

        write(*,*) '   outm:      Output renumbered input file'
        write(*,*) '   outm cont: Output contact slide-line data'
        write(*,*) '   outm defo: Output renumbered deformed mesh'
        write(*,*) '   outm doma: Output mesh for domains'
        write(*,*) '   outm elem: Output mesh for single element'
        write(*,*) '   outm bina: Output renumbered binarynput file'

      elseif(pcomp(name,'renu',4)) then

        write(*,*) '   renu: Output list of old/renumbered nodes'

      elseif(pcomp(name,'show',4)) then

        write(*,*) '   show: show current solution parameters'
        write(*,*) '   show cont: show user contact types and variables'
        write(*,*) '   show dict: show dictionary of array allocation'
        write(*,*) '   show elem: show user element types'
        write(*,*) '   show mate: show material type use'
        write(*,*) '   show part: show partition data.'

      elseif(pcomp(name,'scre',4)) then

        write(*,*) '   screen,<on,off>: Set plot to screen on/off'

      elseif(pcomp(name,'comm',4)) then

        write(*,*) '   comment <message>:to screen when in batch mode'

      elseif(pcomp(name,'outp',4)) then

        write(*,*) '   output <array>:Output sparse matrix format'
        write(*,*) '          <array>  [tang,utan,lmas,mass,cmas,umas'
        write(*,*) '                    damp,cdam,udam,dr,form]'
        write(*,*) '         filename = same as array specified'

      elseif(pcomp(name,'list',4)) then

        write(*,*) '   list,,i: i is list number'
        write(*,*) '            add list of nodes to output'

      elseif(pcomp(name,'tplo',4)) then

        write(*,*) '   tplot,,<interval> - time history plots'
        write(*,*) '   options: disp n1 n2 x y z  n1 = dof'
        write(*,*) '            velo n1 n2 x y z  n2 = node'
        write(*,*) '            acce n1 n2 x y z'
        write(*,*) '            reac n1 n2 x y z'
        write(*,*) '            stre n1 n2 x y z  n1 = component'
        write(*,*) '            hist n1 n2 x y z  n2 = element'
        write(*,*) '            elem n1 n2 x y z'
        write(*,*) '            user n1 n2 x y z'
        write(*,*) '            cont n1 n2 x y z'
        write(*,*) '            arcl n1 n2'
        write(*,*) '            rsum n1 n2'
        write(*,*) '            sums n1 n2'
        write(*,*) '            ener: Energy and momenta'
        write(*,*) '            show: Adds list to output file'

      elseif(pcomp(name,'para',4)) then

        write(*,*) '   parameter,name,value: set name to value'

      elseif(pcomp(name,'func',4)) then

        write(*,*) '   func,fname: execute function fname'

      elseif(pcomp(name,'part',4)) then

        write(*,*) '   part,,i: Set partition i to active'

      elseif(pcomp(name,'mono',4)) then

        write(*,*) '   mono: Monolithic, no partitions'

      elseif(pcomp(name,'capt',4)) then

        write(*,*) '   caption label: caption for next contout'

      elseif(pcomp(name,'get',3)) then

        write(*,*) '   get pname: get value of parameter pname'

      endif

      write(*,*) ' '

      end

      subroutine pman_plot(name)

      implicit   none

      character  name*4
      logical    pcomp

      if(pcomp(name,'fram',4)) then

        write(*,*) '   fram ifrm: ifrm = 0 for whole screen'
        write(*,*) '              ifrm = 1 for upper-left'
        write(*,*) '              ifrm = 2 for upper-right'
        write(*,*) '              ifrm = 3 for lower-left'
        write(*,*) '              ifrm = 4 for lower-right'
        write(*,*) '              ifrm = 5 for legend Box'

      elseif(pcomp(name,'wipe',4)) then


        write(*,*) '   wipe ifrm: ifrm = 0 for whole screen'
        write(*,*) '              ifrm = 1 for upper-left'
        write(*,*) '              ifrm = 2 for upper-right'
        write(*,*) '              ifrm = 3 for lower-left'
        write(*,*) '              ifrm = 4 for lower-right'
        write(*,*) '              ifrm = 5 for legend Box'

      elseif(pcomp(name,'fact',4)) then

        write(*,*) '   fact value: multiply plot by value'

      elseif(pcomp(name,'cent',4)) then

        write(*,*) '   cent s1 s2: center plot at (s1,s2)'
        write(*,*) '               si = 0.5 is default'

      elseif(pcomp(name,'cart',4)) then

        write(*,*) '   cart: cartesian view'

      elseif(pcomp(name,'line',4)) then

        write(*,*) '   line,value,width: value is type'
        write(*,*) '                     deviced dependent'

      elseif(pcomp(name,'symm',4)) then

        write(*,*) '   symm x1 x2 x3: symmetry reflection'

      elseif(pcomp(name,'cont',4)) then

        write(*,*) '   contour k1 k2 k3: k1 = solution component'
        write(*,*) '                     k2 = # lines (fill if <= 0)'
        write(*,*) '                     k3 = 0: superpose mesh'
        write(*,*) '                     k3 > 0: no mesh'

      elseif(pcomp(name,'cwir',4)) then

        write(*,*) '   cwir k1: k1 = solution component number'


      elseif(pcomp(name,'velo',4)) then

        write(*,*) '   velocity k1 k2 k3: k1 = velocity component'
        write(*,*) '                      k2 = # lines (fill if <= 0)'
        write(*,*) '                      k3 = 0: superpose mesh'
        write(*,*) '                      k3 > 0: no mesh'

      elseif(pcomp(name,'vwir',4)) then

        write(*,*) '   vwir k1: k1 = velocity component number'


      elseif(pcomp(name,'acce',4)) then

        write(*,*) '   accelerate k1 k2 k3: k1 = component'
        write(*,*) '                        k2 = # lines (fill if <= 0)'
        write(*,*) '                        k3 = 0: superpose mesh'
        write(*,*) '                        k3 > 0: no mesh'

      elseif(pcomp(name,'awir',4)) then

        write(*,*) '   awir k1: k1 = acceleration component number'

      elseif(pcomp(name,'ndat',4)) then

        write(*,*) '   ndata xxxx no: xxxx = [disp stre pstre]'
        write(*,*) '                  no = filenumber'

      elseif(pcomp(name,'outl',4)) then

        write(*,*) '   outline: outline of mesh region'

      elseif(pcomp(name,'load',4)) then

        write(*,*) '   load k1 k2: Plot load vector'
        write(*,*) '               k1 > 0 tip on node'
        write(*,*) '               k1 < 0 tail on node'
        write(*,*) '               k2 scale factor'

      elseif(pcomp(name,'rfor',4)) then

        write(*,*) '   rfor k1 k2: radial follower load'
        write(*,*) '               k1 > 0 tip on node'
        write(*,*) '               k1 < 0 tail on node'
        write(*,*) '               k2 scale factor'

      elseif(pcomp(name,'disp',4)) then

        write(*,*) '   disp k1 k2: displacement vector'
        write(*,*) '               k1 > 0 tip on node'
        write(*,*) '               k1 < 0 tail on node'
        write(*,*) '               k2 scale factor'

      elseif(pcomp(name,'mesh',4)) then

        write(*,*) '   mesh k1: mesh of current material number'
        write(*,*) '            k1 < 0 alters line color'

      elseif(pcomp(name,'stre',4)) then

        write(*,*) '   stress k1 k2 k3: k1 = component to contour'
        write(*,*) '                    k2 = # lines (fill if <,= 0)'
        write(*,*) '                    k3 = 0: superpose mesh'
        write(*,*) '                    k3 > 0: no mesh'

      elseif(pcomp(name,'stra',4)) then

        write(*,*) '   strain k1 k2 k3: k1 = component to contour'
        write(*,*) '                    k2 = # lines (fill if <,= 0)'
        write(*,*) '                    k3 = 0: superpose mesh'
        write(*,*) '                    k3 > 0: no mesh'

      elseif(pcomp(name,'flux',4)) then

        write(*,*) '   flux k1 k2 k3: k1 = component to contour'
        write(*,*) '                  k2 = # lines (fill if <,= 0)'
        write(*,*) '                  k3 = 0: superpose mesh'
        write(*,*) '                  k3 > 0: no mesh'

      elseif(pcomp(name,'hist',4)) then

        write(*,*) '   history k1 k2 k3: k1 = component to contour'
        write(*,*) '                     k2 = # lines (fill if <,= 0)'
        write(*,*) '                     k3 = 0: superpose mesh'
        write(*,*) '                     k3 > 0: no mesh'

      elseif(pcomp(name,'pstr',4)) then

        write(*,*) '   pstress k1 k2 k3: k1 = component to contour'
        write(*,*) '                     k2 = # lines (fill if <,= 0)'
        write(*,*) '                     k3 = 0: superpose mesh'
        write(*,*) '                     k3 > 0: no mesh'

      elseif(pcomp(name,'estr',4)) then

        write(*,*) '   estress k1 k2 k3: k1 = component to contour'
        write(*,*) '                     k2 = # lines (fill if <,= 0)'
        write(*,*) '                     k3 = 0: superpose mesh'
        write(*,*) '                     k3 > 0: no mesh'

      elseif(pcomp(name,'swir',4)) then

        write(*,*) '   swire k1: k1 = component to contour'

      elseif(pcomp(name,'pwir',4)) then

        write(*,*) '   pwire k1: k1 = component to contour'

      elseif(pcomp(name,'ewir',4)) then

        write(*,*) '   ewire k1: k1 = component to contour'

      elseif(pcomp(name,'node',4)) then

        write(*,*) '   node k1 k2; Plot nodes k1 to k2 & numbers'
        write(*,*) '           k1 = 0 show all nodes & numbers'
        write(*,*) '           k1 < 0 show all nodes, no numbers'

      elseif(pcomp(name,'boun',4)) then

        write(*,*) '   boun <k1>: Plot nodes k1 to k2 & numbers'
        write(*,*) '              k1 = 0 restraints up to dof = 3'
        write(*,*) '              k1 > 0 restraints of dof k1'

      elseif(pcomp(name,'elem',4)) then

        write(*,*) '   element k1 k2: Plot elmts k1 to k2 numbers'
        write(*,*) '              k1 = k2 = 0, plot all'

      elseif(pcomp(name,'zoom',4)) then

        write(*,*) '   zoom k1 k2: Set window between nodes k1 & k2'

      elseif(pcomp(name,'colo',4)) then

        write(*,*) '   color k1 k2: set color'
        write(*,*) '         k1 <  0  greyscale postscript'
        write(*,*) '         k1 >= 0  color postscript'
        write(*,*) '         k2 =  0  standard color order'
        write(*,*) '         k2 != 0  reversed color order'

      elseif(pcomp(name,'fill',4)) then

        write(*,*) '   fill k1: fill current material in color k1'

      elseif(pcomp(name,'text',4)) then

        write(*,*) '   text k1 x y: Put text at (x,y)'
        write(*,*) '               (0<x<1.4; 0<y<1) in color k1'

      elseif(pcomp(name,'size',4)) then

        write(*,*) '   size k1: set text size to k1'
        write(*,*) '            device dependent'

      elseif(pcomp(name,'cvar',4)) then

        write(*,*) '   cvar k1 k2 k3: contact variable k1 from'
        write(*,*) '                  history vector k2 for pair k3'
        write(*,*) '                  k1 = 0 -> def = 1'
        write(*,*) '                  k2 = 0 -> def = 2, plot from CH2'
        write(*,*) '                  k3 = 0 -> plot for all pairs'

      elseif(pcomp(name,'eigv',4)) then

        write(*,*) '   eigv k1 k2 k3: eigenvector k1 & material k2'
        write(*,*) '                  k2 = 0 for all materials'
        write(*,*) '                  k3 = dof to contour'

      elseif(pcomp(name,'bord',4)) then

        write(*,*) '   border k1: plot border in color k1'

      elseif(pcomp(name,'scal',4)) then

        write(*,*) '   scale cs: Set plot scale to cs'

      elseif(pcomp(name,'axis',4)) then

        write(*,*) '   axis x y: plot axis at screen coords (x,y)'

      elseif(pcomp(name,'pers',4)) then

        write(*,*) '   pers k1: perspective view'
        write(*,*) '            k1 = 0: at default location'
        write(*,*) '            k1 > 0: input perspective data'

      elseif(pcomp(name,'show',4)) then

        write(*,*) '   show: display current parameters'

      elseif(pcomp(name,'hide',4)) then

        write(*,*) '   hide: hide mesh not visible'

      elseif(pcomp(name,'defo',4)) then

        write(*,*) '   deform scale resize escale: deform plots'
        write(*,*) '          scale: displ. multiplier (default=1)'
        write(*,*) '          resize > 0: do not rescale plot'
        write(*,*) '          escale: eigen multiplier (default=1)'

      elseif(pcomp(name,'unde',4)) then

        write(*,*) '   undeform,,resize escale: deform plots'
        write(*,*) '          resize > 0: do not rescale plot'
        write(*,*) '          escale: eigen multiplier (default=1)'

      elseif(pcomp(name,'post',4)) then

        write(*,*) '   post k1: PostScrip outputs'
        write(*,*) '            k1 = 0: portrait'
        write(*,*) '            k1 = 1: landscape'
        write(*,*) '            give before and after plot'

      elseif(pcomp(name,'reac',4)) then

        write(*,*) '   reaction k1 k2: nodal reactions'
        write(*,*) '                   k1 > 0 tip at node'
        write(*,*) '                   k2 = scale value'

      elseif(pcomp(name,'eige',4)) then

        write(*,*) '   eige k1 c2: last element eigenvector k1'
        write(*,*) '               c2 = color'

      elseif(pcomp(name,'mate',4)) then

        write(*,*) '   mate k1: restrict plots to material k1'

      elseif(pcomp(name,'back',4)) then

        write(*,*) '   back k1: set postscript backgound color'
        write(*,*) '            k1 > 0: white'
        write(*,*) '            k1 = 0: black'

      elseif(pcomp(name,'clip',4)) then

        write(*,*) '   clip xdir xmin xmax: clip between xmin & xmax'
        write(*,*) '                        in direction xdir'

      elseif(pcomp(name,'titl',4)) then

        write(*,*) '   title : place problem title on plot'

      elseif(pcomp(name,'mark',4)) then

        write(*,*) '   mark k1: place max/min on contour plot'
        write(*,*) '            k1 > 0: turns off marks'

      elseif(pcomp(name,'caption',4)) then

        write(*,*) '   caption text: define caption for next contour'

      elseif(pcomp(name,'pick',4)) then

        write(*,*) '   pick: use mouse to pick plot window'

      elseif(pcomp(name,'pbou',4)) then

        write(*,*) '   pbou idir: use mouse to set boundary conditions'
        write(*,*) '              in direction idir'

      elseif(pcomp(name,'pfor',4)) then

        write(*,*) '   pfor idir value: use mouse to set node force'
        write(*,*) '                    in direction idir with value'

      elseif(pcomp(name,'pnod',4)) then

        write(*,*) '   pnod: pick a node to get number'

      elseif(pcomp(name,'quad',4)) then

        write(*,*) '   quad: set plot quadrant'

      elseif(pcomp(name,'real',4)) then

        write(*,*) '   real: set component real (complex)'

      elseif(pcomp(name,'imag',4)) then

        write(*,*) '   imag: set component imaginary (complex)'

      elseif(pcomp(name,'eyes',4)) then

        write(*,*) '   eyes: set perspective from eyes region'

      elseif(pcomp(name,'dofs',4)) then

        write(*,*) '   dofs k1 k2 k3: set ki > 0 for dof to plot'

      elseif(pcomp(name,'prof',4)) then

        write(*,*) '   profile k1: upper profile if k1 = 0'
        write(*,*) '               total profile if k1 > 0'

      elseif(pcomp(name,'prax',4)) then

        write(*,*) '   prax,k1,k2,k3: principal stress directions'
        write(*,*) '                  k1 principal stress component'
        write(*,*) '                  k2 < 0 negative compnents only'
        write(*,*) '                  k2 > 0 positive compnents only'
        write(*,*) '                  k2 = 0 all compnents'
        write(*,*) '                  k3 > 0 optional colorEshift'
        write(*,*) '                  k3 < 0 scaling factor'

      elseif(pcomp(name,'pair',4)) then

        write(*,*) '   pair k1 k2 k3: contact pairs k1 to k2'
        write(*,*) '                  k3 = 0 -> slave & master surf'
        write(*,*) '                  k3 = 1 -> slave surf'
        write(*,*) '                  k3 = 2 -> master surf'

      elseif(pcomp(name,'clea',4)) then

        write(*,*) '   clear: wipe screen area'

      elseif(pcomp(name,'dplo',4)) then

        write(*,*) '   dplot k1: use mouse to pick line for displ. plot'

      elseif(pcomp(name,'vplo',4)) then

        write(*,*) '   vplot k1: use mouse to pick line for veloc. plot'

      elseif(pcomp(name,'aplo',4)) then

        write(*,*) '   aplot k1: use mouse to pick line for accel. plot'

      elseif(pcomp(name,'splo',4)) then

        write(*,*) '   splot k1: use mouse to pick line for stress plot'

      elseif(pcomp(name,'prom',4)) then

        write(*,*) '   prompt <on/off>: set gaphics prompts on/off'

      elseif(pcomp(name,'defa',4)) then

        write(*,*) '   dafault <on/off>: set gaphics defaults on/off'

      elseif(pcomp(name,'labe',4)) then

        write(*,*) '   label <on/off>: turn on/off contour label area'

      elseif(pcomp(name,'exno',4)) then

        write(*,*) '   exno k1: Plot exterior nodes'
        write(*,*) '            k1 = 0: nodes & numbers'
        write(*,*) '            k1 < 0: nodes no numbers'

      elseif(pcomp(name,'wind',4)) then

        write(*,*) '   window k1: Set window to k1'

      elseif(pcomp(name,'logo',4)) then

        write(*,*) '   logo: place feap logo'

      elseif(pcomp(name,'time',4)) then

        write(*,*) '   time <on/off>: turn time label on/off'

      elseif(pcomp(name,'rang',4)) then

        write(*,*) '   range pmin pmax: set range of contours between'
        write(*,*) '                    pmin and pmax'

      elseif(pcomp(name,'nora',4)) then

        write(*,*) '   norang: turn off range setting'

      elseif(pcomp(name,'rect',4)) then

        write(*,*) '   rectangular: use Cartesian components'

      elseif(pcomp(name,'cyli',4)) then

        write(*,*) '   cyli x0 y0 z0: use cylindridal components'
        write(*,*) '                  origin at: (x0 y0 z0)'

      elseif(pcomp(name,'sphe',4)) then

        write(*,*) '   sphe x0 y0 z0: use spherical components'
        write(*,*) '                  origin at: (x0 y0 z0)'

      elseif(pcomp(name,'full',4)) then

        write(*,*) '   full: use full screen (Windows only)'

      elseif(pcomp(name,'regi',4)) then

        write(*,*) '   region <n1 n2>: plot regions n1 to n2 only'
        write(*,*) '                   All regions if n1 = n2 = 0'

      elseif(pcomp(name,'norm',4)) then

        write(*,*) '   normal,,n2: plot normals of frame elements'
        write(*,*) '               n2 scales size of vectors'

      elseif(pcomp(name,'inte',4)) then

        write(*,*) '   interval <int>: plot contours at int time steps'

      elseif(pcomp(name,'swee',4)) then

        write(*,*) '   sweep <ang,inc>: 3d axisymmetric mesh of 2d'
        write(*,*) '                    problem. ang = angle in degrees'
        write(*,*) '                    inc = number sweep increments.'
        write(*,*) '                    Plot in perspective view'

      elseif(pcomp(name,'edef',4)) then

        write(*,*) '   edeform c1: set eigenvector scale factor to c1'

      elseif(pcomp(name,'elpl',4)) then

        write(*,*) '   elplot ne: plot mesh of element ne only'

      elseif(pcomp(name,'stri',4)) then

        write(*,*) '   string,,k1: plot text at bottom of plot'
        write(*,*) '               k1 device dependent text size'
        write(*,*) '               N.B. text input after command'

      elseif(pcomp(name,'acti',4)) then

        write(*,*) '   activate k1 k2 k3: Activate regions k1 to k2 in'
        write(*,*) '                      k3 increments'

      elseif(pcomp(name,'deac',4)) then

        write(*,*) '   deactivate k1 k2 k3: Deactivate regions k1 to k2'
        write(*,*) '                        in k3 increments'

      elseif(pcomp(name,'jpeg',4)) then

        write(*,*) '   jpeg: write a JPEG file of screen'

      endif

      write(*,*) ' '

      end
