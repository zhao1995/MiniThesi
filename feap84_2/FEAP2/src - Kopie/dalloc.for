      subroutine deallocall
c-----------------------------------------------------------------------
c.... Purpose: Deallocation of all arrays
c-----------------------------------------------------------------------
      USE adap
      USE arcl
      USE aunr
      USE back1
      USE dirdat
      USE dyndat
      USE edgdat
      USE eig1
      USE epspu
      USE errin3
      USE ext1
      USE ext2
      USE fdata !??
      USE fe2dat
      USE fodata
      USE hdata
      USE hdatam
      USE hidden
      USE isbfgs
      USE isbfgs1
      USE iscsr
      USE isgmr
      USE isogeo
      USE ispgmr
      USE ispcg
      USE isprec
      USE jinteg
      USE maclg
      USE mate
      USE mdata
      USE mdat2
      USE mxasz
      USE ndata
      USE nolink
      USE pardi
      USE pcrit
      USE pdata3
      USE pdata7
      USE pindex
      USE plotdrp
      USE pnodn
      USE qload
      USE rsum
      USE slid1
      USE slid2
      USE slid3
      USE slid4
      USE slu
      USE smpak
      USE strnam
      USE subdt
      USE sumdt
      USE uneig
      USE working
      USE ximp

      USE doalloc
      implicit none
      logical tfl, ldummy

c.... adap 
      call rdealloc(adaerrot, tflb)

      call rdealloc(adanx, tflb)
      call rdealloc(adanf, tflb)
      call rdealloc(adanu, tflb)


c.... arcl 
      call rdealloc(arclm1     ,arcf)
      call rdealloc(arclm2     ,arcf)

c.... aunr 
      call rdealloc(adr       ,anr)

c.... back1
      call rdealloc(ustore    ,ldummy)

c.... dirdat
      call rdealloc(basea     ,tfl)

c.... dyndat
      call rdealloc(dynrea    ,fldyn)

c.... edgdat 
      call idealloc(edge1     ,tfl)
      call idealloc(edge2     ,tfl)
      call idealloc(edge3     ,tfl)
      call idealloc(edge4     ,tfl)

c.... eig1
      call rdealloc(eigk1     ,eigflg)
      call rdealloc(eigk2     ,eigflg)

c.... epspu
      call rdealloc(epspg    ,fepspg)

c.... errin3
      call rdealloc(e_ome    ,ldummy)

c.... ext1
      call rdealloc(extkh    ,extflg)
      call rdealloc(extkc    ,extflg)
      call rdealloc(extkd    ,extflg)
      call rdealloc(extke    ,extflg)
      call rdealloc(extkdh   ,extflg)
      call rdealloc(extkz1   ,extflg)
      call rdealloc(extkz2   ,extflg)

c.... fe2dat
      call rdealloc( mfe2a    ,ldummy)
      call rdealloc(mfe2g1    ,ldummy)
      call rdealloc(mfe2g2    ,ldummy)
      call rdealloc(mfe2g3    ,ldummy)
      call rdealloc(mfe2ix    ,ldummy)
      call rdealloc(mfe2ht    ,ldummy)
      call rdealloc(mfe2tau   ,ldummy)

c.... fodata
      call rdealloc(aifour1   ,fouflg)
      call rdealloc(aifour2   ,fouflg)

c.... hdata
      call rdealloc(gh1       ,ldummy)
      call rdealloc(gh2       ,ldummy)
      call rdealloc(gh3       ,ldummy)

c.... hdatam
      call rdealloc(eh1       ,ldummy)
      call rdealloc(eh2       ,ldummy)
      call rdealloc(eh3       ,ldummy)

c...  hidden    
      call rdealloc(dist      ,ldummy)
      call idealloc(idis      ,ldummy)
      call idealloc(mia       ,ldummy)

c.... isbfgs
      call rdealloc(bfgsbo    ,ldummy)
      call rdealloc(bfgsbd    ,ldummy)
      call rdealloc(bfgsbv    ,ldummy)
      call rdealloc(bfgsbw    ,ldummy)
      call rdealloc(bfgsbt    ,ldummy)

c.... isbfgs1
      call rdealloc(bfgsbs    ,ldummy)

c.... iscsr
      call idealloc(csrja     ,ldummy)
      call idealloc(csrka     ,ldummy)

c.... isgmr
      call rdealloc(rmgmrx    ,ldummy)
      call rdealloc(rmgmrss   ,ldummy)
      call rdealloc(rmgmrhh   ,ldummy)
      call rdealloc(rmgmrrs   ,ldummy)
      call rdealloc(rmgmrc    ,ldummy)
      call rdealloc(rmgmrs    ,ldummy)

c.... isogeo
      call idealloc(AInmpq    ,ldummy)
      call idealloc(AIninc    ,ldummy)
      call idealloc(AInien    ,ldummy)
      call idealloc(AInipa    ,ldummy)
      call idealloc(AInstre   ,ldummy)
      call rdealloc(AInkv1    ,ldummy)
      call rdealloc(AInkv2    ,ldummy)

c.... ispgmr
      call rdealloc(rmpgmrx   ,ldummy)
      call rdealloc(rmpgmrvv  ,ldummy)
      call rdealloc(rmpgmrrs  ,ldummy)
      call rdealloc(rmpgmrc   ,ldummy)
      call rdealloc(rmpgmrs   ,ldummy)

c.... ispcg
      call rdealloc(amcgz     ,ldummy)
      call rdealloc(amcgzz    ,ldummy)
      call rdealloc(amcgr     ,ldummy)
      call rdealloc(amcgrr    ,ldummy)
      call rdealloc(amcgp     ,ldummy)
      call rdealloc(amcgpp    ,ldummy)
      call rdealloc(amcgx     ,ldummy)

c.... isprec
      call rdealloc(rmpcalu   ,ldummy)
      call rdealloc(rmpcwl    ,ldummy)
      call rdealloc(rmpcwu    ,ldummy)

      call idealloc(impcjlu   ,ldummy)
      call idealloc(impcju    ,ldummy)
      call idealloc(impcjwl   ,ldummy)
      call idealloc(impcjwu   ,ldummy)
      call idealloc(impcjr    ,ldummy)
      call idealloc(impclevs  ,ldummy)

c.... jinteg
      call idealloc(ajint     ,jflgu)

c.... maclg
      call rdealloc(macl      ,tfl)

c.... mate
      call idealloc(matenew   ,flmat)
      call idealloc(mateorg   ,flmat)

c.... mdata
      call rdealloc(edis      ,tfl)
      call rdealloc(ecor      ,tfl)
      call rdealloc(etem      ,tfl)
      call rdealloc(epve      ,tfl)
      call rdealloc(ekma      ,tfl)
      call rdealloc(edma      ,tfl)
      call rdealloc(coor      ,tfl)
      call rdealloc(gloa      ,tfl)
      call rdealloc(glo0      ,tfl)
      call rdealloc(gu        ,tfl)
      call rdealloc(gtem      ,tfl)

      call idealloc(eeqn      ,tfl)
      call idealloc(nmat      ,tfl)
      call idealloc(psid      ,tfl)
      call idealloc(econ      ,tfl)
      call idealloc(jdt12     ,tfl)

c.... mdat2
      call rdealloc(aang      ,tfl)
      call rdealloc(bang      ,tfl)

c.... mxasz
      call idealloc(optin     ,tfl)

c.... ndata
      call rdealloc(gstiff    ,fl(3))
      fl(4) = .true.
      call rdealloc(dampm     ,tfl)
      call rdealloc(massm     ,fl(5))
      fl(6) = .true.
      call rdealloc(trans     ,tfl)

c.... nolink
      call idealloc(link1     ,tfl)
      call idealloc(link2     ,tfl)
      call idealloc(link3     ,tfl)

c.... pardi
      call rdealloc(drpar     ,lrhp)

c.... pcrit
      call rdealloc(crit     ,ldummy) !cww?? clfl?

c.... pdata7
      call idealloc(aipma     ,tfl)

c.... pindex
      call idealloc(apost     ,lindex)

c.... plotdrp
      call rdealloc(drp       ,fdrp)

c.... pnodn
      call idealloc(gtie      ,tfl)
      call idealloc(tecon     ,tfl)
      call rdealloc(parvvel   ,flparv)
      call rdealloc(parvacce  ,flparv)

c.... qload
      call rdealloc(aqloa     ,ldummy)

c...  rsum    
      call idealloc(irpt      ,ldummy)

c.... slid1
      call idealloc(cl00      ,ldummy)
      call idealloc(cl01      ,ldummy)
      call idealloc(cl02      ,ldummy)
      call idealloc(cl04      ,ldummy)
      call idealloc(cl05      ,ldummy)
      call idealloc(cl06      ,ldummy)
      call idealloc(cl07      ,ldummy)

      call rdealloc(cl03      ,ldummy)
      call rdealloc(cl08      ,ldummy)
      call rdealloc(cl09      ,ldummy)
      call rdealloc(cl10      ,ldummy)
      call rdealloc(cl11      ,ldummy)
      call rdealloc(cl12      ,ldummy)

c.... slid2
      call idealloc(cl22      ,ldummy)
      call idealloc(cl23      ,ldummy)

      call rdealloc(cl13      ,ldummy)
      call rdealloc(cl14      ,ldummy)
      call rdealloc(cl15      ,ldummy)
      call rdealloc(cl16      ,ldummy)
      call rdealloc(cl17      ,ldummy)
      call rdealloc(cl18      ,ldummy)
      call rdealloc(cl19      ,ldummy)
      call rdealloc(cl20      ,ldummy)
      call rdealloc(cl21      ,ldummy)
      call rdealloc(cl24      ,ldummy)
      call rdealloc(cl25      ,ldummy)

c.... slid3
      call idealloc(cl26      ,ldummy)
      call idealloc(cl27      ,ldummy)
      call idealloc(cl28      ,ldummy)
      call idealloc(cl29      ,ldummy)
      call idealloc(cl31      ,ldummy)

      call rdealloc(cl30      ,ldummy)

c.... slid4
      call rdealloc(cl33      ,ldummy)
      call rdealloc(cl34      ,ldummy)
      call rdealloc(cl35      ,ldummy)

c.... slu
      call rdealloc(drslu     ,ldummy)
      call rdealloc(diagslu   ,ldummy)

c.... smpak
      call rdealloc(smsidr1     ,ldummy)
      call rdealloc(smsirsp     ,ldummy)

      call idealloc(smsperm     ,ldummy)
      call idealloc(smsperm1    ,ldummy)
      call idealloc(smsisp      ,ldummy)

c.... strnam
      call rdealloc(strea     ,plfl)

c.... subdt
      call rdealloc(eigv      ,ldummy)
      call rdealloc(eigd      ,ldummy)
      call rdealloc(eigia     ,ldummy)

c.... sumdt
      call rdealloc(summ      ,flsum)

c.... uneig
      call rdealloc(aeigr     ,ldummy)
      call rdealloc(aeigi     ,ldummy)
      call rdealloc(aeigv     ,ldummy)
      call rdealloc(eigmma     ,ldummy)

c.... working
      call rdealloc(plo       ,tfl)
      call rdealloc(dr        ,tfl)

c.... ximp
      call rdealloc(uimp      ,flimp)

      return
      end
