module src_global

Implicit none
Include 'variablen.incf'

Integer, parameter :: mpnt = 32
Real*8 :: u(mpnt),v(mpnt),wt(mpnt),xkr(mpnt),dkr(mpnt)
Complex*16 ::pjzhx(mrank1,mrank1,mpnt),&
             pjzjx(mrank1,mrank1,mpnt),&
             rbjz(mrank1,mpnt),&
             rbhx(mrank1,mpnt), &
             refrc, eps, tmt(m2rank,m2rank,mrank1), pvz, phy
Real*8 :: rbjx(mrank1,mpnt), wvnm
Integer :: kxysym
REAL :: scshp, deq, dmx, wcnxl, wcnzc, wcnlm, &
          temp, scmix, refre, refim
REAL :: caupp(2), &
    cdtyp(2),&
    caavr(2),&
    cadev(2),&
        calow(2)
INTEGER :: canum(2), ndr
REAL ::  extm(4,4,mdr),&
     bmum(4,4,mdr), &
         stks(14,mdr)
REAL :: cdprb(mca2,mdr),&        !Canting parameters
    cdpth(mca2,mdr),&
        cdpab(mca2,mdr),&
    cdpdd(mca2,mdr),&
    cdpab2(mca2,mdr),&
        cdpdd2(mca2,mdr),&
    cdpabd(mca2,mdr)
Integer ::  cdpnum(mdr), n2rank, n2mode
Complex*16 :: tmat(m2rank,m2rank)
!Complex*16 :: aaa(m2rank,m2rank),&
!              b(m2rank,m2rank)
Save pjzhx, pjzjx, rbjz, rbhx, rbjx, wvnm, refrc, eps, kxysym, tmat,  &
  u, v, wt, xkr, dkr, scshp, deq, dmx, wcnxl, wcnzc, wcnlm, refre, refim, &
       caupp, cdtyp, caavr, cadev, calow, temp, scmix,  cdprb, cdpth, &
        cdpab,&
    cdpdd,&
    cdpab2,&
        cdpdd2,&
    cdpabd,ndr,  &
        cdpnum, extm, bmum, stks,  tmt, n2rank, n2mode, pvz, phy

End Module src_global
