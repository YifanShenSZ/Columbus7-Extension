cmcscf2.f
cmcscf part=2 of 10.  workspace allocation routine
cversion=5.6 last modified: 16-jan.2001
c
c
c***********************************************************************
c
c   the following computer programs contain work performed
c   partially or completely by the argonne national laboratory
c   theoretical chemistry group under the auspices of the office
c   of basic energy sciences, division of chemical sciences,
c   u.s. department of energy, under contract w-31-109-eng-38.
c
c   these programs may not be (re)distributed without the
c   written consent of the argonne theoretical chemistry group.
c
c   since these programs are under development, correct results
c   are not guaranteed.
c
c***********************************************************************
c  things to do in this program:
c  replace orthos.
c  replace orthot.
c  use 1-index transformations for large-scale 2-nd order convergence.
c  optionally eliminate mcuft program by inlining the walk generation.
c  replace flags(*) with integer array, allow keywords in input.
c  replace /csym/ with long version to eliminate equivalences.
c  indent do-loops and if blocks.
c
      subroutine driver( core, lcore, mem1, ifirst )
c
c  this program performs mcscf wave function optimization.
c
c  written by: ron shepard
c              theoretical chemistry group
c              chemistry division
c              argonne national laboratory
c              argonne, il 60439
c              internet: shepard@tcg.anl.gov
c              bitnet: shepard at anlchm
c
c  version log:
c  20-sep-96  state averaging extended for different symetries (-md)
c  17-jul-96  for flags(22) dynamicaly split the integral
c             transformation (-md)
c  17-jul-96  for flags 22 or 24 no sif2anl step necessary (-md)
c  11-july-96  ANL modifications (sif routines for density matrices
c              and phase factors for resolved orbitals)
c              according to:
c  14-jul-93  sifs i/o routines installed in wmcd1f() and wmcd2f()
c             -gfg, -hl
c  13-jun-96  direct computation of selected closed-shell contributions
c             to B*r added (-md)
c  13-jun-96  direct integral transformation added (-md)
c  13-jun-96  minimum restriction on number of act. orb. removed (-md)
c  13-jun-96  M*s in PSCI direct (-md)
c  13-jun-96  state averaging added (tk/md)
c  13-jun-96  optimized out-of-core integral transformation added (-tk/m
c  11-jul-94  label 'morbl' written to restart file. ahhc
c  04-may-92  eric stahlberg's MRPT changes merged. -rls
c  24-oct-91  doub(*) computation fixed in rddrt(). -rls/j.tilson
c  28-sep-91  hacked version to work with 4.1.1.0a2 COLUMBUS. -rls
c  08-sep-90  mcd*fl unit numbers initialized. -rls/tk
c  17-jan-90  wnorm equivalence bug fixed in readin. -rls/eas
c  06-dec-89  status='scratch' removed from all open statements. -rls
c  25-mar-89  sun, titan, alliant port. cmdc changes. -rls
c  07-jan-89  call wtstv1() from wmcd1f(). -rls
c  15-nov-88  change density file structure. -rls
c  01-sep-88  split wrmode(). -rls
c  11-jul-88  unicos version (rls).
c  04-feb-88  fps integral header change, nmotx in tran routines (rls).
c  15-jun-87  write density matrices on file for gradient
c  04-apr-86  iwrdim bug fixed in seg (rls).
c  04-apr-86  da2cr added to uexp (rls).
c  19-mar-86  fps master conversion (rls).
c  17-mar-86  vax/cray master version written (rls).
c  11-oct-85  flags(20) maximum overlap psci vector (rls).
c  23-jul-85  namelist input version (rls).
c  04-jun-85  cray version changes (rls).
c  25-apr-85  mctran files modified and added (rls).
c  05-apr-85  iterative hmc(*) diagonalization added (rls).
c  15-mar-85  inactive-orbital blocks added (rls).
c  19-oct-84  symmetry-dependent drt and ft added (rls).
c  19-apr-84  first alpha-test version 2.00 written by ron shepard
c             and developed in collaboration with i. shavitt from
c             jan-83 to apr-84.
c
       use tranlibrary, only: daclrd
       implicit none
c  ## parameter & common section
c
c michal2{
c
c Silmar 19/09/02 - common blocks and variables for cosmo
c
       integer maxorb,norb,maxdens,freemem
       real*8 Vsolv
       parameter(maxorb=1023)
       common/solvent/Vsolv(maxorb*(maxorb+1)/2),
     . norb,maxdens
       real*8 Vnenuc
       common/Vne3/Vnenuc
c
       integer cosmocalc
      common/cosmo2/cosmocalc
c
c Silmar - files, common block and variables  needed for cosmo-module
c
      integer maxat,natoms,maxnps,maxrep
      parameter(maxat=350, maxnps = 2000,maxrep=8)
      integer mmtype(maxat),typecalc
      real*8 xyz(3,maxat),xyzpoints(3),sumq1,sumq2,sumq
      real*8 phi(maxnps)
      real*8 qcos(maxnps),qcoszero(maxnps),symf
      real*8 ediel,elast,xfinal(maxrep),yfinal(maxrep),zfinal(maxrep)
      real*8 rtolcosmo
      real*8 phizero(maxnps),dcos(3,maxat)
      real*8 edielcorr
      character*26 atext,string
      REAL*8, allocatable :: a1mat(:),a2mat(:,:),a3mat(:)
      REAL*8, allocatable :: cosurf(:,:)

c
      character*3 symgrp
      integer nps,npspher
      integer dbglvl,isymf,np
c
      include 'coss.inc'
c michal2}
c
*@ifndef direct
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
*@endif
c
      character*8 program,version
      parameter (program='MCSCF')
      parameter (version='5.5')
c
      integer nbfmxp,ninfomx,nenrgymx,nmapmx
      parameter (nbfmxp=1023,ninfomx=6,nenrgymx=5,nmapmx=2)
c
      real*8    zero, one, two
      parameter(zero=0d0,one=1d0,two=2d0)
c
c    #  maximal No. of different DRTs
      integer maxnst
      parameter (maxnst=8)
c
c    #  maximal No. of states in one DRT
      integer mxavst
      parameter (mxavst = 50)
c
c    #  maximal total number of states
      integer maxstat
      parameter (maxstat = maxnst*mxavst)
c
      integer nbukmx
      parameter (nbukmx=255)
c
      integer navst,nst,navst_max
      real*8 heig,wavst
      common /avstat/ wavst(maxnst,mxavst), heig(maxnst,mxavst),
     & navst(maxnst),nst,navst_max
c
      integer iunits
      common/cfiles/iunits(55)
      integer nlist, nin, aoints, aoints2, mocoef, nslist,
     & nrestp, ntemp, scrtfl, mcd1fl, mcd2fl, mos, hdiagf,
     & info2, mrptfl,tmin,mocout
      equivalence (iunits(1),nlist)
      equivalence (iunits(2),nin)
      equivalence (iunits(3),aoints)
      equivalence (iunits(4),aoints2)
      equivalence (iunits(7),mocoef)
      equivalence (iunits(9),nslist)
      equivalence (iunits(12),nrestp)
      equivalence (iunits(13),ntemp)
      equivalence (iunits(14),scrtfl)
      equivalence (iunits(15),mcd1fl)
      equivalence (iunits(16),mcd2fl)
      equivalence (iunits(19),mos)
      equivalence (iunits(20),hdiagf)
      equivalence (iunits(21),info2)
      equivalence (iunits(22),mrptfl)
      equivalence (iunits(23),tmin)
      equivalence (iunits(24),mocout)
c
      integer nrow_f,ncsf_f,nwalk_f,ssym
      common/drtf/ nrow_f(maxnst),ncsf_f(maxnst),
     & nwalk_f(maxnst),ssym(maxnst)
c
c     cpt2(ist+1) = cpt2(ist) + forbyt(nsym*nrow_f(ist))
c     cpt3(ist+1) = cpt3(ist) + forbyt(nwalk_f(ist))
c     cpt4(ist+1) = cpt4(ist) + atebyt(ncsf_f(ist)*navst(ist))
c     cpt2_tot = cpt2_tot + nsym*nrow_f(ist)
c     cpt3_tot = cpt3_tot + nwalk_f(ist)
c     cpt4_tot = cpt4_tot + ncsf_f(ist)*navst_(ist)
c
      integer cpt2,cpt2_tot,cpt3,cpt3_tot,cpt4,cpt4_tot
      integer ncsf_max
      common/counter/cpt2(maxnst+1),cpt3(maxnst+1),cpt4(maxnst+1),
     & ncsf_max,cpt2_tot,cpt3_tot,cpt4_tot
c
      integer msize,csize
      common/size/ msize(maxnst),csize(4,maxnst)
c michal2{
c
c Silmar - 23/09/02 - The nbpsy array is transfered to a
c common block in case of cosmo calculation. With this procedure
c the adresses of WORK array in driverDALTON are calculated using this
c information
c
      common/cosmoblock/nbfpsy
c
c michal2}
c
      integer nxy, mult, nsym, nsymm, nxtot
      common/csymb/nxy(8,42),mult(8,8),nsym,nxtot(24)
      integer nbpsy(8),nbfpsy(8),nsb(8),nnsb(8),n2sb(8),nbft,nntb, n2tb
      equivalence (nxy(1,1),nbpsy(1))
      equivalence (nxy(1,2),nsb(1))
      equivalence (nxy(1,3),nnsb(1))
      equivalence (nxy(1,4),n2sb(1))
      equivalence (nxtot(1),nbft)
      equivalence (nxtot(2),nntb)
      equivalence (nxtot(3),n2tb)
      integer nmpsy(8),nsm(8),nnsm(8),n2sm(8),nmot,nntm,n2tm
      equivalence (nxy(1,5),nmpsy(1))
      equivalence (nxy(1,6),nsm(1))
      equivalence (nxy(1,7),nnsm(1))
      equivalence (nxy(1,8),n2sm(1))
      equivalence (nxtot(4),nmot)
      equivalence (nxtot(5),nntm)
      equivalence (nxtot(6),n2tm)
      integer ndpsy(8),nsd(8),nnsd(8),n2sd(8),ndot,nntd,n2td
      equivalence (nxy(1,9),ndpsy(1))
      equivalence (nxy(1,10),nsd(1))
      equivalence (nxy(1,11),nnsd(1))
      equivalence (nxy(1,12),n2sd(1))
      equivalence (nxtot(7),ndot)
      equivalence (nxtot(8),nntd)
      equivalence (nxtot(9),n2td)
      integer napsy(8),nsa(8),nnsa(8),n2sa(8),nact,nnta,n2ta
      equivalence (nxy(1,13),napsy(1))
      equivalence (nxy(1,14),nsa(1))
      equivalence (nxy(1,15),nnsa(1))
      equivalence (nxy(1,16),n2sa(1))
      equivalence (nxtot(10),nact)
      equivalence (nxtot(11),nnta)
      equivalence (nxtot(12),n2ta)
      integer nvpsy(8),nsv(8),nnsv(8),n2sv(8),nvrt,nntv,n2tv
      equivalence (nxy(1,17),nvpsy(1))
      equivalence (nxy(1,18),nsv(1))
      equivalence (nxy(1,19),nnsv(1))
      equivalence (nxy(1,20),n2sv(1))
      equivalence (nxtot(13),nvrt)
      equivalence (nxtot(14),nntv)
      equivalence (nxtot(15),n2tv)
      integer tsymdd(8),tsmdd2(8),tsymad(8),tsymaa(8)
      integer tsymvd(8),tsymva(8),tsymvv(8),tsmvv2(8)
      equivalence (nxy(1,21),tsymdd(1))
      equivalence (nxy(1,22),tsmdd2(1))
      equivalence (nxy(1,23),tsymad(1))
      equivalence (nxy(1,24),tsymaa(1))
      equivalence (nxy(1,25),tsymvd(1))
      equivalence (nxy(1,26),tsymva(1))
      equivalence (nxy(1,27),tsymvv(1))
      equivalence (nxy(1,28),tsmvv2(1))
      integer nadt,nvdt,nvat
      equivalence (tsymad(1),nadt)
      equivalence (tsymvd(1),nvdt)
      equivalence (tsymva(1),nvat)
      integer nsad(8),nsvd(8),nsva(8)
      equivalence (nxy(1,29),nsad(1))
      equivalence (nxy(1,30),nsvd(1))
      equivalence (nxy(1,31),nsva(1))
      integer nsbm(8),nbmt,ndimd,ndima,ndimv
      equivalence (nxy(1,32),nsbm(1))
      equivalence (nxtot(16),nbmt)
      equivalence (nxtot(17),ndimd)
      equivalence (nxtot(18),ndima)
      equivalence (nxtot(19),ndimv)
      integer irb1(8),irb2(8)
      equivalence (nxy(1,33),irb1(1))
      equivalence (nxy(1,34),irb2(1))
      integer irm1(8),irm2(8)
      equivalence (nxy(1,35),irm1(1))
      equivalence (nxy(1,36),irm2(1))
      integer ird1(8),ird2(8)
      equivalence (nxy(1,37),ird1(1))
      equivalence (nxy(1,38),ird2(1))
      integer ira1(8),ira2(8)
      equivalence (nxy(1,39),ira1(1))
      equivalence (nxy(1,40),ira2(1))
      integer irv1(8),irv2(8)
      equivalence (nxy(1,41),irv1(1))
      equivalence (nxy(1,42),irv2(1))
c
        integer nmotx
       parameter (nmotx=1023)
c
      integer iorbx, ix0
      common/corbc/iorbx(nmotx,8),ix0(3)
      integer nndx(nmotx),symb(nmotx),symm(nmotx),symx(nmotx)
      integer invx(nmotx),orbidx(nmotx),orbtyp(nmotx),fcimsk(nmotx)
      equivalence (iorbx(1,1),nndx(1))
      equivalence (iorbx(1,2),symb(1))
      equivalence (iorbx(1,3),symm(1))
      equivalence (iorbx(1,4),symx(1))
      equivalence (iorbx(1,5),invx(1))
      equivalence (iorbx(1,6),orbidx(1))
      equivalence (iorbx(1,7),orbtyp(1))
      equivalence (iorbx(1,8),fcimsk(1))
c
      integer addpt, numint, bukpt, intoff, szh
      common/caddb/addpt(24),numint(12),bukpt(12),intoff(12),szh(15)
      integer nbuk
      equivalence (bukpt(12),nbuk)
c
      integer hbci
      common/chbci/hbci(4,3)
c
      integer hblkw, nhbw, htape
      common/chbw/hblkw(15*maxstat),nhbw,htape
c
      integer stape, lenbfs, h2bpt
      common/cbufs/stape,lenbfs,h2bpt
c
*@ifdef direct
*      include 'param.h'
*      integer ifock,ipq
*      common /fockij/ ifock(ndi4+1)
*      common /pairij / ipq(ndi4+1)
*@endif
c
      integer perm1, ind1, inds1, npass1
      common/ccone/perm1(2,2),ind1(2),inds1,npass1
c
      integer perm, ind, inds, npassp
      common/cctwo/perm(4,4),ind(4),inds(4),npassp
c
      integer ldamin, ldamax, ldainc
      common/clda/ldamin,ldamax,ldainc
c
      integer lenm2e,n2embf
      integer  lena1e, n1eabf,lenbft,
     &        lend2e,n2edbf,lend1e,n1edbf,ifmt1,ifmt2
      common/cfinfo/lena1e,n1eabf,lenbft,
     &  lend2e,n2edbf,lend1e,n1edbf,ifmt1,ifmt2
c
c  lena1e=length of buffers for basis function 1-e integrals.
c  n1eabf=number of basis function 1-e integrals in each buffer.
c  lenbft=length of the formula tape buffers.
c  lenbfs=length of sorted integral buffers.
c
      integer niter,nstate,ncol,ncoupl,noldv,noldhv,nunitv, nciitr,
     &        mxvadd,nvcimx,nmiter,nvrsmx
      common /inputc/niter,nstate(8),ncol(8),ncoupl,
     &  noldv(8),noldhv(8),nunitv,nciitr,mxvadd,nvcimx,nmiter,nvrsmx
c
      integer nflag
      parameter (nflag=50)
      logical flags
      common/lflags/flags(nflag)
c
      real*8 tol,rtolci
      common/toler/tol(12),rtolci(mxavst)
c
      integer lvlprt
      common/prtout/lvlprt
c
      integer lcinf
      parameter(lcinf=6)
c
      real*8 coninf
      common/cconi/coninf(lcinf)
      real*8 emc,demc,wnorm,knorm,apxde,repnuc
      equivalence (coninf(1),emc)
      equivalence (coninf(2),repnuc)
      equivalence (coninf(3),demc)
      equivalence (coninf(4),wnorm)
      equivalence (coninf(5),knorm)
      equivalence (coninf(6),apxde)
c
      integer occl
      common /ocupy/ occl(nmotx)
c
c
      integer nunit2,len2,irec2
      common/da2/nunit2,len2,irec2
c
      real*8 sifsce
      external sifsce
      character*4 slabelao(nbfmxp)
      character*8 bfnlabao(nbfmxp)
      character*8 bfnlabmo(nbfmxp)
      integer info,nmap,ietype,imtype,map,nenrgy
      real*8 energy
      common /forsif/ info(ninfomx),nenrgy,nmap,
     . ietype(nenrgymx), energy(nenrgymx),
     . imtype(nmapmx), map(nbfmxp*nmapmx),
     . slabelao, bfnlabao
c
      integer eventao
      common/timerao/ eventao
c
c  ##  integer section
c
      integer dlt2,lumorb,sewd
      integer avcors, avc2is, avchbc, apt, avchd, avcsk, aoint2
      integer cpt(50)
      integer dpt,dimsao,dimh1
      integer event1, event2, event3
      integer i,ierr,isym,ist, iter, itert, idim, ii, ij, impt, ivpt,
     &        ifirst,idim1,idim2,ivcmax,icount
      external ivcmax
      integer j
      integer kcoret
      integer lenbuf_f(maxnst),lcore,lcoremain, lenbuk,lenbuf
      integer mapmm(nmotx), mxcsao, mem1
      integer ntitle,ninfo, ndoub, nelt, nela, nlevel, nnact, naar,
     &        nblktr, nipbuk, ndimw, nrowmax, numpbf,nbf_ca, nd_at,
     &        nipsy_mx,nipsy(8)
      integer req2, reqhbc
      integer smult
      integer vpt,infoloc(6)
      integer isymunit,itmp(4),typconv(0:6),itm
      integer nbt,idummy,luonel,itmp2(8),itmpcnt(2*nmotx)
cfp
      real*8 au2ev
      parameter (au2ev = 2.72113961d01)
      integer orbmxpri
      parameter (orbmxpri = 16)
      character*40 cfmt
      parameter (cfmt='(1p3d25.15)')
      integer nxx,adpt1,sdpt1,sdpt2,cnospt,occspt,jstrt,jlast
      integer filerr, syserr
      real*8  deltaab, enev
      integer dst,ndens,densind,ndbra,brai,brinv
      common/dens/dst,ndens,densind(2,mxavst),ndbra,brai(mxavst),
     & brinv(mxavst)
c
c  ## real*8 section
c
      real*8 core(lcore)
      real*8 ecore0,emcold
      real*8 hcore
      real*8 emin,score
      real*8 rtitle(24) 
c
c  ##  logical section
c
      logical lopen,op
      logical qsort,qcoupl,qconvg
      integer getfreeunit
c
c  ## character section
c
      character*16 convrg
      character*60 fname
      character*4 slabel(20)
      character*80 title(50)
      character*3 mtype(nmotx)
      character*3 tmplabel(8)
      integer nvalid,irc
      integer mcirc,mcnsym,mcnbpsy(8)
     . ,mctmp,mcmap(nbfmxp),mcone,mcz,mcopt
c
c  ##  external section
c
      integer  atebyt,  forbyt
      external atebyt,  forbyt
c
c  ##  data section
      data typconv /1,2,4,6,9,12,16/
c
c
c  this has to be done in a more automatic and flexible
c  manner ( think for example of a change of the default record
c  length) the best way to do is probably to use the sifcfg routines
c  to calculate the appropriate numbers - then there is only ONE
c  parameter to change in order to make a proper adjustment
c
c  n2ebf_1=number of mo 2-e integrals in each buffer.(if nbft.le.255).
c  n2ebf_2=number of mo 2-e integrals in each buffer.(if nbft.gt.255).
      integer  n2ebf_1,  n2ebf_2
      data  n2ebf_1/2730/
      data  n2ebf_2/2047/
c
c
c--------------------------------------------------------------
c
cmd
*@ifdef crayctss
*CC  create the dropfile.
*      call link('unit5=tty,unit6=unit5//')
*@endif
c
      fname = 'mcscfls'
      call trnfln( 1, fname )
      open(unit=nlist,file=fname,status='unknown')
c
      call ibummr(nlist)
      call timer(' ',0,event1,nlist)
c
      fname = 'mcscfin'
      call trnfln( 1, fname )
      open(unit=nin,file=fname,status='old')
c
      fname = 'mcscfsm'
      call trnfln( 1, fname )
      open(unit=nslist,file=fname,status='unknown')
c

      call headwr(nlist,program,version)
6010  write(nlist,7000)
7000  format(' This program allows the',
     &  ' csf mixing coefficient and orbital expansion coefficient'
     &  /' optimization using the graphical unitary group approach',
     &  ' and the exponential'/,' operator mcscf method.')
      write(nlist,7010)
7010  format
     & (' references:',t15,'r. shepard and j. simons,
     & '' int. j. quantum chem. symp. 14, 211 (1980).'
     & /t15,'r. shepard, i. shavitt, and j. simons,',
     & ' j. chem. phys. 76, 543 (1982).')
      write(nlist,7020)
7020  format
     & (t15,'r. shepard in "ab initio methods in quantum chemistry',
     & ' ii" advances in chemical'/t19,
     & 'physics 69, edited by k. p. lawley (wiley, new york, ',
     & '1987) pp. 63-200.')
      write(nlist,7030)
7030  format
     &  (' Original autor: Ron Shepard, ANL'
     &  /,' Later revisions: Michal Dallos, University Vienna')
c
      call who2c( 'MCSCF', nlist )
      call headwr(nlist,'MCSCF','5.4.0.2 ')
c
      write(nlist,7040)lcore,
     & real(lcore)/(128.d+00*1024.0d+00)
7040  format
     & (' Workspace allocation information:'/,
     & 3x,i13,1x,'of real*8 words (',f8.2,
     & ' MB) of work space has been allocated.')
c
c  initialize local variables.
c
      repnuc=zero
      ecore0=zero
      emc=zero
      demc=zero
      emcold=zero
      knorm=zero
      wnorm=zero
      apxde=zero
      call izero_wr(8,nbpsy,1)
c
c  read header information from the integral file:
c   1:buffer
c
      cpt(1)=1
c
cmd
c     first read only the flags
      call readin(.true.)

c  molcas stuff

      if (flags(26)) then
        call link_it()
        call getenvinit()
        call fioinit()
        call NameRun('RUNFILE')
      endif 

c
      write(nlist,
     &'(//,1x,3(1h*),''  Integral file informations  '',3(1h*)//)')
c
      if (flags(22).or.flags(24))then
c     write(nlist,*)'reading aoints'
c     fname = 'aoints'
c     call trnfln( 1, fname )
c     open(unit=aoints,file=fname,status='old',form='unformatted')
c     call rdhead(aoints,ntitle,title,repnuc,nbft,nsym,nbpsy,ierr)

c  get rid of this file by first invoking kora_md (statistics run)
c  at the beginning of the driver routine and afterwards just
c  accessing the appropriate arrays, i.e.
c     nsym and nlamda can be found in common block /irreps/
c      and correspond to nsym and nbpsy
c     repnuc can be found in common block /rep/ or if
c     this common block is not yet initialized can be calculated
c     via call nucrep(xyz,charge,repnuc,natoms)
c      with xyz,charge,natoms in common /infoa/
c

      write(nlist,*)'reading repnuc'
      fname = 'info2'
      call trnfln( 1, fname )
      open(unit=info2,file=fname,status='old',form='formatted')
      read(info2,*)nsym
      if (nsym.gt.8) call bummer('uncorect No. of irreps =',nsym,faterr)
      read(info2,*)(nbpsy(i),i=1,nsym)
      read(info2,*)repnuc
      close(unit=info2)
      nbft = 0
      do isym = 1,nsym
         nbft = nbft + nbpsy(isym)
      enddo
cmd   write(nlist,6120)(title(i),i=1,ntitle)
      write(nlist,6130)repnuc,nbft,nsym,(i,nbpsy(i),i=1,nsym)
cmd   write(nslist,6120)(title(i),i=1,ntitle)
      write(nslist,6130)repnuc,nbft,nsym,(i,nbpsy(i),i=1,nsym)
6120  format(/' integral file header information:'/(1x,a))
6130  format(' repulsion energy:',1pe25.15,' hartree.'/
     &  ' total number of basis functions:',i4/
     &  ' total number of irreps:',i2/' (isym:nbpsy)',8(i3,':',i3))
c  get rid of this stuff
cmd
c  lenm2e=length of buffers for mo 2-e integrals.
c  n2embf=number of mo 2-e integrals in each buffer.

c     if(nbft.le.255) then
c        lenm2e = le2e_1
c        n2embf=n2ebf_1
c        write(nlist,*)'No of nbft is less than 255'
c        write(nlist,*)'lenm2e,n2embf',lenm2e,n2embf
c     else
c        lenm2e = le2e_2
c        n2embf=n2ebf_2
c        write(nlist,*)'No of nbft is greater than 255'
c        write(nlist,*)'lenm2e,n2embf',lenm2e,n2embf
c     endif
c use the following commands
c      input:
c      lrecal: physical record length
c      nbft : number of basis functions
c      returns reclength for 1e-integrals (l1rec)
c      returns  number of integrals per buffer (n1max)
c     call sifcfg(1,4096,nbft,0,ifmt1,l1rec,n1max,ierr)
c     the same for the 2e-stuff
      call sifcfg(2,4096,nbft,0,ifmt2,lenm2e,n2embf,ierr)

cmd
      else
c  molcas
      if (flags(26)) then
         ntitle=1
         title(1)='SEWARD INTEGRALS'
         call get_dscalar('PotNuc',repnuc)
         nenrgy=1
         energy(nenrgy)=repnuc
         ietype(nenrgy)=0
         call mcget_iscalar('nsym',mcnsym)
         nsym=mcnsym
         call mcget_iarray('nbas',mcnbpsy,mcnsym)
         mctmp=24 
         call mcget_carray('irreps',tmplabel,mctmp)
     
         nbft=0
         do i=1,nsym
            slabelao(i)(2:4)=tmplabel(i)
            slabelao(i)(1:1)=' '
            nbpsy(i)=mcnbpsy(i)
            nbft=nbft+nbpsy(i)
         enddo
         nmap=2
         if (nbft.gt.nbfmxp) then
           call bummer('mcscf2.f: seward nbft=',nbft,2)
         elseif (nmap.gt.nmapmx) then
           call bummer ('mcscf2.f: seward nmap=',nmap,2)
         endif
         mctmp=8*nbft
         call mcget_carray('unique basis names',bfnlabao,mctmp)
c
c       we cannot obtain currently proper bfn-to-center and
c       orbital type map from the RUNFILE infos
c       instead use the $Project.SymInfo file (molcas.SymInfo)
c
         isymunit=-1
         do i=80,90
         inquire(unit=i,opened=lopen)
         if (.not.lopen) then
           isymunit=i
           exit
         endif
         enddo
         if (isymunit.eq.-1) then
          call bummer('could not find valid unitno (isymunit)',0,2)
         endif
         open(unit=isymunit,file='molcas.SymInfo')
         read(isymunit,*)
         read(isymunit,*)
         imtype(1)=3
         imtype(2)=4
         do i=1,nbft
          read(isymunit,*,err=91,end=92) (itmp(j),j=1,4)
          map(i)=itmp(2)
          map(i+nbft) = typconv(itmp(3)) 
c #of funct, unique centre, L, M , # of sym.ad.functions , Phases
         enddo
         goto 93
 91      call bummer('molcas.SymInfo: error encountered',0,2)
 92      call bummer('molcas.SymInfo: eof encountered',0,2)
 93      continue

         write(nlist,'(/1x,''Integral file header information:'')')
         write(nlist,6110) (title(i),i=1,ntitle)
         repnuc = sifsce(nenrgy,energy,ietype)
         write(nlist,6100) 'total ao core energy =',repnuc
         write(nslist,6100) 'total ao core energy =',repnuc
         write(nlist,
     & '(//,3x,6(1h*),''  Basis set information:  '',6(1h*)/)')
         write(nlist,'('' Number of irreps:'',t37,i1)')nsym
         write(nlist,'('' Total number of basis functions:'',t35,i3)')
     & nbft
        write(nlist,'(/,1x,''irrep no.'',t23,8(i3,2x))')(i,i=1,nsym)
        write(nlist,'(1x,''irrep label'',t23,8(a4,1x))')
     & (slabelao(i),i=1,nsym)
        write(nlist,'(1x,''no. of bas.fcions.'',t23,8(i3,2x))')
     & (nbpsy(i),i=1,nsym)
         if (lvlprt.ge.1) then
         write(nlist,6110) 'input basis function labels, i:bfnlab(i)='
         write(nlist,6080) (i,bfnlabao(i),i=1,nbft)
         endif
        write(nlist,*) 'map-vector 1 , imtype=',imtype(1)
        write(nlist,'(20(i3,1x))') (map(i),i=1,nbft)
        write(nlist,*) 'map-vector 2 , imtype=',imtype(2)
        write(nlist,'(20(i3,1x))') (map(i+nbft),i=1,nbft)

      elseif (flags(27)) then 
c       DALTON2 integrals 
         ntitle=1
         title(1)='DALTON2 INTEGRALS'
c
c        take symmetry info from infofl
c
         luonel=getfreeunit()
          open(unit=luonel,file='infofl',form='formatted')
            read(luonel,*)
            read(luonel,*)
            read(luonel,*) nsym 
            read(luonel,97) (slabelao(i)(1:4),i=1,nsym)
 97         format(8(1x,a4))
            rewind(luonel)
            read(luonel,*)
            read(luonel,*) (nbpsy(i),i=1,nsym)
         close(luonel)
        CALL GPOPEN(LUONEL,'AOONEINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
        READ(LUONEL) RTITLE,itmp(1),(itmp2(I),I = 1,nsym),repnuc
        if (lvlprt.ge.0) then
        write(6,185) repnuc,itmp(1),(itmp2(i),i=1,nsym)
 185    format('DALTON2 INFO: repnuc=',f12.6,' nsym=',i4,
     .        ' nbpsy(*)=',8i4)
        endif
        nbft=0
        do i=1,nsym
          if (nbpsy(i).ne.itmp2(i))
     .     call bummer('inconsistent AOONEINT and infofl file',0,2)
           nbft=nbft+nbpsy(i)
        enddo
         nenrgy=1
         energy(nenrgy)=repnuc
         ietype(nenrgy)=0

         nmap=2
         if (nbft.gt.nbfmxp) then
           call bummer('mcscf2.f: dalton nbft=',nbft,2)
         elseif (nmap.gt.nmapmx) then
           call bummer ('mcscf2.f: dalton nmap=',nmap,2)
         endif
         CALL MOLLAB('SYMINPUT',LUONEL,6)
         READ (LUONEL) nbt,(itmpcnt(I),I=1,2*nbft)
         IF (NBT .NE. nbft) THEN
           call bummer('inconsistent data within SYMINPUT',0,2)
         END IF
         DO I=1,nbft
            write(bfnlabao(i),187) itmpcnt(i),itmpcnt(i+nbft)
            do j=1,8
               if (bfnlabao(i)(j:j).eq.' ') bfnlabao(i)(j:j)='_'
            enddo
  187       format(a4,a4)
         enddo
        CALL GPCLOSE(LUONEL,'KEEP')

         isymunit=getfreeunit()
         open(unit=isymunit,file='SYMINFO')
         read(isymunit,*)
         read(isymunit,*)
         imtype(1)=3
         imtype(2)=4
         do i=1,nbft
          read(isymunit,*,err=191,end=192) (itmp(j),j=1,4)
          map(i)=itmp(2)
          map(i+nbft) = typconv(itmp(3))
c #of funct, unique centre, L, M , # of sym.ad.functions , Phases
         enddo
         goto 193
 191      call bummer('SYMINFO: error encountered',0,2)
 192      call bummer('SYMINFO: eof encountered',0,2)
 193      continue

         write(nlist,'(/1x,''Integral file header information:'')')
         write(nlist,6110) (title(i),i=1,ntitle)
         repnuc = sifsce(nenrgy,energy,ietype)
         write(nlist,6100) 'total ao core energy =',repnuc
         write(nslist,6100) 'total ao core energy =',repnuc
         write(nlist,
     & '(//,3x,6(1h*),''  Basis set information:  '',6(1h*)/)')
         write(nlist,'('' Number of irreps:'',t37,i1)')nsym
         write(nlist,'('' Total number of basis functions:'',t35,i3)')
     & nbft
        write(nlist,'(/,1x,''irrep no.'',t23,8(i3,2x))')(i,i=1,nsym)
        write(nlist,'(1x,''irrep label'',t23,8(a4,1x))')
     & (slabelao(i),i=1,nsym)
        write(nlist,'(1x,''no. of bas.fcions.'',t23,8(i3,2x))')
     & (nbpsy(i),i=1,nsym)
         if (lvlprt.ge.1) then
         write(nlist,6110) 'input basis function labels, i:bfnlab(i)='
         write(nlist,6080) (i,bfnlabao(i),i=1,nbft)
         endif
        write(nlist,*) 'map-vector 1 , imtype=',imtype(1)
        write(nlist,'(20(i3,1x))') (map(i),i=1,nbft)
        write(nlist,*) 'map-vector 2 , imtype=',imtype(2)
        write(nlist,'(20(i3,1x))') (map(i+nbft),i=1,nbft)

      else 

        fname='aoints'
        call trnfln (1,fname)
        open(unit=aoints,file=fname,status='old',form='unformatted')
        inquire (unit=aoints,name=fname)
        write(nlist,*) 'input integral file : ' // fname
c       read header info

        call sifrh1(aoints,ntitle,nsym,nbft,
     .  ninfo,nenrgy,nmap,ierr)
        if (ierr.ne.0) then
        call bummer('mcscf2.f: from sifrh1, ierr=',ierr,2)
        elseif (ntitle.gt.50) then
        call bummer('mcscf2.f: ntitle=',ntitle,2)
        elseif (nbft.gt.nbfmxp) then
        call bummer('mcscf2.f: nbft=',nbft,2)
        elseif (nenrgy.gt.nenrgymx) then
       call bummer ('mcscf2.f: nenrgy=',nenrgy,2)
        elseif (nmap.gt.nmapmx) then
        call bummer ('mcscf2.f: nmap=',nmap,2)
        endif
c    read additional header info
        call sifrh2(
     .  aoints,ntitle,nsym,nbft,ninfo,nenrgy,nmap,title,
     .  nbpsy,slabelao,info,bfnlabao,ietype,energy,imtype,
     .  map,ierr)

        if (ierr.ne.0)
     .  call bummer('mcscf2.f: from sifrh2, ierr=',ierr,2)
c       if ((info(2).ne.4096).or.(info(4).ne.4096)) then
c       write(nlist,*) 'AO integral program did not write ',
c    .  'integral files with record length 4096'
c       call bummer('AO integral record length <> 4096',
c    .  info(2),2)
c       endif

        write(nlist,'(/1x,''Integral file header information:'')')
        write(nlist,6110) (title(i),i=1,ntitle)
        write(nlist,'(/'' Core type energy values:'')')
        call sifpre(nlist,nenrgy,energy,ietype)
        repnuc = sifsce(nenrgy,energy,ietype)
        write(nlist,6100) 'total ao core energy =',repnuc
        write(nslist,6100) 'total ao core energy =',repnuc
        write(nlist,
     &  '(//,3x,6(1h*),''  Basis set information:  '',6(1h*)/)')
        write(nlist,'('' Number of irreps:'',t37,i1)')nsym
        write(nlist,'('' Total number of basis functions:'',t35,i3)')
     &  nbft
        write(nlist,'(/,1x,''irrep no.'',t23,8(i3,2x))')(i,i=1,nsym)
        write(nlist,'(1x,''irrep label'',t23,8(a4,1x))')
     & (slabelao(i),i=1,nsym)
        write(nlist,'(1x,''no. of bas.fcions.'',t23,8(i3,2x))')
     & (nbpsy(i),i=1,nsym)
         if (lvlprt.ge.1) then
         write(nlist,6050) 'info(*) =', (info(i),i=1,ninfo)
         write(nlist,6110) 'input basis function labels, i:bfnlab(i)='
         write(nlist,6080) (i,bfnlabao(i),i=1,nbft)
         endif

       endif
c
6050  format(1x,a,8i5)
6080  format(10(i4,':',a8))
6100  format(1x,a,f15.9)
6110  format(1x,a)

      if (flags(26) .or. flags(27)) then
         call sifcfg(1,4096,nbft,0,ifmt1,info(2),info(3),ierr)
         call sifcfg(2,4096,nbft,0,info(6),info(4),info(5),ierr)
         write(nlist,*) 'using sifcfg values:',(info(i),i=2,6)
      else
c  generate ifmt parameters locally
         call sifcfg(1,4096,nbft,0,ifmt1,infoloc(2),infoloc(3),ierr)
         call sifcfg(2,4096,nbft,0,infoloc(6),infoloc(4),infoloc(5),
     .                  ierr)
c  check for consistency with info parameter from existing aoints file
c  argos uses ibvtyp=1, i.e.
c  call sifcfg(2,4096,nbft,1,ifmt2,infoloc(4),infoloc(5),ierr)
           do i=2,6
              if (info(i).ne.infoloc(i)) then
               write(nlist,*)'inconsistent sifs parameter'
               write(nlist,'(a,3i6)')' i,info(i),infoloc(i):'
     $ ,i,info(i),infoloc(i)
               call bummer('inconsistent sifs parameter ',i,0)
               endif
           enddo

        call sifo2f(aoints,aoints2,'aoints2',info,aoint2,ierr)
c       aoint2 is the unit number of the 2e ao integral file
        if (ierr.ne.0)
     .  call bummer('from sifo2f: ierr=',ierr,2)
      endif


      lenm2e=info(4)
      n2embf=info(5)
c     use the same record lengths for all integral files
      lena1e=info(2)
      n1eabf=info(3)
      lend2e=info(4)
      n2edbf=info(5)
      lend1e=info(2)
      n1edbf=info(3)
c
      endif
c
c  read in the user input:
c
      call readin(.false.)
      close(unit=nin)
c
c     open drt file
c
      ncsf_max=1
      cpt2(1) = 0
      cpt3(1) = 0
      cpt4(1) = 0
      cpt2_tot = 0
      cpt3_tot = 0
      cpt4_tot = 0
      nrowmax = 0
c
      write(nlist,
     &'(//,3x,6(1h*),''  DRT info section  '',6(1h*)/)')
c
        lenbuf=0
        do ist=1,nst
c
        write(nlist,'(/'' Informations for the DRT no. '',i2)')ist
        ntitle=ntitle+1
c
c  read the drt information from the drtfile:
        call rdhdrt(
     &   title(ntitle), ndoub,       nact,          nrow_f(ist),
     &   ssym(ist),     smult,       nelt,          nela,
     &   nwalk_f(ist),  ncsf_f(ist), lenbuf_f(ist), nmotx,
     &   slabel,        ist)
c
        if (ncsf_f(ist).eq.0) call bummer(
     &   'driver: no csf in the symmetry No ',ist,faterr)
c
        cpt2(ist+1) = cpt2(ist) + forbyt(nsym*nrow_f(ist))
        cpt3(ist+1) = cpt3(ist) + forbyt(nwalk_f(ist))
        cpt4(ist+1) = cpt4(ist) + atebyt(ncsf_f(ist)*navst(ist))
        cpt2_tot = cpt2_tot + forbyt(nsym*nrow_f(ist))
        cpt3_tot = cpt3_tot + forbyt(nwalk_f(ist))
        cpt4_tot = cpt4_tot + atebyt(ncsf_f(ist)*navst(ist))
        if (nrowmax.lt.nrow_f(ist)) nrowmax=nrow_f(ist)
        if (ncsf_max.lt.ncsf_f(ist)) ncsf_max=ncsf_f(ist)
c
        if (lenbuf_f(ist).gt.lenbuf) lenbuf=lenbuf_f(ist)
c
       enddo ! end of: do ist=1,nst
c
c  allocate space for drt arrays:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:syml,8:a,9:b,10:l,
c   11:y,12:nj,13:njs,14:doub,15:modrt_1,16:syml_1,17:buf
c
      nlevel=nact+1
      nnact=nact*(nact+1)/2
      cpt(2)=cpt(1)+cpt2_tot
      cpt(3)=cpt(2)+cpt2_tot
      cpt(4)=cpt(3)+cpt2_tot
      cpt(5)=cpt(4)+cpt3_tot
      cpt(6)=cpt(5)+cpt3_tot
      cpt(7)=cpt(6)+forbyt(nlevel)
      cpt(8)=cpt(7)+forbyt(nlevel)
      cpt(9)=cpt(8)+forbyt(nrowmax)
      cpt(10)=cpt(9)+forbyt(nrowmax)
      cpt(11)=cpt(10)+forbyt(4*nrowmax)
      cpt(12)=cpt(11)+forbyt(4*nsym*nrowmax)
      cpt(13)=cpt(12)+forbyt(nlevel)
      cpt(14)=cpt(13)+forbyt(nlevel)
      cpt(15)=cpt(14)+forbyt(ndoub)
      cpt(16)=cpt(15)+forbyt(nlevel)
      cpt(17)=cpt(16)+forbyt(nlevel)
      cpt(18)=cpt(17)+forbyt(lenbuf)

c
      do ist=1,nst
c
        do j=1,nmot
        orbtyp(j) = 3
        enddo
        call rddrt(ndoub,nrow_f(ist),nwalk_f(ist),ncsf_f(ist),
     &   lenbuf_f(ist),core(cpt(17)),core(cpt(14)),
     &   core(cpt(6)),core(cpt(7)),
     &   core(cpt(12)),core(cpt(13)),core(cpt(8)),core(cpt(9)),
     &   core(cpt(10)),core(cpt(11)),core(cpt(1)+cpt2(ist)),
     &   core(cpt(2)+cpt2(ist)), core(cpt(3)+cpt2(ist)),
     &   core(cpt(5)+cpt3(ist)),core(cpt(4)+cpt3(ist)),occl,
     &   flags(9),core(cpt(15)),core(cpt(16)),ist)
c
        if (navst(ist).gt.ncsf_f(ist)) then
         write(nlist,7050)
     &    ist,ncsf_f(ist),navst(ist),ist
         call bummer('Incorect navst value',navst(ist),faterr)
        endif
7050     format(/' In the DRT no. ',i1,' there are only ',
     & i3,' CSFs, so ',
     & i1,' states are not allowed to be averaged!',
     & /' Change the navst(',i1,') value or the MCDRT!')
c
      enddo ! end of : do ist=1,nst
c
      write(nlist,*)
      call timer('initialization',2,event1,nlist)
c
c     ... for a 'scf' calculation
      if(nact.eq.0) then
         do ist = 1,nst
         ncsf_f(ist) = 1
         enddo
      endif
c
c  using information read from aoints, aoints2 file, drt file,
c  and user input,
c  calculate remaining symmetry array pointers in /csymb/ and
c  the orbital pointer arrays in /corbc/, and determine the
c  addressing array space requirements.
c
      call add1
c
      do 30 i=1,nmotx
         mapmm(i)=i
30    continue
c
c  determine the non-redundant rotations in the active orbital space:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:iwop
c
      cpt(8)=cpt(7)+forbyt(nnact)
      call faar(naar,core(cpt(7)),core(cpt(8)))
c
c  calculate addressing arrays and remaining symmetry arrays:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add
c
cmb   navst
      call add2(naar,core(cpt(8)))
      write(nlist,'(/'' Number of active-double rotations: '',t40,i6)')
     & nadt
      write(nlist,'('' Number of active-active rotations: '',t40,i6)')
     & naar
      write(nlist,'('' Number of double-virtual rotations: '',t40,i6)')
     & nvdt
      write(nlist,'('' Number of active-virtual rotations: '',t40,i6)')
     & nvat
      write(nslist,'(/'' Number of active-double rotations: '',t40,i6)')
     & nadt
      write(nslist,'('' Number of active-active rotations: '',t40,i6)')
     & naar
      write(nslist,'('' Number of double-virtual rotations: '',t40,i6)')
     & nvdt
      write(nslist,'('' Number of active-virtual rotations: '',t40,i6)')
     & nvat
c
c  allocate space for 1-e integral arrays and orbital coefficients.
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:
c
      cpt(9)=cpt(8)+forbyt(addpt(24))
      cpt(10)=cpt(9)+atebyt(nntb)
      cpt(11)=cpt(10)+atebyt(nntb)
      cpt(12)=cpt(11)+atebyt(nbmt)
      cpt(13)=cpt(12)+atebyt(nntm)
c
c  check to make sure enough workspace is available for the calculation:
c  allocate space for integral sort (avcors), for second-half of the
c  mo integral sort (avc2is), and for the hessian block construction
c  (avchbc).
c  *** assume 255 buckets of maximal length ***
c
c   13:lastb(nbuk),14:buk(lenbuk),15:ibuk(npbuk),16:bufs(lenbfs),
c   17:core(avc2is)

       call getlenbfs(lenbfs,numint,nlist) 
c
      cpt(14)=cpt(13)+forbyt(nbukmx)
      avcors=lcore-cpt(14)+1
       if (flags(22)) then
        req2=ldamax+forbyt(numpbf(ldamax)+1)+lenbfs+nvdt
        avcors = avcors-nvdt
       else
        req2=ldamax+forbyt(numpbf(ldamax)+1)+lenbfs
       endif
      avc2is=avcors-req2

c
c   13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,30:qvv,32:scr
c
cmb   please mind - hvect dimensioning changed...
cmd   second change in hmcvec dimensions
c     ncsf --> ncsf_max
c
      navst_max = navst(1)
      do ist=2,nst
          if (navst_max.lt.navst(ist)) navst_max = navst(ist)
          enddo
      reqhbc=ncsf_max*navst_max+nntd+nadt+nnta+nvdt+
     &  nvat+nntv+nnta+numint(1)+n2td+nadt+n2ta+nvdt+nvat+
     &  nadt+nnta+nvat+n2tv
      avchbc=avcors-reqhbc
c
c
c  determine the number of ft reads, integral sorting parameters,
c  and hessian block construction parameters for this iteration:

      call blkasn(avchbc,avc2is,flags(2))
c         ...return means that everything fits as allocated.
      if (nbuk.gt.nbukmx) then
       write(nlist,'(''internal error in MOSORT space allocation'')')
       write(nlist,'(''contact the program administrator'')')
       call bummer(' nbuk no large, nbuk=',nbuk,faterr)
      endif ! (nbuk.gt.nbukmx)

c
c  open the scratch integral file and sort the ao 2-e integrals:
c
c     fname = 'mcscr1'
c     call trnfln( 1, fname )
c     open(unit=scrtfl,file=fname,status='unknown',
c    & form='unformatted')
c
c  allocate core space for 2-e transformation.  on virtual machines
c  such as the vax, it may be useful to restrict this space to force
c  large blocks to be transformed using out-of-core routines.
      kcoret=lcore
      if(.not.(flags(18).and.(niter.eq.1)))then
         if (flags(22).or.flags(24)) then
c    # initialize the timer for the direct AO integral evalulation
*@ifdef direct
*           call timer(' ',1,eventao,6)
*           call timer(' ',5,eventao,6)
*           call mosave(core(cpt(11)),slabel)
*           call kora_md(1,1,mxcsao,nbf_ca)
*@else
            call bummer('DIRECT MCSCF IS NOT AVAILABLE',0,faterr)
*@endif
c
         else
c
c      ...sort the ao 2-e integrals.
         if (.not. (flags(26) .or. flags(27))) then
           call sifsk1(aoint2,info,ierr)
           if (ierr.ne.0) call bummer('from sifsk1, ierr=',ierr,2)
         endif
c
c
         endif
      endif
c
c  allocate space and read in the basis function 1-e integrals:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:buf,14:val,15:ilab
c
      cpt(14)=cpt(13)+atebyt(lena1e)
      cpt(15)=cpt(14)+atebyt(n1eabf)

      if (flags(22).or.flags(24)) then
*@ifdef direct
*        write(nlist,
*     &    '(/,"s-ao and h-ao matrices are calculated direct")')
*        dimsao=cpt(11)-cpt(10)
*        dimh1=cpt(10)-cpt(9)
*        call drive_oneint(core(cpt(10)),core(cpt(9)),atebyt(nntb),
*     &  mxcsao)
*@endif
      else

c   read s integrals
c   this routines returns  cpt9 : s(*)  cpt10: h(*)
c    check for sufficient space

       if (flags(26)) then
          mcone=1
          mcz=0
          mcopt=6
          call rdonexx
     .    (mcirc,mcopt,'Mltpl  0',mcone,core(cpt(10)),mcone,nvalid,
     .    'ONEINT')
           irc=mcirc
          if (irc.ne.0 .or. nvalid.gt.nntb) 
     .    call bummer('reading seward overlap failed',irc,2)
          call rdonexx
     .     (mcirc,mcopt,'oneham  ',mcone,core(cpt(9)),mcone,nvalid,
     .     'ONEINT')
           irc=mcirc
           if (irc.ne.0.or. nvalid.gt.nntb) 
     .     call bummer('reading seward oneham failed',irc,2)
           hcore=0.0d0
           score=0.0d0
       elseif (flags(27)) then
      CALL GPOPEN(LUONEL,'AOONEINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
        hcore=0.0d0
        score=0.0d0
      call getdalton1('OVERLAP ',LUONEL,core(cpt(10)),nntb,0,0,0,0)                    
      call getdalton1('ONEHAMIL',LUONEL,core(cpt(9)),nntb,0,0,0,0)
      CALL GPCLOSE(LUONEL,'KEEP')
       else
          rewind(aoints)
          if (cpt(13)-cpt(9).lt.3*nntb)
     . call bummer ('insufficient space cpt13-cpt9=',cpt(13)-cpt(9),2)
          call sifrsh(aoints,info,core(cpt(13)),core(cpt(14)),
     .             core(cpt(15)),nsym,nbpsy,mapmm,nntb,
     .             core(cpt(9)),score,hcore,symb,ierr)
          if (ierr.ne.0)
     .    call bummer('from sifrsh: ierr=',ierr,2)
c core part of veff
          repnuc=repnuc +hcore




c   change the order of s and h array
c      s:9 h:10 scr:11
          call dcopy_wr(nntb,core(cpt(9)),1,
     .       core(cpt(11)),1)
c      s:9 h:10 s:11
         call dcopy_wr(nntb,core(cpt(10)),1,core(cpt(9)),1)
         call dcopy_wr(nntb,core(cpt(11)),1,core(cpt(10)),1)
c   that's it.
         call sifsk1(aoint2,info,ierr)
         if (ierr.ne.0)
     .    call bummer ('from sifsk1 ierr=',ierr,2)
       endif
       endif

      if(flags(4))call plblks('s-ao matrix',
     &  core(cpt(10)),nsym,nbpsy,' bfn',1,nlist)
      if(flags(4))call plblks('h-ao matrix',
     &  core(cpt(9)),nsym,nbpsy,' bfn',1,nlist)
c
c
c
c  open the remaining sequential files required by the program:
c
      fname = 'restart'
      call trnfln( 1, fname )
      open(unit=nrestp,file=fname,status='unknown',
     &  form='unformatted')
c
ct    fname = 'mcscr2'
ct    open(unit=ntemp,file=fname,status='unknown',form='unformatted')
cmd
c
      fname = 'mcscr3'
      call trnfln( 1, fname )
      open(unit=stape,file=fname,status='unknown',form='unformatted')
c
      fname = 'mchess'
      call trnfln( 1, fname )
      open(unit=htape,file=fname,status='unknown',form='unformatted')
c
      if ((.not.flags(11)).and.flags(12)) then
      fname = 'hdiagf'
      open(unit=hdiagf,file=fname,status='unknown',form='unformatted')
      endif
c
c  get the initial orbital coefficients:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:t,14:ht,15:u,16:sx
c
      cpt(14)=cpt(13)+atebyt(nbmt)
      cpt(15)=cpt(14)+atebyt(nntm)
      cpt(16)=cpt(15)+atebyt(n2tm)
      call inorbc(core(cpt(11)),core(cpt(9)),core(cpt(10)),
     &  core(cpt(13)),core(cpt(14)),core(cpt(15)),core(cpt(16)))
cmd
c     initialization of the direct 2-e evaliation
c      (means statistic run fir KORA)

      if (flags(22).or.flags(24)) then
*@ifdef direct
*         call mosave(core(cpt(11)),slabel)
*         call kora_md(1,0,mxcsao,nbf_ca)
*         idim1 = ivcmax(8,nbpsy)
*         idim2 = ivcmax(8,nmpsy)
*         idim = max(idim1*idim1,idim2*idim2,ipq(nbf_ca+1),
*     &    ifock(nbf_ca+1))
*@endif
      endif
cmd
c
c  begin the mcscf iterative procedure:
c
      write(nlist,*)
      write(nslist,*)
      write(nslist,6231)
      convrg='*not conv.*'
      qconvg=.false.
      iter=0
c
c michal2{
c
      if(cosmocalc.ne.0) then
c
c Cosmo initialization
c
c****** Statements for cosmo calculation ********
c Silmar - initialization of Vsolv at the beginning
                 dbglvl = 0
      call wzero(maxdens,Vsolv,1)
      allocate(cosurf(3,2*maxnps),stat=ierr)
       if (ierr.ne.0)
     .  call bummer('allocating cosurf failed',0,2)
      allocate(a1mat(maxnps*(maxnps+1)/2),stat=ierr)
       if (ierr.ne.0)
     .  call bummer('allocating a1mat failed',0,2)
      allocate(a2mat(maxnps,maxnps),stat=ierr)
       if (ierr.ne.0)
     .  call bummer('allocating a2mat failed',0,2)
      allocate(a3mat(maxnps*(maxnps+1)/2),stat=ierr)
       if (ierr.ne.0)
     .  call bummer('allocating a3mat failed',0,2)
c
        Call cosmoinitial(mmtype,natoms,xyz,symgrp,nsymm)
c
         call cosmo(1,xyz,mmtype,natoms,nsymm,cosurf,nps,npspher
     . ,phi,qcos,ediel,elast,a1mat,a2mat,a3mat,dcos)
         call cosmo(2,xyz,mmtype,natoms,nsymm,cosurf,nps,npspher
     . ,phi,qcos,ediel,elast,a1mat,a2mat,a3mat,dcos)
c
           if(dbglvl.eq.1) then
         write(nlist,*)
         write(nlist,*)'nps,nsymm,npspher = ',nps,nsymm,npspher
         write(nlist,*)'cosurf = '
         write(nlist,*)((cosurf(i,j),i=1,3),j=1,nps)
         write(nlist,*)
         write(nlist,*)'Is the cosmo initialization correct???'
c
           endif ! dbglvl
c
      endif ! for cosmocalc.ne.0
c
c michal2}
c
      do 1000 itert=1,niter
      iter=itert
c  for vdisk
      qsort=iter.eq.1
c     qsort=.true.
      write(nlist,6200)iter
6200  format(/t15,' starting mcscf iteration...',i4)
      call timer(' ',2,event1,nlist)
      call timer(' ',1,event2,nlist)
c
c  determine if orbital-state coupling is included this iteration:
c
cmd   qcoupl=flags(2).and.flags(11).and.
cmd  &       (iter.ge.ncoupl).and.(wnorm.le.tol(9))
      qcoupl=flags(2).and.
     &       (iter.ge.ncoupl).and.(wnorm.le.tol(9))
      if(qcoupl)then
          write(nlist,6230)' '
      else
          write(nlist,6230)' not '
      endif
c
c  orthonormalize the orbital coefficients:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:scr
c
      call sorthn(core(cpt(11)),core(cpt(10)),core(cpt(13)))
      if(flags(4))call prblks('orthonormalized orbitals',
     &  core(cpt(11)),nsym,nbpsy,nmpsy,' bfn','  mo',1,nlist)
c
c michal2{
c
c
      if(cosmocalc.ne.0.and.iter.ge.2) then
c
c Silmar - free memory to be used by driverDALTON
c
       freemem = lcore -(cpt(13)-1)
c
       do i = 1,8
c
       nbfpsy(i)=nbpsy(i)
c
       enddo
c
       write(nlist,*)'************************************************'
       write(nlist,*)'**Call for electrostatic potential calculation**'
       write(nlist,*)'************************************************'
c
c
      call driverDALTON(core(cpt(13)),cosurf,phi,qcos,nps
     . ,freemem,0,0)
c
c             write(nlist,*)
c             write(nlist,*)'cosurf_and_phi'
c          do i = 1,nps
c      write(nlist,11239) cosurf(1,i),cosurf(2,i),cosurf(3,i)
c     . ,phi(i)
c          enddo
c             write(nlist,*)'end_of_phi'
c             write(nlist,*)
c
       write(nlist,*)
      write(nlist,*)'**************************************************'
      write(nlist,*)'***Cosmo calculation for ground state*************'
      write(nlist,*)'**************************************************'
       write(nlist,*)
c
         call cosmo(3,xyz,mmtype,natoms,nsymm,cosurf,nps,npspher
     . ,phi,qcos,ediel,elast,a1mat,a2mat,a3mat,dcos)
c
             write(nlist,*)
             write(nlist,*)'ediel after cosmo(3', ediel
             write(nlist,*)
c             write(nlist,*)'qcos after cosmo(3 = '
c             write(nlist,*)(qcos(i),i=1,nps)
c             write(nlist,*)
c
            sumq1=0.0d0
            sumq2=0.0d0
            sumq=0.0d0
c
             write(nlist,*)
             write(nlist,*)'cosurf_and_qcos_from cosmo(3'
              do i=1,nps
c         write(nlist,11239) cosurf(1,i),cosurf(2,i),cosurf(3,i)
c     . ,qcos(i)
                 if(qcos(i).lt.0) then
              sumq1 = sumq1 + qcos(i)
                 else
              sumq2 = sumq2 + qcos(i)
                 endif
              enddo
              sumq = sumq1 + sumq2
            write(nlist,*)'end_of_qcos'
            write(nlist,*)
            write(nlist,*)'Sum of negative charges = ',sumq1
            write(nlist,*)'Sum of positive charges = ',sumq2
            write(nlist,*)'Total sum = ',sumq
            write(nlist,*)
c
c This copying is to keep the qcos array non-modified
c
         call dcopy_wr(nps,qcos,1,qcoszero,1)
c
c         write(nlist,*)'cosurf = '
c              do i=1,nps
c         write(nlist,11239) cosurf(1,i),cosurf(2,i),cosurf(3,i),
c     . qcoszero(i)
c              enddo
c
11239    format(4f16.12)
c
c Call for solvent modified integrals
c
c        call wzero(nps,qcos,1)
c
c The charges have to be multiplied by
c fepsi before they are used to calculate the solvent modified one-
c electron integrals, also for ground-state calculation!!
c
c
             write(nlist,*) 'fepsi = ', fepsi
c
         call dscal_wr(nps,fepsi,qcoszero,1)
c
      write(nlist,*)'*************************************************'
      write(nlist,*)'*Call for solvent-modified integrals calculation*'
      write(nlist,*)'*************************************************'
c
      call driverDALTON(core(cpt(13)),cosurf,phi,qcoszero,nps
     . ,freemem,1,0)
c
c The solvent modified nuc. rep. energy is transfered to
c the nucrep variable
c
         repnuc = Vnenuc
c
      endif! for cosmocalc.ne.0
c
c michal2}
c
c  transform 1-e h(*):
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,9:add,10:h1ao,10:sao,
c   11:c,12:h1mo,13:scr
c
c michal2{
c
           if(iter.ge.2.and.cosmocalc.ne.0) then
c
c      call plblks('Vsolv matrix',
c     &  Vsolv,nsym,nmpsy,'  ao',1,nlist)
c
      call tran1(core(cpt(11)),Vsolv,core(cpt(12)),
     &  core(cpt(13)))
c
c       call plblks('h-mo matrix',
c     &  core(cpt(12)),nsym,nmpsy,'  mo',1,nlist)
c
           else !for cosmocalc.ne.0
c
c michal2}
c
      call tran1(core(cpt(11)),core(cpt(9)),core(cpt(12)),
     &  core(cpt(13)))
c
c michal2{
c
           endif !for cosmocalc.ne.0
c
c michal2}
c
      write(nlist,'(/," *** Starting integral transformation ***")')
      if(flags(4))call plblks('h-mo matrix',
     &  core(cpt(12)),nsym,nmpsy,'  mo',1,nlist)
c
c  transform the 2-e integrals:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:scr
c
      if(.not.flags(18))then
         call timer(' ',1,event3,nlist)
         if (flags(22).or.flags(24)) then
*@ifdef direct
*           call mosave(core(cpt(11)),slabel)
*           call kora_md(0,0,mxcsao,nbf_ca)
*@endif
c
         else

          sewd=0
          lumorb=0
          dlt2=0
          if (flags(26))  sewd=1
          if (flags(27))  dlt2=1
          if (flags(29))  lumorb=1
          lcoremain=lcore

c
c      sort the MO integrals 
c      additional parameters  nbuk,nvdt,cpt(35),cpt(14),cpt(13)
        if (flags(22)) then
          cpt(35)=cpt(13)+forbyt(nbuk)
          cpt(14)=cpt(35)+nvdt
c     initialize the bdg2e
c     call wzero(nvdt,core(cpt(35)),1)
        else
         cpt(14) = cpt(13)+forbyt(nbuk)
        endif
c
c  sort the transformed integrals:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,35:bdg2e,14:scratch

c


          call driver_wrap( core(cpt(14)),core(cpt(14)),lcore-cpt(14),  
     .                    nlist,
     .             nbft,nvpsy,nsym,nbpsy,sewd,dlt2,lumorb,core(cpt(11)),
     .                    core(cpt(10)),qcoupl,qsort,qconvg,lcoremain,
     .                    flags(31),
     .                    core(cpt(13)),core(cpt(8)),core(cpt(12)),
     .                    core(cpt(35)),lenbuk,nipbuk,ecore0 )
c         avchbc=avchbc-lcore+lcoremain
c         avc2is=avc2is-lcore+lcoremain
c         avcors=avcors-lcore+lcoremain
         endif
         call timer('2-e transformation',4,event3,nlist)
      else
c       restart case, single iteration (density) 
c       skip MO trafo 
        if (flags(22)) then
        call bummer(' prior mosort, skip MO trafo, disable flag 22',0,0)
          cpt(35)=cpt(13)+forbyt(nbuk)
          cpt(14)=cpt(35)+nvdt
c     initialize the bdg2e
c     call wzero(nvdt,core(cpt(35)),1)
        else
         cpt(14) = cpt(13)+forbyt(nbuk)
        endif
        flags(22)=.false.
      endif
      flags(18)=.false.
c
      avcors=lcore-cpt(14)+1
c
c
      call mosort(
     & core(cpt(14)),  core(cpt(13)),   avcors,    avc2is,
     & core(cpt(8)),   core(cpt(12)),   ntemp,     lenm2e,
     & n2embf,         lenbuk,          nipbuk,    ecore0,
     & core(cpt(35)))
c
c  transpose 1-e hamiltonian matrix from orbital-packed to
c  sub-block packed.
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,35:bdg2e,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv
c
c     change in hmcvec dimension ncsf --> cpt4tot
      cpt(15)=cpt(14)+cpt4_tot
      cpt(16)=cpt(15)+atebyt(n2td)
      cpt(17)=cpt(16)+atebyt(nadt)
      cpt(18)=cpt(17)+atebyt(n2ta)
      cpt(19)=cpt(18)+atebyt(nvdt)
      cpt(20)=cpt(19)+atebyt(nvat)
      cpt(21)=cpt(20)+atebyt(tsymvv(1))
c
cmd
      if(flags(22)) then
c     dimensions for scratch fields used only in u_2e
*@ifdef direct
*         cpt(22)=cpt(21)+ atebyt(idim)
*         cpt(23)=cpt(22)+ atebyt(idim)
*         cpt(24)=cpt(23)+ atebyt(idim)
*         cpt(25)=cpt(24)+ atebyt(idim)
*         call u_2e(core(cpt(20)),core(cpt(19)),core(cpt(11)),
*     &    core(cpt(21)),core(cpt(22)),core(cpt(23)),core(cpt(24)),
*     &    idim,mxcsao,core(cpt(25)),lcore)
*@endif
      else
         call wzero(cpt(21)-cpt(19),core(cpt(19)),1)
      endif
cmd
      dpt=ix0(1)+1
      apt=ix0(2)+1
      vpt=ix0(3)+1
      call tltpsp(invx(dpt),invx(apt),invx(vpt),
     &  core(cpt(12)),core(cpt(15)),core(cpt(16)),core(cpt(17)),
     &  core(cpt(18)),core(cpt(19)),core(cpt(20)))

      if(flags(4).and.(ndot.ne.0))call plblks('modified h-mo matrix',
     &  core(cpt(12)),nsym,nmpsy,'  mo',1,nlist)
c
c  solve for the appropriate h(mc) eigenvector:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,35:bdg2e,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:scr
c
c     avchd=lcoremain-cpt(21)+1
      avchd=lcore-cpt(21)+1
      emc = 0.d+00
      rewind iunits(20)
      do ist = 1,nst
      ncol(ist) = min(ncol(ist),ncsf_f(ist))
      call hmcvec(
     & core(cpt(21)), avchd,        ncsf_f(ist),    repnuc,
     & ecore0,        lenbft,       core(cpt(13)),  lenbuk,
     & nipbuk,        core(cpt(6)), core(cpt(17)),  core(cpt(8)),
     &  core(cpt(1)+cpt2(ist)),   core(cpt(2)+cpt2(ist)),
     &  core(cpt(3)+cpt2(ist)),   core(cpt(4)+cpt3(ist)),
     &  core(cpt(5)+cpt3(ist)),   wnorm,  emc, core(cpt(14)+cpt4(ist)),
     &  ist)
      enddo
      demc=emcold-emc
      emcold=emc
c
c  construct density matrices, fock matrices, and hessian blocks:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,35:bdg2e,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:scr
c
      cpt(22)=cpt(21)+atebyt(nnta)
      cpt(23)=cpt(22)+atebyt(numint(1))
      cpt(24)=cpt(23)+atebyt(n2td)
      cpt(25)=cpt(24)+atebyt(nadt)
      cpt(26)=cpt(25)+atebyt(n2ta)
      cpt(27)=cpt(26)+atebyt(nvdt)
      cpt(28)=cpt(27)+atebyt(nvat)
      cpt(29)=cpt(28)+atebyt(nadt)
      cpt(30)=cpt(29)+atebyt(nnta)
      cpt(31)=cpt(30)+atebyt(nvat)
      cpt(32)=cpt(31)+atebyt(n2tv)
c
      call hbcon(avchbc,core(cpt(32)),
     &  core(cpt(21)),core(cpt(22)),core(cpt(13)),core(cpt(8)),
     &  lenbft,lenbuk,nipbuk,core(cpt(14)),core(cpt(7)),
     &  naar,core(cpt(1)),core(cpt(2)),
     &  core(cpt(3)),core(cpt(4)),core(cpt(5)),core(cpt(6)),
     &  core(cpt(15)),core(cpt(16)),core(cpt(17)),core(cpt(18)),
     &  core(cpt(19)),core(cpt(20)),
     &  core(cpt(23)),core(cpt(24)),core(cpt(25)),core(cpt(26)),
     &  core(cpt(27)),core(cpt(28)),core(cpt(29)),core(cpt(30)),
     &  core(cpt(31)))
c
c  construct the gradient vector from the fock matrices.
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,35:bdg2e,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w
c
      ndimw=nadt+naar+nvdt+nvat
      call wcon(core(cpt(24)),core(cpt(28)),core(cpt(25)),core(cpt(26)),
     &  core(cpt(27)),core(cpt(7)),nadt,naar,nvdt,nvat,ndimw,
     &  core(cpt(32)),wnorm)
c
c  solve for the orbital corrections:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,35:bdg2e,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:scr
c
c
      cpt(33)=cpt(32)+atebyt(ndimw)
      cpt(34)=cpt(33)+atebyt(ndimw)
c     avcsk=lcoremain-cpt(34)+1
      avcsk=lcore-cpt(34)+1

cmb   in contrast to the original, when state averaging, solvek must
cmb   get all eigenvalues, in order to construct properly the multi-
cmb   ple "m"-blocks of the hessian -heig instead of emc...
cmb   original ncsf passed to solvek, too - not only ncfx

cmd
      if (flags(25)) then
      write(nlist,'(/,"Total memory alocated:",i15)') lcore
      write(nlist,'("Free memory before entering SOLVEK:",i15)') avcsk
      write(nlist,'("xbar  :",i15)') cpt(2)-cpt(1)
      write(nlist,'("xp    :",i15)') cpt(3)-cpt(2)
      write(nlist,'("z     :",i15)') cpt(4)-cpt(3)
      write(nlist,'("r     :",i15)') cpt(5)-cpt(4)
      write(nlist,'("ind   :",i15)') cpt(6)-cpt(5)
      write(nlist,'("modrt :",i15)') cpt(7)-cpt(6)
      write(nlist,'("iorder:",i15)') cpt(8)-cpt(7)
      write(nlist,'("add   :",i15)') cpt(9)-cpt(8)
      write(nlist,'("h1ao  :",i15)') cpt(10)-cpt(9)
      write(nlist,'("sao   :",i15)') cpt(11)-cpt(10)
      write(nlist,'("c     :",i15)') cpt(12)-cpt(11)
      write(nlist,'("h1mo  :",i15)') cpt(13)-cpt(12)
      if (flags(22)) then
         write(nlist,'("lastb :",i15)') cpt(35)-cpt(13)
         write(nlist,'("bdg2e :",i15)') cpt(14)-cpt(35)
         write(nlist,'("hvec  :",i15)') cpt(15)-cpt(14)
      else
         write(nlist,'("lastb :",i15)') cpt(14)-cpt(13)
         write(nlist,'("bdg2e :",i15)') 0
         write(nlist,'("hvec  :",i15)') cpt(15)-cpt(14)
      endif
      write(nlist,'("udd   :",i15)') cpt(16)-cpt(15)
      write(nlist,'("uad   :",i15)') cpt(17)-cpt(16)
      write(nlist,'("uaa   :",i15)') cpt(18)-cpt(17)
      write(nlist,'("uvd   :",i15)') cpt(19)-cpt(18)
      write(nlist,'("uva   :",i15)') cpt(20)-cpt(19)
      write(nlist,'("uvv   :",i15)') cpt(21)-cpt(20)
      write(nlist,'("d1    :",i15)') cpt(22)-cpt(21)
      write(nlist,'("d2    :",i15)') cpt(23)-cpt(22)
      write(nlist,'("fdd   :",i15)') cpt(24)-cpt(23)
      write(nlist,'("fad   :",i15)') cpt(25)-cpt(24)
      write(nlist,'("faa   :",i15)') cpt(26)-cpt(25)
      write(nlist,'("fvd   :",i15)') cpt(27)-cpt(26)
      write(nlist,'("fva   :",i15)') cpt(28)-cpt(27)
      write(nlist,'("qad   :",i15)') cpt(29)-cpt(28)
      write(nlist,'("qaa   :",i15)') cpt(30)-cpt(29)
      write(nlist,'("qva   :",i15)') cpt(31)-cpt(30)
      write(nlist,'("qvv   :",i15)') cpt(32)-cpt(31)
      write(nlist,'("w     :",i15)') cpt(33)-cpt(32)
      write(nlist,'("k     :",i15,/)') cpt(34)-cpt(33)
      endif
cmd
      call solvek(
     & avcsk,          core(cpt(34)),   nadt,           naar,
     & nvdt,           nvat,            ndimw,          core(cpt(7)),
     & core(cpt(8)),   wnorm,           core(cpt(32)),  core(cpt(14)),
     & core(cpt(33)),  knorm,           apxde,          core(cpt(11)),
     & core(cpt(21)),  core(cpt(35)),   core(cpt(27)),  core(cpt(30)),
     & core(cpt(24)),  core(cpt(28)),   core(cpt(23)),  core(cpt(31)),
     & idim,           mxcsao,          lenbft,         core(cpt(6)),
     & core(cpt(2)),   core(cpt(3)),    core(cpt(1)),   core(cpt(4)),
     & core(cpt(5)),   lenbuk,          nipbuk,         core(cpt(13)),
     & core(cpt(17)),  core(cpt(16)),   core(cpt(18)),  core(cpt(19)),
     & qcoupl)
c
c     call da2cr
      call daclrd(5)
c
c  transform the orbitals for the next iteration:
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cn,38:sa,39:sx,40:ia,41:ib
c
      cpt(35)=cpt(34)+atebyt(nbmt)
      cpt(36)=cpt(35)+atebyt(nbmt)
      cpt(37)=cpt(36)+atebyt(nmot)
      cpt(38)=cpt(37)+atebyt(nbmt)
      cpt(39)=cpt(38)+atebyt(n2tm)
      cpt(40)=cpt(39)+atebyt(nmot)
      cpt(41)=cpt(40)+forbyt(nmot)
      cpt(42)=cpt(41)+forbyt(nmot)
      cpt(43)=cpt(42)+atebyt(nmot)
      cpt(44)=cpt(43)+atebyt(nmot)
c     if(cpt(44)-1 .gt. lcoremain)then
      if(cpt(44)-1 .gt. lcore)then
         write(nlist,6260)'motran',(cpt(i),i=1,44)
         call bummer('mcscf: allocation error. lcore=',lcore,faterr)
      endif
      call motran(core(cpt(11)),core(cpt(33)),core(cpt(23)),
     &  core(cpt(25)),core(cpt(29)),core(cpt(31)),core(cpt(21)),
     &  core(cpt(7)),naar,
     &  core(cpt(34)),core(cpt(35)),core(cpt(36)),core(cpt(37)),
     &  core(cpt(38)),core(cpt(39)),core(cpt(40)),core(cpt(41)),
     &  core(cpt(42)),core(cpt(43)) )
c
c  switch orbital coefficients in c(*) and cn(*) for the next
c  mcscf iteration.
c
      call dswap_wr(nbmt,core(cpt(37)),1,core(cpt(11)),1)
c
c  save iteration information on the restart file.
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:sa,39:sx,40:ia,41:ib
c
      qconvg=abs(demc).le.tol(1).and.wnorm.le.tol(2).and.knorm.le.tol(3)
     &  .and.apxde.le.tol(4)
      call restrt(
     &  nsym,           lcinf,          nmot,           nbmt,
     &  nmpsy,          nbpsy,          ntitle,         qconvg,
     &  title,          coninf,         core(cpt(34)),  core(cpt(35)),
     &  core(cpt(36)),  core(cpt(11)),  core(cpt(14)),  core(cpt(37)))
c
c  end of mcscf iteration:
c
      call timer('mcscf iteration',4,event2,nlist)
c
c  check convergence of mcscf iterations:
c
      if(qconvg)then
         write(nlist,6210)
     &    'all mcscf convergence criteria are satisfied.'
         convrg='*converged*'
         go to 1010
      else
         write(nlist,6210)
     &    'not all mcscf convergence criteria are satisfied.'
      endif
      write(nlist,6220)iter,emc,demc,wnorm,knorm,apxde,convrg
      write(nslist,6225)iter,emc,demc,wnorm,knorm,apxde,qcoupl,convrg
c
1000  continue
1010  continue

          qconvg=.true.
c  darf nicht core (cpt(13): lcore) als scratch benutzen
c  ueberschreibt MOs! falls qconvg=true wird keine AO-MO Trafo
c  durchgefuehrt; alternativ muessen die noch weiter benoetigten
c  core-Bereiche nach cpt(13) nach vorne verschoben werden und
c  mit core(top) state core(cpt(13)) aufgerufen werden.
c  die letzten beiden Zeilen dummy Argumente duerfen nicht verwendet werden.

          call driver_wrap( core(cpt(13)),core(cpt(13)),lcore-cpt(13),  
     .                    nlist,
     .            nbft,nvpsy,nsym,nbpsy,sewd,dlt2,lumorb,core(cpt(11)),
     .                    core(cpt(10)),qcoupl,qsort,qconvg,lcoremain,
     .            flags(31),
     .            core(cpt(13)),core(cpt(8)),core(cpt(12)),
     .            core(cpt(35)),lenbuk,nipbuk,ecore0 )


      
c
c michal2{
c
         if(cosmocalc.ne.0) then
c
c    Energy of last iteration, to be used by cosmo
c
         if(nsymm.gt.0) then
         call replicate(nsymm,1.d0,1.d0,1.d0,np,xfinal,yfinal,zfinal)
         isymf = np
         symf = dble(isymf)
         endif
c
           if(symf.eq.0) symf=1.0d0
c
         elast = emc-symf*ediel
c
         write(nlist,*)'symf = ', symf
         write(nlist,*)'elast = ', elast
         write(nlist,*)'ediel = ', ediel
c
       write(nlist,*)
      write(nlist,*)'**************************************************'
      write(nlist,*)'***Final cosmo calculation for ground state*******'
      write(nlist,*)'**************************************************'
       write(nlist,*)
c
c Final electrostatic potential calculation
c
       freemem = lcore -(cpt(38)-1)
c
       call driverDALTON(core(cpt(38)),cosurf(1,nps+1),phi,
     .qcos,npspher,freemem,0,0)
c
c             write(nlist,*)
c             write(nlist,*)'cosurf_and_pot_on_the_outer_cavity'
c
c          do i = 1,npspher
c      write(nlist,11239) cosurf(1,i+nps),cosurf(2,i+nps),
c     .  cosurf(3,i+nps),phi(i)
c          enddo
c             write(nlist,*)'end_of_cosurf_and_pot_on_the_outer_cavity'
             write(nlist,*)
       call cosmo(4,xyz,mmtype,natoms,nsymm,cosurf,nps,npspher
     . ,phi,qcos,ediel,elast,a1mat,a2mat,a3mat,dcos)
c
         endif !for cosmocalc.ne.0
c
c michal2}
c
c  write out convergence info for the final iteration:
c
      write(nlist,6210)'final mcscf convergence values:'
      write(nlist,6220)iter,emc,demc,wnorm,knorm,apxde,convrg
      write(nslist,6210)'final mcscf convergence values:'
      write(nslist,6225)iter,emc,demc,wnorm,knorm,apxde,qcoupl,convrg
c
c    find smallest total energy
c
      do ii=1,nst
       do ij=1,navst(ii)
        emin=min(emin,heig(ii,ij))
       enddo
      enddo

c
c    print out relative to lowest state (from 2.15.2.6)  
c
      emin=100.0d0
      do ii=1,nst
       do ij=1,navst(ii)
        emin=min(emin,heig(ii,ij))
       enddo
      enddo


      write(nslist,6221)
      write(nlist,6221)
      do ii=1,nst
       do ij=1,navst(ii)
         enev = (heig(ii,ij) - emin)*au2ev
         write(nslist,6222) ii,ij,wavst(ii,ij),heig(ii,ij),enev
         write(nlist,6222) ii,ij,wavst(ii,ij),heig(ii,ij),enev
       enddo
      enddo
      write(nlist,6223)
      write(nslist,6223)
6221  format(1x//,1x//3x,
     .'---------Individual total energies for all states:----------')
6222  format(3x,'DRT #',i1,' state #',i2,' wt ',f5.3,
     .    ' total energy= ',f18.9,', rel. (eV)=',f11.6)
6223  format(3x,15('----')//)


c
c  print the last set of orbitals.
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:sa,39:sx,40:ia,41:ib
c
      cpt(39)=cpt(38)+atebyt(nmot)
      if(cpt(39)-1 .gt. lcore)then
        write(nlist,6260)'mowrit call:',(cpt(i),i=1,39)
        call bummer('mcscf: allocation error. lcore=',lcore,faterr)
      endif

1112    format (' MO ',i4)
      do i=1,nbft
          write(bfnlabmo(i),1112) i
      enddo
cfp: write out the mocoef_mc and nocoef_mc file here, rather than through calling the mofmt.x program
c
c assign dummy occupations to MCSCF orbitals
c    docc -> 2, act -> 1, virt -> 0
      do i = 1,nmot
         if (orbtyp(i).eq.1) then
            core(cpt(38)+i-1) = two
         elseif (orbtyp(i).eq.2) then
            core(cpt(38)+i-1) = one
         elseif (orbtyp(i).eq.3) then
            core(cpt(38)+i-1) = zero        
         else
            core(cpt(38)+i-1) = -1d0
         endif
      enddo
      
      ntitle=4
      title(1)='MO-coefficients from mcscf.x'
      title(2)=' with dummy occupation 1.0 for active orbitals'
      
      write(title(3),6100) 'total ao core energy =',repnuc
      write(title(4),1113) emc
1113  format ('MCSCF energy =',f18.9)
      
      open(unit=mocout,file='mocoef_mc',status='unknown',
     &    form='formatted')
      call mowrit(mocout,10,filerr,syserr,ntitle,title,cfmt,nsym,
     &     nbpsy,nmpsy,slabelao,core(cpt(37)))
      call mowrit(mocout,20,filerr,syserr,ntitle,title,cfmt,nsym,
     &     nbpsy,nmpsy,slabelao,core(cpt(37)))
      call mowrit(mocout,40,filerr,syserr,ntitle,title,cfmt,nsym,
     &     nbpsy,nmpsy,slabelao,core(cpt(38)))
         close(unit=mocout)

      open (unit=mocout,file='mocoef_mc.lumorb',form='formatted',
     .      status='unknown')
      call wrtmof_molcas(mocout, nsym,nmpsy,nbpsy,core(cpt(37)),'mos',
     .      'mocoef_mc')
      call wrtmof_molcas(mocout, nsym,nmpsy,nbpsy,core(cpt(38)),'occ',
     .      'mocoef_mc')
      close (mocout)


      
c
      ntitle=3
      title(1)='NO-coefficients from mcscf.x'
      write(title(2),6100) 'total ao core energy =',repnuc
      write(title(3),1113) emc
      
      open(unit=mocout,file='nocoef_mc',status='unknown',
     &    form='formatted')
      call mowrit(mocout,10,filerr,syserr,ntitle,title,cfmt,nsym,
     &     nbpsy,nmpsy,slabelao,core(cpt(35)))
      call mowrit(mocout,20,filerr,syserr,ntitle,title,cfmt,nsym,
     &     nbpsy,nmpsy,slabelao,core(cpt(35)))
      call mowrit(mocout,40,filerr,syserr,ntitle,title,cfmt,nsym,
     &     nbpsy,nmpsy,slabelao,core(cpt(36)))
         close(unit=mocout)
      open (unit=mocout,file='nocoef_mc.lumorb',form='formatted',
     .status='unknown')
      call wrtmof_molcas(mocout, nsym,nmpsy,nbpsy,core(cpt(35)),'mos',
     .   'nocoef_mc')
      call wrtmof_molcas(mocout, nsym,nmpsy,nbpsy,core(cpt(36)),'occ',
     .   'nocoef_mc')
      close (mocout)
c
cfp: print out MOs only if there are not too many
      if ((nbft.le.orbmxpri).or.flags(32))then
        if (flags(22).or.flags(24)) then
          call prblks('mcscf orbitals of the final iteration,',
     &  core(cpt(37)),nsym,nbpsy,nmpsy,' bfn','  mo',1,nlist)
        else
          call prsbkc('mcscf orbitals of the final iteration,',
     .    core(cpt(37)),nsym,slabelao,
     .    nbpsy,nmpsy,bfnlabao,bfnlabmo, 1,nlist)
        endif
      else
        write(nlist,*)'MO-coefficient print-out skipped (no flag 32)'
        write(nlist,*)'They may be found in the MOCOEF directory.'
      endif
c
      if(flags(17))then
c        ...print the natural orbitals and occupations.
         impt=cpt(35)
         ivpt=cpt(36)
         icount=1
         do 1050 isym=1,nsym
            if(nmpsy(isym).ne.0)then
c              write(nlist,'(/10x,a,i3,x3,a)')
c    &    'natural orbitals of the final iteration,block',isym,slabelao(isym)
               write(nlist,1123)'natural orbitals of the final
     & iteration,block', isym,slabelao(isym)
1123  format(/10x,a,i3,3x,' - ',a)
              if ((nbft.le.orbmxpri).or.flags(32))then
                call prvbkc(' ',core(impt),core(ivpt),
     &          nbpsy(isym),nbpsy(isym),nmpsy(isym),
     &          bfnlabao(icount),bfnlabmo,' occ(*)=',1,nlist)
              else
                jlast = 0
                do jstrt = 1, nmpsy(isym), 8
                  jlast = min( nmpsy(isym), jlast+8 )
         write(nlist,1120) ( bfnlabmo(j), j = jstrt, jlast )
         write(nlist,1121) ' occ(*)=', (core(ivpt+j-1), j=jstrt,jlast)
c         write(nlist,*)
                enddo
              endif
              impt=impt+nbpsy(isym)*nmpsy(isym)
              ivpt=ivpt+nmpsy(isym)
              icount = icount + nbpsy(isym)
            endif
1050     continue
      endif
1120   format(8x,8(6x,a8,1x))
1121   format(1x,a8,8f15.8)
c
      if(flags(21))then
         call timer(' ',1,event3,nlist)
ctm
c already set
ctm
c  write out d1(*), f(*), and q(*) from the last iteration
c  to the file mcd1fl.
c
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:buffer,39:valbuf,40:labbuf,41:d1,42:fmc,43:qmc,44:mapmm,
c   45: mapim
c
         cpt(39)=cpt(38)+atebyt(lend1e)
         cpt(40)=cpt(39)+atebyt(n1edbf)
         cpt(41)=cpt(40)+forbyt(2*n1edbf)
         cpt(42)=cpt(41)+atebyt(nntm)
         cpt(43)=cpt(42)+atebyt(nntm)
         cpt(44)=cpt(43)+atebyt(nntm)
         cpt(45)=cpt(44)+forbyt(nmot)
         cpt(46)=cpt(44)+forbyt(nmot)
         cpt(47)=cpt(46)+forbyt(nmot)
         cpt(48)=cpt(47)+forbyt(nmot)
         cpt(49)=cpt(48)+atebyt(nmot)
         if(cpt(49)-1 .gt. lcore)then
            write(nlist,6260)'wmcd1f call:',(cpt(i),i=1,49)
            call bummer('mcscf: allocation error. lcore=',lcore,faterr)
         endif

c
         fname = 'mcd1fl'
         call trnfln( 1, fname )
         open(unit=mcd1fl,file=fname,status='unknown',
     &    form='unformatted')
c
         call wmcd1f(nlist,mcd1fl,demc,knorm,wnorm,apxde,repnuc,
     &    orbtyp,
     &    core(cpt(38)),core(cpt(39)),core(cpt(40)),ntitle,title,
     &    nsym,nmot,nmpsy,ndpsy,napsy,nvpsy,
     &    emc,invx(ix0(1)+1),invx(ix0(2)+1),invx(ix0(3)+1),
     &    ird1,ird2,ira1,ira2,irv1,irv2,irm1,irm2,
     &    nsa,n2sa,nsm,nnsm,nntm,nsd,n2sd,nsv,n2sv,fcimsk,
     &    nndx,core(cpt(21)),core(cpt(25)),core(cpt(23)),
     &    core(cpt(29)),core(cpt(28)),core(cpt(30)),core(cpt(31)),
     &    core(cpt(41)),core(cpt(42)),core(cpt(43)),core(cpt(44)),
     &    core(cpt(45)),core(cpt(46)),core(cpt(47)),core(cpt(48)))
              close(unit=mcd1fl)
c
c  write out the 2-particle density matrix to mcd2fl.
c
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:buffer,39:valbuf,40:labbuf
c
         cpt(39)=cpt(38)+atebyt(lend2e)
         cpt(40)=cpt(39)+atebyt(n2edbf)
         cpt(41)=cpt(40)+forbyt(4*n2edbf)
c
         fname = 'mcd2fl'
         call trnfln( 1, fname )
c
c - note sifo2f and sifc2f are called in wmcd2f
c
        call wmcd2f(nlist,mcd2fl,
     &    core(cpt(38)),core(cpt(39)),core(cpt(40)),nact,
     &    invx(ix0(2)+1),symx(ix0(2)+1),mult,ira1,
     &    ira2,ndot,invx(ix0(1)+1),symx(ix0(1)+1),nsm,
     &    core(cpt(21)),core(cpt(22)),fname,1d0)
         close(unit=mcd2fl)

c
         call timer('writing the mc density files required',4,
     &    event3,nlist)
      endif
      if(flags(30) )then
      write(nlist,*)'Computing the requested mcscf (transition)
     & density matrices (flag 30)'
      write(nlist,*)'Reading mcdenin ...'
      call rdtmin
      write(nlist,*)'Number of density matrices (ndens):',ndens
      if(ndens.eq.0)then
          call bummer('flag(30) but no density specified',0,wrnerr)
          goto 1053
      endif
      write(nlist,*)'Number of unique bra states (ndbra):',ndbra
c  compute the transition density matrices
c
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:sd1s,39:ad1s,40:sd2s,41:tden,42:atden
c
c sd1s
      cpt(39)=cpt(38)+atebyt(nnta*ndens)
c ad1s -> actually smaller (zero diagonal now explicitely stored)
      cpt(40)=cpt(39)+atebyt(nnta*ndens)
c sd2s
      cpt(41)=cpt(40)+atebyt(numint(1)*ndens)
c tden
      cpt(42)=cpt(41)+atebyt(ncsf_f(dst)*ndbra*3)
c atden: the antisymmetric transition density and density are computed
c   in any case to facilitate the addressing
c   this could be changed but
c   the overhead for atden is only a third of tden
c   and only 1-particle density matrices are computed
      cpt(43)=cpt(42)+atebyt(ncsf_f(dst)*ndbra)
      if(cpt(43)-1 .gt. lcore)then
        write(nlist,6260)'rdft_grd call:',(cpt(i),i=1,43)
        call bummer('mcscf: allocation error. lcore=',lcore,faterr)
      endif
c       write(nlist,*)'starting rdft_grd; lenbft, lcore',lenbft,lcore
c
         call rdft_grd(core(cpt(43)),lenbft,
     & core(cpt(14)),core(cpt(7)),
     & ncsf_f(dst),core(cpt(1)),core(cpt(2)),core(cpt(3)),core(cpt(4)),
     & core(cpt(5)),core(cpt(6)),
     & core(cpt(8)),
     & navst(dst),core(cpt(38)),
     & core(cpt(39)),core(cpt(40)),core(cpt(41)),core(cpt(42)))
c
6249  format('mcsd',I1,'fl.drt',I1,'.st',I2.2)
6250  format('mc',A1,'d',I1,'fl.drt',I1,'.st',I2.2,'-st',I2.2)
6251  format('symm. mcscf dens. mat. (mcscf.x), DRT ',I1,', state '
     &,I2.2)
6252  format('symm. mcscf trans. dens. mat. (mcscf.x), DRT ',I1,
     &',state ',I2.2,' - DRT ',I1,', state ',I2.2)
6253  format('antisymm. mcscf trans. dens. mat. (mcscf.x), DRT ',I1,
     &',state ',I2.2,' - DRT ',I1,', state ',I2.2)
c
       do nxx=1,ndens
cfp: pointers to the current density matrix
         sdpt1=cpt(38) + atebyt(nnta*(nxx-1))
         sdpt2=cpt(40) + atebyt(numint(1)*(nxx-1))
         adpt1=cpt(39) + atebyt(nnta*(nxx-1))
c
         if (densind(1,nxx).eq.densind(2,nxx))then
            deltaab = 1d0
         else
            deltaab = 0d0
         endif
c
c - note sifo2f and sifc2f are called in wmcd2f
c
c write the symmetric and antisymmetric one-particle (transition) density matrix
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:sd1s,39:ad1s,40:sd2s,41:buffer,42:valbuf,43:labbuf,44:d1,45:mapmm,
c   46:mapim,47:mapar,48:mvres,49:bfnlab
c
         cpt(42)=cpt(41)+atebyt(lend1e)
         cpt(43)=cpt(42)+atebyt(n1edbf)
         cpt(44)=cpt(43)+forbyt(2*n1edbf)
         cpt(45)=cpt(44)+atebyt(nntm)
         cpt(46)=cpt(45)+atebyt(nmot)
         cpt(47)=cpt(46)+atebyt(nmot)
         cpt(48)=cpt(47)+forbyt(nmot)
         cpt(49)=cpt(48)+forbyt(nmot)
         cpt(50)=cpt(49)+forbyt(nmot)
         if(cpt(50)-1 .gt. lcore)then
            write(nlist,6260)'wmcd1f_grd call:',(cpt(i),i=1,50)
            call bummer('mcscf: allocation error. lcore=',lcore,faterr)
         endif
c
c
c -> write antisymm. density to the same file?
c symmetric density
c
         if (densind(1,nxx).eq.densind(2,nxx))then
            write(fname,6249)1,dst,densind(1,nxx)
            write(title(1),6251)dst,densind(1,nxx)
         else
            write(fname,6250)'s',1,dst,densind(1,nxx),densind(2,nxx)
            write(title(1),6252)dst,densind(1,nxx),dst,densind(2,nxx)
         endif
c
         call trnfln( 1, fname )
c         write(nlist,*)'fname d1',fname
         open(unit=mcd1fl,file=fname,status='unknown',
     &    form='unformatted')
c
         call wmcd1f_grd(nlist,mcd1fl,heig(dst,densind(1,nxx)),
     &    heig(dst,densind(2,nxx)),repnuc,
     &    orbtyp,
     &    core(cpt(41)),core(cpt(42)),core(cpt(43)),1,title,
     &    nsym,nmot,nmpsy,ndpsy,napsy,nvpsy,
     &    emc,invx(ix0(1)+1),invx(ix0(2)+1),invx(ix0(3)+1),
     &    ird1,ird2,ira1,ira2,irv1,irv2,irm1,irm2,
     &    nsa,n2sa,nsm,nnsm,nntm,nsd,n2sd,nsv,n2sv,fcimsk,
     &    nndx,core(sdpt1),
     &    core(cpt(44)),core(cpt(45)),core(cpt(46)),core(cpt(47)),
     &    core(cpt(48)),core(cpt(49)),deltaab,.True.)
          close(unit=mcd1fl)          
c
c antisymmetric density
       if(densind(1,nxx).ne.densind(2,nxx))then
         write(fname,6250)'a',1,dst,densind(1,nxx),densind(2,nxx)
         write(title(1),6253)dst,densind(1,nxx),dst,densind(2,nxx)
c
         call trnfln( 1, fname )
c         write(nlist,*)'fname d1',fname
         open(unit=mcd1fl,file=fname,status='unknown',
     &    form='unformatted')
c
         call wmcd1f_grd(nlist,mcd1fl,heig(dst,densind(1,nxx)),
     &    heig(dst,densind(2,nxx)),repnuc,
     &    orbtyp,
     &    core(cpt(41)),core(cpt(42)),core(cpt(43)),1,title,
     &    nsym,nmot,nmpsy,ndpsy,napsy,nvpsy,
     &    emc,invx(ix0(1)+1),invx(ix0(2)+1),invx(ix0(3)+1),
     &    ird1,ird2,ira1,ira2,irv1,irv2,irm1,irm2,
     &    nsa,n2sa,nsm,nnsm,nntm,nsd,n2sd,nsv,n2sv,fcimsk,
     &    nndx,core(adpt1),
     &    core(cpt(44)),core(cpt(45)),core(cpt(46)),core(cpt(47)),
     &    core(cpt(48)),core(cpt(49)),deltaab,.False.)
          close(unit=mcd1fl)          
       endif
c
c  write symmetric two-particle (transition) density matrices
c
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:sd1s,39:ad1s,40:sd2s,41:buffer,42:valbuf,43:labbuf
c
         cpt(42)=cpt(41)+atebyt(lend2e)
         cpt(43)=cpt(42)+atebyt(n2edbf)
         cpt(44)=cpt(43)+forbyt(4*n2edbf)
         if(cpt(44)-1 .gt. lcore)then
            write(nlist,6260)'wmcd2f call:',(cpt(i),i=1,44)
            call bummer('mcscf: allocation error. lcore=',lcore,faterr)
         endif
c
         if (densind(1,nxx).eq.densind(2,nxx))then
            write(fname,6249)2,dst,densind(1,nxx)
         else
            write(fname,6250)'s',2,dst,densind(1,nxx),densind(2,nxx)
         endif
         call trnfln( 1, fname )
c        write(nlist,*)'fname d2',fname
c
          call wmcd2f(nlist,mcd2fl,
     &    core(cpt(41)),core(cpt(42)),core(cpt(43)),nact,
     &    invx(ix0(2)+1),symx(ix0(2)+1),mult,ira1,
     &    ira2,ndot,invx(ix0(1)+1),symx(ix0(1)+1),nsm,
     &    core(sdpt1),core(sdpt2),fname,deltaab)
         close(unit=mcd2fl)
c         call timer('writing the mc transition-density
c     &     files required',4,event3,nlist)
       enddo
c
c  compute state specific NOs
c
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:sd1s,39:ad1s,40:cnos,41:occs,42:sa,43:sx
c
      cpt(41)=cpt(40)+ndens*atebyt(nbmt)
      cpt(42)=cpt(41)+ndens*atebyt(nmot)
      cpt(43)=cpt(42)+atebyt(n2tm)
      cpt(44)=cpt(43)+atebyt(nmot)
c
      if(cpt(44)-1 .gt. lcore)then
         write(nlist,6260)'state spec. NOs',(cpt(i),i=1,44)
         call bummer('mcscf: allocation error. lcore=',lcore,faterr)
      endif
c
6263  format('nocoef_mc.drt',I1,'.st',I2.2)
6264  format('nocoef_mc',A1,'.drt',I1,'.st',I2.2,'-st',I2.2)
6265  format('DRT ',I1,', state ',I2.2)
6266  format('DRT ',I1,',state ',I2.2,' - DRT ',I1,', state ',I2.2)
6272  format(/10x,'state spec. NOs: DRT',I2,', State',I3)
6273  format(/10x,'NOs of ',a,' trans. dens. mat.: DRT',I2,', State',I3,
     &  ' - DRT ',I2,', state ',I3)
1114  format ('state energy =',f18.9)
c     
      do nxx=1,ndens
        cnospt = cpt(40) + atebyt(nbmt*(nxx-1))
        occspt = cpt(41) + atebyt(nmot*(nxx-1))
        sdpt1  = cpt(38) + atebyt(nnta*(nxx-1))
c
        if (densind(1,nxx).ne.densind(2,nxx)) cycle
cfp: Compute NTOs in this case?

        ntitle=3
        write(fname,6263)dst,densind(1,nxx)
        title(1)='state specific NOs'
        write(title(2),6265)dst,densind(1,nxx)
        write(title(3),1114)heig(dst,densind(1,nxx))
        deltaab = 1d0
        write(nlist,6272)dst,densind(1,nxx)
c
         open(unit=mocout,file=fname,status='unknown',
     &    form='formatted')
      call stspnos(core(cpt(37)),core(sdpt1),core(cnospt),
     &  core(occspt),title,ntitle,core(cpt(42)),core(cpt(43)),deltaab)
         close(unit=mocout)
c
       if(flags(17))then
c        ...print the state specific natural orbitals and occupations.
         impt=cpt(40) + atebyt(nbmt*(nxx-1))
         ivpt=cpt(41) + atebyt(nmot*(nxx-1))
         icount=1
c write(fname,6250)'s',2,dst,densind(1,nxx),densind(2,nxx)

         do 1052 isym=1,nsym
            if(nmpsy(isym).ne.0)then
c
            write(nlist,'(/10x,a,i3)')
     &          'block',isym
c
              if ((nbft.le.orbmxpri).or.flags(32))then
                call prvbkc(' ',core(impt),core(ivpt),
     &          nbpsy(isym),nbpsy(isym),nmpsy(isym),
     &          bfnlabao(icount),bfnlabmo,' occ(*)=',1,nlist)
              else
                jlast = 0
                do jstrt = 1, nmpsy(isym), 8
                  jlast = min( nmpsy(isym), jlast+8 )
         write(nlist,1120) ( bfnlabmo(j), j = jstrt, jlast )
         write(nlist,1121) ' occ(*)=', (core(ivpt+j-1), j=jstrt,jlast)
c         write(nlist,*)
                enddo
              endif
              impt=impt+nbpsy(isym)*nmpsy(isym)
              ivpt=ivpt+nmpsy(isym)
              icount = icount + nbpsy(isym)
            endif
1052     continue
c endif flags(17)
      endif
c do nxx=1,ndens
      enddo

c
c endif flags(30)
1053  continue
      endif
c
c   1:xbar,2:xp,3:z,4:r,5:ind,6:modrt,7:iorder,8:add,9:h1ao,10:sao,
c   11:c,12:h1mo,13:lastb,14:hvec,15:udd,16:uad,17:uaa,18:uvd,
c   19:uva,20:uvv,21:d1,22:d2,23:fdd,24:fad,25:faa,26:fvd,27:fva,
c   28:qad,29:qaa,30:qva,31:qvv,32:w,33:k,34:cr,35:cno,36:occ,
c   37:cl,38:sd1s,39:ad1s,40:cnos,41:occs,42:cl,43:mnl,44:ms,45:scr1,46:scr2

c   37:cl,38:mnl,39:ms,40:scr1,41:scr2
c
c  ##  calculate mulliken population analysis
c
      cpt(44) = cpt(43) + forbyt(nbft)
      cpt(45) = cpt(44) + forbyt(nbft)
       if(cpt(44)-1 .gt. lcore)then
        write(nlist,6260)'mulpop_map call:',(cpt(i),i=1,44)
        call bummer('mcscf: allocation error. lcore=',lcore,faterr)
       endif
c
      call mulpop_map(
     & bfnlabao,        imtype,    map,     core(cpt(43)),
     & core(cpt(44)),   nmap,      nd_at,   mtype,
     & nbft,            nsym,      nlist, flags(26).or.flags(27))
c
      do i=1,8
      nipsy(i)=ndpsy(i)+napsy(i)
      enddo
c
      nipsy_mx=nipsy(1)
      do isym=2,nsym
      if (nipsy(isym).gt.nipsy_mx) nipsy_mx=nipsy(isym)
      enddo
      cpt(46)=cpt(45) + atebyt(9*nd_at)
      cpt(47)=cpt(46) + atebyt(9*nd_at*nipsy_mx)
       if(cpt(47)-1 .gt. lcore)then
        write(nlist,6260)'mulpop call:',(cpt(i),i=1,47)
        call bummer('mcscf: allocation error. lcore=',lcore,faterr)
       endif
c
      call mulpop(
     & core(cpt(36)), core(cpt(35)), core(cpt(45)), slabelao,
     & core(cpt(43)), core(cpt(44)), mtype,         nmpsy,
     & nd_at,         core(cpt(46)), core(cpt(10)), nlist,
     & nbft,          nmpsy,         nsym)
c
cfp: Mulliken populations for individual states
c
      if (flags(30)) then
        do nxx=1,ndens
          cnospt = cpt(40) + atebyt(nbmt*(nxx-1))
          occspt = cpt(41) + atebyt(nmot*(nxx-1))
          if (densind(1,nxx).eq.densind(2,nxx))then
            write(nlist,*)'Mulliken population for:'
            write(nlist,6265)dst,densind(1,nxx)
          else
            write(nlist,*)'Off-diagonal Mulliken population for:'
            write(nlist,6266)dst,densind(1,nxx),dst,densind(2,nxx)
          endif


          call mulpop(
     & core(occspt), core(cnospt), core(cpt(45)), slabelao,
     & core(cpt(43)), core(cpt(44)), mtype,         nmpsy,
     & nd_at,         core(cpt(46)), core(cpt(10)), nlist,
     & nbft,          nmpsy,         nsym)

c do nxx=1,ndens
        enddo
      endif
c
c  close and delete the scratch files...
c
*@ifndef noclose
      inquire(unit=scrtfl,opened=op,iostat=itm)
      if (op .and. itm.eq.0)
     .close(unit=scrtfl,status='delete')
      inquire(unit=ntemp,opened=op,iostat=itm)
      if (op .and. itm.eq.0)
     .close(unit=ntemp)
      inquire(unit=stape,opened=op,iostat=itm)
      if (op .and. itm.eq.0)
     .close(unit=stape,status='delete')
      inquire(unit=htape,opened=op,iostat=itm)
      if (op .and. itm.eq.0) then
      if (flags(22)) then
      close(unit=htape,status='delete')
      else
      close(unit=htape)
      endif
      endif
cmd
      if ((.not.flags(11)).and.flags(12)) then
      inquire(unit=htape,opened=op,iostat=itm)
      if (op .and. itm.eq.0)
     & close(unit=hdiagf)
      endif
*@endif
c
      call timer('mcscf',7,event1,nlist)
      return
c
6210  format(/1x,a)
6220  format(' iter= ',i4,' emc=',f18.10,' demc=',1pe11.4,
     &  ' wnorm=',e11.4,' knorm=',e11.4,' apxde=',e11.4,4x,a)
6225  format(i5,f18.10,1pe11.3,e11.3,e11.3,e11.3,2x,l1,3x,a)
6231  format(' iter     emc (average)    demc       wnorm      knorm
     &      apxde  qcoupl')
6230  format(/' orbital-state coupling will',a,
     &  'be calculated this iteration.')
6260  format(/' uexps allocation error ',a,' cpt(*)='/(1x,10i8))
c
      end
c
c****************************************************************
c
      block data
c
      implicit integer(a-z)
c     avepar.h
c
c    #  maximal No. of different DRTs
      integer maxnst
      parameter (maxnst=8)
c
c    #  maximal No. of states in one DRT
      integer mxavst
      parameter (mxavst = 50)
c
c    #  maximal total number of states
      integer maxstat
      parameter (maxstat = maxnst*mxavst)
c
c
      common/csymb/nxy(8,42),mult(8,8),nsym,nxtot(24)
c
      common/ccone/perm1(2,2),ind1(2),inds1,npass1
c
      common/cctwo/perm(4,4),ind(4),inds(4),npassp
c
      common/chbw/hblkw(15*maxstat),nhbw,htape
c
      common/cbufs/stape,lenbfs,h2bpt
c
c
      common/da2/nunit2,len2,irec2
c
c
      common/clda/ldamin,ldamax,ldainc
c
      integer  lena1e, n1eabf,lenbft,
     &        lend2e,n2edbf,lend1e,n1edbf,ifmt1,ifmt2
      common/cfinfo/lena1e,n1eabf,lenbft,
     &  lend2e,n2edbf,lend1e,n1edbf,ifmt1,ifmt2

c
c  standard 8-by-8 symmetry multiplication table.
c
      data mult/
     & 1,2,3,4,5,6,7,8,
     & 2,1,4,3,6,5,8,7,
     & 3,4,1,2,7,8,5,6,
     & 4,3,2,1,8,7,6,5,
     & 5,6,7,8,1,2,3,4,
     & 6,5,8,7,2,1,4,3,
     & 7,8,5,6,3,4,1,2,
     & 8,7,6,5,4,3,2,1/
c
c  permutation array for 1-particle density matrix indices.
c
      data perm1/1,2, 2,1/
c
c  permutation array for 2-particle density matrix indices.
c
      data perm/1,2,3,4, 2,1,3,4, 3,4,1,2, 4,3,1,2/
c
c  file unit numbers.
c
c
c  buffer length information:
c  lend1e=length of buffers for the mcscf d1 file
c  n1edbf=number of d1 elements in each buffer
c  lend2e=length of buffers for the mcscf d2 file
c  n2edbf=number of d2 elements in each buffer
c  lena1e=length of buffers for basis function 1-e integrals.
c  n1eabf=number of basis function 1-e integrals in each buffer.
c  lenbft=length of the formula tape buffers.
c  lenbfs=length of sorted integral buffers.
c
ctm
c default values (max. 255 basis functions)
      data lend1e/4096/, n1edbf/3272/
      data lend2e/4096/, n2edbf/2730/
      data lena1e/4096/, n1eabf/3272/
ctm
      data lenbft/2047/
      data lenbfs/65537/
c
c  ldamin=minimum da record length for integral sorting.
c  ldamax=maximum da record length.
c  ldainc=size increment for da record length
c
      data ldamin/127/, ldamax/4095/, ldainc/64/
c
      end

      subroutine driver_wrap(core,icore,lcore,nlist,nmot,
     .                       nexo,nsym,nbpsy,sewd,dlt2,lum,c,sao,
     .                   qcoupla,qsorta,qconvga,lcoremain,flag31a,
     .                   lastb,add,h1,bdg2e,lenbuk,nipbuk,ecore0)

      use tranlibrary
      implicit none
      real*8 core(*),c(*),sao(*) 
      integer icore(*),lcore,nlist,nmot,nexo(*),nsym
      integer nbpsy(*),sewd,dlt2,lum
      integer i,me,lcoremain,dimc,dimsao
      logical qcoupl,qsort,qconvg,flag31
      logical qcoupla,qsorta,qconvga,flag31a
      integer ga_nodeid,dtoi,itod
      integer lastb(*),add(*), lenbuk,nipbuk
      real*8  h1(*),ecore0,bdg2e(*)

c
c YifanShenSZ modification 1
      logical:: flag_YifanShenSZ
c end modification 1
c

      integer nflag
      parameter (nflag=50)
      logical flags
      common/lflags/flags(nflag)
*@ifdef int64
       INTEGER, PARAMETER :: INTLEN=1
*@else
*      INTEGER, PARAMETER :: INTLEN=2
*@endif
        dtoi(i) = intlen*(i-1) + 1
        itod(i) = (i-1)/intlen + 1



       qcoupl=qcoupla
       qconvg=qconvga
       qsort=qsorta
       flag31=flag31a
*@ifdef parallel
*       if (ga_nodeid().gt.0 .or. qconvg) goto 900
*@else 
       if (qconvg) return
*@endif
c
c     # prnopt == option%TRANOP(1)
c     #        = 0 : minimum print.
c     #        = 1 : orbital and transformation block information.
c     #        = 2 : print mo coefficients.
c     #        = 3 : indexing arrays and one-electron integral arrays.
c     #        = 4 : print final transformed integrals.
      if (flags(32))then
        options%tranop(1)=2           
      else
        options%tranop(1)=1
      endif 

c
c     # chkopt == option%TRANOP(2)
c     #        > 0 : program stops after determining which integral
c     #              blocks will be transformed in-core and out-of-core.
      options%tranop(2)=0           
c
c     # ortopt == option%TRANOP(3)
c     #         = 0 : stop if orbitals are not orthonormal.
c     #         = 1 : continue if orbitals are not orthonormal.
c     #         = 2 : schmidt orthonormalize the orbitals.
      options%tranop(3)=0           
c
c     # denopt == option%TRANOP(6)
c     #        = 0 : program expects to transform 1 and 2 -electron
c     #              integrals.
c     #        = 1 : program expects to transform 1 and 2 -particle
c     #              density matrices.
      options%tranop(6)=0           
c
c     # outlab == option%TRANOP(7)
c     #        = 0 : default output orbital labels.
c     #              for now, use tout:NNN, eventually use labels read
c     #              from the mocoef file.
c     #        = 1 : slab:NNN output labels. slab=clabs%SLABEL(isym) and
c     #              NNN is the symmetry-reduced orbital number.
c     #        = 2 : tout:NNN output labels. NNN is the global
c     #              output orbital number.
      options%tranop(7)=0           
c     #        NGRAD  (seward interface)
      options%tranop(8)=0           
c     #  force out of core 
*@ifdef parallel
*      options%tranop(9)=1
*@else
      options%tranop(9)=0
*@endif
c
c     # default threshhold value for small integrals.
      acracy%thresh=5.0d-12

c YifanShenSZ modification 2
      inquire(file='mcscfin_YifanShenSZ',exist=flag_YifanShenSZ)
      if(flag_YifanShenSZ) then
          open(unit=99,file='mcscfin_YifanShenSZ',status='old')
              read(99,*)acracy%thresh
          close(99)
      end if
c end modification 2

c
c     # default all integrals transformed
c     # possible values 0 : 0-ext integrals, only
c     #                 1 : 0,1-ext integrals 
c     #                 2 : 0,1,2-ext integrals 
c     #                 3 : 0,1,2,3-ext integrals 
c     #                 4 : 0,1,2,3,4-ext integrals 
       inputc%nextint=2
c
c     # output record length parameters.
c     # -1 = use sifs default values as defined by sifcfg().
c     #  0 = copy the input file parameters.
c     #  n = use the value n subject to any additional restrictions
c     #      imposed by sifcfg().
      clrc%lrc1mx=-1        
      clrc%lrc2mx=-1          

c
c     # sequential scratch file record length.
      clrc%lrcscr=65000
c
c     # initialize the entire cmapmm%MAPIN(*), cmapmm%FREEZE(*), and cmapmm%MAPOUT(*) arrays.
c     # there is no drt information on frozen core /frozen virtual in the mcdrtfl
c     # mo-to-level is constrained to active orbitals plus implicit ordering 
c     #           docc-active-virt 
      do 10 i = 1, nmot
         cmapmm%mapin(i)  = i
         cmapmm%mapout(i) = i
         cmapmm%freeze(i) = i
10    continue
c
c     # initialize inputc%NSYMAO, inputc%NAOPSY, and cmapmm%MAPIN(1) so that the defaults
c     # will be assigned later.
      inputc%nsymao=nsym
      inputc%naopsy(1:nsym)=nbpsy(1:nsym)
c
c     # initialize inputc%NSYMMO, inputc%NMOPSY(1), and cmapmm%MAPOUT(1), so that the
c     # defaults will be assigned later.
      inputc%nsymmo=-1  
c     inputc%nmopsy(1:nsym)=nbpsy(1:nsym)
      inputc%nmopsy(1)=-1
      csymb%NMPSY(1:nsym) = nbpsy(1:nsym) 
c     # set the default inputc%FSPLIT so that 1 output file will be created.
      inputc%fsplit=1     
c
c     # set the maximum space to use for in-core and out-of-core
c     # transformations.  these are otherwise-undocumented debugging
c     # variables. -rls
      inputc%szsgmx=2**35  
      inputc%incrmx=2**35 
      inputc%outcmx=2**35 
c     # default values for DA record length
c      cinfda%ldamax=16384
c
c     Ist ausser fuer C1 keine gute Wahl
cfp   adjust ldamax to the basis set size
      cinfda%ldainc=64     
      cinfda%ldamax=max(nmot*nmot/2+cinfda%ldainc+1,65000)

      cinfda%ldamin=127    
      csymb%nexo(1:nsym)=nexo(1:nsym)
      molcas%seward=sewd   
      molcas%lumorb=lum    
      molcas%dalton2=dlt2    
c     use vdisk for sorted integral file 
*@ifdef parallel
*       options%fsort=2 
*c     if (flag31) options%fsort=2 
*@else
      if (flag31) options%fsort=1 
*@endif
c
c
c     # adjust cinfda%LDAMAX so that it satisfies the equation
c     #    cinfda%LDAMAX = cinfda%LDAMIN + m * cinfda%LDAINC
c     # for some integer m, and limited by the input value of cinfda%LDAMAX.
c     # this fixes a logic error in subroutine seg() that can occur
c     # for large basis sets. 20-feb-02 -rls
c
*@ifndef parallel 
      cinfda%LDAMAX = cinfda%LDAMIN + ((cinfda%LDAMAX - cinfda%LDAMIN)
     .                  /cinfda%LDAINC) * cinfda%LDAINC
*@endif 

c
c     # write out the variables for this job.
      if (options%tranop(1).gt.0) then
      write(nlist,6030) options%tranop(1:3),options%tranop(6),
     &  cmapmm%MAPIN(1), inputc%NSYMAO, inputc%NAOPSY(1), 
     .  cmapmm%FREEZE(1),  
     &  cmapmm%MAPOUT(1), inputc%NSYMMO, inputc%NMOPSY(1), 
     .  inputc%FSPLIT, options%tranop(7),molcas%SEWARD,  
     &  molcas%LUMORB,molcas%DALTON2,inputc%nextint
6030  format(/' module tranlib input parameters:'// 
     &  ' prnopt    =',i6,', chkopt    =',i6,
     .  ',ortopt    =',i6, ', denopt    =',i6/ 
     &  ' mapin(1 ) =',i6,', nsymao    =',i6,
     .  ', naopsy(1) =',i6, ', freeze(1) =',i6/ 
     &  ' mapout(1) =',i6,', nsymmo    =',i6,
     .  ', nmopsy(1) =',i6, ', fsplit    =',i6/ 
     &  ' outlab    =',i6,', seward    =',i6,
     .  ', lumorb    =',i6, ', DALTON2   =',i6/ 
     &  ' nextint   =',i6)
c
      write(nlist,6040) cinfda%LDAMIN, cinfda%LDAMAX, cinfda%LDAINC
6040  format(' LDAMIN    =',i6,', LDAMAX    =',i6,', LDAINC    =',i6)
c
      write(nlist,6050) clrc%LRC1MX, clrc%LRC2MX, clrc%LRCSCR
6050  format(' LRC1MX    =',i6,', LRC2MX    =',i6,', LRCSCR    =',i6)
c
      write(nlist,6060) acracy%THRESH
6060  format(/' THRESH    =',1pe12.4,'  [cutoff threshold]')
      endif 
 900  dimc=0
      dimsao=0
      do i=1,nsym
        dimc=dimc+nbpsy(i)*nbpsy(i)
        dimsao=dimsao+nbpsy(i)*(nbpsy(i)+1)/2
      enddo
        
*@ifdef parallel
*      continue
*      me=ga_nodeid()
*      call brdcast1 (qconvg,qsort,dimc,dimsao)
*      if (me.eq.0) then
*        call brdcast2 (c,dimc)
*      else 
*        call brdcast2 (core,dimc)
*        lcore=lcore-dimc 
*      endif
*      call ga_sync()
*      if (qconvg) return
**
*c     if (qconvg)  then 
*c        call daclrd(2)
*c        return
*c     endif
*c
*c     broadcast all data needed to initialize driver_tran
*c     from process 0 to all other processes including qconv 
*c
*@else 
       me=0 
*@endif
      if (me.eq.0) then 
      call driver_tran( core,icore, lcore, nlist,c,qsort,me,lcoremain)
      else
      call driver_tran( core(dimc+1),icore(dtoi(dimc+1)), lcore, 6,
     .                 core,qsort, me,lcoremain)
      endif
*@ifdef parallel
*      if (ga_nodeid().gt.0 .and. (.not. qconvg) ) goto 900
*@endif
      return
      end 
*
      integer function getfreeunit()
c
c     # dynamically allocate unit numbers
c
      implicit none
      integer iunit,i
      logical op
      iunit=-1
      do i=7,99
       inquire(unit=i,opened=op)
       if (.not. op) then
          iunit=i
          exit
       endif
      enddo
      if (iunit.eq.-1) then
       call bummer('getfreeunit: all units in use',99,2)
      endif
       getfreeunit=iunit
      return
      end 
*
        subroutine drvk2(l1,l2,l3)
         return
         end subroutine drvk2
*
         subroutine cmpctr(l1,l2,l3)
         return
         end subroutine cmpctr
*
         subroutine cmpcts(l1,l2,l3)
         return
         end subroutine cmpcts
*
         subroutine crapso(l1,l2,l3)
         return
         end subroutine crapso
         subroutine dorys(l3)
         return
         end subroutine dorys
*

      subroutine getlenbfs(lenbfsdef,numint,nlist)

c  perform the second-half bin sort of 2-e integral types 6-11.
c
c  written by ron shepard.
c
      implicit integer(a-z)
      integer numint(12),nlist 
c
      common/csymb/nxy(8,42),mult(8,8),nsym,nxtot(24)
      integer nvpsy(8)
      equivalence (nxy(1,17),nvpsy(1))
      equivalence (nxtot(7),ndot)
      equivalence (nxtot(10),nact)
      equivalence (nxtot(13),nvrt)
      integer tsymdd(8),tsmdd2(8),tsymad(8),tsymvd(8)
      integer tsymva(8),tsymvv(8),tsmvv2(8)
      equivalence (nxy(1,21),tsymdd(1))
      equivalence (nxy(1,22),tsmdd2(1))
      equivalence (nxy(1,23),tsymad(1))
      equivalence (nxy(1,25),tsymvd(1))
      equivalence (nxy(1,26),tsymva(1))
      equivalence (nxy(1,27),tsymvv(1))
      equivalence (nxy(1,28),tsmvv2(1))
c
c
      integer nmotx
      parameter (nmotx=1023)
c
      common/corbc/iorbx(nmotx,8),ix0(3)
      integer symx(nmotx)
      equivalence (iorbx(1,4),symx(1))
      equivalence (ix0(1),ix0d),(ix0(2),ix0a)
c
c  prepare for sort: initialize /cis2/
c  add0  = address offset for placing integrals in core.
c  ibpt  = beginning of current buffer space.
c  ic0   = running offset within core(*).
c  icmax = address of the last integral in core(*).
c  ic0mx = minimum of icmax and the address of last buffer
c          element in core(*), ic0mx=min(icmax,ibpt+lenbfs-1).
c
       size_stape=0
       lenbfs=0
      
c
c  type-6 integrals: off6(pq)+adddd(ij), off6(pq)+adddd2(j,i)
c
      if(numint(6).eq.0)go to 700
      do 620 p=1,nact
      psym=symx(ix0a+p)
      do 610 q=1,p
      pqsym=mult(symx(ix0a+q),psym)
      nterm=tsymdd(pqsym)
      lenbfs=max(lenbfs,nterm)
      nterm=tsmdd2(pqsym)
      lenbfs=max(lenbfs,nterm)
610   continue      
620   continue
700   continue
c
c  type-7 integrals: off7(pq)+addvv(ab)
c
      if(numint(7).eq.0)go to 800
      do 720 p=1,nact
      psym=symx(ix0a+p)
      do 710 q=1,p
      pqsym=mult(symx(ix0a+q),psym)
      nterm=tsymvv(pqsym)
      lenbfs=max(lenbfs,nterm)
710   continue            
720   continue
c
800   continue
c
c  type-8 integrals: off8(pq)+addvv2(b,a)
c
      if(numint(8).eq.0)go to 900
      do 820 p=1,nact
      psym=symx(ix0a+p)
      do 810 q=1,p
      pqsym=mult(symx(ix0a+q),psym)
      nterm=tsmvv2(pqsym)
      lenbfs=max(lenbfs,nterm)
810   continue       
820   continue
c
900   continue
c
c  type-9 integrals: off9(ai)+addva(b,p)
c
      if(numint(9).eq.0)go to 1000
      nai=tsymvd(1)
      do 910 ai=1,nai
      nterm=tsymva(1)
      lenbfs=max(lenbfs,nterm)
910   continue     
1000  continue
c
c  type-10 integrals: off10(ai)+addad(j,p)
c
c     if(bukpt(11).eq.0)go to 1100
      if(numint(10).eq.0)go to 1100
      nai=tsymvd(1)
      do 1010 ai=1,nai
      nterm=tsymad(1)
      lenbfs=max(lenbfs,nterm)
1010  continue        
1100  continue
c
c  type-11 integrals: off11(ij)+addvv3(b,a)
c
      if(numint(11).eq.0)go to 1200
      do 1120 i=1,ndot
      na=nvpsy(symx(ix0d+i))
      do 1110 j=1,i
      nb=nvpsy(symx(ix0d+j))
      nterm=na*nb
      lenbfs=max(lenbfs,nterm)
1110  continue      
1120  continue
1200  continue

      write(nlist,*) 'lenbfsdef=',lenbfsdef,' lenbfs=',lenbfs 
      if (lenbfsdef.lt.lenbfs) then
        write(nlist,*) 'default value of lenbfs too small '
        write(nlist,*) 'increasing lenbfs from',lenbfsdef,' to ', lenbfs
        lenbfsdef=lenbfs
      endif
      write(nlist,*) ' number of integrals per class 1:11 (cf adda '
1201  format(1x,a,i2,a,i12)
      write(nlist,1201) 'class ',1,' (pq|rs):         #',numint(1)
      write(nlist,1201) 'class ',2,' (pq|ri):         #',numint(2)
      write(nlist,1201) 'class ',3,' (pq|ia):         #',numint(3)
      write(nlist,1201) 'class ',4,' (pi|qa):         #',numint(4)
      write(nlist,1201) 'class ',5,' (pq|ra):         #',numint(5)
      write(nlist,1201) 'class ',6,' (pq|ij)/(pi|qj): #',numint(6)
      write(nlist,1201) 'class ',7,' (pq|ab):         #',numint(7)
      write(nlist,1201) 'class ',8,' (pa|qb):         #',numint(8)
      write(nlist,1201) 'class ',9,' p(bp,ai)         #',numint(9)
      write(nlist,1201) 'class ',10,'p(ai,jp):        #',numint(10)
      write(nlist,1201) 'class ',11,'p(ai,bj):        #',numint(11)
      return
      end
