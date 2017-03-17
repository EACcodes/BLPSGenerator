!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv
!! NAME
!! scfcv
!!
!! FUNCTION
!! Self-consistent-field convergence.
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! ground state and optionally to compute forces and energy.
!! This routine is called to compute forces for given atomic
!! positions or else to do non-SCF band structures.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, GMR, AR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points which can fit in memory;
!!   |    set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |    for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  fatvshift=factor to multiply dtset%atvshift
!!  iapp=indicates the eventual suffix to be appended to the generic
!!     output root
!!     if 0 : no suffix to be appended (called directly from gstate)
!!     if positive : append "_TIM//iapp" (called from move or brdmin)
!!     if -1 : append "_TIM0" (called from brdmin)
!!     if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)= # atoms of each type.
!!  ndtpawuj=size of dtpawuj
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and
!!     related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for
!!     each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real
!!     spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points
!!     and spins
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions; if mkmem>=nkpt, these are kept in a disk file.
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtpawuj(ndtpawuj)= data used for the automatic determination of U
!!     (relevant only for PAW+U) calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for
!!     the electron-positron annihilation
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet
!!     finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible
!!     zone data
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!     for the "fine" grid (see NOTES below)
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2)
!!     at each k point
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic
!!     occupancies
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic
!!     translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!     forces and its components, the stress tensor) of a ground-state
!!     computation (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  scf_history <type(scf_history_type)>=arrays obtained from previous
!!     SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic
!!     energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic
!!     coordinates at output, current xred is transferred to xred_old
!!
!! NOTES
!! It is worth to explain THE USE OF FFT GRIDS:
!! ============================================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      scfcv_new
!!
!! CHILDREN
!!      ab6_mixing_deallocate,ab6_mixing_new,ab6_mixing_use_disk_cache
!!      abi_etsf_init,afterscfloop,berryphase_new,check_kxc,chkdilatmx
!!      chkpawovlp,cprj_alloc,cprj_free,ctocprj,destroy_paw_an,destroy_paw_ij
!!      energies_init,energy,etotfor,expibr,extraprho,first_rec,fourdp,fresid
!!      getcut,getmpw,getng,getph,hdr_update,init_metricrec,init_paw_an
!!      init_paw_ij,initmpi_seq,initylmg,int2char4,ioarr,kpgio,leave_test
!!      metric,newrho,newvtr,nhatgrid,nullify_paw_an,nullify_paw_ij,odamix
!!      out_geometry_xml,out_resultsgs_xml,outscfcv,pawdenpot,pawdij
!!      pawfgrtab_free,pawfgrtab_init,pawmknhat,pawtwdij,pawuj_red,prc_mem_free
!!      prtefield,prtene,rhohxc,rhotov,scprqt,setnoccmmp,setrhoijpbe0,setsym
!!      setup_positron,setvtr,sphereboundary,status,symdij,symzat,timab,vtorho
!!      vtorhorec,vtorhotf,wrtout,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine scfcv(atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,&
&  dtset,ecore,eigen,electronpositron,fatvshift,hdr,iapp,indsym,&
&  initialized,irrzon,kg,mcg,mpi_enreg,my_natom,nattyp,ndtpawuj,nfftf,npwarr,occ,&
&  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,&
&  pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
 use m_xmpi
 use m_profiling
 use m_wffile
 use m_paw_toolbox
 use m_rec
 use m_ab6_mixing
 use m_errors
 use m_efield
 use mod_prc_memory
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use m_results_gs ,      only : results_gs_type
 use m_scf_history,      only : scf_history_type
 use m_energies,         only : energies_type, energies_init
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_paw_dmft,         only : paw_dmft_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_pawcprj,          only : cprj_type, cprj_alloc, cprj_free
 use m_header,           only : hdr_update

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfcv'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_27_toolbox_oop
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_57_iovars
 use interfaces_61_ionetcdf
 use interfaces_62_iowfdenpot
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_68_recursion
 use interfaces_68_rsprc
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => scfcv
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iapp,mcg,my_natom,ndtpawuj,pwind_alloc
 integer,intent(inout) :: initialized,nfftf
 real(dp),intent(in) :: cpus,ecore,fatvshift
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(inout) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(recursion_type),intent(inout) :: rec_set
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),intent(inout) :: scf_history
 type(wffile_type),intent(inout) :: wffnew,wffnow
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in) :: symrec(3,3,dtset%nsym)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), pointer :: taug(:,:),taur(:,:)
 real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: xred(3,dtset%natom)
 real(dp), intent(inout) :: xred_old(3,dtset%natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
 type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type), intent(inout) :: paw_dmft

!Local variables -------------------------
!scalars
 integer,parameter :: level=110,response=0
 integer :: accessfil,afford,bantot,choice
 integer :: computed_forces,cplex,ctocprj_choice,dbl_nnsclo,dielop,dielstrt,dimdmat
 integer :: fformr,forces_needed,errid
!integer :: dtset_iprcel
 integer :: iatom,ider,idir,ierr,iexit,ii,iidir,ikpt,impose_dmat,denpot
 integer :: initialized0,iorder_cprj,ipert,ipositron,ir,isave,iscf10,ispden
 integer :: ispmix,istep,istep_mix,itypat,izero,lmax_diel,lpawumax,mcprj,me
 integer :: mgfftdiel,mgfftf,moved_atm_inside,moved_rhor,my_nspinor,n1xccc
 integer :: n3xccc,ncpgr,nfftdiel,nfftmix,ngrvdw,nhatgrdim,nk3xc,nkxc
 integer :: npawmix,npwdiel,nstep,nzlmopt,optberry,optcut,optene,optgr0
 integer :: optgr1,optgr2,option,optrad,optres,optxc,prtden,prtkden,prtfor,prtxml,quit
 integer :: quit_sum,rdwr,rdwrpaw,spaceComm,spaceComm_fft,stress_needed,unit_out
 integer :: usecprj,usexcnhat,useylmgr,v_size
 real(dp) :: boxcut,compch_fft,compch_sph,deltae,diecut,diffor,ecut
 real(dp) :: ecutf,ecutsus,edum,elast,etotal,fermie,gsqcut
 real(dp) :: maxfor,res2,residm,ucvol,val_max
 real(dp) :: val_min,vxcavg,vxcavg_dum
 character(len=4) :: tag
 character(len=1500) :: message
 character(len=fnlen) :: fildata,kgnam
 type(MPI_type) :: mpi_enreg_diel
 type(energies_type) :: energies
 type(ab6_mixing_object) :: mix
 logical :: VERBOSE=.FALSE.
 logical :: recompute_cprj=.false.
!arrays
 integer :: ngfft(18),ngfftdiel(18),ngfftf(18),ngfftmix(18),npwarr_diel(1)
 integer :: npwtot_diel(1)
 integer,allocatable :: dimcprj(:),gbound_diel(:,:),irrzondiel(:,:,:),kg_diel(:,:)
 integer,allocatable :: lmn_size(:),l_size_atm(:)
 integer,allocatable :: indsym_dum(:,:,:),symrec_dum(:,:,:)
 real(dp) :: dielar(7),dphase(3),dummy2(6),favg(3),gmet(3,3),gprimd(3,3)
 real(dp) :: kpt_diel(3),pel(3),pel_cg(3),pelev(3),pion(3),ptot(3),red_ptot(3) !!REC
 real(dp) :: rhodum(1),rmet(3,3),strsxc(6),strten(6),tollist(12)
 real(dp) :: tsec(2),vnew_mean(dtset%nspden),vres_mean(dtset%nspden)   !!HONG
 real(dp) :: efield_old_cart(3),  efield_test_cart(3), ptot_cart(3)   !!HONG
 real(dp) :: red_efield1(3),red_efield2(3),red_efieldbar_lc(3),red_efield2_old(3)
! red_efield1(3),red_efield2(3) is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009)
! red_efield1(3) for fixed ebar calculation, red_efield2(3) for fixed reduced d calculation, in mixed BC
! red_efieldbar_lc(3) is local reduced electric field, defined by Eq.(28) of Nat. Phys. suppl. (2009)
! pbar(3) and dbar(3) are reduced polarization and displacement field, defined by Eq.(27) and (29) Nat. Phys. suppl. (2009)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable :: dielinv(:,:,:,:,:),dtn_pc(:,:)
 real(dp),allocatable :: fcart(:,:),forold(:,:),fred(:,:),gresid(:,:)
 real(dp),allocatable :: grewtn(:,:),grhf(:,:),grnl(:),grvdw(:,:),grxc(:,:)
 real(dp),allocatable :: kxc(:,:),nhat(:,:),nhatgr(:,:,:),nvresid(:,:)
 real(dp),allocatable :: ph1d(:,:),ph1ddiel(:,:),ph1df(:,:)
 real(dp),allocatable :: phnonsdiel(:,:,:),shiftvector(:)
 real(dp),allocatable :: susmat(:,:,:,:,:),synlgr(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),vxctau(:,:,:),workr(:,:),xccc3d(:),ylmdiel(:,:)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:)
 type(cprj_type),allocatable :: cprj(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)

 !==========chenh ===============================
 character(len=500) :: msg,structure
 integer, dimension(8) :: dateTime
 integer  :: i, j, k, myicount, c_nsym, maxcount,tmp_int
 integer  :: n1, n2, n3, eEvalcount_c,io
 logical  :: tigger_mode,ex
 real(dp) :: histW(4),net_value

 ! parameters to restart BFGS if max line search is reached
 logical  :: restart=.false.
 integer  :: nline,nBFGS,maxLine=5
 real(dp) :: last_x(dtset%nfft), &
             last_x2(dtset%nfft), &
             last_ke,last_fe, &
             pre_rhor(dtset%nfft), &
             pre_rhor2(dtset%nfft)

 ! mpi vars
 integer :: spaceworld, mpi_err

 ! control flags
 Integer :: bLoadDen,bLoadVr,bAllowSym,bSpin

 !BFGS...
 integer,parameter :: mmax=10

 character(len=60) :: task, BFGS_csave

 integer  ::   BFGS_nbd(dtset%nfft),& 
               BFGS_iwa(3*dtset%nfft), &
               BFGS_iwa2(6*dtset%nfft), &
               BFGS_isave(44), & 
               BFGS_print=-1
 logical  ::   BFGS_lsave(4)

 real(dp) ::   BFGS_l(dtset%nfft),&
               BFGS_l2(2*dtset%nfft),&
               BFGS_u(dtset%nfft), & 
               BFGS_u2(2*dtset%nfft), &
               BFGS_dsave(29), &
               BFGS_wa(2*mmax*dtset%nfft+4*dtset%nfft+12*mmax*mmax+12*mmax), &
               BFGS_wa2(4*mmax*dtset%nfft+8*dtset%nfft+12*mmax*mmax+12*mmax)

 real(dp) ::   cgW, & 
               cgG(dtset%nfft), &                ! gradient for total/spin up electron
               cgG2(dtset%nfft), &               ! gradient for spin down electrons in spin polarized case
               cgGCombine(2*dtset%nfft), &       ! double sized array for gradient for both spin up/down electrons
               ref_rhor(dtset%nfft), &
               ref_rhor2(dtset%nfft), &          ! the similar array for spin down electron, if needed in spinpol case
               My_vtrial(2*dtset%nfft), &        ! the vtrial array double-sized in case of spin polarized
               eold1_c,eold2_c,wtol,ektol, & 
               tmp1,tmp2,tmp3, & 
               vtmp(dtset%nfft,2)

 real(dp) :: ke_firstOrder

! penalty funciton related
 real(dp) :: penFun, penLambda, lapvfg(dtset%nfft), lapvfg2(dtset%nfft)

! MPI related 
 integer :: nproc


! *********************************************************************

!DEBUG
!write(std_out,*) '->scfcv: enter'
!write(std_out,*) '++++++++++++++'
!write(std_out,*) 'iapp=',iapp
!write(std_out,*) 'ndtpawuj=',ndtpawuj
!write(std_out,*) 'pwind_alloc=',pwind_alloc
!write(std_out,*) 'initialized=',initialized
!write(std_out,*) 'nfftf=',nfftf
!write(std_out,*) 'cpus=',cpus
!write(std_out,*) 'ecore=',ecore
!write(std_out,*) 'fatvshift=',fatvshift
!write(std_out,*) 'atindx=',atindx
!write(std_out,*) 'indsym=',indsym
!write(std_out,*) 'irrzon=',irrzon
!!  write(std_out,*) 'kg=',kg
!write(std_out,*) 'nattyp=',nattyp
!write(std_out,*) 'npwarr=',npwarr
!write(std_out,*) 'pwind=',pwind
!write(std_out,*) 'symrec=',symrec
!!  write(std_out,*) 'cg=',cg
!write(std_out,*) 'eigen=',eigen
!write(std_out,*) 'occ=', occ
!write(std_out,*) 'phnons=',phnons
!write(std_out,*) 'pwnsfac=',pwnsfac
!write(std_out,*) 'rprimd='
!do ii=1,3
!write(std_out,*) rprimd(:,ii)
!end do
!write(std_out,*) 'resid=',resid
!write(std_out,*) 'xred='
!do ii=1,dtset%natom
!write(std_out,*) xred(:,ii)
!end do
!write(std_out,*) 'xred_old='
!do ii=1,dtset%natom
!write(std_out,*) xred_old(:,ii)
!end do
!write(std_out,*) 'ylm=',ylm
!write(std_out,*) 'ylmgr=',ylmgr

 DBG_ENTER("COLL")

 call timab(238,1,tsec)
 call timab(54,1,tsec)

!Structured debugging if prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' scfcv : enter '
   call wrtout(std_out,message,'COLL')
 end if

!######################################################################
!Initializations - Memory allocations
!----------------------------------------------------------------------

 call status(0,dtfil%filstat,iexit,level,'allocate/init ')

!MPI communicators
 spaceComm=mpi_enreg%comm_cell
 me=xcomm_rank(spaceComm)
 spaceComm_fft=mpi_enreg%comm_fft

!Save some variables from dataset definition
 nstep=dtset%nstep
!dtset_iprcel = dtset%iprcel
 ecut=dtset%ecut
 ecutf=ecut
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 iscf10=mod(dtset%iscf,10)
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff
 dielstrt=0

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Some variables need to be initialized/nullified at start
 nullify(grhor,lrhor,elfr)
 quit=0 ; dbl_nnsclo=0 ;
 dielop=0 ; strsxc=zero
 deltae=zero ; elast=zero ;
 results_gs%residm=zero;results_gs%res2=zero
 results_gs%deltae=zero;results_gs%diffor=zero
 call energies_init(energies)
 if (dtset%positron/=0.and.initialized/=0) then
   energies%e0_electronpositron =results_gs%energies%e0_electronpositron
   energies%e_electronpositron  =results_gs%energies%e_electronpositron
   energies%edc_electronpositron=results_gs%energies%edc_electronpositron
   maxfor=zero
 end if
 if (dtset%nstep==0) energies%e_fermie=results_gs%energies%e_fermie
 energies%e_corepsp = ecore / ucvol
 fermie=energies%e_fermie
 isave=0 !initial index of density protection file
 optres=merge(0,1,dtset%iscf<10)
 usexcnhat=0;usecprj=0
 initialized0=initialized
 ipert=0;idir=0;cplex=1
 istep_mix=1
 ipositron=electronpositron_calctype(electronpositron)
 unit_out=0;if (dtset%prtvol >= 10) unit_out=ab_out

!Stresses and forces flags
 forces_needed=0;prtfor=0
 if ((dtset%optforces==1.or.dtset%ionmov==4.or.abs(tollist(3))>tiny(0._dp))) then
   if (dtset%iscf>0.and.nstep>0) forces_needed=1
   if (nstep==0) forces_needed=2
   prtfor=1
 else if (dtset%iscf>0.and.dtset%optforces==2) then
   forces_needed=2
 end if

 stress_needed=0
 if (dtset%optstress>0.and.dtset%iscf>0.and.dtset%prtstm==0.and. (nstep>0.or.dtfil%ireadwf==1)) stress_needed=1

!This is only needed for the tddft routine, and does not
!correspond to the intented use of results_gs (should be only
!for output of scfcv
 etotal  =results_gs%etotal

!Get FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
   mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
 else
   mgfftf=dtset%mgfft;ngfftf(:)=ngfft(:)
 end if

!We create Output files when required
 if (dtset%accesswff == IO_MODE_ETSF) then
#if defined HAVE_TRIO_ETSF_IO
!  Compute this lmn_size stuff
   ABI_ALLOCATE(lmn_size,(psps%npsp))
   if(psps%usepaw==1) then
     lmn_size(:) = pawtab(1:psps%npsp)%lmn_size
   else
     lmn_size(:) = psps%lmnmax
   end if

!  Create an ETSF file for each required files

   if (dtset%prtden /= 0) then ! Case of density.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_den, 1, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtelf /= 0) then ! Case of electron localization function
     call abi_etsf_init(dtset,dtfil%fnameabo_app_elf, 1, .true., lmn_size, psps, wvl%wfs)
     if (dtset%nspden==2) then
!      Case of spin-dependent electron localization function
       call abi_etsf_init(dtset,dtfil%fnameabo_app_elf_up, 1, .true., lmn_size, psps, wvl%wfs)
       call abi_etsf_init(dtset,dtfil%fnameabo_app_elf_down, 1, .true., lmn_size, psps, wvl%wfs)
     end if
   end if

   if (dtset%prtgden /= 0) then ! Case of gradient of electron density.
     call abi_etsf_init(dtset, dtfil%fnameabo_app_gden1, 1, .true., lmn_size, psps, wvl%wfs)
     call abi_etsf_init(dtset, dtfil%fnameabo_app_gden2, 1, .true., lmn_size, psps, wvl%wfs)
     call abi_etsf_init(dtset, dtfil%fnameabo_app_gden3, 1, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtkden /= 0) then ! Case of kinetic energy density.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_kden, 1, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtlden /= 0) then ! Case of Laplacian of electron density.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_lden, 1, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtwf == 1) then ! Case of wavefunctions.
     call abi_etsf_init(dtset, dtfil%fnameabo_app_wfk, 2, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtvxc /= 0) then ! Case of Exchange-correlation potential.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_vxc, 24, .true., lmn_size, psps, wvl%wfs)
   end if
!  FIXME: append other possibilities.
!  * ETSFIO cases of only correlation or only exchange are not abinit
!  options
!  * the fix for VHA,VHXC,POT,STM is dirty: they are flagged as
!  exchange correlation pot files, except stm, which is flagged density.
!  START dirty treatment

   if (dtset%prtvha /= 0) then ! Case of Hartree potential.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_vha, 24, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtvhxc /= 0) then ! Case of Hartree+XC potential.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_vhxc, 24, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtpot /= 0) then ! Case of total potential.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_pot, 24, .true., lmn_size, psps, wvl%wfs)
   end if

   if (dtset%prtstm /= 0) then ! Case of STM output.
     call abi_etsf_init(dtset,dtfil%fnameabo_app_stm, 1, .true., lmn_size, psps, wvl%wfs)
   end if
   ABI_DEALLOCATE(lmn_size)
#endif
 end if

!Entering a scfcv loop, printing data to XML file if required.
 prtxml=0;if (me==0.and.dtset%prtxml==1) prtxml=1
 if (prtxml == 1) then
!  scfcv() will handle a scf loop, so we output the scfcv markup.
   write(ab_xml_out, "(A)") '    <scfcvLoop>'
   write(ab_xml_out, "(A)") '      <initialConditions>'
!  We output the geometry of the dataset given in argument.
!  xred and rprimd are given independently since dtset only
!  stores original and final values.
   call out_geometry_XML(dtset, 4, dtset%natom, rprimd, xred)
   write(ab_xml_out, "(A)") '      </initialConditions>'
 end if

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor and res2 are here initialized to 0)
 choice=1 ; diffor=zero ; res2=zero
 ABI_ALLOCATE(fcart,(3,dtset%natom))
 ABI_ALLOCATE(fred,(3,dtset%natom))
 fred(:,:)=zero
 fcart(:,:)=results_gs%fcart(:,:) ! This is a side effect ...
!results_gs should not be used as input of scfcv
!HERE IS PRINTED THE FIRST LINE OF SCFCV

 call scprqt(choice,cpus,deltae,diffor,dtset,&
& eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
& dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
& maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
& occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
& psps%usepaw,vxcavg,dtset%wtk,xred)

!Various allocations (potentials, gradients, ...)
 ABI_ALLOCATE(forold,(3,dtset%natom))
 ABI_ALLOCATE(grnl,(3*dtset%natom))
 ABI_ALLOCATE(gresid,(3,dtset%natom))
 ABI_ALLOCATE(grewtn,(3,dtset%natom))
 ABI_ALLOCATE(grxc,(3,dtset%natom))
 ABI_ALLOCATE(synlgr,(3,dtset%natom))
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
 ABI_ALLOCATE(vhartr,(nfftf))
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vpsp,(nfftf))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vxctau,(nfftf,dtset%nspden*dtset%usekden,4))
 ngrvdw=0;if (dtset%vdw_xc==5) ngrvdw=dtset%natom
 ABI_ALLOCATE(grvdw,(3,ngrvdw))
 forold(:,:)=zero ; gresid(:,:)=zero ; pel(:)=zero
 vtrial(:,:)=zero; vxc(:,:)=zero
 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))

!Allocations/initializations for PAW only
 lpawumax=-1

 if(psps%usepaw==1) then
!  Variables/arrays related to the fine FFT grid
   ABI_ALLOCATE(nhat,(nfftf,dtset%nspden*psps%usepaw))
   if (nstep==0) nhat=zero
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(my_natom))
   if (my_natom>0) then
     ABI_ALLOCATE(l_size_atm,(my_natom))
     do iatom=1,my_natom
       itypat=pawrhoij(iatom)%itypat
       l_size_atm(iatom)=pawtab(itypat)%lcut_size
     end do
     call pawfgrtab_init(pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat,&
&     mpi_atmtab=mpi_enreg%my_atmtab,mpi_comm_atom=mpi_enreg%comm_atom)
     ABI_DEALLOCATE(l_size_atm)
   end if
   compch_fft=-1.d5
   usexcnhat=maxval(pawtab(:)%usexcnhat)
   if (usexcnhat==0.and.dtset%ionmov==4.and.dtset%iscf<10) then
     message = ' You cannot simultaneously use ionmov=4 and such a PAW psp file !'
     MSG_ERROR(message)
   end if

!  Variables/arrays related to the PAW spheres
   ABI_DATATYPE_ALLOCATE(paw_ij,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_an,(my_natom))
   call nullify_paw_an(paw_an)
   call nullify_paw_ij(paw_ij)
   option=0
   if (dtset%iscf==22) option=1
   call init_paw_an(dtset%natom,dtset%ntypat,0,dtset%nspden,cplex,&
&   dtset%pawxcdev,dtset%typat,pawang,pawtab,paw_an,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call init_paw_ij(paw_ij,cplex,dtset%nspinor,&
&   dtset%nsppol,dtset%nspden,dtset%pawspnorb,dtset%natom,dtset%ntypat,&
&   dtset%typat,pawtab,has_dij=1,has_dijso=1,has_dijhat=option,&
&   has_pawu_occ=1,has_exexch_pot=1,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   compch_sph=-1.d5
   ABI_ALLOCATE(dimcprj,(dtset%natom))
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
     dimcprj(iatom)=pawtab(itypat)%lmn_size
     if (pawtab(itypat)%usepawu>0) lpawumax=max(pawtab(itypat)%lpawu,lpawumax)
   end do

 end if ! PAW

!WVL - since wavelets change the size of the box, dont need dilatmax.
 if (dtset%usewvl == 0) then
!  Check that the possible change of unit cell size has not lead to a too large increase
   call chkdilatmx(dtset%dilatmx,rprimd,dtset%rprimd_orig(1:3,1:3,1))
 end if

!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (dtset%iscf>=0) then
   dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
   dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
   dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
   dielar(7)=dtset%diemix;if (dtset%iscf>=10) dielar(7)=dtset%diemixmag
   ABI_ALLOCATE(nvresid,(nfftf,dtset%nspden))
   if (nstep==0) nvresid=zero
   ABI_ALLOCATE(dtn_pc,(3,dtset%natom))
!  The next arrays are needed if iscf==5 and ionmov==4,
!  but for the time being, they are always allocated
   ABI_ALLOCATE(grhf,(3,dtset%natom))
!  Additional allocation for mixing within PAW
   npawmix=0
   if(psps%usepaw==1) then
     do iatom=1,my_natom
       itypat=pawrhoij(iatom)%itypat
       pawrhoij(iatom)%use_rhoijres=1
       sz1=pawrhoij(iatom)%cplex*pawtab(itypat)%lmn2_size
       sz2=pawrhoij(iatom)%nspden
       ABI_ALLOCATE(pawrhoij(iatom)%rhoijres,(sz1,sz2))
       do ispden=1,pawrhoij(iatom)%nspden
         pawrhoij(iatom)%rhoijres(:,ispden)=zero
       end do
       ABI_ALLOCATE(pawrhoij(iatom)%kpawmix,(pawtab(itypat)%lmnmix_sz))
       pawrhoij(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       pawrhoij(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+pawrhoij(iatom)%nspden*pawtab(itypat)%lmnmix_sz*pawrhoij(iatom)%cplex
     end do
   end if
   if (dtset%iscf > 0) then
     denpot = AB6_MIXING_POTENTIAL
     if (dtset%iscf > 10) denpot = AB6_MIXING_DENSITY
     if (psps%usepaw==1.and.dtset%pawmixdg==0) then
       ispmix=AB6_MIXING_FOURRIER_SPACE;nfftmix=dtset%nfft;ngfftmix(:)=ngfft(:)
     else
       ispmix=AB6_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
     end if
     call ab6_mixing_new(mix, iscf10, denpot, ispmix, &
&     nfftmix, dtset%nspden, npawmix, errid, message, dtset%npulayit)
     if (errid /= AB6_NO_ERROR) then
       MSG_ERROR(message)
     end if
     if (dtset%mffmem == 0) then
       call ab6_mixing_use_disk_cache(mix, dtfil%fnametmp_fft)
     end if
   end if
 end if ! iscf>0

!Here, allocate arrays for computation of susceptibility and dielectric matrix or for TDDFT
 if( (nstep>0 .and. dtset%iscf>=0) .or. dtset%iscf==-1 ) then !MF
!  The dielectric stuff is performed in sequential mode; set mpi_enreg_diel accordingly
   call initmpi_seq(mpi_enreg_diel)

!  Here, for TDDFT, artificially set iprcel . Also set a variable to reduce the memory needs.
   afford=1
   if(dtset%iscf==-1) then
!    dtset%iprcel=21
     afford=0
   end if

!  First compute dimensions
   if(dtset%iprcel>=21 .or. dtset%iscf==-1)then
!    With dielop=1, the matrices will be computed when istep=dielstrt
!    With dielop=2, the matrices will be computed when istep=dielstrt and 1
     dielop=1
     if(dtset%iprcel>=41)dielop=2
     if((dtset%iprcel >= 71).and.(dtset%iprcel<=79)) dielop=0 !RSkerker preconditioner do not need the susceptibility matrix
!    Immediate computation of dielectric matrix
     dielstrt=1
!    Or delayed computation
     if(modulo(dtset%iprcel,100)>21 .and. modulo(dtset%iprcel,100)<=29)dielstrt=modulo(dtset%iprcel,100)-20
     if(modulo(dtset%iprcel,100)>31 .and. modulo(dtset%iprcel,100)<=39)dielstrt=modulo(dtset%iprcel,100)-30
     if(modulo(dtset%iprcel,100)>41 .and. modulo(dtset%iprcel,100)<=49)dielstrt=modulo(dtset%iprcel,100)-40
     if(modulo(dtset%iprcel,100)>51 .and. modulo(dtset%iprcel,100)<=59)dielstrt=modulo(dtset%iprcel,100)-50
     if(modulo(dtset%iprcel,100)>61 .and. modulo(dtset%iprcel,100)<=69)dielstrt=modulo(dtset%iprcel,100)-60
!    Get diecut, and the fft grid to be used for the susceptibility computation
     diecut=abs(dtset%diecut)
     if( dtset%diecut<0.0_dp )then
       ecutsus=ecut
     else
       ecutsus= ( sqrt(ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
     end if
!    Impose sequential calculation
     ngfftdiel(1:3)=0 ; ngfftdiel(7)=100 ; ngfftdiel(9)=0; ngfftdiel(8)=dtset%ngfft(8);ngfftdiel(10:18)=0
     if(dtset%iscf==-1)ngfftdiel(7)=102


     call getng(dtset%boxcutmin,ecutsus,gmet,mpi_enreg_diel%me_fft,mgfftdiel,nfftdiel,ngfftdiel,&
&     mpi_enreg_diel%nproc_fft,dtset%nsym,1,mpi_enreg_diel%paral_kgb,dtset%symrel,&
&     use_gpu_cuda=dtset%use_gpu_cuda)

!    Compute the size of the dielectric matrix
     kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
     call getmpw(diecut,dtset%exchn2n3d,gmet,(/1/),kpt_diel,mpi_enreg_diel,npwdiel,1)
     lmax_diel=0
     if (psps%usepaw==1) then
       do ii=1,dtset%ntypat
         lmax_diel=max(lmax_diel,pawtab(ii)%lcut_size)
       end do
     end if
   else
     npwdiel=1
     mgfftdiel=1
     nfftdiel=1
     lmax_diel=0
   end if

!  Now, performs allocation
   ABI_ALLOCATE(dielinv,(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
   ABI_ALLOCATE(susmat,(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
   ABI_ALLOCATE(kg_diel,(3,npwdiel))
   ABI_ALLOCATE(gbound_diel,(2*mgfftdiel+8,2))
   ABI_ALLOCATE(irrzondiel,(nfftdiel**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
   ABI_ALLOCATE(phnonsdiel,(2,nfftdiel**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
   ABI_ALLOCATE(ph1ddiel,(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw))
   ABI_ALLOCATE(ylmdiel,(npwdiel,lmax_diel**2))
!  Then, compute the values of different arrays
   if(dielop>=1)then
     call status(0,dtfil%filstat,iexit,level,'kpgio(sus)    ')
!    Note : kgnam is dummy, npwarr_diel is dummy, npwtot_diel is dummy
!    This kpgio call for going from the suscep FFT grid to the diel sphere
     npwarr_diel(1)=npwdiel

     call kpgio(diecut,dtset%exchn2n3d,gmet,(/1/),kg_diel,kgnam,&
&     kpt_diel,1,(/1/),1,'COLL',mpi_enreg_diel,npwdiel,&
&     npwarr_diel,npwtot_diel,dtset%nsppol,tmp_unit)
     call sphereboundary(gbound_diel,1,kg_diel,mgfftdiel,npwdiel)
     if (dtset%nsym>1 .and. dtset%iscf>=0 ) then
!      Should replace this initialization of irrzondiel and phnonsdiel through setsym by a direct call to irrzg
       ABI_ALLOCATE(indsym_dum,(4,dtset%nsym,dtset%natom))
       ABI_ALLOCATE(symrec_dum,(3,3,dtset%nsym))
       call setsym(indsym_dum,irrzondiel,dtset%iscf,dtset%natom,&
&       nfftdiel,ngfftdiel,dtset%nspden,dtset%nsppol,dtset%nsym,phnonsdiel,&
&       dtset%symafm,symrec_dum,dtset%symrel,dtset%tnons,dtset%typat,xred)
       ABI_DEALLOCATE(indsym_dum)
       ABI_DEALLOCATE(symrec_dum)
     end if
     if (psps%usepaw==1) then
       call getph(atindx,dtset%natom,ngfftdiel(1),ngfftdiel(2),&
&       ngfftdiel(3),ph1ddiel,xred)
       call initylmg(gprimd,kg_diel,kpt_diel,1,mpi_enreg_diel,&
&       lmax_diel,npwdiel,dtset%nband,1,npwarr_diel,dtset%nsppol,0,&
&       rprimd,tmp_unit,tmp_unit,ylmdiel,rhodum)
     end if
   end if

 else
   npwdiel=1
   mgfftdiel=1
   nfftdiel=1
 end if

 nkxc=0
!TDDFT - For a first coding
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==0)nkxc=0 !MF no xc kernel
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==1)nkxc=0 !MF no xc kern, but (later) RPA ok
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==2)nkxc=1 !MF LDA xc kernel + (later) RPA
 if (dtset%iscf==-1 .and. dtset%nspden==1) nkxc=2
 if (dtset%iscf==-1 .and. dtset%nspden==2) nkxc=3
!Eventually need kxc-LDA when susceptibility matrix has to be computed
 if (dtset%iscf>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79)) nkxc=min(2*dtset%nspden-1,3)
!Eventually need kxc-LDA for residual forces (when density mixing is selected)
 if (dtset%iscf>=10.and.dtset%usewvl==0.and.forces_needed>0 .and. &
& abs(dtset%iprcch)>=1.and.abs(dtset%iprcch)<=6.and.abs(dtset%iprcch)/=5) then
   if (dtset%xclevel==1.or.dtset%iprcch>=0) nkxc=min(2*dtset%nspden-1,3)
   if (dtset%xclevel==2.and.dtset%nspden==2.and.dtset%iprcch<0) nkxc=23
 end if
 if (nkxc>0) then
   call check_kxc(dtset%ixc,dtset%optdriver)
 end if
 ABI_ALLOCATE(kxc,(nfftf,nkxc))

!This flag will be set to 1 just before an eventual change of atomic
!positions inside the iteration, and set to zero when the consequences
!of this change are taken into account.
 moved_atm_inside=0
!This flag will be set to 1 if the forces are computed inside the iteration.
 computed_forces=0

 if(dtset%wfoptalg==2)then
   ABI_ALLOCATE(shiftvector,((dtset%mband+2)*dtset%nkpt))
   val_min=-1.0_dp
   val_max=zero
 else
   ABI_ALLOCATE(shiftvector,(1))
 end if

 call status(0,dtfil%filstat,iexit,level,'berryphase    ')

!!PAW+DMFT: allocate structured datatype paw_dmft if dtset%usedmft=1
!call init_sc_dmft(dtset%dmftbandi,dtset%dmftbandf,dtset%mband,dtset%nkpt,&
!&  dtset%nsppol,dtset%usedmft,paw_dmft,dtset%usedmft)
!call print_sc_dmft(paw_dmft)

!Electric field initializations: initialize pel_cg(:) and p_ion(:)
!HONG  finite displacement fields with PAW needs to test !
 if (dtset%berryopt == 4 .or. abs(dtset%berryopt) == 5 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7 .or. &
& dtset%berryopt ==14 .or. dtset%berryopt ==16 .or. dtset%berryopt ==17) then   !!HONG
!  finite fields with PAW needs cprj
   mcprj=0
   if (psps%usepaw==1) then
     usecprj=1;my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
     mcprj=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
     ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
     ncpgr = 3 ! so that gradients wrt atom position may be included, may add strain later
     call cprj_alloc(cprj,ncpgr,dimcprj)
     iatom=0 ; iorder_cprj=1 ! retain ordering of input list
     ctocprj_choice = 2 ! compute cprj and derivatives wrt atomic position
     useylmgr = 0
!    all arguments to ctocprj are defined already except ph1d, do that here
     call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
     call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,iatom,idir,iorder_cprj,&
&     dtset%istwfk,kg,dtset%kptns,dtset%mband,mcg,mcprj,dtset%mgfft,dtset%mkmem,&
&     mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,dtset%nband,&
&     dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,&
&     dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,&
&     dtset%typat,ucvol,dtfil%unpaw,dtfil%unkg,dtfil%unylm,useylmgr,wffnow,&
&     xred,ylm,ylmgr)
   end if
   optberry=1     ! compute polarization only
   pel_cg(:) = zero;pelev=zero
   call berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,gprimd,hdr,psps%indlmn,kg,&
&   psps%lmnmax,dtset%mband,mcg,mcprj,dtset%mkmem,mpi_enreg,dtset%mpw,my_natom,&
&   dtset%natom,npwarr,dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawrhoij,pawtab,&
&   pel_cg,pelev,pion,ptot,red_ptot,pwind,&
&   pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&   unit_out,usecprj,psps%usepaw,wffnow,xred,psps%ziontypat)

   dtefield%red_ptot1(:)=red_ptot(:)

   efield_old_cart(:)=dtset%efield(:)   !!HONG

   if ( dtset%berryopt ==16 .or. dtset%berryopt ==17) then   !!HONG
     do iidir=1,3
       red_efield2(iidir)=zero
       red_efield2_old(iidir)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,iidir))
     end do
   end if

!  write(std_out,*)pel_cg(1),pel_cg(2),pel_cg(3)
!  write(std_out,*)pelev(1),pelev(2),pelev(3)
!  deallocate cprj
   if (psps%usepaw==1) then
!    pel_cg(:) = pel_cg(:) + pelev(:) ! electronic polarization includes on-site term
!    above line suppressed because in the PAW case, pel_cg already includes all on-site
!    terms and pelev should not be added in additionally. We are computing pelev separately for
!    reporting purposes only.
!    13 June 2012 J Zwanziger
     usecprj=0
     call cprj_free(cprj)
     ABI_DATATYPE_DEALLOCATE(cprj)
   end if
 end if

 if (dtset%iscf==22) energies%h0=zero

 call timab(54,2,tsec)

!##################################################################
!PERFORM ELECTRONIC ITERATIONS
!##################################################################

!Offer option of computing total energy with existing
!wavefunctions when nstep<=0, else do nstep iterations
!Note that for non-self-consistent calculations, this loop will be exited
!after the first call to vtorho
!Pass through the first routines even when nstep==0

!! berryopt == 14   !!!!!!!! ebar field
!Finite displacement field - update efield   !!HONG

 if (dtset%berryopt == 14 .and. quit /=1) then
!  ! Convert polarization to cartesian coords

   ptot_cart(:)=zero
   do iidir = 1,3
     ptot_cart(iidir)=rprimd(iidir,1)*red_ptot(1) + rprimd(iidir,2)*red_ptot(2) + &
&     rprimd(iidir,3)*red_ptot(3)
   end do
   ptot_cart(:)=ptot_cart(:)/ucvol

!  save this value in order to print the final value of real electric field, comparing with the desired red_fieldbar
   dtefield%efield2(:)=dtset%efield(:)

   do iidir=1,3
     red_efieldbar_lc(iidir) = dot_product(dtset%efield(:),rprimd(:,iidir))
     dtefield%efield_dot(iidir) = dot_product(dtset%efield(:),rprimd(:,iidir))
   end do

!  !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
   write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar-field:'

   call wrtout(std_out,message,'COLL')
   call prtefield(dtset,dtefield,std_out,rprimd)

   if(dtset%prtvol>=10)then
     call wrtout(ab_out,message,'COLL')
     call prtefield(dtset,dtefield,ab_out,rprimd)
   end if

!  updating E field
   do iidir =1,3   ! desired E field
     efield_test_cart(iidir)=gprimd(iidir,1)*dtset%red_efieldbar(1) + &
&     gprimd(iidir,2)*dtset%red_efieldbar(2)+gprimd(iidir,3)*dtset%red_efieldbar(3)
   end do

!  if not convergence well, need to add some code here to make sure efield_test_cart(:) not change much
   dtset%efield(:) = efield_test_cart(:)

 end if  ! berryopt ==14


!!!==========================================================================
!!                  start SCF loop
!!!===========================================================================

 do istep=1,max(1,nstep)
!  istep_mix=istep                  !!HONG

   call timab(240,1,tsec)
   if (moved_atm_inside==1 .or. istep==1) then
!    ##############################################################
!    The following steps are done once for a given set of atomic
!    coordinates or for the nstep=1 case
!    --------------------------------------------------------------

!    Eventually symmetrize atomic coordinates over space group elements:
     call symzat(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

     if (dtset%usewvl == 0) then
!      Get cut-off for g-vectors
       if (psps%usepaw==1) call wrtout(std_out,' FFT (fine) grid used in SCF cycle:','COLL')
       call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

!      Compute structure factor phases and large sphere cut-off (gsqcut):
       call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)

       if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
         call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
       else
         ph1df(:,:)=ph1d(:,:)
       end if
     end if

!    Initialization of atomic data for PAW
     if (psps%usepaw==1) then
!      Check for non-overlapping spheres
       call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,xred)

!      Identify parts of the rectangular grid where the density has to be calculated
       optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
       if (forces_needed==1.or.(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)) then
         optgr1=dtset%pawstgylm;if (stress_needed==1) optrad=1; if (dtset%pawprtwf==1) optrad=1
       end if

       if (abs(dtset%berryopt) == 5) optrad = 1 ! need pawfgrtab%rfgd for mag field

       call status(istep,dtfil%filstat,iexit,level,'call nhatgrid ')
       call nhatgrid(atindx1,gmet,my_natom,dtset%natom,&
&       nattyp,ngfftf,psps%ntypat,optcut,optgr0,optgr1,optgr2,optrad,&
&       pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       mpi_comm_fft=spaceComm_fft)

!      magnetic field specific initialization
       if (abs(dtset%berryopt) == 5) then
         call expibr(dtefield,gprimd,my_natom,dtset%natom,pawfgrtab,xred,&
&         mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
       end if
     end if

!    If we are inside SCF cycle or inside dynamics over ions,
!    we have to translate the density of previous iteration
     moved_rhor=0

     if (initialized/=0.and.dtset%usewvl == 0.and.ipositron/=1.and. &
&     (abs(dtset%iprcch)==2.or.abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)) then
       moved_rhor=1
       if (abs(dtset%iprcch)==2) then
         option=2
         ABI_ALLOCATE(workr,(nfftf,dtset%nspden))
         call status(istep,dtfil%filstat,iexit,level,'call fresid   ')
         call fresid(dtset,gresid,mpi_enreg,nfftf,ngfftf,&
&         psps%ntypat,option,pawtab,rhor,rprimd,&
&         ucvol,workr,xred,xred_old,psps%znuclpsp)
         rhor=workr
         ABI_DEALLOCATE(workr)
       else if (abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6) then
         call status(istep,dtfil%filstat,iexit,level,'call extraprho')
         scf_history%icall=scf_history%icall+1
         call extraprho(atindx,atindx1,cg,dtset,gmet,gprimd,gsqcut,&
&         scf_history%icall,kg,mcg,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&         my_natom,nattyp,nfftf,ngfftf,npwarr,psps%ntypat,pawrhoij,pawtab,&
&         ph1df,psps,psps%qgrid_vl,rhor,rprimd,scf_history,ucvol,&
&         psps%usepaw,xred,xred_old,ylm,psps%ziontypat,psps%znuclpsp)
       end if
       call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     end if

   end if ! moved_atm_inside==1 .or. istep==1

!  Initialize/update data in the electron-positron case
   if (dtset%positron<0.or.(dtset%positron>0.and.istep==1)) then
     call setup_positron(atindx,atindx1,cg,dtfil,dtset,ecore,eigen,&
&     etotal,electronpositron,energies,forces_needed,fred,gmet,gprimd,&
&     grewtn,grvdw,gsqcut,hdr,initialized0,indsym,istep,istep_mix,kg,&
&     kxc,maxfor,mcg,mgfftf,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfftf,ngrvdw,nhat,&
&     nkxc,npwarr,nvresid,occ,optres,paw_ij,pawang,pawfgr,pawfgrtab,&
&     pawrhoij,pawtab,ph1df,ph1d,psps,rhog,rhor,rprimd,&
&     stress_needed,strsxc,symrec,ucvol,vhartr,vpsp,vxc,&
&     wffnow,xccc3d,xred,ylm,ylmgr)
     ipositron=electronpositron_calctype(electronpositron)
   end if

   if ((moved_atm_inside==1 .or. istep==1).or. (dtset%positron<0.and.istep_mix==1)) then
!    PAW only: we sometimes have to compute compensation density
!    and eventually add it to density from WFs
     nhatgrdim=0
     if (psps%usepaw==1.and.(dtset%positron>=0.or.ipositron/=1) &
&     .and.((usexcnhat==0) &
&     .or.(dtset%xclevel==2.and.(dtfil%ireadwf/=0.or.dtfil%ireadden/=0.or.initialized/=0)) &
&     .or.(dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0))) then
       call timab(558,1,tsec)
       nhatgrdim=0;if (dtset%xclevel==2) nhatgrdim=usexcnhat*dtset%pawnhatxc
       ider=2*nhatgrdim;izero=0
       if (nhatgrdim>0)   then
         ABI_ALLOCATE(nhatgr,(cplex*nfftf,dtset%nspden,3*nhatgrdim))
       end if
       call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
&       nfftf,ngfftf,nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,&
&       nhatgr,nhat,pawrhoij,pawrhoij,pawtab,k0,rprimd,ucvol,xred,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       mpi_comm_fft=spaceComm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0)
       if (dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0) then
         rhor(:,:)=rhor(:,:)+nhat(:,:)
         call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       end if
       call timab(558,2,tsec)
     end if

!    The following steps have been gathered in the setvtr routine:
!    - get Ewald energy and Ewald forces
!    - compute local ionic pseudopotential vpsp
!    - eventually compute 3D core electron density xccc3d
!    - eventually compute vxc and vhartr
!    - set up vtrial

     call status(istep,dtfil%filstat,iexit,level,'call setvtr   ')
     if (dtset%usewvl == 0) then
       optene = 4 * optres
       if(dtset%iscf==-3)optene=4
     else
!      We need the Hartree energy for the wavefunctions mixing
       optene = 1
     end if

     call setvtr(atindx1,dtset,energies,gmet,gprimd,grewtn,grvdw,gsqcut,&
&     istep,kxc,mgfftf,moved_atm_inside,moved_rhor,mpi_enreg,&
&     nattyp,nfftf,ngfftf,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,psps%ntypat,&
&     n1xccc,n3xccc,optene,pawtab,ph1df,psps,rhog,rhor,rmet,rprimd,&
&     strsxc,ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,&
&     xccc3d,xred,electronpositron=electronpositron,&
&     taug=taug,taur=taur,vxctau=vxctau)

!    DEBUG
!    write(std_out,*)' scfcv : after setvtr, energies%e_hartree=',energies%e_hartree
!    ENDDEBUG

     if (nhatgrdim>0.and.nstep>0)  then
       ABI_DEALLOCATE(nhatgr)
     end if

!    Recursion Initialisation
     if(dtset%userec==1 .and. istep==1)  then
       rec_set%quitrec = 0
!      --At any step calculate the metric
       call Init_MetricRec(rec_set%inf,rec_set%nl%nlpsp,rmet,ucvol,rprimd,xred,&
&       dtset%ngfft(1:3),dtset%natom,rec_set%debug)
       if(initialized==0)  call first_rec(dtset,psps,rec_set)
     end if

!    End the condition of atomic position change or istep==1
   end if
   call timab(240,2,tsec)
   call timab(241,1,tsec)

!  ######################################################################
!  The following steps are done at every iteration
!  ----------------------------------------------------------------------

!  PAW: Compute energies and potentials in the augmentation regions (spheres)
!  Compute pseudopotential strengths (Dij quantities)
   if (psps%usepaw==1)then
!    "on-site" energies, potentials, densities computation
     nzlmopt=0;if (istep_mix==2.and.dtset%pawnzlm>0) nzlmopt=-1
     if (istep_mix>2) nzlmopt=dtset%pawnzlm
     option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
     do iatom=1,my_natom
       itypat=paw_ij(iatom)%itypat
       v_size=paw_an(iatom)%lm_size;if (dtset%pawxcdev==0) v_size=paw_an(iatom)%angl_size
       paw_ij(iatom)%has_dij=1;paw_ij(iatom)%has_dijhartree=1
       paw_an(iatom)%has_vxc=1
       ABI_ALLOCATE(paw_ij(iatom)%dijhartree,(pawtab(itypat)%lmn2_size))
!      already allocated in init_paw_an
       if (associated(paw_an(iatom)%vxc1))  then
         ABI_DEALLOCATE(paw_an(iatom)%vxc1)
       end if
       ABI_ALLOCATE(paw_an(iatom)%vxc1 ,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
       if (associated(paw_an(iatom)%vxct1))  then
         ABI_DEALLOCATE(paw_an(iatom)%vxct1)
       end if
       ABI_ALLOCATE(paw_an(iatom)%vxct1,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
       if (pawtab(itypat)%useexexch>0) then
         if (associated(paw_an(iatom)%vxc_ex))  then
           ABI_DEALLOCATE(paw_an(iatom)%vxc_ex)
         end if
         ABI_ALLOCATE(paw_an(iatom)%vxc_ex,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
       end if
     end do

!    Local exact exch.: impose occ. matrix if required
     if (dtset%useexexch>0) then
       call setrhoijpbe0(dtset,psps%indlmn,initialized0,istep,istep_mix,psps%lmnmax,&
&       spaceComm,my_natom,dtset%natom,dtset%ntypat,pawrhoij,pawtab,dtset%typat,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if

     call status(istep,dtfil%filstat,iexit,level,'call pawdenpot')

!    Computation of on-site densities/potentials/energies
     call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,dtset%ixc,my_natom,dtset%natom,&
&     dtset%nspden,psps%ntypat,nzlmopt,option,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,&
&     pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,psps%znuclpsp,&
&     mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     electronpositron=electronpositron)

!    PAW+U: impose density matrix if required
     if (dtset%usepawu>0.and.(ipositron/=1)) then
       impose_dmat=0
       if ((istep<=abs(dtset%usedmatpu)).and.(dtset%usedmatpu<0.or.initialized0==0)) impose_dmat=1
       if (impose_dmat==1.or.dtset%dmatudiag/=0) then
         dimdmat=0;if (impose_dmat==1) dimdmat=2*lpawumax+1
         call setnoccmmp(0,dimdmat,&
&         dtset%dmatpawu(1:dimdmat,1:dimdmat,1:dtset%nsppol*dtset%nspinor,1:dtset%natpawu*impose_dmat),&
&         dtset%dmatudiag,impose_dmat,indsym,my_natom,dtset%natom,dtset%natpawu,&
&         dtset%nspinor,dtset%nsppol,dtset%nsym,dtset%ntypat,paw_ij,pawang,dtset%pawprtvol,&
&         pawrhoij,pawtab,dtset%spinat,dtset%symafm,dtset%typat,0,dtset%usepawu,&
&         mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!        Reinitalize mixing if PAW+U and occupation matrix now allowed to change
!        For experimental purpose...
         if ((dtset%userib==1234).and.(istep==abs(dtset%usedmatpu)+1).and.(dtset%usedmatpu<0.or.initialized0==0)) istep_mix=1
       end if
     end if
!    Dij computation
     call status(istep,dtfil%filstat,iexit,level,'call pawdij   ')
     call pawdij(cplex,dtset,dtset%enunit,fatvshift,gprimd,ipert,&
&     my_natom,dtset%natom,nfftf,ngfftf,dtset%nspden,psps%ntypat,&
&     paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,&
&     pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,k0,ucvol,&
&     vtrial,vxc,xred,&
&     mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     mpi_comm_fft=spaceComm_fft,electronpositron=electronpositron)

!    if berryopt = \pm 5, need the phase twisted dij terms
     if(abs(dtset%berryopt)==5) then
       call pawtwdij(dtefield,gprimd,my_natom,dtset%natom,nfftf,dtset%nspden,psps%ntypat,&
&       paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,vtrial,vxc,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if
     do iatom=1,my_natom
       itypat=paw_ij(iatom)%itypat
       ABI_DEALLOCATE(paw_ij(iatom)%dijhartree)
       paw_an(iatom)%has_vxc=0
       paw_ij(iatom)%has_dijhartree=0
       ABI_DEALLOCATE(paw_an(iatom)%vxc1)
       ABI_DEALLOCATE(paw_an(iatom)%vxct1)
       if (pawtab(itypat)%useexexch>0)  then
         ABI_DEALLOCATE(paw_an(iatom)%vxc_ex)
       end if
     end do
     call status(istep,dtfil%filstat,iexit,level,'call symdij   ')
     call symdij(gprimd,psps%indlmn,indsym,ipert,psps%lmnmax,my_natom,dtset%natom,dtset%nsym,&
&     psps%ntypat,0,paw_ij,pawang,dtset%pawprtvol,rprimd,dtset%symafm,symrec,&
&     mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     if (dtset%iscf==22) then
       call symdij(gprimd,psps%indlmn,indsym,ipert,psps%lmnmax,my_natom,dtset%natom,dtset%nsym,&
&       psps%ntypat,1,paw_ij,pawang,dtset%pawprtvol,rprimd,dtset%symafm,symrec,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if

   end if

!  Write out occupancies to dtpawuj-dataset
   if (dtset%usepawu>0.and.dtset%macro_uj>0.and.istep>1.and.ipositron/=1) then
     call pawuj_red(dtset,dtpawuj,fatvshift,my_natom,dtset%natom,dtset%ntypat,&
     paw_ij,pawrad,pawtab,ndtpawuj,&
&     mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   end if

   call timab(241,2,tsec)

!  No need to continue and call vtorho, when nstep==0
   if(nstep==0)exit

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------
   call timab(56,1,tsec)
   call status(istep,dtfil%filstat,iexit,level,'loop istep    ')

   if(dtset%iscf>=0)then
     write(message, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,message,'COLL')
   end if

!  DEBUG XG 120518 Please do not remove : apparently needed to avoid segfault on ibm6
!  jmb 2012   write(std_out, '(a)' )' scfcv : 1 '
!  ENDDEBUG

!  The next flag says whether the xred have to be changed in the current iteration
   moved_atm_inside=0
   if(dtset%ionmov==4 .and. mod(iapp,2)/=1 .and. dtset%iscf>=0 )moved_atm_inside=1
   if(dtset%ionmov==5 .and. iapp/=1 .and. istep==1 .and. dtset%iscf>=0)moved_atm_inside=1

!  The next flag says whether the forces have to be computed in the current iteration
   computed_forces=0
   if ((dtset%optforces==1 .and. dtset%usewvl == 0).or.(moved_atm_inside==1)) computed_forces=1
   if (abs(tollist(3))>tiny(0._dp)) computed_forces=1
   if (dtset%iscf<0) computed_forces=0
   if ((istep==1).and.(dtset%optforces/=1)) then
     if (moved_atm_inside==1) then
       write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&       ' scfcv : WARNING -',ch10,&
&       '  Although the computation of forces during electronic iterations',ch10,&
&       '  was not required by user, it is done (required by the',ch10,&
&       '  choice of ionmov input parameter).'
       call wrtout(std_out,message,'COLL')
     end if
     if (abs(tollist(3))+abs(tollist(7))>tiny(0._dp)) then
       write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&       ' scfcv : WARNING -',ch10,&
&       '  Although the computation of forces during electronic iterations',ch10,&
&       '  was not required by user, it is done (required by the',ch10,&
&       '  "toldff" or "tolrff" tolerance criteria).'
       call wrtout(std_out,message,'COLL')
     end if
   end if
   if ((istep==1).and.(dtset%optforces==1).and. dtset%usewvl == 1) then
     write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&     ' scfcv : WARNING -',ch10,&
&     '  Although the computation of forces during electronic iterations',ch10,&
&     '  was required by user, it has been disable since the tolerence',ch10,&
&     '  is not on forces (force computation is expensive in wavelets).'
     call wrtout(std_out,message,'COLL')
   end if

   call timab(56,2,tsec)

!  ######################################################################
!  Compute the density rho from the trial potential
!  ----------------------------------------------------------------------

   call timab(242,1,tsec)
!  Compute the density from the trial potential
   if (dtset%tfkinfunc==0) then

   ! ======================================================================================
   ! =============================== chenh trap the code ==================================
   ! ============================= Steven add in spin case ================================
   ! ======================================================================================

   !call xcomm_world(mpi_enreg,spaceworld)
   spaceworld = spaceComm

   ! initialize parameters:
   My_vtrial = 0.d0
   bAllowsym = 0;
   bLoadVr   = 1;
   bSpin     = 1;    !Default spin number is 1, namely spin unpolarized case

   ! BFGS ....
   task = 'START'
   BFGS_nbd = 0 
   nline    = 0
   nBFGS    = 0

   !------------------
   myicount = 0;
   cgW      = 0.d0 ;
   maxcount = 200;
   ektol = -1.d0 ;   !control the loop of vtorhor, to converge more, make it smaller
   n1=ngfft(1); 
   n2=ngfft(2); 
   n3=ngfft(3);
   histW = 0.d0
   
   open(unit=3333,access='sequential',action='write',status='replace', file='diff',form='formatted');
   close (3333);
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MAIN LOOP START  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
   
   
   ! load parameters for job ----------------------------------
   call load_param(bLoadDen,bLoadVr,maxcount,bAllowsym,bSpin,wtol,structure,ektol,penLambda)
   
   ! test if we are running in tigger mode
   ! in tigger mode, the code will monitor a file called "message.out"
   ! and reads refden.in when it is triggered
   if (bSpin==2) then
     call c_wrtlog('The calculation is spin-polarized!')
   endif

   if (structure(1:11)=='TIGGER_MODE') then 
     tigger_mode = .true.
     write(message,'(a)')'=== THE CODE RUNS IN TIGGER_MODE ==='
   else
     tigger_mode = .false.
   endif

5555 continue 

   ! if in tigger_mode, we wait for the message.out file 
   !====================================================
   if (tigger_mode) then 
     task      = 'START'
     nBFGS     = 0       ! iter counter for BFGS
     maxLine   = 5       ! max line search limit
     nline     = 0       ! iter counter for line search during BFGS
     myicount  = 0       ! total f/g iter
     istep_mix = 1
!!     vtrial   = 0.d0
!!     last_x   = vtrial(:,1)
     histW    = 0.d0
     cgW      = 0.d0
     cgG      = 0.d0
     penfun   = 0.d0
     call c_wrtlog('')
     call c_wrtlog('=============================================')
     call c_wrtlog('[tigger_mode] waiting for file: message.out')
     call c_wrtlog('=============================================')
     call c_wrtlog('')

     ! check the existence of the message.out
     !========================================
     do
       inquire(file='message.out', exist=ex)
       if ( ex .eqv. .false. )  then
         call sleep(2)
         cycle
       else
         exit
       endif
     enddo
     write(6,'(a,a)') '[tigger_mode] found file message.out'

     ! read in message.out file 
     !==========================
     call xbarrier_mpi(spaceworld) ! force synchronization
     do
       open(unit=111,file='message.out',action='read',iostat=io)
       read(111,*,iostat=io)message,message
       close(111)
       if (io/=0 .or. message(1:3)/='end') then
         call sleep(1)
       else
         call c_wrtlog('[tigger_mode] get the message.out file, move on.')
         exit
       endif
     enddo

     ! notify uopt to delete message.out file 
     ! NOTE:
     !  This is a water-proof way to communicate with the UOPT
     !==========================================================
     call xbarrier_mpi(spaceworld) ! force synchronization
     if (mpi_enreg%me==0) then
       open(unit=111,file='please_remove_message',action='write',iostat=io)
       write(111,*)'remove message'
       call flush(111)
       close(111)
       ! wait until message.out is removed by UOPT
       do 
         inquire(file='message.out', exist=ex)
         if ( ex )  then
           call sleep(1)
           cycle
         else
           exit
         endif
       enddo
       ! remove our signal: please_remove_message file
       open(unit=111,file='please_remove_message',action='read',iostat=io)
       close(111,status='delete')
     endif
     call xbarrier_mpi(spaceworld) ! force synchronization

   endif 
   !!=====================
   !! End of tigger_mode
   !!=====================



   ! WRITE INFORMATIONS
   !=======================
   call c_wrtlog(' ')
   call c_wrtlog('WELCOMEMSG_BF')
   call c_wrtlog(' ')
   call date_and_time(VALUES=dateTime)    !today
   write(msg, '(A,i2,A,i2,A,i4,A,i2.2,A,i2.2,A,i2.2,A,i3.3)') "Start at: ", &
               dateTime(3),"/",  dateTime(2),"/",  dateTime(1)," : ",             &
               dateTime(5), ":", dateTime(6), ":", dateTime(7), ".", dateTime(8)
   call c_wrtlog(msg);
   CALL C_WRTLOG('====== PARAMETERS =======');
   write (message, '(A,I10,A,I5,A,I5,A,I5)')  &
    'nfft=',dtset%nfft,' n1=',n1,' n2=',n2,' n3=',n3
   call c_wrtlog(message)
   write (message, '(A,es12.4,a,es12.4)') & 
    'ektol = ',ektol,' wtol=',wtol
   call c_wrtlog(message)
   write (message, '(A,E15.7,A,9I2,a,I2,a,E15.7,a,E15.7)') & 
    'ecut=',dtset%ecut,' nkpt=',dtset%kptrlatt, ' occ=', dtset%occopt, & 
    ' tsmear=', dtset%tsmear, ' ucvol=', ucvol
   call c_wrtlog(message);
   write(message,*) 'penalty coef = ',penLambda
   call c_wrtlog(message)
   
   ! check the 'nsym' flag in the ABINIT input file
   if (bAllowSym == 0) then
     if (dtset%nsym/=1) then
       call c_wrtlog(">>> You want to the symmtry or not, bAllowsym=0 but nsym!=1, confusing... stop.")
       EXIT
     end if
   else
     call c_wrtlog('>>> You set the bAllowSym = 1')
     if (dtset%nsym /= 1) then
       call c_wrtlog('>>> And nsym>1, You have switched symmetry ON.')
     else
       call c_wrtlog('>>> The dtset%nsym=1 but you set bAllowSym/=0 Confusing. STOP')
       EXIT
     end if
   end if

   ! READ IN THE REFERENCE DENSITY
   if (bLoadDen==0 .and. bloadvr==1 .and. dtset%nsym==1) then
       call c_wrtlog("You want to do Vr/structure_factor? You need to commnet off nsym=1 in the input file!")
       EXIT
   end if
   if (bLoadDen == 1) then
     open(unit=1111, access="sequential", action="read", &
       file='refden.in', form="formatted", status="old");
     read(1111,*) ref_rhor;
     close(1111);
     tmp1 = 0.d0
     do i=1,dtset%nfft
        tmp1 = tmp1 + ref_rhor(i)
     end do
     tmp1 = tmp1*ucvol/real(dtset%nfft, kind=dp)
     write(msg,'(a,2es12.4,a,f)') 'min/max(ref_rhor1)=',minval(ref_rhor),maxval(ref_rhor),' sum=',tmp1
     call c_wrtlog(msg)
     call c_wrtlog('Loaded refden.in file.');
     
     ! For spin-polarized case, load spin down density too
     if (bSpin == 2) then
       open(unit=1113, access="sequential", action="read", &
         file='refden2.in', form="formatted", status="old");
       read(1113,*) ref_rhor2;
       close(1113);
       tmp1 = 0.d0
       do i=1,dtset%nfft
         tmp1 = tmp1 + ref_rhor2(i)
       end do
       tmp1 = tmp1*ucvol/real(dtset%nfft, kind=dp)
       write(msg,'(a,2es12.4,a,f)') 'min/max(ref_rhor2)=',minval(ref_rhor2),maxval(ref_rhor2),' sum=',tmp1
       call c_wrtlog(msg)
       call c_wrtlog('Loaded refden2.in file.');
     endif

   else
     call c_wrtlog('not load the refden.in, be cautious!')
   end if

   ! PREPARE VTRIAL
   !==================
   if (bLoadVr == 1) then
     open (unit=1111, access="sequential", action="read", &
           file='refpot.in', form="formatted", status="old");
     read (1111, *) vtrial(:,1);
     close(1111);
     
     if(bSpin==2) then
       open (unit=1112, access="sequential", action="read", &
             file='refpot2.in', form="formatted", status="old");
       read (1112, *) vtrial(:,2);
       close(1112);
     endif
    

     ! To do dividing structure factor thing
     call c_wrtlog("Load in refpot.in and set it to vtrial")
     call c_wrtlog("have symtrized the Vtrial")
     call c_wrtlog("Start to make atomic Vg ...")
     
     ! debug
     call c_wrtlog('xred:')
     do i=1, dtset%natom
        write (msg,*) xred(:,i)
        call c_wrtlog(msg)
     end do
     ! end debug
     
     if(bSpin==1) then
       call symrhg(1,gprimd,irrzon,mpi_enreg,dtset%nfft,&
                   dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%ngfft,&
                   dtset%nspden,dtset%nsppol,dtset%nsym,dtset%paral_kgb,phnons, &
                   rhog,vtrial,rprimd,dtset%symafm,dtset%symrel)
       call c_DoAtomV(rhog,vtrial,psps,nattyp,dtset%mgfft,gmet,&
                      gsqcut,mpi_enreg,dtset,ph1d,ucvol,structure)

     else
       call c_wrtlog("Divding S(q) does not support spin=2 calculation!")
       call c_wrtlog("Set spin to 1 and load two spin channel potentials seperately!")
       call c_wrtlog("Just do it twice!")
     endif

     EXIT
   else 
     ! vtrial is not loaded 
     ! prepare init vtrial based on refden
     !=======================================
     if (bLoadDen == 0) then
       call c_wrtlog('error! You set both bLoadpot and bLoadDen == 0! stop!')
       EXIT
     end if
     if (restart .eqv. .false. .and. tigger_mode .eqv. .false.) then  ! mohan changed "==" to "eqv"
       if(bSpin==1) then
         rhor(:,1) = ref_rhor(:);
       else
         rhor(:,1) = ref_rhor(:) + ref_rhor2(:)! strange Abinit seems store total density on the first dimension
         rhor(:,2) = ref_rhor(:)               ! and spin up density on the second dimension 
       endif

       call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,dtset%nfft,ngfft,dtset%paral_kgb,0);
       ! ABINIT v6 setvtr subroutine
       call setvtr(atindx1,dtset,energies,gmet,gprimd,grewtn,grvdw,gsqcut,&
&      istep,kxc,mgfftf,moved_atm_inside,moved_rhor,mpi_enreg,&
&      nattyp,nfftf,ngfftf,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,psps%ntypat,&
&      n1xccc,n3xccc,optene,pawtab,ph1df,psps,rhog,rhor,rmet,rprimd,&
&      strsxc,ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,&
&      xccc3d,xred,electronpositron=electronpositron,&
&      taug=taug,taur=taur,vxctau=vxctau)

       call c_wrtlog("set vtrial = vhart + vxc + vpsp.")
       ! ==== if spin-polarized, prepare My_vtrial
       if(bSpin==2) then
         do i=1,dtset%nfft
           My_vtrial(i) = vtrial(i,1)
           My_vtrial(i+dtset%nfft) = vtrial(i,2)
         enddo
       endif
       ! ===========================================
     end if
   endif

   call c_wrtlog('Inverting the Kohn-Sham equation (v6.0)....')

   !--------------------------------------------------------------------------------------
   !----------------------- main loop started , by chen huang-----------------------
   !--------------------------------------------------------------------------------------
   DO  
       ! BEGIN CHEN'S LOOP... I KNOW THIS IS A INFINITE LOOP
       if (task(1:5) == 'START') goto  2012

       call c_wrtlog('')
       write(message, '(a,i4,a)')  & 
         '------------------- TRIED POTENTIAL ', MYICOUNT,' TIMES --------------------';
       !================ output vtrial info ============================
       call c_wrtlog(message)
       write(message,'(a,2es14.6,a,es14.6)')'min/max(vtrial1):', & 
          minval(vtrial(:,1)),maxval(vtrial(:,1)),' AMP(centered)=',maxval(vtrial(:,1))-minval(vtrial(:,1))
       call c_wrtlog(message)
       
       if(bSpin==2) then
         write(message,'(a,2es14.6,a,es14.6)')'min/max(vtrial2):', &
             minval(vtrial(:,2)),maxval(vtrial(:,2)),' AMP(centered)=',maxval(vtrial(:,2))-minval(vtrial(:,2))
         call c_wrtlog(message)
       endif
       !================================================================

       eEvalcount_c = 0;
       eold1_c = 0.0d0
       eold2_c = 0.0d0
       ektol   = 1e-7
       write (msg,'(a,es12.4)')'tol=',ektol
       call c_wrtlog(msg)
 
       ! vtrial -> rhor loop
       !======================
       if(bSpin==1) then
         call c_wrtlog('         ek             enl            eeig            wfres        rho_resid')
       else
         call c_wrtlog('         ek             enl            eeig            wfres        rho_resid1            rho_resid2')
       endif


       do while ( eEvalcount_c <= 100 ) 

         call xbarrier_mpi(spaceworld) ! force synchronization
         
         !=================================
         ! NOTE: iscf>10 optres = 1 (den mixing) 
         !       iscf<10 optres = 0 (pot mixing)
         !       optres=0: the new value of the density is computed in place of the input value
         !              1: only the density residual is computed ; the input density is kept
         optres       = 0
!         dtset%tolwfr = 1e-5
!         dtset%nline  = 10
!         dtset%nnsclo = 2
!         dtset%iscf = -3
        
         !=============== record the old density=============
         if(bSpin==1) then
           pre_rhor = rhor(:,1)
         else
           pre_rhor = rhor(:,2)
           pre_rhor2 = rhor(:,1) - rhor(:,2)
         endif
         !===================================================


!         call status(istep,dtfil%filstat,iexit,level,'call vtorho   ')
         call vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&
&         dielop,dielstrt,dphase,dtefield,dtfil,dtset,&
&         eigen,electronpositron,energies,etotal,gbound_diel,&
&         gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
!&         istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
&         istep_mix,istep_mix,kg,kg_diel,kxc,lmax_diel,  mcg,  mgfftdiel,mpi_enreg,&
&         my_natom,   dtset%natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&         npwarr,npwdiel,res2,   psps%ntypat,nvresid,occ,&
&         computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
&         pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
&         pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,rmet,rprimd,&
&         susmat,symrec,taug,taur,ucvol,wffnew,wffnow,vtrial,wvl,xred,&
&         ylm,ylmgr,ylmdiel,vxctau=vxctau)
!         call status(istep,dtfil%filstat,iexit,level,'after vtorho  ')
         
         istep_mix = istep_mix + 1
         
         !======================= output energies ===================
         if(bSpin==1) then
           write(msg,'(5(es16.8))') & 
           energies%e_kinetic,& 
           energies%e_nonlocalpsp, & 
           energies%e_eigenvalues, & 
           residm, &
           maxval(abs(rhor(:,1)-pre_rhor))
         else
           write(msg,'(6(es16.8))') &
             energies%e_kinetic,&
             energies%e_nonlocalpsp, &
             energies%e_eigenvalues, &
             residm, &
             maxval(abs(rhor(:,2)-pre_rhor)), &
             maxval(abs(rhor(:,1)-rhor(:,2)-pre_rhor2))
         endif
         call c_wrtlog(msg)
         !=============================================================

         !===================== see if ek converged ===================
         eEvalcount_c = eEvalcount_c + 1;
         if (eEvalcount_c>=4) then
           if (dabs(energies%e_kinetic-eold1_c)<ektol .and.  &
               dabs(energies%e_kinetic-eold2_c)<ektol ) then
             call c_wrtlog('=== converged to ektol ===')
             exit
           else
             ! restore 1st one
             eold2_c = eold1_c
             eold1_c = energies%e_kinetic
           endif
         else
           eold2_c = eold1_c
           eold1_c = energies%e_kinetic
         end if         
         !=============================================================
       end do   ! end of vtrial->rhor Loop
       !=======================
    
       ! compute XC and hatree energy
       vtmp = vtrial 
       call rhotov(dtset,energies,gprimd,gsqcut,   istep,   kxc,mpi_enreg,nfftf,ngfftf, &
&        nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,&
&        optene,optres,optxc,&
&        rhog,rhor,rprimd,strsxc,ucvol,psps%usepaw,usexcnhat,&
&        vhartr,vnew_mean,vpsp,vres_mean,res2,vtmp,vxcavg,vxc,   wvl,  xccc3d,  xred, &
&        electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau)
      
       ! compute total energy 
       energies%e_entropy = - energies%entropy * dtset%tsmear ! entropy energy
       etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&               energies%e_localpsp + energies%e_corepsp + &
&               energies%e_entropy +  energies%e_ewald + energies%e_nonlocalpsp
       
       !=========== update My_vtrial array if in spin-polarized case============
       if(bSpin==2) then
         do i=1,dtset%nfft
           My_vtrial(i) = vtrial(i,1)
           My_vtrial(i+dtset%nfft) = vtrial(i,2)
         Enddo
       endif
       !========================================================================

       !=============== output rho info =====================
       if(bSpin==1) then
         write(message,'(a,2es12.4,a,f)') & 
           'max/min(rhor):',minval(rhor),maxval(rhor), & 
           ' sum:',sum(rhor)*ucvol/dble(dtset%nfft)
         call c_wrtlog(message)
       else
         write(message,'(a,2es12.4,a,f)') &
           'max/min(rhor1):',minval(rhor(:,2)),maxval(rhor(:,2)), &
           ' sum:',sum(rhor(:,2))*ucvol/dble(dtset%nfft)
         call c_wrtlog(message)

         write(message,'(a,2es12.4,a,f12.4)') & !mohan updated 2016-05-07
           'max/min(rhor2):',minval(rhor(:,1)-rhor(:,2)),maxval(rhor(:,1)-rhor(:,2)), &
             ' sum:',sum(rhor(:,1)-rhor(:,2))*ucvol/dble(dtset%nfft)
             call c_wrtlog(message)
       endif
       !======================================================
       !=================== output energies ==================
       write(message,'(a,es16.8)')'etotal     : ', etotal
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'e_entropy  : ', energies%e_entropy
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'kinetic    : ', energies%e_kinetic
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'hartree    : ', energies%e_hartree
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'e_xc       : ', energies%e_xc
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'e_localpsp : ', energies%e_localpsp
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'e_nonlocal : ', energies%e_nonlocalpsp
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'e_corepsp  : ', energies%e_corepsp
       call c_wrtlog(message)
       write(message,'(a,es16.8)')'e_ewald    : ', energies%e_ewald
       call c_wrtlog(message)
       !========================================================

       !===================== get G =========
       CALL C_WRTLOG('=== GET G ===');
       if(bSpin==1) then
         cgG = ref_rhor - rhor(:,1);                 !assume spin number = 1
       else
         cgG = ref_rhor - rhor(:,2);
         cgG2 = ref_rhor2 - (rhor(:,1)-rhor(:,2))
       endif

       write(message,'(a,3es16.8,a)')'(scfcv) AverNorm/min/max(cgG):',Sum(abs(cgG))/size(cgG), &
                                                                      minval(cgG),maxval(cgG),' (no penalty)'
       call c_wrtlog(message)
       if(bSpin==2) then
         write(message,'(a,3es16.8,a)')'(scfcv) AverNorm/min/max(cgG2):',Sum(abs(cgG2))/size(cgG2), &
                                                                         minval(cgG2),maxval(cgG2),' (no penalty)'
         call c_wrtlog(message)
       endif
       !======================================================
       !=================== calculate penalty ================
       call laplacian(gprimd,mpi_enreg,nfftf,1,ngfftf,dtset%paral_kgb,rdfuncr=vtrial(:,1),laplacerdfuncr=lapvfg)
       if(bSpin==2) then
         call laplacian(gprimd,mpi_enreg,nfftf,1,ngfftf,dtset%paral_kgb,rdfuncr=vtrial(:,2),laplacerdfuncr=lapvfg2)
       endif
       cgG = cgG - 2._DP * penLambda * lapvfg
       if(bSpin==2) cgG2 = cgG2 - 2._DP * penLambda * lapvfg2
       write(message,'(a,3es16.8,a)')'(scfcv) AverNorm/min/max(cgG):',Sum(abs(cgG))/size(cgG), &
                                                                      minval(cgG),maxval(cgG),' (with penalty)'
       call c_wrtlog(message)
       if(bSpin==2) then
         write(message,'(a,3es16.8,a)')'(scfcv) AverNorm/min/max(cgG2):',Sum(abs(cgG2))/size(cgG2), &
                                                                         minval(cgG2),maxval(cgG2),' (with penalty)'
         call c_wrtlog(message)
         !============= also update the cgGCombine array ==========
         do i=1,dtset%nfft
           cgGCombine(i) = cgG(i)
           cgGCombine(i+dtset%nfft) = cgG2(i)
         enddo
         !==========================================================
       endif
       !========================================================
       
       if(bSpin==1) then
         write(message,'(a,3es12.4)') & 
           '(scfcv) AverNorm/min/max(density difference):', & 
           Sum(abs(ref_rhor-rhor(:,1)))/size(ref_rhor), minval(ref_rhor-rhor(:,1)),maxval(ref_rhor-rhor(:,1))
         call c_wrtlog(message)
       else
         write(message,'(a,3es12.4)') &
           '(scfcv) AverNorm/min/max(spin up density difference):', &
           Sum(abs(ref_rhor-rhor(:,2)))/size(ref_rhor), minval(ref_rhor-rhor(:,2)),maxval(ref_rhor-rhor(:,2))
         call c_wrtlog(message)
         write(message,'(a,3es12.4)') &
           '(scfcv) AverNorm/min/max(spin dn density difference):', &
           Sum(abs(ref_rhor2-(rhor(:,1)-rhor(:,2))))/size(ref_rhor2), &
           minval(ref_rhor2-(rhor(:,1)-rhor(:,2))),maxval(ref_rhor2-(rhor(:,1)-rhor(:,2)))
         call c_wrtlog(message)
       endif

       !!
       !! IMPORTANT NOTES:
       !!
       !!   W MUST contains the entropy term, if we are using smearing
       !!   In smearing, we are populating electrons (Mermin ensemble DFT)
       !!   this makes replace the total energy functional with the free energy.
       !!
       !!   Moreover, in the derivation of the dW[V]/dV in Wu-Yang OEP 
       !!  
       !!     dW/dV = \int{ d(Ts+entropy)/dV * drho/dV } + (rho-rho0)
       !!
       !!   The entropy term makes sure that the first term in the integral 
       !!   to be the chemical potential, therefore the first term is finally zero.
       !!
       !============================get W ========================
       CALL C_WRTLOG('=== GET W ===');
       if(bSpin==1) then
         cgW    = -(  energies%e_kinetic     + &            
                    energies%e_nonlocalpsp + &        
                    energies%e_entropy     + & 
                    dot_product(vtrial(:,1),rhor(:,1)-ref_rhor)*ucvol/dble(nfftf) & 
                   )
         penfun = -penLambda*dot_product(vtrial(:,1),lapvfg)*ucvol/dble(nfftf)
       else
         cgW    = -(  energies%e_kinetic     + &
                      energies%e_nonlocalpsp + &
                      energies%e_entropy     + &
                      dot_product(vtrial(:,1),rhor(:,2)-ref_rhor)*ucvol/dble(nfftf) + &
                      dot_product(vtrial(:,2),rhor(:,1)-rhor(:,2)-ref_rhor2)*ucvol/dble(nfftf) &
                   )

         penfun = -penLambda*dot_product(vtrial(:,1),lapvfg)*ucvol/dble(nfftf) &
                  -penLambda*dot_product(vtrial(:,2),lapvfg2)*ucvol/dble(nfftf)
       endif
       cgW    = cgW + penfun
       !=============================================================

!! NOTE:
!!   First order correction to calculate Ts[ref_rhor]
!!   based on what we actually obtained => Ts[rhor].
!!   In my potential formalism for embedding theory paper, I disscussed in detail
!!   how to make this first order correction to the sum of Ts + NLPS + entropy energies
!!   Here I obtain the Ts by subtracting the NLPS and entropy energies from cgW
       ke_firstOrder = - cgW + penfun - energies%e_nonlocalpsp - energies%e_entropy
       write(message,'(a,es16.8)')'first-order corrected KE=',ke_firstOrder
       call c_wrtlog(message)

       !---------------------------------- 
       !-------------- BFGS --------------
       !----------------------------------
2012   continue        
       CALL C_WRTLOG('=== PERFORM   BFGS ===')
       BFGS_print = 99
       if(bSpin==1) then
         call setulb(dtset%nfft,mmax,vtrial(:,1),BFGS_l,BFGS_u,BFGS_nbd,cgW,cgG,  & 
                   0._DP,0._DP,BFGS_wa,BFGS_iwa, & 
                   task,BFGS_print,BFGS_csave,BFGS_lsave,BFGS_isave,BFGS_dsave)
       else
         call setulb(2*dtset%nfft,mmax,My_vtrial,BFGS_l2,BFGS_u2,BFGS_nbd,cgW,cgGCombine,  &
                   0._DP,0._DP,BFGS_wa2,BFGS_iwa2, &
                   task,BFGS_print,BFGS_csave,BFGS_lsave,BFGS_isave,BFGS_dsave)
       endif

       !================= set My_vtrial back to vtrial(:,:) ===================
       !================= and cgWCombine back to cgW and cgW2 =================
       

       if(bSpin==2) then
         do i=1,dtset%nfft
           vtrial(i,1) = My_vtrial(i)
           vtrial(i,2) = My_vtrial(i+dtset%nfft)
           cgG(i) = cgGCombine(i)
           cgG2(i) = cgGCombine(i+dtset%nfft)
         enddo
       endif
       !=======================================================================
!! debug        
!!       vtrial = vtrial - sum(vtrial)/dble(nfftf)
!! end of debug       

       CALL C_WRTLOG('=== CHECK THE FLAG ===')
       IF ( task(1:5) .EQ. 'NEW_X')  THEN 
         CALL c_wrtlog('TASK = NEW_X ')
         nline   = 0
         nBFGS   = nBFGS + 1
         last_x  = vtrial(:,1)        ! back up current vtrial
         if(bSpin==2) then
           last_x2 = vtrial(:,2)
         endif
         last_ke = ke_firstOrder      ! back up last ke
         last_fe = energies%e_fermie  ! back up fermi 
         write(message,'(a,I4,a,g22.14,a,ES12.4,a)') &
          '(scfcv) iter:',nBFGS,' cgW(with penalty)=',cgW,' penalty=', penfun,' (new iter)'
         call c_wrtlog(message)
         if (myicount+1 > maxcount) then
           call c_wrtlog('max bfgs iter reached. write out density and potential.')
           exit
         endif ! mycount > maxcount
         ! store current cgW value
         histW(1) = histW(2)
         histW(2) = histW(3)
         histW(3) = histW(4)
         histW(4) = cgW
         if (myicount > 5) then
           if ( abs(histW(4)-histW(3))<wtol .and.  & 
                abs(histW(3)-histW(2))<wtol .and.  & 
                abs(histW(2)-histW(1))<wtol) then 
             write(message,'(a,es12.4,a)') & 
               '(converged) In the last three BFGS iterations, cgW changes less than ', wtol,' hartree'
             call c_wrtlog(message)
             write(message,'(a,es20.12,a)') & 
               '(converged) final etotal = ',etotal,' hartree'
             call c_wrtlog(message)
             exit  ! jump out of BFGS loop
           endif
         endif
         goto 2012
       
       elseif (task(1:2) .EQ. 'FG') then
         ! NEED MORE LINE SEARCH
         if (nline>0) then 
           write(message,'(a,I3,a,g22.14,a,ES12.4,a)') &
            '(scfcv) nline:',nline,' cgW(with penalty)=',cgW,' penalty=', penfun,' '
           call c_wrtlog(message)
           call c_wrtlog('TASK = FG, continue line search')
         endif
!         if (nline>=5) then 
!           write(message,'(a)') & 
!             '(scfcv) max line search (5) reached , restart BFGS.'
!           call c_wrtlog(message)
!           task = 'START'
!           restart = .true. 
!           vtrial(:,1) = last_x
!           myicount = 0
!           nline    = 0 
!           goto 5555
!         endif
         nline = nline + 1
       else 
         call c_wrtlog('what is Task?, I should not get to this point, EXIT anyway, task=>')
         call c_wrtlog(task)
         exit
       endif 
       ! NEED MORE BFGS 
       myicount = myicount +1
       call c_wrtlog('=== NEXT BFGS ===');
       call load_param(bLoadDen,bLoadVr,maxcount,bAllowsym,bSpin,wtol,structure,ektol,penLambda)
   END DO  
   !------------------------------------
   ! end of BFGS loop
   !------------------------------------

!! debug 
!!  tigger_mode = .true.
!!  goto 5555
!! end of debug
   
!!=======================
!! OUTPUT DATA
!! ONLY MASTER NODE 
!!=======================       
   if (tigger_mode) then 

     call xbarrier_mpi(spaceworld) ! force synchronization

     if (mpi_enreg%me==0) then
       open(file='oep_results.dat',action='write',unit=111,form='unformatted')
       write(111)ke_firstOrder
       write(111)energies%e_nonlocalpsp
       write(111)energies%e_entropy
       write(111)-vtrial(:,1)+energies%e_fermie

       call c_wrtlog('[tigger_mode] tspot=(-vtrial+fermie) dumpped to ts.pot')
       close(111)
     
       ! WRITE MESSAGE FOR UOPT 
       open(file='message.oep',action='write',unit=111)
       write(111,'(a)')'done   end'
       close(111)
       call c_wrtlog('[tigger_mode] message.oep written')
       call c_wrtlog('[tigger_mode] go back to wait for next signal')
     endif  ! master node

     call xbarrier_mpi(spaceworld) ! force synchronization
     goto 5555  ! go back to the top,
                ! wait for the tigger again

   endif  ! TIGGER MODE

!!======================================
!! END OF TIGGER MODE, ROLL BACK NOW
!!======================================
   
    call c_wrtlog('write final POT DEN VHA VXC files...')

    ! xc and hartree potential ...
    call rhotov(dtset,energies,gprimd,gsqcut,   istep,  kxc,mpi_enreg,nfftf,ngfftf, &
&     nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,&
&     optene,optres,optxc,&
&     rhog,rhor,rprimd,strsxc,ucvol,psps%usepaw,usexcnhat,&
&     vhartr,vnew_mean,vpsp,vres_mean,res2,vtrial,vxcavg,vxc,   wvl,   xccc3d,  xred, &
&     electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau)

    ! output POT DEN files
    call outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,  dtefield,  dtfil,&
&    dtset,ecut,eigen,electronpositron,elfr,etotal,energies%e_fermie,&
&    gmet,gprimd,grhor,hdr,istep_mix,kg,lrhor,dtset%mband,  mcg,mcprj,  dtset%mgfft,&
&    dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,  my_natom,  dtset%natom,nattyp,&
&    nfftf,ngfftf,nhat,dtset%nkpt,npwarr,dtset%nspden,&
&    dtset%nsppol,dtset%nsym,psps%ntypat,n3xccc,occ,  paw_dmft,  pawang,pawfgr,pawfgrtab,&
&    pawrad,pawrhoij,pawtab,paw_an,paw_ij,dtset%prtvol,psps,&
&    rhor,rprimd,taur,ucvol,usecprj,wffnow,vhartr,vtrial,vxc,  wvl%den,  xccc3d,xred)

   open(unit=199,access="sequential", action="write",file='res_vion.dat', form="formatted", status="replace")   
   do i=1,dtset%nfft
     write(199, *)vtrial(i,1)-vxc(i,1)-vhartr(i)
   enddo 
   close(199)

   if(bSpin==2) then
     open(unit=299,access="sequential", action="write",file='res_vion2.dat', form="formatted", status="replace")   
     do i=1,dtset%nfft
       write(299, *)vtrial(i,2)-vxc(i,2)-vhartr(i)
     enddo 
     close(299)
   endif

   call c_wrtlog('')
   call c_wrtlog('-------------------------------------------------------')
   write (message, *) 'Kinetic Energy = ', energies%e_kinetic
   call c_wrtlog(message)
   write (message, *) 'Hartree Energy = ', energies%e_hartree
   call c_wrtlog(message)
   write (message, *) 'XC Energy      = ', energies%e_xc
   call c_wrtlog(message)
   write (message, *) 'Entropy        = ', energies%e_entropy
   call c_wrtlog('-------------------------------------------------------')

   call c_wrtlog('END AT:')
   call date_and_time(VALUES=dateTime)    !today
   write(msg, '(A,i2,A,i2,A,i4,A,i2.2,A,i2.2,A,i2.2,A,i3.3)') "", &
               dateTime(3),"/",  dateTime(2),"/",  dateTime(1)," : ",             &
               dateTime(5), ":", dateTime(6), ":", dateTime(7), ".", dateTime(8)
   call c_wrtlog(msg)
   exit

   ! ===================================================================================
   ! ======================================== END OF MY CODE ===========================
   ! ===================================================================================
     
   !  if(VERBOSE)then
   !    call wrtout(std_out,'*. Compute the density from the trial potential (vtorho)',"COLL")
   !  end if

   !  call status(istep,dtfil%filstat,iexit,level,'call vtorho   ')
   !  call vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&
!&     dielop,dielstrt,dphase,dtefield,dtfil,dtset,&
!&     eigen,electronpositron,energies,etotal,gbound_diel,&
!&     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
!&     istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
!&     my_natom,dtset%natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
!&     npwarr,npwdiel,res2,psps%ntypat,nvresid,occ,&
!&     computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
!&     pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
!&     pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,rmet,rprimd,&
!&     susmat,symrec,taug,taur,ucvol,wffnew,wffnow,vtrial,wvl,xred,&
!&     ylm,ylmgr,ylmdiel,vxctau=vxctau)
!     call status(istep,dtfil%filstat,iexit,level,'after vtorho  ')

   else if (dtset%tfkinfunc==1) then
     call wrtout(std_out,'WARNING : THOMAS FERMI',"COLL")
     call vtorhotf(dtfil,dtset,energies%e_kinetic,energies%e_nonlocalpsp,&
&     energies%entropy,energies%e_fermie,gprimd,grnl,irrzon,mpi_enreg,&
&     dtset%natom,nfftf,dtset%nspden,dtset%nsppol,dtset%nsym,phnons,&
&     rhog,rhor,rprimd,ucvol,vtrial)
     residm=zero
     energies%e_eigenvalues=zero
   end if

!  Recursion method
   if(dtset%userec==1)then
     call vtorhorec(dtset,&
&     energies%e_kinetic,energies%e_nonlocalpsp,energies%entropy,energies%e_eigenvalues,energies%e_fermie,&
&     grnl,initialized,irrzon,nfftf,phnons,&
&     rhog,rhor,vtrial,rec_set,istep-nstep,rprimd,gprimd)
     residm=zero
   end if

   if(dtset%wfoptalg==2)then
     do ikpt=1,dtset%nkpt
       shiftvector(1+(ikpt-1)*(dtset%mband+2))=val_min
       shiftvector(2+(ikpt-1)*(dtset%mband+2):ikpt*(dtset%mband+2)-1)=&
&       eigen((ikpt-1)*dtset%mband+1:ikpt*dtset%mband)
       shiftvector(ikpt*(dtset%mband+2))=val_max
     end do
   end if

   call timab(242,2,tsec)

!  ######################################################################
!  Skip out of step loop if non-SCF (completed)
!  ----------------------------------------------------------------------

!  Indeed, nstep loops have been done inside vtorho
   if (dtset%iscf<0) exit

!  ######################################################################
!  In case of density mixing or wavelet handling, compute the total energy
!  ----------------------------------------------------------------------
   call timab(60,1,tsec)
   if (dtset%iscf>=10 .or. dtset%usewvl == 1) then
     optene = 1 ! use double counting scheme (default)
     if (dtset%usewvl == 1 .and. dtset%iscf == 0) then
       optene = 0 ! use direct scheme for computation of energy
     else if (dtset%iscf == 22) then
       optene = -1
     end if

!    if the mixing is the ODA mixing, compute energy and new density here
     if (dtset%iscf==22) then
       call odamix(deltae,dtset,&
&       elast,energies,etotal,gprimd,gsqcut,kxc,dtefield%mag_cart,mpi_enreg,&
&       my_natom,nfftf,ngfftf,nhat,nkxc,psps%ntypat,nvresid,n3xccc,optres,&
&       paw_ij,paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
&       red_ptot,psps,rhog,rhor,rprimd,strsxc,taug,taur,ucvol,psps%usepaw,&
&       usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred)
     end if
!    If the density mixing is required, compute the total energy here
     call etotfor(atindx1,deltae,diffor,dtset,&
&     elast,electronpositron,energies,&
&     etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,grvdw,&
&     grxc,gsqcut,indsym,kxc,dtefield%mag_cart,maxfor,mgfftf,mpi_enreg,my_natom,&
&     nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,psps%ntypat,nvresid,n1xccc,n3xccc,&
&     optene,computed_forces,optres,pawang,pawfgrtab,pawrhoij,pawtab,&
&     ph1df,red_ptot,psps,rhog,rhor,rprimd,symrec,synlgr,&
&     psps%usepaw,vhartr,vpsp,vxc,wvl%descr,wvl%den,xccc3d,xred)
   end if
   call timab(60,2,tsec)

!  ######################################################################
!  In case of density mixing, check the exit criterion
!  ----------------------------------------------------------------------
   if (dtset%iscf>=10 .or. (dtset%usewvl == 1 .And. dtset%iscf > 0)) then
!    Check exit criteria
     call timab(52,1,tsec)
     call status(istep,dtfil%filstat,iexit,level,'call scprqt   ')
     choice=2
     if(paw_dmft%use_dmft==1) call prtene(dtset,energies,std_out,psps%usepaw)
     call scprqt(choice,cpus,deltae,diffor,dtset,&
&     eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
&     dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
&     maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
&     occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
&     psps%usepaw,vxcavg,dtset%wtk,xred,electronpositron=electronpositron)
!    Exit criteria for the recursion method
     if(dtset%userec==1.and.rec_set%quitrec==2)quit=1

     if (istep==nstep) quit=1

!    If criteria in scprqt say to quit, then exit the loop over istep.

     quit_sum=quit
     call xsum_mpi(quit_sum,spaceComm,ierr)
     if (quit_sum > 0) quit=1


     call timab(52,2,tsec)

     if (quit==1) exit
   end if

!  ######################################################################
!  Mix the total density (if required)
!  ----------------------------------------------------------------------
   call timab(68,1,tsec)

   if (dtset%iscf>=10 .and.dtset%iscf/=22.and. dtset%usewvl == 0) then
!    If LDA dielectric matrix is used for preconditionning, has to update here Kxc
     if (nkxc>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79) &
&     .and.((istep==1.or.istep==dielstrt).or.(dtset%iprcel>=100))) then
       optxc=10

!      to be adjusted for the call to rhohxc
       nk3xc=1
       if(dtset%usekden/=0)then
         call rhohxc(dtset,edum,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,0,nkxc,nk3xc,dtset%nspden,&
&         n3xccc,optxc,rhog,rhor,rprimd,dummy2,0,vhartr,vxc,vxcavg_dum,&
&         xccc3d,taug=taug,taur=taur,vxctau=vxctau)
       else
         call rhohxc(dtset,edum,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,0,nkxc,nk3xc,dtset%nspden,&
&         n3xccc,optxc,rhog,rhor,rprimd,dummy2,0,vhartr,vxc,vxcavg_dum,&
&         xccc3d,taug=taug,taur=taur)
       end if
     end if

     call status(istep,dtfil%filstat,iexit,level,'call newrho   ')
     call newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,&
&     dtset,etotal,fcart,pawfgr%fintocoa,&
&     gmet,grhf,gsqcut,initialized,ispmix,istep_mix,kg_diel,kxc,&
&     mgfftf,mix,pawfgr%coatofin,moved_atm_inside,mpi_enreg,my_natom,nattyp,nfftf,&
&     nfftmix,ngfftf,ngfftmix,nkxc,npawmix,npwdiel,nvresid,psps%ntypat,&
&     n1xccc,pawrhoij,pawtab,ph1df,psps,rhog,rhor,&
&     rprimd,susmat,psps%usepaw,vtrial,wvl%descr,wvl%den,xred)
   end if   ! iscf>=10

   call timab(68,2,tsec)

!  ######################################################################
!  Additional computation in case of an electric field or electric displacement field
!  ----------------------------------------------------------------------

!  In case of an electric field calculation, need polarization
!  to compute electric enthalpy instead of energy.

!  In case of an electric displacement field calculation, need polarization
!  to compute internal energy instead of energy.

   call timab(239,1,tsec)

!  if (psps%usepaw==0.and.dtset%berryopt == 4) then
   if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7 .or.  &
&   dtset%berryopt == 14 .or. dtset%berryopt ==16 .or. dtset%berryopt == 17) then  !!HONG
!    When using symmetry, it is costly to update polarization from changes in Zak
!    phases. It is better to call berryphase here.
!    Update polarization by adding increment from the SCF step
!    pel_cg(1:3) = pel_cg(1:3) + dtefield%sdeg*dphase(1:3)/two_pi
!    ptot(1:3) = pel_cg(1:3) + pion(1:3)
!    write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
!    &    ' scfcv: Polarization from accumulated change in Berry phase:',ch10,&
!    &    ' (reduced coordinates, a. u., without correcting for branch cuts)',ch10,&
!    &    '     Electronic: ', (pel_cg(ii), ii = 1, 3), ch10,&
!    &    '     Ionic:      ', (pion(ii), ii = 1, 3), ch10, &
!    &    '     Total:      ', (ptot(ii), ii = 1, 3)
!    call wrtout(std_out,message,'COLL')
!    berrys phase with PAW needs cprj
     mcprj=0
     if (psps%usepaw==1) then
       usecprj=1;my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
       mcprj=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
       ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
       ncpgr = 3
       ctocprj_choice = 2 ! compute cprj and derivatives wrt atomic position
       call cprj_alloc(cprj,ncpgr,dimcprj)
       iatom=0 ; iorder_cprj=1 ! cprj are not ordered
       call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,&
&       iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&       dtset%mband,mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&       dtset%mpw,dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,&
&       dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,&
&       dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&       ucvol,dtfil%unpaw,dtfil%unkg,dtfil%unylm,0,wffnow,xred,ylm,ylmgr)
     end if
     call berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,gprimd,hdr,psps%indlmn,kg,&
&     psps%lmnmax,dtset%mband,mcg,mcprj,dtset%mkmem,mpi_enreg,dtset%mpw,my_natom,dtset%natom,npwarr,&
&     dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawrhoij,pawtab,pel_cg,pelev,pion,ptot,red_ptot,pwind,&  !!REC
&    pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&     unit_out,usecprj,psps%usepaw,wffnow,xred,psps%ziontypat)

     dtefield%red_ptot1(:)=red_ptot(:)
!    !REC start take ptot fom berrypahse
!    red_ptot(:) = pel_cg(:) + pion(:) + pelev(:)
!    above line suppressed because in the PAW case, pel_cg already includes all on-site
!    terms and pelev should not be added in additionally. We are computing pelev separately for
!    reporting purposes only.
!    13 June 2012 J Zwanziger
!    !     red_ptot(:) = pel_cg(:) + pion(:)
!    above line is suppressed because red_ptot is already
!    ! calculated in berryphase_new.F90 and remapped in [-1,1] (and other operation...).
!    ! 25 June 2012 Jiawang Hong

     if (dtset%prtvol >= 10)then
       write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
&       ' scfcv: New value of the polarization:',ch10,&
&       ' (reduced coordinates, a. u.)',ch10,&
&       '     Electronic berry phase:       ', (pel_cg(ii), ii = 1, 3)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
       if(psps%usepaw==1) then
         write(message,'(a,3(e16.9,2x))')&
&         '     ...includes PAW on-site term: ', (pelev(ii), ii = 1, 3)
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end if
       write(message,'(a,3(e16.9,2x),a,a,3(e16.9,2x))')&
&       '     Ionic:                        ', (pion(ii), ii = 1, 3), ch10, &
&       '     Total:                        ', (red_ptot(ii), ii = 1, 3) !!REC
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
     if (psps%usepaw==1) then
!      pel_cg(:) = pel_cg(:) + pelev(:) ! electronic polarization includes on-site term
!      above line suppressed because in the PAW case, pel_cg already includes all on-site
!      terms and pelev should not be added in additionally. We are computing pelev separately for
!      reporting purposes only.
!      13 June 2012 J Zwanziger
       usecprj=0
       call cprj_free(cprj)
       ABI_DATATYPE_DEALLOCATE(cprj)
     end if

!    ! Convert polarization to cartesian coords
     ptot_cart(:)=zero
     do iidir = 1,3
       ptot_cart(iidir)=rprimd(iidir,1)*red_ptot(1) + rprimd(iidir,2)*red_ptot(2) + &
&       rprimd(iidir,3)*red_ptot(3)
     end do
     ptot_cart(:)=ptot_cart(:)/ucvol


!    !!

!    !===================================================================================================
!    !                                       OUTPUT  for fixed E
!    !===================================================================================================

     if (dtset%berryopt == 4) then

!      !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
       write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced E-field:'
       call wrtout(std_out,message,'COLL')
       call prtefield(dtset,dtefield,std_out,rprimd)
       if(dtset%prtvol>=10)then
         call wrtout(ab_out,message,'COLL')
         call prtefield(dtset,dtefield,ab_out,rprimd)
       end if
     end if

!    =====================================================================================
!    !                                      fixed D calculation
!    !====================================================================================
     if (dtset%berryopt == 6) then
       if (istep > 1) then
!        ! update efield taking damping into account dfield is in cartesian in dtset structure (contains input value)
!        ! same goes for efield - update the dtset%efield value
         efield_test_cart(:)=dtset%ddamp*(dtset%dfield(:)-4.0d0*pi*ptot_cart(:))+&
&         (1.0d0-dtset%ddamp)*efield_old_cart(:)

!        ! test whether change in efield in any direction exceed maxestep, if so, set the
!        ! change to maxestep instead   ! need optimized !
         do idir = 1,3

           if (dabs(efield_test_cart(idir)-efield_old_cart(idir)) > dabs(dtset%maxestep)) then

             write(std_out,'(a,a,i5)') "JH - ","  E-field component:",idir
             write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(idir), &
&             ",    E(n-1)=",efield_old_cart(idir), ",    E(n)-E(n-1)=", efield_test_cart(idir)-efield_old_cart(idir), &
&             ",    maxestep=",dtset%maxestep


             if (efield_test_cart(idir) > efield_old_cart(idir)) then
               efield_test_cart(idir) = efield_old_cart(idir) + dabs(dtset%maxestep)
             else
               efield_test_cart(idir) = efield_old_cart(idir) - dabs(dtset%maxestep)
             end if
           end if
         end do


         dtset%efield(:) = efield_test_cart(:)

!        !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
         write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced D-field  - updating E-field:'
         call wrtout(std_out,message,'COLL')
         call prtefield(dtset,dtefield,std_out,rprimd)
         if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
         end if

!        ! need to update dtset%efield_dot(:) with new value
         dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
         dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
         dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

       else

         write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced D-field  - Pre E-field:'
         call wrtout(std_out,message,'COLL')
         call prtefield(dtset,dtefield,std_out,rprimd)
         if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
         end if

       end if  ! istep >1

       efield_old_cart(:)=dtset%efield(:)
     end if  ! berryopt ==6

!    !===================================================================================================
!    !                                      fixed reduced d calculation
!    !===================================================================================================
     if (dtset%berryopt == 16) then

       if (istep > 1) then
!        ! update efield taking damping into account reduced red_dfield
!        red_efield2 is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009)

         red_efield2(:)=dtset%ddamp*(dtset%red_dfield(:)-red_ptot(:))+ (1.0d0-dtset%ddamp)*red_efield2_old(:)

!        to calculate unreduced E
         efield_test_cart(:)=(4*pi/ucvol)*(rprimd(:,1)*red_efield2(1)+rprimd(:,2)*red_efield2(2)+rprimd(:,3)*red_efield2(3))

!        ! test whether change in efield in any direction exceed maxestep, if so, set the
!        ! change to maxestep instead   ! need optimized !
         do idir = 1,3
           if (dabs(efield_test_cart(idir)-efield_old_cart(idir)) > dabs(dtset%maxestep)) then

             write(std_out,'(a,a,i5)') "JH - ","  E-field component:",idir
             write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(idir), &
&             ",    E(n-1)=",efield_old_cart(idir), ",    E(n)-E(n-1)=", efield_test_cart(idir)-efield_old_cart(idir), &
&             ",    maxestep=",dtset%maxestep

             if (efield_test_cart(idir) > efield_old_cart(idir)) then
               efield_test_cart(idir) = efield_old_cart(idir) + dabs(dtset%maxestep)
             else
               efield_test_cart(idir) = efield_old_cart(idir) - dabs(dtset%maxestep)
             end if
           end if
         end do

         dtset%efield(:) = efield_test_cart(:)

!        !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
         write(message,'(a,a)')   ch10, 'scfcv: Constant reduced d-field  - updating E-field:'
         call wrtout(std_out,message,'COLL')
         call prtefield(dtset,dtefield,std_out,rprimd)
         if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
         end if

!        ! need to update dtset%efield_dot(:) with new value
!        ! This needs to be deleted  when efield_dot is deleted
         dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
         dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
         dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

       else

         write(message,'(a,a)')   ch10, 'scfcv: Constant reduced d-field  - Pre E-field:'
         call wrtout(std_out,message,'COLL')
         call prtefield(dtset,dtefield,std_out,rprimd)
         if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
         end if

       end if  ! istep > 1

       efield_old_cart(:)=dtset%efield(:)
       red_efield2_old(:)=red_efield2(:)
     end if  ! berryopt ==16


!    !===================================================================================================
!    !                                      fixed reduced d and ebar calculation (mixed BC)
!    !===================================================================================================
     if (dtset%berryopt == 17) then

       if (istep > 1) then
!        ! update efield taking damping into account reduced red_dfield
!        red_efield1 and red_efield2 is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009)
!        red_efield1 for fixed ebar, red_efield2 for fixed d calculation

!        save this value in order to print the final value of real electric field, comparing with the desired red_fieldbar
         dtefield%efield2(:)=dtset%efield(:)

!        write(*,'(a,3i4)') "jfielddir=", (dtset%jfielddir(ii),ii=1,3)

         do idir=1,3
           if (dtset%jfielddir(idir) ==2 ) then    ! direction under fixed d
!            red_efieldbar_lc(idir) = dot_product(dtset%efield(:),rprimd(:,idir))
             dtset%red_efieldbar(idir) = dot_product(dtset%efield(:),rprimd(:,idir)) !  update ebar which is not fixed
             dtefield%efield_dot(idir) = dot_product(dtset%efield(:),rprimd(:,idir))
             red_efield2(idir)=dtset%ddamp*(dtset%red_dfield(idir) - red_ptot(idir)) +  &
&             (1.0d0-dtset%ddamp)*red_efield2_old(idir)         ! d(idir) is fixed, update e(idir)  may need ddamping here

!            write(message,'(a,a,i5,a,i5)')   ch10, 'direction  ', idir,'   for fixed d, value is (2)  ', dtset%jfielddir(idir)
!            call wrtout(ab_out,message,'COLL')
!            call wrtout(std_out,message,'COLL')

           else if (dtset%jfielddir(idir) ==1 ) then   ! direction under fixed ebar
             red_efield2(idir)= (ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir)) !  update e which is not fixed
             dtset%red_dfield(idir)=red_ptot(idir) +  (ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir))  ! update d

!            write(message,'(a,a,i5,a,i5)')   ch10, 'direction  ', idir,'   for fixed ebar, value is (1)  ', dtset%jfielddir(idir)
!            call wrtout(ab_out,message,'COLL')
!            call wrtout(std_out,message,'COLL')

           end if
         end do

         do idir=1,3
           red_efield1(idir)  =(ucvol/(4*pi))*dot_product(dtset%red_efieldbar(:),gmet(:,idir))
         end do


         dtset%red_efield(:)=(red_efield1(:) + red_efield2(:))/2.0d0 ! average reduced efield, one is from fixed ebar part, the other is from fixed d part. This may need to be optimized !!

         write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Reduced efield from fixed ebar:', ch10, &
&         '       e:  ', (red_efield1(ii),ii=1,3), ch10

!        call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')

         write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Reduced efield from fixed d:', ch10, &
&         '       e:  ', (red_efield2(ii),ii=1,3), ch10

!        call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')

         write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Average reduced efield:', ch10, &
&         '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10

!        call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')

!        to calculate unreduced E
         do idir=1,3
           efield_test_cart(idir)  = (4*pi/ucvol)* dot_product(dtset%red_efield(:),rprimd(:,idir))
         end do

!        ! test whether change in efield in any direction exceed maxestep, if so, set the
!        ! change to maxestep instead   ! need optimized !
         do idir = 1,3
           if (dabs(efield_test_cart(idir)-efield_old_cart(idir)) > dabs(dtset%maxestep)) then

             write(std_out,'(a,a,i5)') "JH - ","  E-field component:",idir
             write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(idir), &
&             ",    E(n-1)=",efield_old_cart(idir), ",    E(n)-E(n-1)=", efield_test_cart(idir)-efield_old_cart(idir), &
&             ",    maxestep=",dtset%maxestep

             if (efield_test_cart(idir) > efield_old_cart(idir)) then
               efield_test_cart(idir) = efield_old_cart(idir) + dabs(dtset%maxestep)
             else
               efield_test_cart(idir) = efield_old_cart(idir) - dabs(dtset%maxestep)
             end if
           end if
         end do

         dtset%efield(:) = efield_test_cart(:)

!        !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
         write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar and d-field  - updating E-field:'
         call wrtout(std_out,message,'COLL')
         call prtefield(dtset,dtefield,std_out,rprimd)
         if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
         end if


!        ! need to update dtset%efield_dot(:) with new value
!        ! This needs to be deleted  when efield_dot is deleted
         dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
         dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
         dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

       else

         write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar and d-field  - Pre E-field:'
         call wrtout(std_out,message,'COLL')
         call prtefield(dtset,dtefield,std_out,rprimd)
         if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
         end if

       end if  ! istep > 1

       efield_old_cart(:)=dtset%efield(:)
       red_efield2_old(:)=red_efield2(:)

     end if  ! berryopt ==17


   end if       ! berryopt 4/6/7/14/16/17
   call timab(239,2,tsec)

!  ######################################################################
!  Compute the new potential from the trial density
!  ----------------------------------------------------------------------

   call timab(243,1,tsec)
   if(VERBOSE)then
     call wrtout(std_out,'*. Compute the new potential from the trial density',"COLL")
   end if

!  Set XC computation flag
   optxc=1
   if (nkxc>0) then
     if (dtset%nfreqsus>0) optxc=2
     if (dtset%iscf<0) optxc=2
     if (modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79).and. &
&     dtset%iscf<10.and. &
&     (dtset%iprcel>=100.or.istep==1.or.istep==dielstrt)) optxc=2
     if (dtset%iscf>=10.and.dtset%iprcch/=0.and.abs(dtset%iprcch)/=5) optxc=2
     if (optxc==2.and.dtset%xclevel==2.and.nkxc==3-2*mod(dtset%nspden,2)) optxc=12
   end if

   if (dtset%iscf/=22) then
!    PAW: eventually recompute compensation density (and gradients)
     nhatgrdim=0
     if (psps%usepaw==1) then
       ider=-1;if (dtset%iscf>=10.and.((dtset%xclevel==2.and.dtset%pawnhatxc>0).or.usexcnhat==0)) ider=0
       if (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0) ider=ider+2
       if (ipositron==1) ider=-1
       if (ider>0) then
         nhatgrdim=1
         ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3))
       end if
       if (ider>=0) then
         call timab(558,1,tsec)
         izero=0
         call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,nfftf,ngfftf,&
&         nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,nhatgr,nhat,&
&         pawrhoij,pawrhoij,pawtab,k0,rprimd,ucvol,xred,&
&         mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&         mpi_comm_fft=spaceComm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0)
         call timab(558,2,tsec)
       end if
     end if

!    Compute new potential from the trial density
     call status(istep,dtfil%filstat,iexit,level,'call rhotov')
     optene=2*optres;if(psps%usepaw==1) optene=2
     call rhotov(dtset,energies,gprimd,gsqcut,istep,kxc,mpi_enreg,nfftf,ngfftf, &
&     nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,&
&     optene,optres,optxc,&
&     rhog,rhor,rprimd,strsxc,ucvol,psps%usepaw,usexcnhat,&
&     vhartr,vnew_mean,vpsp,vres_mean,res2,vtrial,vxcavg,vxc,wvl,xccc3d,xred,&
&     electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau)
   end if

   call timab(243,2,tsec)
   call timab(60,1,tsec)

!  This is inside the loop, its not equivalent to the line 1621
   if(moved_atm_inside==1) xred_old(:,:)=xred(:,:)

   if (dtset%iscf<10 .and. dtset%usewvl == 0) then

     if(VERBOSE)then
       call wrtout(std_out,'*. Check exit criteria in case of potential mixing',"COLL")
     end if

!    If the potential mixing is required, compute the total energy here
!    PAW: has to compute here spherical terms
     if (psps%usepaw==1) then
       nzlmopt=0;if (istep_mix==1.and.dtset%pawnzlm>0) nzlmopt=-1
       if (istep_mix>1) nzlmopt=dtset%pawnzlm
       option=2
       do iatom=1,my_natom
         ABI_ALLOCATE(paw_ij(iatom)%dijhartree,(paw_ij(iatom)%lmn2_size))
         paw_ij(iatom)%has_dijhartree=1
       end do
       call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,&
&       dtset%ixc,my_natom,dtset%natom,dtset%nspden,&
&       psps%ntypat,nzlmopt,option,paw_an,paw_an,&
&       paw_ij,pawang,dtset%pawprtvol,pawrad,pawrhoij,dtset%pawspnorb,&
&       pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,psps%znuclpsp,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       electronpositron=electronpositron)
       do iatom=1,my_natom
         ABI_DEALLOCATE(paw_ij(iatom)%dijhartree)
         paw_ij(iatom)%has_dijhartree=0
       end do
     end if
     call status(istep,dtfil%filstat,iexit,level,'call etotfor  ')
     call etotfor(atindx1,deltae,diffor,dtset,&
&     elast,electronpositron,energies,&
&     etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,grvdw,&
&     grxc,gsqcut,indsym,kxc,dtefield%mag_cart,maxfor,mgfftf,mpi_enreg,my_natom,&
&     nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,dtset%ntypat,nvresid,n1xccc, &
&     n3xccc,0,computed_forces,optres,pawang,pawfgrtab,pawrhoij,&
&     pawtab,ph1df,red_ptot,psps,rhog,rhor,rprimd,symrec,synlgr,&
&     psps%usepaw,vhartr,vpsp,vxc,wvl%descr,wvl%den,xccc3d,xred)
   end if
   call timab(60,2,tsec)

!  ######################################################################
!  Check exit criteria in case of potential mixing or wavelet handling
!  ----------------------------------------------------------------------
   if (dtset%iscf<10 .or. (dtset%usewvl == 1 .and. dtset%iscf == 0)) then
!    Check exit criteria
     call timab(52,1,tsec)
     call status(istep,dtfil%filstat,iexit,level,'call scprqt   ')
     choice=2
     call scprqt(choice,cpus,deltae,diffor,dtset,&
&     eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
&     dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
&     maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
&     occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
&     psps%usepaw,vxcavg,dtset%wtk,xred,electronpositron=electronpositron)
     if (istep==nstep.and.psps%usepaw==1) quit=1
     call timab(52,2,tsec)

     call timab(244,1,tsec)
!    exit criteria for the recursion
     if(dtset%userec==1 .and. rec_set%quitrec==2) quit=1

!    If criteria in scprqt say to quit, then exit the loop over istep.
     quit_sum=quit
     call xsum_mpi(quit_sum,spaceComm,ierr)
     if (quit_sum > 0) quit=1

     if (quit==1) then
       do ispden=1,dtset%nspden
         vtrial(:,ispden)=vtrial(:,ispden)+nvresid(:,ispden)+vres_mean(ispden)
       end do
       call timab(244,2,tsec) ! Due to the exit instruction, two timab calls are needed
       exit ! exit the loop over istep
     end if
     call timab(244,2,tsec) ! Due to the exit instruction, two timab calls are needed
   end if

!  ######################################################################
!  Mix the potential (if required) - Check exit criteria
!  ----------------------------------------------------------------------

   call timab(245,1,tsec)
   if (dtset%iscf<10 .and. dtset%usewvl == 0) then

     if(VERBOSE)then
       call wrtout(std_out,'*. Mix the potential (if required) - Check exit criteria',"COLL")
     end if

!    Precondition the residual and forces, then determine the new vtrial
!    (Warning: the (H)xc potential may have been subtracted from vtrial)
     call status(istep,dtfil%filstat,iexit,level,'call newvtr   ')
     call newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&
&     dtn_pc,dtset,energies%e_fermie,etotal,fcart,pawfgr%fintocoa,&
&     gmet,grhf,gsqcut,initialized,ispmix,&
&     istep_mix,kg_diel,kxc,mgfftf,mix,pawfgr%coatofin,&
&     moved_atm_inside,mpi_enreg,my_natom,nattyp,nfftf,nfftmix,&
&     nhat,nhatgr,nhatgrdim,ngfftf,ngfftmix,nkxc,npawmix,npwdiel,&
&     nstep,psps%ntypat,n1xccc,optres,optxc,&
&     pawrhoij,ph1df,psps,rhor,rprimd,susmat,psps%usepaw,&
&     vhartr,vnew_mean,vpsp,nvresid,vtrial,vxc,xred,atindx1,cg,deltae,&
&     dtfil,energies%e_eigenvalues,eigen,energies%e_kinetic,&
&     energies%e_nonlocalpsp,kg,mcg,nfftf,ngfftf,npwarr,n3xccc,occ,optene,&
&     pawfgr,pawtab,resid,rhog,usexcnhat,wffnow,wvl,ylm,xccc3d)
   end if   ! iscf<10

!  ######################################################################
!  END MINIMIZATION ITERATIONS
!  ######################################################################

   if(VERBOSE)then
     call wrtout(std_out,'*. END MINIMIZATION ITERATIONS',"COLL")
   end if

!  The initialisation of the gstate run should be done when this point is reached
   initialized=1

!  This is to save the density for restart.
   if (xmpi_paral==0.or.me==0) then
     prtden=dtset%prtden
     prtkden=dtset%prtkden
     if(prtden<0.or.prtkden<0) then
!      Update the content of the header (evolving variables)
       bantot=hdr%bantot
       if (dtset%positron==0) then
         call hdr_update(bantot,etotal,energies%e_fermie,hdr,dtset%natom,residm,rprimd,occ,&
&         pawrhoij,psps%usepaw,xred,&
&         mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
       else
         call hdr_update(bantot,electronpositron%e0,energies%e_fermie,hdr,dtset%natom,residm,rprimd,occ,&
&         pawrhoij,psps%usepaw,xred,&
&         mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
       end if
     end if

     if (prtden<0) then
       if (mod(istep-1,abs(prtden))==0) then
         isave=isave+1
         call status(0,dtfil%filstat,iexit,level,'call ioarr-den')
         rdwr=2 ; fformr=52 ; rdwrpaw=0
         call int2char4(mod(isave,2),tag)
         fildata=trim(dtfil%fnametmp_app_den)//'_'//tag
         accessfil = 0
         call ioarr(accessfil,rhor, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&         nfftf,pawrhoij,rdwr,rdwrpaw,wvl%den)
       end if
     end if
     if (prtkden<0) then
       if (mod(istep-1,abs(prtkden))==0) then
         isave=isave+1
         call status(0,dtfil%filstat,iexit,level,'call ioarr-kden')
         rdwr=2 ; fformr=52 ; rdwrpaw=0
         call int2char4(mod(isave,2),tag)
         fildata=trim(dtfil%fnametmp_app_kden)//'_'//tag
         accessfil = 0
         call ioarr(accessfil,taur, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&         nfftf,pawrhoij,rdwr,rdwrpaw,wvl%den)
       end if
     end if
   end if

   if (nhatgrdim>0)  then
     ABI_DEALLOCATE(nhatgr)
   end if

   istep_mix=istep_mix+1
   if (ipositron/=0) electronpositron%istep_scf=electronpositron%istep_scf+1

   call timab(245,2,tsec)

 end do ! istep

 call timab(246,1,tsec)

 if (dtset%iscf > 0) then
   call ab6_mixing_deallocate(mix)
 end if

 if (quit==1.and.nstep==1) initialized=1

!######################################################################
!Case nstep==0: compute energy based on incoming wf
!----------------------------------------------------------------------

 if(nstep==0) then
   optene=2*psps%usepaw+optres
   energies%entropy=results_gs%energies%entropy  !MT20070219: entropy is not recomputed in routine energy
   call energy(atindx,atindx1,cg,compch_fft,dtfil,dtset,electronpositron,&
&   energies,eigen,etotal,gsqcut,indsym,irrzon,kg,mcg,mpi_enreg,my_natom,nattyp,&
&   nfftf,ngfftf,nhat,nhatgr,nhatgrdim,npwarr,n3xccc,&
&   occ,optene,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,&
&   phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,taug,taur,&
&   usexcnhat,vhartr,vtrial,vpsp,vxc,wffnow,wvl%wfs,wvl%den,xccc3d,xred,ylm,vxctau=vxctau)
   if (nhatgrdim>0)  then
     ABI_DEALLOCATE(nhatgr)
   end if
 end if ! nstep==0

!######################################################################
!Additional steps after SC iterations, including force, stress, polarization calculation
!----------------------------------------------------------------------

 if (dtset%userec==1) then
   call prtene(dtset,energies,ab_out,psps%usepaw)
   call prtene(dtset,energies,std_out,psps%usepaw)
 end if
 call status(0,dtfil%filstat,iexit,level,'endloop istep ')

!PAW: in some cases, need to recompute <p_lmn|Cnk> projected WF:
!should be output from vtorho (but is "type-sorted" in vtorho, not here)...
 recompute_cprj = psps%usepaw    ==1 .and. &
& (dtset%prtwant  ==2  .or. &
& dtset%prtwant  ==3  .or. &
& dtset%prtnabla > 0  .or. &
& dtset%prtdos   ==3  .or. &
& dtset%berryopt /=0  .or. &
& dtset%kssform  ==3  .or. &
& dtset%pawfatbnd> 0  .or. &
& dtset%pawprtwf > 0)

 if (recompute_cprj) then
   usecprj=1;my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
   mcprj=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
   ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
   if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7 .or. &
&   dtset%berryopt == 14 .or. dtset%berryopt ==16 .or. dtset%berryopt == 17) then ! finite efield needs gradients for forces
     ncpgr = 3
     ctocprj_choice = 2
   else
     ncpgr = 0
     ctocprj_choice = 1
   end if
   call cprj_alloc(cprj,ncpgr,dimcprj)
   iatom=0 ; iorder_cprj=1 ! cprj are not ordered
   call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,&
&   iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&   dtset%mband,mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&   dtset%mpw,dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,&
&   dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&   ucvol,dtfil%unpaw,dtfil%unkg,dtfil%unylm,0,wffnow,xred,ylm,ylmgr)
 else
   mcprj=0;usecprj=0
 end if

 call timab(246,2,tsec)
 call timab(247,1,tsec)

!SHOULD CLEAN THE ARGS OF THIS ROUTINE
 call status(0,dtfil%filstat,iexit,level,'afterscfloop  ')
 call afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&
& deltae,diffor,dtefield,dtfil,dtset,eigen,electronpositron,elfr,&
& energies,etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grhor,grvdw,&
& grxc,gsqcut,hdr,indsym,irrzon,istep,kg,kxc,lrhor,maxfor,mcg,mcprj,mgfftf,&
& moved_atm_inside,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfft,ngfftf,ngrvdw,nhat,&
& nkxc,npwarr,nvresid,occ,optres,optxc,paw_an,paw_ij,pawang,pawfgr,&
& pawfgrtab,pawrhoij,pawtab,pel,pel_cg,ph1d,ph1df,phnons,pion,prtfor,&
& prtxml,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&
& rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,taug,&
& taur,tollist,usecprj,usexcnhat,vhartr,vpsp,vxc,vxcavg,wffnow,wvl,&
& xccc3d,xred,ylm,ylmgr)

!Before leaving the present routine, save the current value of xred.
 xred_old(:,:)=xred(:,:)

 call timab(247,2,tsec)

!######################################################################
!All calculations in scfcv are finished. Printing section
!----------------------------------------------------------------------

 call timab(248,1,tsec)

 call status(istep,dtfil%filstat,iexit,level,'call outscfcv ')

 call outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtefield,dtfil,&
& dtset,ecut,eigen,electronpositron,elfr,etotal,energies%e_fermie,&
& gmet,gprimd,grhor,hdr,istep_mix,kg,lrhor,dtset%mband,mcg,mcprj,dtset%mgfft,&
& dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,my_natom,dtset%natom,nattyp,&
& nfftf,ngfftf,nhat,dtset%nkpt,npwarr,dtset%nspden,&
& dtset%nsppol,dtset%nsym,psps%ntypat,n3xccc,occ,paw_dmft,pawang,pawfgr,pawfgrtab,&
& pawrad,pawrhoij,pawtab,paw_an,paw_ij,dtset%prtvol,psps,&
& rhor,rprimd,taur,ucvol,usecprj,wffnow,vhartr,vtrial,vxc,wvl%den,xccc3d,xred)

 if(associated(elfr))then
   ABI_DEALLOCATE(elfr)
   nullify(elfr)
 end if

 if(associated(grhor))then
   ABI_DEALLOCATE(grhor)
   nullify(grhor)
 end if

 if(associated(lrhor))then
   ABI_DEALLOCATE(lrhor)
   nullify(lrhor)
 end if

 if(dtset%prtkden/=0 .or. dtset%prtelf/=0)then
   ABI_DEALLOCATE(taur)
 end if

 call timab(248,2,tsec)

 call timab(61,1,tsec)
 call leave_test()
 call timab(61,2,tsec)

 call timab(249,1,tsec)

!Transfer eigenvalues and occupation computed by BigDFT in afterscfloop to eigen.
#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 1) then
   if (dtset%nsppol == 1) then
     eigen = wvl%wfs%ks%orbs%eval
     occ = wvl%wfs%ks%orbs%occup
   else
     eigen(1:wvl%wfs%ks%orbs%norbu) = wvl%wfs%ks%orbs%eval(1:wvl%wfs%ks%orbs%norbu)
     eigen(dtset%mband + 1:dtset%mband + wvl%wfs%ks%orbs%norbd) = &
&     wvl%wfs%ks%orbs%eval(wvl%wfs%ks%orbs%norbu + 1:wvl%wfs%ks%orbs%norb)
     occ(1:wvl%wfs%ks%orbs%norbu) = wvl%wfs%ks%orbs%occup(1:wvl%wfs%ks%orbs%norbu)
     occ(dtset%mband + 1:dtset%mband + wvl%wfs%ks%orbs%norbd) = &
&     wvl%wfs%ks%orbs%occup(wvl%wfs%ks%orbs%norbu + 1:wvl%wfs%ks%orbs%norb)
   end if
 end if
#endif

!Debugging : print the different parts of rhor, as well as vxc
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(dtset%prtvol==-level)then
   write(message,'(a)') '   ir     vxc(ir)     rhor(ir)     '
   call wrtout(std_out,message,'COLL')
   do ir=1,nfftf
     if(ir<=11 .or. mod(ir,301)==0 )then
       write(message,'(i5,a,2es13.6)')ir,' ',vxc(ir,1),rhor(ir,1)
       call wrtout(std_out,message,'COLL')
       if(dtset%nspden==2)then
         write(message,'(a,2es13.6)')'      ',vxc(ir,2),rhor(ir,2)
         call wrtout(std_out,message,'COLL')
       end if
     end if
   end do
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(dtset%prtvol==-level)then
   write(message,'(a,i0,a)')' scfcv : exit. prtvol=-',level,', debugging mode => stop '
   MSG_ERROR(message)
 end if

!######################################################################
!Deallocate memory and save results
!----------------------------------------------------------------------
 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

 call prc_mem_free()

 if (recompute_cprj) then
   usecprj=0
   call cprj_free(cprj)
   ABI_DATATYPE_DEALLOCATE(cprj)
 end if

 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(fred)
 ABI_DEALLOCATE(forold)
 ABI_DEALLOCATE(grnl)
 ABI_DEALLOCATE(gresid)
 ABI_DEALLOCATE(grewtn)
 ABI_DEALLOCATE(grvdw)
 ABI_DEALLOCATE(grxc)
 ABI_DEALLOCATE(synlgr)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(vhartr)
 ABI_DEALLOCATE(vtrial)
 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(vxc)
 ABI_DEALLOCATE(vxctau)
 ABI_DEALLOCATE(xccc3d)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(shiftvector)
 if (dtset%iscf>=0) then
   ABI_DEALLOCATE(dtn_pc)
   ABI_DEALLOCATE(grhf)
   ABI_DEALLOCATE(nvresid)
 end if

 if((nstep>0.and.dtset%iscf>0).or.dtset%iscf==-1) then
   ABI_DEALLOCATE(dielinv)
   ABI_DEALLOCATE(gbound_diel)
   ABI_DEALLOCATE(irrzondiel)
   ABI_DEALLOCATE(kg_diel)
   ABI_DEALLOCATE(phnonsdiel)
   ABI_DEALLOCATE(susmat)
   ABI_DEALLOCATE(ph1ddiel)
   ABI_DEALLOCATE(ylmdiel)
 end if

 if (psps%usepaw==1) then
   if (dtset%iscf>0) then
     do iatom=1,my_natom
       pawrhoij(iatom)%lmnmix_sz=0
       pawrhoij(iatom)%use_rhoijres=0
       ABI_DEALLOCATE(pawrhoij(iatom)%kpawmix)
       ABI_DEALLOCATE(pawrhoij(iatom)%rhoijres)
     end do
   end if
   ABI_DEALLOCATE(nhat)
   call pawfgrtab_free(pawfgrtab)
   call destroy_paw_an(paw_an)
   call destroy_paw_ij(paw_ij)
   ABI_DATATYPE_DEALLOCATE(pawfgrtab)
   ABI_DATATYPE_DEALLOCATE(paw_an)
   ABI_DATATYPE_DEALLOCATE(paw_ij)
   ABI_DEALLOCATE(dimcprj)
 end if

!Restore some variables in the dtset
!Here, for TDDFT, iprcel was artificially set.
!if(dtset%iscf==-1) then
!dtset%iprcel = dtset_iprcel
!end if

 if (prtxml == 1) then
   write(ab_xml_out, "(A)") '      <finalConditions>'
!  We output the final result given in results_gs
   call out_resultsgs_XML(dtset, 4, results_gs, psps%usepaw)
   write(ab_xml_out, "(A)") '      </finalConditions>'
   write(ab_xml_out, "(A)") '    </scfcvLoop>'
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(249,2,tsec)
 call timab(238,2,tsec)

 DBG_EXIT("COLL")

end subroutine scfcv
!!***

! new subroutine for inverting KS (created by Chen)
subroutine load_param(bLoadDen,bLoadVr,maxcount,bAllowsym,bSpin,wtol,structure,ektol,penLambda)
   implicit none
   character(len=500) :: structure,msg
   integer :: bLoadDen, bLoadVr, maxcount, bAllowsym, bSpin, itmp1
   real*8  :: wtol,ektol,penLambda,tmp1

   ! load parameters for job ----------------------------------
   open (unit=1111, access="sequential", action="read", file='param.in', form="formatted", status="old");
   read (1111,*) bLoadDen,msg     ! 1: code will load density file, 0: will not load
   read (1111,*) bLoadVr,msg      ! 1: code will load potential file, 0: code will not load potential file
   read (1111,*) maxcount,msg     ! Max loop for our CG optimizer to find V_{bulk}, usually set to 200
   read (1111,*) tmp1,msg         ! 1: We do nonlocal inverting on the density
   read (1111,*) itmp1,msg        ! 1: code will do V_{bulk} to V_{atom} job, 0: code will just do invert KS job
   read (1111,*) bAllowsym,msg    ! 1: code allow symmetry (for bVatom=1 case), 0: code does not allow symmetry (for bVatom=0 case)
   read (1111,*) wtol,msg         ! if cgW changes less than wtol, program converges
   read (1111,*) structure,msg    ! FCC or BCC or SC or DIA 
   read (1111,*) ektol,msg        ! if ek changes less than ektol, SCF converges
   read (1111,*) penLambda,msg    ! coeff for the penalty function
   read (1111,*) bSpin, msg       ! the number of Spin, if b=1 (default) it is spin-unpolarized; b=2 for spinpolarized
   close (1111)

   return
end subroutine load_param

subroutine c_DoAtomV(vg, vtrial,psps,nattyp,mgfft,gmet, gsqcut, mpi_enreg, dtset, ph1d, ucvol,structure)
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
implicit none

type(pseudopotential_type) :: psps
type(dataset_type) :: dtset
type(MPI_type) :: mpi_enreg
integer :: nattyp(psps%ntypat)
real(dp) :: gmet(3,3)
integer :: mgfft
real(dp) :: gsqcut
real(dp) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
real(dp) :: ucvol
real(dp) :: vtrial(dtset%nfft,1);
real(dp) :: vg(2,dtset%nfft), vr(dtset%nfft)

character(len=500) :: structure
real(dp) :: vag(dtset%nfft,2)  ! vag(:,1):|G| vag(:,2) Vg

  call c_wrtlog("welcome to do the Vg/structure_factor function. Good luck ....")
  vr(:) = vtrial(:,1);
  
  call c_GetAtomVg(gmet,gsqcut,mgfft,mpi_enreg,dtset%natom,nattyp, &
  &                dtset%nfft,dtset%ngfft,psps%ntypat,1,ph1d,ucvol,vg,vtrial(:,1),vag,structure)
  
  call c_wrtlog('--- Done the Vg_Atom')
stop

end subroutine


subroutine c_GetAtomVg(gmet,gsqcut,mgfft,mpi_enreg,natom,nattyp,nfft,ngfft,ntypat,option,ph1d,ucvol,vg,vr,c_vag,structure)
 
  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use defs_wvltypes
  use defs_parameters
  use defs_rectypes

 implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: mgfft,natom,nfft,ntypat,option
  real(dp),intent(in) :: gsqcut,ucvol
  type(MPI_type),intent(inout) :: mpi_enreg
!arrays
  integer,intent(in) :: nattyp(ntypat),ngfft(18)
  real(dp),intent(in) :: gmet(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 

!Local variables-------------------------------
!scalars
  integer,parameter :: im=2,re=1
  integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ierr,ig1,ig2,ig3,ii,ir,isign,itypat
  integer :: jj,me_fft,me_g0,n1,n2,n3,nproc_fft,nri,old_paral_level,shift1
  real*8 :: dqm1,ee,ff,gmag,gsq,gsquar,ph1,ph12i,ph12r,ph1i,ph1r,ph2,ph2i,ph2r
  real*8 :: x1,x2,x3,y1,y2,y3, phi, phr, sfi,sfr
  character(len=5000) :: message,structure,msg

! *************************************************************************

! chen code in
  real(dp) :: vg1(2,nfft),vg(2,nfft),c_vag(nfft,2),vr(nfft)  !1st poistion store |G|, sencond position store Vg
  real(dp) :: c_vagacc(nfft,3), tmpi, tmpr
  integer :: num1, num2, num3, iq, counter
  logical :: bFirstRun

!Define G^2 based on G space metric gmet.
  gsq(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
  & dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
  & dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
  
!Real and imaginary parts of phase--statment functions:
! phr(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
! phi(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3

!-----
  counter = 0
  num1=0; num2=0; num3=0
  c_vagacc(:,:)=0.0d0
  c_vag(:,:) = 0.0d0

!Keep track of total time spent in mklocl
  n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
  me_fft=ngfft(11)
  nproc_fft=ngfft(10)

!Zero out array to permit accumulation over atom types below:
  if (option/=1 .or. ntypat/=1) then
    call c_wrtlog("error, option!=1 or ntypat/=1, stop in C_GetVAtom subroutine ")
    stop
  end if

  bFirstRun = .true. 

  call c_wrtlog("")
  write (message,*) 'gsqcut=>', gsqcut
  call c_wrtlog(message)
  call c_wrtlog("")

  open(unit=1990,access='sequential',action='write',status='replace', file='c_vag.out',form='formatted');

  ia1=1
! ia1,ia2 sets range of loop over atoms:
  ia2=ia1+nattyp(1)-1  ! there is only 1 type of atom

  write (message, *) '# of atoms in the primitive cell=', ia2 
  call c_wrtlog(message)  

  ii=0  ! index of our outpur array c_vag(ii)

  ! start the loop in G-space
  do i3=1,n3   
    if (i3>n3/2+1) then
      ig3 = i3 - (n3+1)
    else if(i3 == n3/2+1) then
      cycle ! reach the bounary, skip the rest of the code
    else
      ig3=i3-1
    end if
    
    do i2=1,n2
      if (i2>n2/2+1) then
        ig2 = i2 - (n2+1)
      else if(i2 == n2/2+1) then
        cycle  ! reach the boundary, skip the rest of the code
      else
        ig2=i2-1
      end if
    
      do i1=1,n1
        if (i1>n1/2+1) then
          ig1 = i1 - (n1+1)
        else if (i3 == n3/2+1 ) then
          cycle
        else
          ig1=i1-1
        end if
        gsquar=two_pi**2*gsq(ig1,ig2,ig3)  ! ABINIT's grimd is not our convential G, diffes by 2*pi factor!
        iq = (i3-1)*(n2*n1)+(i2-1)*n1+i1

! we drop the |G| that beyond gsqcut
        if (gsquar/two_pi**4 > gsqcut) then
          write (message, *) counter
          call c_wrtlog(message)
          counter = counter + 1
          cycle
        endif

!       Assemble structure factor over all atoms of given type:
        sfr=zero
        sfi=zero
        do ia=ia1,ia2
          x1 = ph1d(1,(ia-1)*(2*n1+1)+ig1+n1+1);
          y1 = ph1d(2,(ia-1)*(2*n1+1)+ig1+n1+1);
          x2 = ph1d(1,(ia-1)*(2*n2+1)+ia2*(2*n1+1)+ n2+ig2+1)
          y2 = ph1d(2,(ia-1)*(2*n2+1)+ia2*(2*n1+1)+ n2+ig2+1)
          x3 = ph1d(1,(ia-1)*(2*n3+1)+ia2*(2*n1+1+2*n2+1)+ n3+ig3+1)
          y3 = ph1d(2,(ia-1)*(2*n3+1)+ia2*(2*n1+1+2*n2+1)+ n3+ig3+1)
          sfr = sfr + phr(x1,y1,x2,y2,x3,y3)
          sfi = sfi - phi(x1,y1,x2,y2,x3,y3)
        end do
    
        if (structure(1:3) == 'FCC' .or. structure(1:3)=='BCC' .or. structure(1:2)=='SC') then
          ! FCC or BCC
          if (bFirstRun==.true.) then
            call c_wrtlog('doing FCC or BCC  structure ... ')
            bFirstRun = .false.
          end if
          if (vg(2,iq) > 1.0d-10) then 
            !   vg is not real
            call c_wrtlog("bulk's V is not real! stop! error!==>")
            write (message, *) vg(1,iq), vg(2,iq)
            call c_wrtlog(message)
            stop
          else
            ii=ii+1
            num1 = num1 +1
            c_vag(ii,1) = dsqrt(gsquar); 
            c_vag(ii,2) = vg(1,iq) / sfr *ucvol
          end if

        else if (structure(1:4)=='ATOM') then
          if (bFirstRun==.true.) then
            call c_wrtlog('---- doing atom...')
            bFirstRun =.false.
          endif
          if (vg(2,iq) > 1.0d-10) then 
          !   vg is not real
            call c_wrtlog("You are doing Atom density but bulk's V is not real! stop! error!")
            write (message, *) vg(1,iq), vg(2,iq)
            call c_wrtlog(message)
            stop
          else
            ii=ii+1
            num1 = num1 +1
            c_vag(ii,1) = dsqrt(gsquar); 
            sfr = 1; 
            c_vag(ii,2) = vg(1,iq) / sfr *ucvol ;
          endif
        else
        ! in fact the following codes should work for all kinds of structure as long as the structure factor
        ! is calcualted correctly into sfr and sfi. This should be the case!
        ! else if (structure(1:3) == 'DIA') then
        ! Other structures
          if (bFirstRun==.true.) then
            write (message, *) 'doing ', structure(1:3), ' structure ... '
            call c_wrtlog(message)
            bFirstRun = .false.
          end if
           
          if (sfi<1.0d-12 .and. sfr>1.0d-12) then
            num1 = num1+1
            ii = ii +1 
            c_vag(ii,1) = dsqrt(gsquar)
            c_vag(ii,2) = vg(1,iq)/sfr*ucvol
        
          else if (sfi>1.0d-12 .and. sfr<1.0d-12) then
            num2 = num2+1
            ii = ii + 1
            c_vag(ii,1) = dsqrt(gsquar)
            c_vag(ii,2) = vg(2,iq) /sfi*ucvol
        
          else if (sfi>1.0d-12 .and. sfr>1.0d-12) then
            tmpr = vg(1,iq)/sfr*ucvol
            tmpi = vg(2,iq)/sfi*ucvol
           
            if (dabs(tmpr-tmpi)>1.0d-12) then   
            ! if this cannot hold which means your density does not obey the symmetry of the structure
            !  or ... you are really in trouble!
              write(message, *)'Problem found at |G|=',gsquar,'G=',ig1,ig2,ig3,ii
              call c_wrtlog(message)
                write(message, *)'tmpr=',tmpr, ' tmpi=', tmpi, ' sfr:', sfr, ' sfi=', sfi
              call c_wrtlog(message)
                write(message, *)'vg()', vg(:,iq)
              call c_wrtlog(message)
              call c_wrtlog("exit")
             stop
            else
              num3 = num3+1
              ii = ii + 1
              c_vag(ii,1) = dsqrt(gsquar); 
              c_vag(ii,2) = tmpr
            end if              
          end if      
        end if
        write (1990, *) c_vag(ii,1), c_vag(ii,2)         
      end do
    end do
  end do
  close(1990)

  ! output
  write(msg, *) '# of G points considered:', num1, '(sfi~0,sfr/=0) + ',char(10), & 
                num2,'(sfr~0,sfi/=0) + ', char(10), & 
        num3,'(both sfi and sfr/=0)', char(10), & 
        num1+num2+num3,  ' <- total number processed.', char(10), &
          'Ratio=', dble(num1+num2+num3)/dble(nfft)
  call c_wrtlog(msg)

end subroutine c_GetAtomVg

function phr(x1,y1,x2,y2,x3,y3)
real*8 x1,x2,x3,y1,y2,y3,phr
 phr=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 return
end function

function phi(x1,y1,x2,y2,x3,y3)
real*8 x1,x2,x3,y1,y2,y3,phi
 phi =(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 return
end function

subroutine c_wrtlog(msg)

  implicit none

 include 'mpif.h'

!#if defined HAVE_CONFIG_H
!#include "config.h"
!#endif

  character(len=*), intent(in) :: msg;
  INTEGER :: me, ierr
  INTEGER,save :: master = 0
  
  me = 0

  !Determine who I am
  call MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)

      if(me/=master) RETURN

      ! only master node write to invKS_log file
      ! open file  ...
      open(unit=2222,access='append',action='write',status='unknown', file='invKS_log',form='formatted');
      
      IF (msg(1:13) .eq. 'WELCOMEMSG_BF') THEN
        WRITE(2222,*)'---------------------------------------------------------'
        write(2222,*)'      Revised version based on ABINIT v7.0.5             '
        write(2222,*)'                                                         '
        write(2222,*)'-    This program inverts Kohn-Sham equation for        -'
        write(2222,*)'-    a given electron density.                          -'
        write(2222,*)'-    A BFGS and parallel version                        -'
        write(2222,*)'-    Author: Chen Huang 2011                            -'
        write(2222,*)'---------------------------------------------------------'
      ELSE
        write (2222, '(A)') trim(msg);
      ENDIF
      

      CLOSE (2222)

end subroutine
