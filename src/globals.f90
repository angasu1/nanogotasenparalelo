Module global_variables   
  use types                             
  !General parameters                   
  real(rk),allocatable::coor(:,:)       
  type :: atom                          
  character(len=2)::nom !Name           
  integer(ik)::atnum !Atomic number     
  real(rk)::ma  !Mass                   
  end type atom                         
  type(atom),allocatable::att(:)        
  character(len=80) :: outdir='../results/'
  character(len=80) :: potfilesdir='../src/potfiles/'
  character(len=80) :: title='notitle'  
  real(rk),parameter:: navog=6.02214129e23_rk
  real(rk),parameter:: ev2har=3.6749325e-2_rk,har2ev=1._rk/ev2har
  real(rk),parameter:: ev2J=1.60217657e-19_rk,J2ev=1._rk/ev2J
  real(rk),parameter:: har2cm=219474.631_rk,cm2har=1._rk/har2cm
  real(rk),parameter:: bo2ar=0.5291772_rk,ar2bo=1._rk/bo2ar
  real(rk),parameter:: m2ar=1.0e10_rk,ar2m=1.0_rk/m2ar
  real(rk),parameter:: amu2intmu=1.0_rk/9.5485336e-3_rk
  real(rk),parameter:: au2fs=2.418884326505e-2_rk,fs2au=1.0_rk/au2fs
  real(rk),parameter:: tol=1.0e-14_rk,big=1.0e10_rk,tolb=1.0e-6_rk
  real(rk),parameter:: amu2kg=1.660538921e-27_rk,kg2amu=1.0_rk/amu2kg
  real(rk),parameter:: au2kg=9.10938291e-31_rk,kg2au=1.0_rk/au2kg
  real(rk),parameter:: amu2au=1.660538921e-27_rk/9.10938291e-31_rk
  real(rk),parameter:: au2amu=1.0_rk/amu2au
  real(rk),parameter:: tau2fs=2.41888432e-2_rk,fs2tau=1.0_rk/tau2fs
  real(rk),parameter:: pi=dacos(-1.0_rk),boltz=8.6173324e-5_rk!in eV/K 
  real(rk),parameter:: ei(3)=(/1.0_rk,0.0_rk,0.0_rk/)
  real(rk),parameter:: ej(3)=(/0.0_rk,1.0_rk,0.0_rk/)
  real(rk),parameter:: ek(3)=(/0.0_rk,0.0_rk,1.0_rk/)
  logical::fixran=.false.               
  integer(ik)::simtyp,moltyp,n_at       
  integer(ik)::n_atxmol                 
  logical::is=.false.,qgau=.false.      
  logical::jac=.false.                  
  logical::hhon=.true.,hmon=.true.,mmon=.true.
  real(rk)::requil=0.0_rk,frotmol=1.0_rk,frothe=1.0_rk
  real(rk)::rfb,rfbmax,rfbmin,drfb,dtau,rmax,rmin,umb_ang,xn(3),initd=0.0_rk
  real(rk)::diffcHe,diffchered,diffccm,A,B,C,dhe,mHe,dmassred,mtot
  real(rk)::factorrotx=1.0_rk,factorroty=1.0_rk,factorrotz=1.0_rk,rmaxpot,rexit=20._rk
  integer(ik):: nw,state,node,nstps,nhe,nmon=1,nruns,norg,lmax,nmesh,istat
  integer(ik)::idum,jtot,npoin,nrad,counter=0
  character (len=80)::molname           
  integer(ik)::pottermsp(20)=1 !Terms of the potential that are going to be considered
  real(rk),allocatable::potweig(:)      
  integer(ik),allocatable::potterms(:) !Terms of the potential that are going to be considered
  integer(ik)                            ::ptyp=1,ntermspot=0       
  integer(ik), parameter                 :: maxLinea=10000
  integer(ik)                            :: nLinea, nColumna = 22 !22, lam20
  real(rk), dimension(:,:), allocatable  :: vlam, plam
  character(4)                           :: temp

end Module global_variables  
