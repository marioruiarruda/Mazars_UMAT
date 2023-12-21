!     ------------------------------------------------------------------
!     |     DATE           PROGRAMMER           DESCRIPTION OF CHANGE  |  
!     |     ====           ==========           =====================  |
!     |  04/09/2022      MARIO RUI ARRUDA           STABILIZATIOBN     |
!     |  IST LISBON                                   ALGORITHM        |
!     ------------------------------------------------------------------

!     #COPYRIGHT 2022 BY MÁRIO RUI TIAGO ARRUDA
!     ALL RIGHTS RESERVED. NO PART OF THIS SUBROUTINE MAY BE REPRODUCED OR 
!     USED IN ANY MANNER WITHOUT WRITTEN PERMISSION OF THE COPYRIGHT OWNER.

!     If using this UMAT in future papers please cite: 
!     https://doi.org/10.1016/j.engfracmech.2021.108129
!     and give credit to the original authors of this UMAT.
  
!     ------------------------------------------------------------------
!     ----------ABAQUS INPUT VARIABLES IN UMAT SUBROUTINE---------------
!     ------------------------------------------------------------------
subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, &
           ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, & 
           nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, &
           npt, layer, kspt, kstep, kinc)

IMPLICIT NONE
  
  
!     ------------------------------------------------------------------
!     -------------------ABAQUS DIMENSION VARIABLES---------------------
!     ------------------------------------------------------------------
 
!     ------------------ABAQUS CHARACTER VARIABLES----------------------
character(len=80) :: cmname 

!     -------------------ABAQUS INTEGER VARIABLES-----------------------  
integer :: ntens,ndi,nshr,nstatv,nprops,noel,npt,kspt,kstep,kinc,nprecd,layer
 
!     -------------------ABAQUS REAL VARIABLES--------------------------
real(kind=8) :: celent,sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt 
 
!     -------------------ABAQUS ARRAY VARIABLES-------------------------
real(kind=8) :: stress(ntens),statev(nstatv),ddsdde(ntens,ntens),&
                ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),&
                predef(1),dpred(1),props(nprops),coords(3),drot(3,3),&
                dfgrd0(3,3),dfgrd1(3,3),time(2)
      
      
!     ------------------------------------------------------------------    
!     -----------------DECLARATION OF VARIABLES-------------------------
!     ------------------------------------------------------------------

integer :: i,j,k,l,m,n,kiter,ktotal

real(kind=8) :: E,G,nu,dt,dc,damage,ds,dto,dco,damageo,dso,eta,coef_eeq,fctm,tauc,&
                alphat,alphac,GfIt,GfIc,GfII,ed0,ed,eeq,edo,sd0,coef_int,damagev,&
                pt,pc,ps,dmaxt,dmaxc,dmaxs,At,Ac,Bt,Bc,dvt,dvc,dvs,dvto,dvco,&
					      dvso,Dcoef,ec1,ec2,fcm,eeq_coef,alphato,alphaco,ddtdeeq,ddcdeeq
                           
!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5

!     ------------------INITIATIONS OF ARRAYS---------------------------
real(kind=8) :: eij(ntens),sij(ntens),sije(ntens),sijeo(ntens),eijo(ntens),eij_e(ntens),&
                sijo(ntens),deij(ntens),dsij(ntens),ddsddee(ntens,ntens),&
                eijt(3),eijc(3),psije(3),peij(3),peijbp(3),psijebp(3),psijebn(3),&
					      eijtbp(3),eijcbp(3),peijo(3),peijbpo(3),eijbp(ntens),dddej(ntens),&
                depdeij(3,ntens),pveij(ntens,ntens),dt_eij(ntens),dc_eij(ntens)

!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3-------------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0, FOUR=4.d0,&
                           HALF=0.5d0,TOL0=1.d-8,mONE=-1.d0
											 						 
!     --------------ELASTIC AND MECHANICAL PROPERTIES-------------------

E=props(1)      ! LONGITUDINAL ELASTIC MODULUS
nu=props(2)     ! POISSON COEFICIENT
fctm=props(3)   ! CRACKING STRESS
GfIt=props(5)   ! FRACTURE ENERGY FOR TENSION
GfIc=props(6)   ! FRACTURE ENERGY FOR COMPRESSION

eta=props(8)    ! VISCOUS REGULARIZATION COEFICIENT
dmaxt=props(9)  ! MAXIMUM TENSION DAMAGE ALLOWED
dmaxc=props(10) ! MAXIMUM COMPRESSION DAMAGE ALLOWED

pt=props(12)    ! RESIDUAL TENSILE STRESS
pc=props(13)    ! RESIDUAL COMPRESSION STRESS

At=props(15)    ! At PARAMETER FOR TENSION    
Ac=props(16)    ! AC PARAMETER FOR COMPRESSION 
Bt=props(17)    ! Bt PARAMETER FOR TENSION   
Bc=props(18)    ! BC PARAMETER FOR COMPRESSION
ec1=props(19)   ! strain ec1 for compression   
ec2=props(20)   ! strain ec2 for compression
fcm=props(21)   ! maximum compression stress

G=E/(TWO*(ONE+nu)) ! SHEAR ELASTIC MODULUS
ed0=fctm/E         ! CRACKING STRAIN

!     ---------------STATE FIELD VARIABLES FOR ABAQUS-------------------
if (kinc == 1) then 
	do i=1, nstatv
		statev(i)=ZERO  ! VARIABLES INITIATION (if not=0, depending on compiler LINUX/WINDOWS)
  end do
end if

if (kinc==1) then
  statev(7)=ed0     ! PROBLEMATIC FOR MORE THAN ONE LOAD STEP
end if

dto=statev(1)      ! TENSION DAMAGE FROM PREVIOUS INCREMENT
dco=statev(2)      ! COMPRESSION DAMAGE FROM PREVIOUS INCREMENT

dvto=statev(4)     ! VISCOUS TENSION DAMAGE FROM PREVIOUS INCREMENT
dvco=statev(5)     ! VISCOUS COMPRESSION DAMAGE FROM PREVIOUS INCREMENT

edo=statev(7)      ! UPDATED CRACKING STRAIN FROM PREVIOUS INCREMENT
alphato=statev(8)  ! PREVIOUS TENSION CONSTANT
alphaco=statev(9)  ! PREVIOUS COMPRESSION CONSTANT
damageo=statev(10) ! PREVIOUS COMPRESSION CONSTANT

dt=dto             ! INITIATION OF TENSION DAMAGE VARIABLE
dc=dco             ! INITIATION OF COMPRESSION DAMAGE VARIABLE
ds=dso             
dvt=dvto           ! INITIATION OF TENSION VISCOUS DAMAGE VARIABLE
dvc=dvco           ! INITIATION OF COMPRESSION VISCOUS DAMAGE VARIABLE
dvs=dvso

ed=edo             ! INITIATION OF THE LOADING FUNCTION
damage=damageo     ! INITIATION OF DAMAGE VARIABLE
alphat=alphato     ! INITIATION OF TENSION CONSTANT
alphac=alphaco     ! INITIATION OF COMPRESSION CONSTANT
coef_eeq=ONE       ! INITIATION OF COMPRESSIVE INTEGER


!     ------------------------------------------------------------------
!     -------FIELD PREDICTOR FOR STRAIN AND EQUIVALENT STRAIN-----------
!     ------------------------------------------------------------------

!     --------------------FINAL STRAIN INCREMENT------------------------
do i=1, ntens
  eij(i)=stran(i)+dstran(i)          ! TOTAL STRAIN IN THE BEGGIND OF THE INCREMENT
  eijo(i)=stran(i)                   ! STRAIN FROM PREVIOUS INCREMENT
	deij(i)=dstran(i)                  ! STRAIN INCREMENT
  sijo(i)=stress(i)                  ! STRESS FROM PREVIOUS INCREMENT
  const4=dtime/(dtime+eta)           ! FROM IMPLICIT TO EXPLICIT WHEN dtime<<eta
  eij_e(i)=stran(i)+const4*dstran(i) ! FOR EULER INTEGRATION (ZERO,HALF,ONE)
end do

!     ------------CALCULATION OF PRINCIPAL EFECTIVE STRAIN--------------
if (ndi==3) then ! FOR 3D SOLID ANALYSIS
    call SPRINC(eij_e,peij,2,ndi,nshr)       ! CALCULATION OF PRINCIPAL STRAIN 2 but STRESS 1 // FROM ABAQUS LIBRARY
    call SPRINC(eijo,peijo,2,ndi,nshr)
    call SPRIND(eij_e,peij,pveij,2,ndi,nshr)  ! pveij(k,i) k -  principal value j - directions
else            ! FOR 2D PLANE ANALYSIS 
    call calc_strain_princ_2D(ntens,nu,eij_e,peij)
    call calc_strain_princ_2D(ntens,nu,eijo,peijo)
end if

!     ----------CALCULATION MACAULAY BRACKET EFECTIVE STRAIN------------
call calc_macaulay_bracket(ndi,peij,peijbp,ONE)
call calc_macaulay_bracket(ndi,peijo,peijbpo,ONE)


!     ------------------------------------------------------------------
!     --------------FIELD PREDICTOR FOR EFECTIVE STRESS-----------------
!     ------------------------------------------------------------------

!     ---------------CALCULATION OF EFECTIVE STIFFNESS------------------
call calc_stiffness(ntens,ndi,ddsddee,Dcoef,E,G,nu,ZERO)
  
!     --------------PRINCIPAL 3D STRAIN TO 3D STRESS--------------------
call calc_psij_k_peij(ntens,ndi,peijo,psije,ddsddee,E,nu) 
  
!     -------------PRINCIPAL EFECTIVE STRESS WITH BRACKET -------------
call calc_macaulay_bracket(ndi,psije,psijebp,ONE)
call calc_macaulay_bracket(ndi,psije,psijebn,mONE)

!     ----------EQUIVALENTE STRAIN CORRECTOR FOR CONFINEMENT ----------
!call calc_confinement_coef(ndi,psije,psijebn,coef_eeq)

!     -----------------TENSION AND COMPRESSION STRAIN -----------------
call calc_signal_strain(ndi,E,nu,eijt,psijebp)
call calc_signal_strain(ndi,E,nu,eijc,psijebn)
  
!     ------POSITIVE TENSION AND COMPRESSION STRAIN WITH BRACKET -------
call calc_macaulay_bracket(ndi,eijt,eijtbp,ONE)
call calc_macaulay_bracket(ndi,eijc,eijcbp,ONE)
  
!     ----------------TENSION AND COMPRESSION PARAMETERS----------------
call calc_alpha_e(eijtbp,eijcbp,ndi,alphat,alphac)
    
statev(8)=alphat
statev(9)=alphac


!     ------------------------------------------------------------------
!     ------------------FAILURE CRITERIA VERIFICATION-------------------
!     ------------------------------------------------------------------

!     ------------------LOADING FUNCTION CALCULATION--------------------
eeq=sqrt(peijbp(1)**2+peijbp(2)**2+peijbp(3)**2)
!eeq=eeq*coef_eeq

!     ------------------------DAMAGE EVOLUTION--------------------------
if ((eeq-ed)>ZERO) then ! ALTERNATLY USE if ((eeq-ed0)>=ZERO)

  call calc_damage_t(GfIt,E,nu,dto,dt,dvto,dvt,ed0,At,Bt,eeq,celent,pt,eta,dtime,ddtdeeq,ed)
  
  call calc_damage_c(GfIc,E,nu,dco,dc,dvco,dvc,ed0,Ac,Bc,ec1,ec2,eeq,celent,pc,eta,fcm,dtime,ddcdeeq,ed)
  
  dt=min(dt,dmaxt)   ! USED FOR IMPROVING CONVERGENCE DURING TENSION DAMAGE 
  dc=min(dc,dmaxc)   ! USED FOR IMPROVING CONVERGENCE DURING COMPRESSION DAMAGE
	statev(1)=dt
	statev(2)=dc
  
  dvt=min(dvt,dmaxt) ! USED FOR IMPROVING CONVERGENCE DURING VISCOUS TENSION DAMAGE 
  dvc=min(dvc,dmaxc) ! USED FOR IMPROVING CONVERGENCE DURING VISCOUS COMPRESSION DAMAGE
	statev(4)=dvt
	statev(5)=dvc
  
end if

statev(7)=max(eeq,ed)
  
!     ---------------------UPDATED VISCOUS DAMAGE-----------------------  
damage=alphat*dt+alphac*dc
damagev=alphat*dvt+alphac*dvc
statev(10)=damage


!     ------------------------------------------------------------------
!     ---------FIELD CORRECTOR FOR STRESS AND SECANT STIFFNESS----------
!     ------------------------------------------------------------------

call calc_stiffness(ntens,ndi,ddsdde,Dcoef,E,G,nu,damagev)

!     ------------------------UPDATED STRESS----------------------------
call calc_stress(ndi,ntens,stress,eij,ddsdde)

!     ---------------------UPDATED STRAIN ENERGY------------------------
call calc_energy(ndi,ntens,sse,stress,sijo,dstran)

return
end




!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||USER SUBROUTINES||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||




!     __________________________________________________________________
!     ___________SUBROUTINE FOR MACAULY BRACKET OPERATION_______________
!     __________________________________________________________________

subroutine calc_macaulay_bracket(ndi,mat,mat_b,signal)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i
integer, intent(in):: ndi

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: mat(3),signal
real(kind=8), intent(inout) :: mat_b(3)

!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2---------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,3
	if (signal>=ZERO) then
		mat_b(i)=abs(max(ZERO,mat(i)))
	else
		mat_b(i)=abs(min(ZERO,mat(i)))
	end if
end do

return
end subroutine


!     __________________________________________________________________
!     _________SUBROUTINE FOR STRAIN STIFFNESS MULTIPLICATION___________
!     __________________________________________________________________

subroutine calc_stress(ndi,ntens,mat_si,mat_ej,mat_Cij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: mat_Cij(ntens,ntens), mat_ej(ntens)
real(kind=8), intent(inout) :: mat_si(ntens)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
  mat_si(i)=ZERO
  do j=1,ndi ! NORMAL STRAIN
	  mat_si(i)=mat_si(i)+mat_Cij(i,j)*mat_ej(j)
  end do
end do
do i=ndi+1,ntens ! SHEAR STRAIN
  mat_si(i)=mat_Cij(i,i)*mat_ej(i) 
end do

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE FOR DEFORMATION ENERGY__________________
!     __________________________________________________________________

subroutine calc_energy(ndi,ntens,energy,sij,sij_old,deij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: sij(ntens), sij_old(ntens), deij(ntens)
real(kind=8), intent(inout) :: energy

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ntens
  energy=energy+HALF*(sij_old(i)+sij(i))*deij(i)
end do

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE STIFFNESS ASSEMBLER____________________
!     __________________________________________________________________

subroutine calc_stiffness(ntens,ndi,Cij,Dvar,E,G,nu,d)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens,ndi

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: E,G,nu,d
real(kind=8), intent(inout) :: Cij(ntens,ntens), Dvar

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

if (ndi==3) then ! FOR 3D SOLID ANALYSIS
  Dvar=E/((ONE+nu)*(ONE-TWO*nu))
  Cij(1,1)=(ONE-d)*(ONE-nu)*Dvar
  Cij(1,2)=(ONE-d)*nu*Dvar
  Cij(1,3)=Cij(1,2)
  Cij(2,1)=Cij(1,2)
  Cij(2,2)=Cij(1,1)
  Cij(2,3)=Cij(1,2)
  Cij(3,1)=Cij(1,2)
  Cij(3,2)=Cij(1,2)
  Cij(3,3)=Cij(1,1)
  Cij(4,4)=(ONE-d)*G
  Cij(5,5)=Cij(4,4)
  Cij(6,6)=Cij(4,4)
  Cij(4,1)=ZERO;Cij(4,2)=ZERO;Cij(4,3)=ZERO;Cij(4,5)=ZERO;Cij(4,6)=ZERO
  Cij(5,1)=ZERO;Cij(5,2)=ZERO;Cij(5,3)=ZERO;Cij(5,4)=ZERO;Cij(5,6)=ZERO
  Cij(6,1)=ZERO;Cij(6,2)=ZERO;Cij(6,3)=ZERO;Cij(6,4)=ZERO;Cij(6,5)=ZERO
else            ! FOR 2D SOLID ANALYSIS
  Dvar=ONE/(ONE-nu**2)
  Cij(1,1)=(ONE-d)*E*Dvar
  Cij(2,2)=Cij(1,1)
  Cij(3,3)=(ONE-d)*G
  Cij(1,2)=(ONE-d)*nu*E*Dvar
  Cij(2,1)=Cij(1,2)
  Cij(1,3)=ZERO;Cij(3,1)=ZERO
  Cij(2,3)=ZERO;Cij(3,2)=ZERO
end if

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE FOR 2D PRINCIPAL VALUES_________________
!     __________________________________________________________________

subroutine calc_strain_princ_2D(ntens,nu,mat_eij,pij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8) :: C,Cm,R
real(kind=8), intent(in)  :: nu,mat_eij(ntens)
real(kind=8), intent(inout) :: pij(3)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0
                           
C=(mat_eij(1)+mat_eij(2))*HALF        ! CENTER OF MOHR CIRCULE
Cm=(mat_eij(1)-mat_eij(2))*HALF       ! AVERAGE TENSOR OF MOHR CIRCULE
R=sqrt(Cm**2+(mat_eij(3)/TWO)**2)     ! RADIUS OF MOHR CIRCULE
pij(1)=C+R
pij(2)=C-R
pij(3)=-nu/(ONE-nu)*(pij(1)+pij(2))

return
end subroutine


!     __________________________________________________________________
!     ____________SUBROUTINE FOR TENSILE DAMAGE EVOLUTIION______________
!     __________________________________________________________________

subroutine calc_damage_t(Gf,E,nu,dto,dt,dvto,dvt,ed0,A,B,eeq,catLc,pt,eta,dtime,dtde,ed)

IMPLICIT NONE

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8) :: eu,ep
real(kind=8), intent(in)  :: Gf,E,nu,dto,dvto,ed0,A,B,eeq,catLc,pt,eta,dtime,ed
real(kind=8), intent(inout) :: dt,dvt,dtde

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

dtde=ZERO
if (Gf<=ZERO) then
    dt=ONE-ed0*(ONE-A)/eeq-A/exp(B*(eeq-ed0))
else
  eu=Gf/(catLc*E*ed0)+ed0           ! STRAIN FOR LINEAR SOFTENING
  ep=ed0-(eu-ed0)*log(pt)           ! STRAIN AT RESIDUAL STRESS FOR EXPONENTIAL SOFTENING
  !ep=eu-pt*(eu-ed0)                ! STRAIN AT RESIDUAL STRESS FOR LINEAR SOFTENING
  if (eeq<ep) then
    dt=ONE-ed0/eeq*exp((ed0-eeq)/(eu-ed0))
    dtde=(ed0/eeq**2)*exp((ed0-eeq)/(eu-ed0))+&
         (ed0/eeq)*ONE/(eu-ed0)*exp((ed0-eeq)/(eu-ed0))
    ! dt=eu*(eeq-ed0)/(eeq*(eu-ed0))         ! LINEAR SOFTENING
    ! dtde=(eu*eeq*(eu-ed0)-eu*(eeq-ed0)*(eu-ed0))/(eeq*(eu-ed0))**2
  else
    dt=ONE-pt*ed0/eeq
    dtde=pt*ed0/eeq**2
  end if
end if

dt=max(dto,dt) 
dvt=eta/(eta+dtime)*dvto+dtime/(eta+dtime)*dt
dvt=max(dvto,dvt)
if (eeq<ed) then
  dtde=ZERO
end if
  
return
end subroutine 


!     __________________________________________________________________
!     _________SUBROUTINE FOR COMPRESSIVE DAMAGE EVOLUTIION_____________
!     __________________________________________________________________

subroutine calc_damage_c(Gf,E,nu,dco,dc,dvco,dvc,ed0,A,B,e1,e2,eeq,catLc,pc,eta,fcm,dtime,dcde,ed)

IMPLICIT NONE

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5,const6,ct_k1,ct_k2

!     ------------------------REAL VARIABLES----------------------------
real(kind=8) :: eu,ep,eeqc,ecb,ek
real(kind=8), intent(in)  :: Gf,E,nu,dco,dvco,ed0,A,B,e1,e2,eeq,catLc,pc,eta,fcm,dtime,ed
real(kind=8), intent(inout) :: dc,dvc,dcde

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

dcde=ZERO
if (Gf<=ZERO) then
  dc=ONE-ed0*(ONE-A)/eeq-A/exp(B*(eeq-ed0))
else
  eu=TWO*Gf/(catLc*fcm)-(e2-e1)                                ! ULTIMATE CRUSHING STRAIN
  const4=-fcm/(eu-e2)                                          ! m VARIABLE FROM SOFTENING EQUATION
  const5=fcm-const4*e2                                         ! b VARIABLE FROM SOFTENING EQUATION b=fcm-m*e2
	ct_k1=fcm/(eu-e2)                                            ! k1 VARIABLE FROM SOFTENING EQUATION
  ct_k2=fcm+ct_k1*e2                                           ! k2 VARIABLE FROM SOFTENING EQUATION
  ep=(pc*fcm-const5)/const4                                    ! PLASTIC STRAIN INITIATION
  eeqc=eeq/(nu*sqrt(TWO))
  if (eeqc<=e1) then
    ecb=eeqc/e1                                                ! eta VARIABLE FROM EUROCODE 2 EQUATION 3.14
    ek=1.05d0*E*e1/fcm                                         ! k VARIABLE FROM EUROCODE 2 EQUATION 3.14
    const6=(ek*ecb-ecb**2)/(ONE+(ek-TWO)*ecb)                  ! EQUATION 3.14 FROM EUROCODE 2
    dc=ONE-const6*fcm/(E*eeqc)
    const1=(ONE+(ek-TWO)*ecb)*E*eeqc
    const2=(ek*ecb-ecb**2)
    const3=((ek-TWO*ecb)*const1-const2*(ek-TWO)*E*eeqc)*fcm/const1**2
    dcde=-(const3*ONE/e1*ONE/(nu*sqrt(TWO)) &
         -const2*fcm/((ONE+(ek-TWO)*ecb)*E*eeqc**2)*ONE/(nu*sqrt(TWO)))
  else if (eeqc>e1 .and. eeqc<=e2) then
    dc=ONE-fcm/(E*eeqc)
    dcde=fcm/(E*eeqc**2)*ONE/(nu*sqrt(TWO))
  else if (eeqc>e2 .and. eeqc<=ep) then
    dc=ONE+ct_k1/E-ct_k2/(E*eeqc)
    dcde=ct_k2/(E*eeqc**2)*ONE/(nu*sqrt(TWO))
  else
    dc=ONE-pc*fcm/(E*eeqc)
    dcde=pc*fcm/(E*eeqc**2)*ONE/(nu*sqrt(TWO))
  end if
end if

dc=max(dco,dc)
dvc=eta/(eta+dtime)*dvco+dtime/(eta+dtime)*dc
dvc=max(dvco,dvc)
if (eeq<ed) then
  dcde=ZERO
end if
  
return
end subroutine 


!     __________________________________________________________________
!     ________SUBROUTINE FOR TENSION AND COMPRESSION PARAMETERS_________
!     __________________________________________________________________

subroutine calc_alpha_e(ebp_t,ebp_c,ndi,at,ac)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi

!     ------------------------REAL VARIABLES----------------------------
real(kind=8) :: et,ec,ev
real(kind=8), intent(in)  :: ebp_t(3),ebp_c(3)
real(kind=8), intent(inout) :: at,ac

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

et=ZERO
ec=ZERO
do i=1, 3
  et=et+ebp_t(i)
  ec=ec+ebp_c(i)
end do
ev=et+ec

at=min(et/ev,ONE)
ac=ONE-at

return
end subroutine


!     __________________________________________________________________
!     __________________SUBROUTINE PRINCIPAL STRAINS____________________
!     __________________________________________________________________

subroutine calc_signal_strain(ndi,E,nu,peij,psij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ndi

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5

!     ------------------------REAL VARIABLES----------------------------
real(kind=8), intent(in)  :: psij(3),E,nu
real(kind=8), intent(inout) :: peij(3)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

peij(1)=(ONE+nu)/E*psij(1)-nu/E*(psij(1)+psij(2)+psij(3))
peij(2)=(ONE+nu)/E*psij(2)-nu/E*(psij(1)+psij(2)+psij(3))
peij(3)=(ONE+nu)/E*psij(3)-nu/E*(psij(1)+psij(2)+psij(3))

return
end subroutine


!     __________________________________________________________________
!     _____________SUBROUTINE PRINCIPAL MULTIPLICATION__________________
!     __________________________________________________________________

subroutine calc_psij_k_peij(ntens,ndi,peij,psij,Cij,E,nu)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5

!     ------------------------REAL VARIABLES----------------------------
real(kind=8), intent(in)  :: E,nu
real(kind=8), intent(inout) :: peij(3),psij(3),Cij(ntens,ntens)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
   psij(i)=ZERO
   do j=1,ndi ! NORMAL STRAIN
     psij(i)=psij(i)+Cij(i,j)*peij(j)
   end do
end do
if (ndi==2) then ! FOR THE 2D PLANE ANALYSIS
  psij(3)=ZERO
end if

return
end subroutine

!     __________________________________________________________________
!     ______________SUBROUTINE FOR CONFINEMENT CORRECTOR________________
!     __________________________________________________________________

subroutine calc_confinement_coef(ndi,pij,psijbn,coef)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5

!     -------------------------REAL VARIABLES--------------------------- 
real(kind=8), intent(in)  :: pij(3),psijbn(3)
real(kind=8), intent(inout) :: coef

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0
                           
if ((pij(1)+pij(2)+pij(3))<ZERO) then
  const1=ZERO
  const2=ZERO
  do i=1, 3
    const1=const1+(psijbn(i))**2
    const2=const2+pij(i)
  end do
  coef=sqrt(const1)/abs(const2)
  coef=min(coef,ONE)
else
  coef=ONE
end if

return
end subroutine