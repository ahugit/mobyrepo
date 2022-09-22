	module params branchtesting branchtesting2 branchtesting3 b3
	after branching
    even more after branching
    issuetesting1branch
    !insert space and keep tabs option in visual studio
    ! testing
	use nrtype 
	use nr 
	use alib, only: one,logit,logitinv,min2pls,min2plsinv,lin2ndim,ndim2lin,sqrtpi,multinom,condmom
	implicit none 
    include 'mpif.h'
    integer :: mysay,mygroup
	integer :: iter,comm,iwritegen 
	!real(dp) :: one=1.0_dp
	!integer(i4b), parameter :: rp = kind(1.0d0)			! kind(1.0) !!!
    real(dp), parameter :: replacement_rate=0.4_dp          !ahu summer18 050318: added replacement rate
    integer(i4b), parameter :: nl=9,ndecile=10
!ahu030622	logical, parameter :: groups=.true.,onlysingles=.true.,onlymales=.false.,onlyfem=.false.,optimize=.true.,chkstep=.false.,condmomcompare=.false.,comparepars=.false.,extramoments=.true.
logical, parameter :: groups=.true.,onlysingles=.false.,onlymales=.false.,onlyfem=.false.,optimize=.false.,chkstep=.false.,condmomcompare=.false.,comparepars=.false.,extramoments=.true. !ahu030622
logical :: nonneg
    logical, parameter :: onthejobsearch=.TRUE. !set in m\ain
	integer(i4b), parameter :: numit=5 !ahumarch1122
    real(dp), dimension(2) :: nonlabinc !=(/ 0.0_dp,0.0_dp /) !(/ 300.0_dp,1100.0_dp /) !ahu summer18 051418: changing it back to parameter and changing dimension to 2 (not educ and educ) !ahu summer18 042318 changing this so it is set at main again
	real(dp), parameter :: eps = 1.0d-6,zero=0.0_dp,epstest=2.0_dp					! could do tiny(.) but that gives a number that is way too small and therefore not appropriate for those places where we check the inequalities or equalities	
	real(dp), parameter :: eps2= 1.0d-6
	integer(i4b), parameter :: nhome=1,nhomep=nl
	logical :: conditional_moments		! can set this in main
	logical :: skriv,yaz,insol,yazmax   
	integer(i4b) :: whereami
	character(len=1), parameter :: runid='r'		! string inserted into output filenames to identify which run !ahu 062413 set this in main instead 
	integer(i4b), parameter :: nco=1,ncop=1
	integer(i4b), parameter :: ntyp=1,ntypp=4   ! types !ahu030622 changed ntypp to 1 (was 4)
	integer(i4b), parameter :: nin  = nco * ntyp * nhome
	integer(i4b), parameter :: ninp = ncop * ntypp * nhomep
    integer(i4b) :: nindex !determined in objf according to groups, it's either nin or ninp
	integer(i4b) :: iwritemom,myhome,mytyp,myco,myindex,mygrank
	logical, parameter :: indsimplexwrite=.false.		! in optimiation with parallel simplex, re-solve model and write moments to file whenever check for,find new best point 
	logical, parameter :: write_welfare=.false.			! write things needed for welfare calculations
	logical, parameter :: grid_usage=.false.			! keep track of how much of each grid is being used in simulations
	logical, parameter :: moremom=.false.
	logical, parameter :: icheck_eqvmvf=.false.,icheck_eqvcvs=.false.,icheck_probs=.false.
	integer(i4b), parameter :: npars    = 93
    character(len=15), dimension(npars) :: parname ! names of each parameter   !ahu 121118 now declkare as global it here instead of getsteps
    real(dp), dimension(npars) :: stepmin,realpartemp,parsforcheck,stepos !ahu 121118
	!character(len=15), dimension(npars) :: parname 
	integer(i4b), parameter :: mna=18,mxa=50   !,mxai=50		!ahu 070713 40	!if you ever change mina from 16 that might be problematic, because of psiddata%home and simdata%home definitions. look in read_data and read simdata for this
    integer(i4b), parameter :: MNAD=MNA-1,MXAD=MXA-1            !ahu jan19 010219
    integer(i4b), parameter :: nh=2,nexp=4,nsimeach=10,neduc=2,nkid=2 !kid is 1 if no kid,2 if yes kid !ahu 0327 changed nsimeach from 10 to 5
	integer(i4b), parameter :: np=3,nz=1 !ag090122 agsept2022 changed nz frmo 1 to 5 !ahu 121818 changed from 3 to 6 !ahu 0327 changed np from 5 to 2
	integer(i4b), parameter :: nqs = (np+2) *  nl
	integer(i4b), parameter :: nq  = (np+2) * (np+2) * nl
	integer(i4b), parameter :: np1=np+1, np2=np+2   !w=np1 is getting laid off in the shock space q and w=np2 is nothing happening. In the state space q, np1 is unemployment and np2 is not in the state space q (I try to check that it is not there, at various points in code)
	integer(i4b), parameter :: nxs = neduc * nexp * nkid
	integer(i4b), parameter :: nx  = nxs * nxs
	integer(i4b), parameter :: ncs = nl+2
	integer(i4b), parameter :: nc  = nl+8
    integer(i4b), parameter :: nepsmove=3, nepskid=2 !ag090122 agsept2022 changed nepsmove frmo 2 to 5 !ahumarch1022 changed nepsmove to 2 from 13
	integer(i4b), parameter :: ndata    = 5233 !5390 !2386   
	integer(i4b), parameter :: ndataobs = 84507 !86873 !41494  
	integer(i4b), parameter :: nsim     = ndata*nsimeach  
	integer(i4b), parameter :: nmom     = 2200 !ahu summer18 050418: changed from 4200 to 498
    integer(i4b) :: calcvar(nmom),calcorr(nmom)
	integer(i4b), parameter :: maxrellength=10
	integer(i4b), parameter :: namelen=60					!if you change this, don't forget to also change a100 in writemoments	
	integer(i4b), parameter :: ma=1,fe=2
    INTEGER(I4B), PARAMETER :: NOCOLLEGE=1,COLLEGE=2
	integer(i4b), parameter, dimension(2) :: agestart=(/ mna,22 /)		!changed this from 18,22 !chanage this back ahu 070312 (/18,22/) !starting age for simulations for each education level
	real(dp), parameter :: mult1=50000.0_dp !ahu jan19 012519
    real(dp), parameter :: multmar=500000.0_dp,mult1c=50000.0_dp  !ahu jan19 012019  !ahu030622 VERY IMPORTANT CHANGE MULTMAR
	real(dp), parameter :: maxhgrid=8.0_dp 
	real(dp), parameter :: tottime=16.0_dp
	real(dp), parameter :: hhours_conv=250.0_dp					! multuiply hours per day by this to get hours per year for hours worked
	real(dp), parameter :: maxh=hhours_conv*maxhgrid				! hhours_conv*hmgrid(ntwork) ! truncation point of labor market hours in data !ahu 071712
	real(dp), parameter, dimension(nh) :: hgrid=(/ 0.0_dp,maxhgrid /) 
	real(dp), parameter, dimension(nh) :: hgrid_ann=hgrid*hhours_conv
	real(dp), parameter :: d1=1.0_dp						! note that now doing calculations as single real(sp) not double real (i modified module1 to allow this)
	real(dp), parameter :: h_parttime=1000.0_dp					! number of hours to be considred part time !ahu 071712
	real(dp), parameter :: h_fulltime=1000.0_dp					! ahu 071712 changed from 1000 too 2000 !ahu 062812 changed this to 1000. 2000. ! number of hours to be considred full time 
        !ahu 021817: note that the employment is decided according to whether hours is more than h_fulltime which is 1000 but annual wages is calculated using wge*h_wmult where h_wmult is 2000. 
        !This is how you did it in the original version. I do the same thing so that the wage numbers are consistent with the previous version. See page 15 last paragraph in the original text.	
    real(dp), parameter :: h_wmult=2000.0_dp                    !what you multiply hourly wages with in order to turn them into annual wages
    real(dp), parameter :: hbar=h_parttime
	real(dp), parameter :: minw=1.0_dp					! lower truncation point of male log wage
	real(dp), parameter :: maxw=150.0_dp                ! upper truncation point of male log wage
	real(dp), parameter :: pen=-99999999.0_dp
	integer(i4b), parameter :: ipen=-99999	
	real(dp), parameter :: init=pen,initi=ipen
    !ahu 122818 changed wmove mar from 100 to 1000
    !ag090122 agsept2022 real(dp), parameter :: wtrans=50.0_dp,wwaged=50.0_dp,wdifww=50.0_dp,wrel=1.0_dp,wmove=10.0_dp,whour=1.0_dp,wwage=10.0_dp,wkid=1.0_dp,wmovemar=1.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    !ag090122 agsept2022: increasing wmovemar to 100 because mar move rates are too low
    real(dp), parameter :: wtrans=50.0_dp,wwaged=50.0_dp,wdifww=50.0_dp,wrel=1.0_dp,wmove=10.0_dp,whour=1.0_dp,wwage=10.0_dp,wkid=1.0_dp,wmovemar=200.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    !real(dp), parameter :: wtrans=10.0_dp,wwaged=1.0_dp,wdifww=1.0_dp,wrel=1.0_dp,wmove=10.0_dp,whour=1.0_dp,wwage=1.0_dp,wkid=1.0_dp,wmovemar=1.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    !real(dp), parameter :: wtrans=50.0_dp,wwaged=50.0_dp,wdifww=1.0_dp,wrel=1.0_dp,wmove=10.0_dp,whour=1.0_dp,wwage=10.0_dp,wkid=1.0_dp,wmovemar=1.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    !real(dp), parameter :: wrel=10.0_dp,wmove=10.0_dp,whour=1.0_dp,wwage=1.0_dp,wkid=1.0_dp,wmovemar=10.0_dp,wmovesin=10.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    !real(dp), parameter :: wrel=1.0_dp,wmove=10.0_dp,whour=1.0_dp,wwage=1.0_dp,wkid=1.0_dp,wmovemar=1000.0_dp,wmovesin=1.0_dp,wwagebymove=10.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    !real(dp), parameter :: wrel=1.0_dp,wmove=1.0_dp,whour=1.0_dp,wwage=1.0_dp,wkid=1.0_dp,wmovemar=1.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp	! weights for moments for married couples. set in objfunc.
    character(len=23), parameter :: datafilename= 'familymigpsid.txt' ! data filename
	!character(len=23), parameter :: initcondfile= 'familymiginit111113.txt' ! filename of initial conditions
	character(len=23) :: momentfile='mom.txt'
    character(len=23) :: momentonlyfile='momonly.txt'
    integer(i4b) :: mominfo(0:5,nmom)
	!model parameters 
	!real(sp):: umdist_condf(np,np),ufdist_condm(np,np)
	real(dp), parameter :: delta=0.96_dp,alf=0.5_dp   !ahumarch1022 delta=0.96_dp changing to 0 to figure out the mumardecrease problem
	real(dp), parameter :: mu_wge(2)=0.0_dp
	real(dp) :: sig_wge(2),mu_mar(ntypp),sig_mar,ro,mu_o , sigo_m,sigo_f
	real(dp) :: uhome(2),alphaed(2,neduc),alphakid(2,nkid),alfhme
	real(dp) :: gam_e,gam_u,cst(ntypp),kcst,ecst,scst,divpenalty,uloc(nl),sig_uloc
	real(dp) :: alf10(nl),alf11,alf12,alf13,alf1t(ntypp)            ! types
	real(dp) :: alf20(nl),alf21,alf22,alf23,alf2t(ntypp)            ! types
	real(dp) :: ptype,pmeet,omega(2),ptypehs(ntypp),ptypecol(ntypp) ! types
	real(dp) :: pkid,psio(12),psil(2),psih(4)
	real(dp) :: popsize(nl)
	integer(i4b) :: distance(nl,nl)
	real(dp) :: wg(np,2),wgt(np),mg(nz,ninp),mgt(nz),best
	type :: initcond
        integer(i4b) :: id
		integer(i4b) :: co				
		integer(i4b) :: sexr			
		integer(i4b) :: hme
		integer(i4b) :: endage		
		integer(i4b) :: edr			
    end type
	type, extends(initcond) :: statevar
        integer(i4b) :: expr
        integer(i4b) :: kidr      
        integer(i4b) :: hhr		! annual hours worked by r (turned into discrete in read_data)
		real(dp) :: logwr,wr	! wr is annual income and logwr is log of annual income. hourly wage (wr_perhour or wsp_perhour) is read from the data (see familymig_2.do to see how it's turned into hourly) and that hourly wage is turned into annual by multiplying it by h_wmult	  
        integer(i4b) :: l		! location		
        integer(i4b) :: rel		! relationship status. -1: not observed, 0: single, 1: married, 2: cohabiting
		integer(i4b) :: rellen	! length of current relationship starting from 1 in first period
        integer(i4b) :: edsp
        integer(i4b) :: expsp
        integer(i4b) :: kidsp
        integer(i4b) :: hhsp	! annual hours worked by spouse (turned into discrete in read_data)
		real(dp) :: logwsp,wsp  ! see explanation for logwr,wr.  
        integer(i4b) :: lsp
        integer(i4b) :: nomiss
        integer(i4b) :: nn,mm,r,typ
	end type
	type :: shock
		real(dp) :: meet 
		real(dp) :: marie
		real(dp) :: meetq 
		real(dp) :: meetx 
		real(dp) :: q
		real(dp) :: x
        real(dp) :: iepsmove
        real(dp) :: typ
	end type	
	type(statevar) :: ones
    type(initcond) :: ones_init
contains




	! get parameters from transformed values. input is free
	! to be any real(sp) resulting parameters are appropriately constrained
	subroutine getpars(par,realpar)
	real(dp), dimension(npars), intent(in)  :: par ! vector of parameters
	real(dp), dimension(npars), intent(out) :: realpar ! vector of parameters
	integer(i4b) :: g,i,j,ed,indust1(ntypp),indust2(ntypp)
    integer(i4b) :: dw0,de,dsex
    real(dp) :: oftemp(3,np1,2,2)

    stepos=0.0_dp
    indust1=0
    indust2=0
    
	realpar=pen 
    parname=''
    j=1
    !ahu jan19 012819: not iterating on ed offers anymore. replacing them with curloc and ofloc offers instead 
	realpar(j)=par(j)               ; parname(j)='emp,cur,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(1)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,cur,m' ; stepos(j)=1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(2)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,of,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp!psio is for the offer function
	psio(3)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,of,m' ; stepos(j)=0.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function !NO MORE OFLOC LAYOFF NONSENSE
	psio(4)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,cur,f' ; stepos(j)=-1.2_dp ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(5)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,cur,f' ; stepos(j)=-1.5_dp ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function
	psio(6)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,of,f' ; stepos(j)=1.0_dp ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function
	psio(7)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,of,f' ; stepos(j)=0.0_dp ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function !NO MORE OFLOC LAYOFF NONSENSE
	psio(8)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,cur,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(9)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,of,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(10)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,cur,f' ; stepos(j)=0.2_dp  ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(11)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,of,f' ; stepos(j)=1.0_dp  ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(12)=realpar(j)	            ; j=j+1
    !print*, 'Here is psio12',j-1

	realpar(j)=par(j)               ; parname(j)='psil(1)' ; stepos(j)=0.0_dp 
	psil(1)=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='alfhme' ; stepos(j)=0.0_dp 
	alfhme=realpar(j)	            ; j=j+1
	realpar(j)=min2pls(par(j))      ; parname(j)='ro'	; stepos(j)=0.5_dp ; if (onlysingles) stepos(j)=0.0_dp !15 !2.0_dp*(1.0_dp/(1.0_dp+exp(-par(j))))-1.0_dp 
	ro=realpar(j)                   ; j=j+1

	realpar(j)=logit(par(j))               ; parname(j)='psih(1)' ; stepos(j)=-1.0_dp !ahu jan19 011719 changing to logit
	psih(1)=realpar(j)	            ; j=j+1
	realpar(j)=logit(par(j))               ; parname(j)='psih(2)' ; stepos(j)=0.0_dp !ahu jan19 011719 changing to logit
	psih(2)=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='psih(3)' ; stepos(j)=0.0_dp !ahu jan19 011519 getting rid of probdown
	psih(3)=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='psih(4)' ; stepos(j)=0.0_dp  !ahu jan19 011519 getting rid of probdown
	psih(4)=realpar(j)	            ; j=j+1

    
    realpar(j)=logit(par(j))        ; parname(j)='pkid' ; stepos(j)=1.0_dp  ; if (onlysingles) stepos(j)=0.0_dp !20
	pkid=realpar(j)                 ; j=j+1
    realpar(j)=logit(par(j))	    ; parname(j)='pmeet' ; stepos(j)=0.5_dp ; if (onlysingles) stepos(j)=0.0_dp !21
	pmeet=realpar(j)                ; j=j+1


    !realpar(j) = mult1c * logit(par(j))              ; parname(j)='uhome(1)' ; stepos(j)=0.2_dp	!mult3*logit(par(2:3)) !22:23
	!uhome(1)=realpar(j)                             ; j=j+1
    !realpar(j) = mult1c * logit(par(j))              ; parname(j)='uhome(2)' ; stepos(j)=0.2_dp	 ; if (onlymales) stepos(j)=0.0_dp !mult3*logit(par(2:3)) !22:23
	!uhome(2)=realpar(j)                             ; j=j+1

    realpar(j) = par(j)             ; parname(j)='uhome(1)' ; stepos(j)=0.5_dp*PAR(J)	 ; if (onlyfem) stepos(j)=0.0_dp !mult3*logit(par(2:3)) !22:23
	uhome(1)=realpar(j)                             ; j=j+1
    realpar(j) = par(j)              ; parname(j)='uhome(2)' ; stepos(j)=0.5_dp*PAR(J)	 ; if (onlymales) stepos(j)=0.0_dp !mult3*logit(par(2:3)) !22:23
	uhome(2)=realpar(j)                             ; j=j+1
    realpar(j)=par(j)            ; parname(j)='ecst'	; stepos(j)=0.5_dp*par(j) !24 !-1.0_dp*mult1c * logit(par(j)) !ahu 112718 changing to only minus from: mult1 * min2pls(par(j))     ! types
    ecst=realpar(j)                                     ; j=j+1               ! types
	realpar(j) = par(j)          ; parname(j)='kcst'	; stepos(j)=0.5_dp*par(j) !25 !-1.0_dp*mult1c * logit(par(j)) !ahu 112718 changing to only minus from: mult1 * min2pls(par(5)) !mult2*logit(par(4:6))	
	kcst=realpar(j)                                     ; j=j+1
	realpar(j) = -1.0_dp*mult1 * logit(par(j))          ; parname(j)='divpenalty'	; stepos(j)=1.5_dp ; if (onlysingles) stepos(j)=0.0_dp !26 !ahu 112718 changing to only minus from: mult1 * min2pls(par(6))                         !ahu summer18 050418: changed from 1000 to 10,000 (mult to mult1)
	divpenalty=realpar(j)                               ; j=j+1
    !print*, 'Here is divpenalty',j-1,divpenalty 

    realpar(j:j+1) = mult1 * logit(par(j:j+1))          ; parname(j)='alphaed(m,ned)' ; parname(j+1)='alphaed(f,ned)'    !27:28   !ahu jan19 012719 changing it yet again back to logit because there is not that much of different in objval between alpha=0 and alpha=-49000    !ahu jan19 012019 changing it back to min2pls  ! noed !ahu 112718 changing to only plus from: mult1*min2pls(par(7:8))   !mult1 * logit(par(7))	
	stepos(j)=0.5_dp  ; if (onlyfem) stepos(j)=0.0_dp ; 	stepos(j+1)=0.5_dp  ; if (onlymales) stepos(j+1)=0.0_dp 
    alphaed(:,1)=realpar(j:j+1)                         ; j=j+2
    !realpar(j:j+1)= mult1 * min2pls(par(j:j+1))           ; parname(j)='alphaed(m,ed)' ; parname(j+1)='alphaed(f,ed)'         !ahu jan19 012019 changing it back to min2pls       ! ed !ahu 112718 changing to only plus from: mult1 * min2pls(par(9:10))	 !mult1 * logit(par(9:10))	
	!stepos(j)=5.0_dp ; 	stepos(j+1)=5.0_dp   ; if (onlymales) stepos(j+1)=0.0_dp 
    realpar(j:j+1)=alphaed(:,1)            ; parname(j)='alphaed(m,ed)' ; parname(j+1)='alphaed(f,ed)'     !ahu jan19 012419 changing this so that there is no alpha by ed no more    !ahu jan19 012019 changing it back to min2pls       ! ed !ahu 112718 changing to only plus from: mult1 * min2pls(par(9:10))	 !mult1 * logit(par(9:10))	
	stepos(j)=0.0_dp  ; if (onlyfem) stepos(j)=0.0_dp ; 	stepos(j+1)=0.0_dp   ; if (onlymales) stepos(j+1)=0.0_dp 
    alphaed(:,2)=realpar(j:j+1)                         ; j=j+2 
    realpar(j:j+1)=mult1 * logit(par(j:j+1))            ; parname(j)='alphakid(m)' ; parname(j+1)='alphakid(f)'          !31:32           !ahu 112718 changing to only plus from: mult1 * min2pls(par(j:j+1))	 !mult1 * logit(par(9:10))	
    stepos(j)=5.0_dp  ; if (onlyfem) stepos(j)=0.0_dp ; 	stepos(j+1)=0.5_dp  ; if (onlymales) stepos(j:j+1)=0.0_dp 
    alphakid(:,2)=realpar(j:j+1)                        ; j=j+2         
    !print*, 'Here is uloc',j
	
    !uloc: 33-41
    !do ed=1,2
        do i=1,nl
            if (i==2) then
			    realpar(j) = 0.0_dp  ; stepos(j)=0.0_dp
			    uloc(i)=0.0_dp
		    else 
			    realpar(j) = par(j) ; stepos(j)=0.5_dp*PAR(J)    !mult1 * min2pls( par(j) )
			    uloc(i)=realpar(j)
		    end if 
            parname(j)='uloc' 
            j=j+1
        end do
	!end do
    !print*, 'Here is alf10',j
	!wage 42: 65
    do i=1,nl
        realpar(j)=par(j) !1.5_dp*min2pls(par(j))+8.5_dp 
		alf10(i)=realpar(j)
        parname(j)='alf10' ; stepos(j)=0.15_dp ; if (i==3) stepos(j)=0.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !dont iterate on alpha10(1) for a second
        j=j+1
    end do 
    !print*, 'Here is alf11 etc',j
    realpar(j)=logit(par(j))                        ; parname(j)='alf11' ; stepos(j)=0.3_dp  ; if (onlyfem) stepos(j)=0.0_dp
	alf11=realpar(j)                                ; j=j+1
    !print*, 'Here is alf12',j	
    realpar(j)=3.0_dp*logit(par(j))                 ; parname(j)='alf12' ; stepos(j)=0.3_dp  ; if (onlyfem) stepos(j)=0.0_dp
    alf12=realpar(j)                                ; j=j+1
    !print*, 'Here is alf13',j	
    realpar(j)=0.0_dp                               ; parname(j)='alf13' ; stepos(j)=0.0_dp   ; if (onlyfem) stepos(j)=0.0_dp  !-1.0_dp*logit(par(j)) 
    alf13=realpar(j)	                            ; j=j+1
    !print*, 'Here is alf20',j
    do i=1,nl
		realpar(j)=par(j) !1.5_dp*min2pls(par(j))+8.5_dp    
		alf20(i)=realpar(j) 
        parname(j)='alf20' ; stepos(j)=0.15_dp  ; if (onlymales) stepos(j)=0.0_dp 
		j=j+1
	end do 
    !print*, 'Here is alf21 etc',j
	realpar(j)=logit(par(j))                        ; parname(j)='alf21' ; stepos(j)=0.3_dp  ; if (onlymales) stepos(j)=0.0_dp 
	alf21=realpar(j)                                ; j=j+1
    !print*, 'Here is alf22',j	    
    realpar(j)=3.0_dp*logit(par(j))                 ; parname(j)='alf22' ; stepos(j)=0.3_dp  ; if (onlymales) stepos(j)=0.0_dp 
	alf22=realpar(j)                                ; j=j+1
    !print*, 'Here is alf23',j	
    realpar(j)=0.0_dp                               ; parname(j)='alf23' ; stepos(j)=0.0_dp  ; if (onlymales) stepos(j)=0.0_dp !-1.0_dp*logit(par(j)) 
	alf23=realpar(j)	                            ; j=j+1
	
    realpar(j:j+1)=logit(par(j:j+1))                ; parname(j:j+1)='sig_wge'	; stepos(j:j+1)=0.7_dp	  ; if (onlyfem) stepos(j)=0.0_dp  ; if (onlymales) stepos(j+1)=0.0_dp !66:67
	sig_wge(1:2)=realpar(j:j+1)                     ; j=j+2
    !sigom and sigof: 68:69
    realpar(j)=par(j)                               ; parname(j)='sigo_m'	; stepos(j)=0.5_dp*PAR(J) ; if (nepsmove==1) stepos(j)=0.0_dp ; if (onlyfem) stepos(j)=0.0_dp
    !print*, "Here it is sigom", j,par(j),realpar(j)
    sigo_m=realpar(j)                                ; j=j+1
    realpar(j)=par(j)                               ; parname(j)='sigo_f'	; stepos(j)=0.5_dp*PAR(J) ; if (nepsmove==1) stepos(j)=0.0_dp ; if (onlymales) stepos(j)=0.0_dp
    !print*, "Here it is sigof", j,par(j),realpar(j)
    sigo_f=realpar(j)                                ; j=j+1

    
    do i=1,ntypp
        if (i==1) then 
            realpar(j)=0.0_dp                           ; parname(j)='ptypehs' ; stepos(j)=0.0_dp
            ptypehs(i)=exp(realpar(j))                  ; indust1(i)=j ; j=j+1
            realpar(j)=0.0_dp                           ; parname(j)='ptypecol'  ; stepos(j)=0.0_dp
            ptypecol(i)=exp(realpar(j))                 ; indust2(i)=j ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf1t'     ; stepos(j)=0.2_dp  ; if (onlyfem) stepos(j)=0.0_dp
            alf1t(i)=realpar(j)                         ; j=j+1
            realpar(j)=par(j)                          ; parname(j)='alf2t'     ; stepos(j)=0.2_dp  ; if (onlymales) stepos(j)=0.0_dp
	        alf2t(i)=realpar(j)                         ; j=j+1
            !realpar(j)= -1.0_dp*mult1c * logit(par(j))   ; parname(j)='cst'       ; stepos(j)=0.5_dp
            !cst(i)=realpar(j)                           ; j=j+1 
            realpar(j)= par(j)                          ; parname(j)='cst'       ; stepos(j)=1.5_dp*par(j) !not iterating on this anymore. see notes. under cost vs. sigo. they are just not sep ident I think. 
            cst(i)=realpar(j)                           ; j=j+1 
            !ahu082822 august2022 print*, 'mumar(1)',j,par(j),multmar, min2pls(par(j)),multmar*min2pls(par(j))
            realpar(j)=multmar * min2pls(par(j))          ; parname(j)='mu_mar'     ; stepos(j)=1.5_dp*par(j)    ; if (onlysingles) stepos(j)=0.0_dp 	    
            mu_mar(i)=realpar(j)                        ; j=j+1      
        else
            realpar(j)=par(j)                            ; parname(j)='ptypehs' ; stepos(j)=1.5_dp*par(j)
            ptypehs(i)=exp(realpar(j))                  ; indust1(i)=j ; j=j+1
            realpar(j)=par(j)                            ; parname(j)='ptypecol'  ; stepos(j)=1.5_dp*par(j)
            ptypecol(i)=exp(realpar(j))                 ; indust2(i)=j ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf1t'     ; stepos(j)=0.2_dp  ; if (onlyfem) stepos(j)=0.0_dp
            alf1t(i)=realpar(j)                         ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf2t'     ; stepos(j)=0.2_dp  ; if (onlymales) stepos(j)=0.0_dp
	        alf2t(i)=realpar(j)                         ; j=j+1
            realpar(j)= par(j)                          ; parname(j)='cst'       ; stepos(j)=1.5_dp*par(j)
            cst(i)=realpar(j)                           ; j=j+1 
            realpar(j)=par(j)                        ; parname(j)='mu_mar'     ; stepos(j)=1.5_dp*par(j)    ; if (onlysingles) stepos(j)=0.0_dp 	    
            mu_mar(i)=realpar(j)                        ; j=j+1      
        end if 
    	!print*, 'Here is cost',cst(i),par(j-1),realpar(j-1)
        !ahu 122818 changed mult1 to multmar 
    end do 
    ptypehs(:)=ptypehs(:)/sum(ptypehs)
    ptypecol(:)=ptypecol(:)/sum(ptypecol)
    do i=1,ntypp
        realpar(indust1(i))=ptypehs(i)
        realpar(indust2(i))=ptypecol(i)
    end do 
 
!ahu jan19 012419: not iterating the below anymore 
psih(2)=0.0_dp
!ECST=-8000.0_DP
!kcst=0.0_dp
alphaed(:,2)=alphaed(:,1)



    !alf11=alf21
    !alf12=alf22
    !alf13=alf23
    !psio(1:4)=psio(5:8)
    !psio(9:10)=psio(11:12)
    
    
    !ahu 121918 make all types the same 
    !    alf1t=0.0_dp
    !    alf2t=0.0_dp
    !    cst=-10000.0_dp
    !    mu_mar=1993.0_dp
    !if (ntypp==1) then 
    !    ptypehs=1.0_dp 
    !    ptypecol=1.0_dp 
    !else 
    !    ptypehs=0.25_dp
    !    ptypecol=0.25_dp
    !end if
    !print*, 'ptypehs', ptypehs
    !print*, 'ptypecol',ptypecol
    !print*, 'alf1t',alf1t
    !print*, 'alf2t',alf2t
    !print*, 'cst',cst
    !print*, 'mu_mar',mu_mar

    

    
    
    !print*, 'ptype', ptypehs, ptypecol
    
    !print*, 'Here is IT ISISISISISI npars', j-1
    !if (j.ne.npars) then
    !    print*, 'j not equal to npars',j,npars
    !    stop
    !end if


    mu_o=0.0_dp
    !sig_o=1000.0_dp
    scst=0.0_dp
    sig_mar=0.0_dp
    
    !alf11=0.0_Dp 
    !alf12=0.0_dp
    !alf13=0.0_dp 
    !sig_wge=0.0001_dp
    
    !if (ntypp==1) then 
    !    ptypehs(1)=1.0_dp    ! types
    !    ptypecol(1)=1.0_dp  ! types
    !else if (ntypp==2) then 
    !    ptypehs(ntypp)=1.0_dp-ptypehs(1)    ! types
    !    ptypecol(ntypp)=1.0_dp-ptypecol(1)  ! types
    !end if     
    !alf1t(1)=0.0_dp                 ! types
    !alf2t(1)=0.0_dp                 ! types
    alphakid(:,1)=0.0_dp

    !cst=0.0_dp
    !ecst=0.0_dp
    !kcst=0.0_dp
    
    
    !***********************
    !ahu 041118 del and remove later:
    !alphaed(2,:)=alphaed(1,:)
    !psio(5:8)=psio(1:4)
    !psio(11:12)=psio(9:10)
    !alphakid(2,:)=alphakid(1,:)
    
    !alf20=alf10
    !alf21=alf11
    !alf22=alf12
    !alf23=alf13
    !uhome(2)=uhome(1)
    !***********************
    
    !psio(1:4)=psio(5:8)
    !psio(9:10)=psio(11:12)
    
    !pkid=0.0_dp
    !alphakid=0.0_dp
    !kcst=0.0_dp 
    !ro=0.0_dp
    !sig_mar=0.0_dp
    !scst=0.0_dp
    
    !alf11=0.0_dp
    !alf21=0.0_dp
    !psio(3:4)=psio(1:2)
    !psio(7:8)=psio(5:6)
    !psio(10)=psio(9)
    !psio(12)=psio(11)
    !uloc(:,2)=uloc(:,1)
    !alphaed(:,2)=alphaed(:,1)
    !alf1t(2)=0.0_dp
    !alf2t(2)=0.0_dp
    !ptypecol=ptypehs
    
    
    !cst=0.0_dp
    !kcst=0.0_dp
    
    !ro=0.0_dp !0.98_dp
!alpha=0.0_dp
!kcst=0.0_dp
!pkid=0.0_dp
!alf1t(2)=alf1t(1)
!alf2t(2)=alf2t(1)
!cst(1)=cst(2)
!uhome=0.0_dp

    !uhome=0.0_dp
    !cst=-150000.0_dp
    !kcst=-150000.0_dp
    !divpenalty=0.0_dp
    !pkid=0.0_dp
    !kcst=0.0_dp
    !alpha(:,2)=alpha(:,1)
    
    
    !sig_o=sig_mar    
    !alpha(1,:)=alpha(2,:) 
    	!if ((.not.optimize).and.(.not.chkstep) ) print*, "sig_mar,mu_mar,npars;,; ", sig_mar,mu_mar,j
    
	!print*, logitinv(alf11),logitinv(alf12),logitinv(alf13)
	!print*, logitinv(alf21),logitinv(alf22),logitinv(alf23)

	!realpar(j)=exp(par(j))
	!sig_uloc=realpar(j) !; j=j+1		
	!realpar(j)=logit(par(j))
	!omega(1)=realpar(j) ; j=j+1
	!realpar(j)=logit(par(j))
	!omega(2)=realpar(j) 
	!realpar(j)=par(j)
	!alf14(ntyp)=realpar(j)  ; j=j+1
	!alf14(1)=0.0_sp
	!realpar(j)=par(j)
	!alf24(ntyp)=realpar(j) ; j=j+1  
	!alf24(1)=0.0_sp
	!realpar(j)=logit(par(j))
	!ptype=realpar(j)	!low type prob 
	!if (j/=npars) then ; print*, "something wrong in getpar! ",j,npars ; stop ; end if
    

    oftemp=0.0_dp
    do dsex=1,2
        do de=1,2            
            do dw0=1,np1
                if ( dw0 <= np ) then 
		            if ( de==1 .and. dsex==1 ) then 
			            oftemp(1:2,dw0,de,dsex)=exp(psio(1:2)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 
		            else if ( de==2 .and. dsex==1 ) then 
			            oftemp(1:2,dw0,de,dsex)=exp(psio(3:4)) 
                    else if ( de==1 .and. dsex==2 ) then 
			            oftemp(1:2,dw0,de,dsex)=exp(psio(5:6)) 
                    else if ( de==2 .and. dsex==2 ) then 
			            oftemp(1:2,dw0,de,dsex)=exp(psio(7:8)) 
                    end if 
		            oftemp(3,dw0,de,dsex)=exp( 0.0_dp )		! nothing happens												
	            else if (dw0 == np1) then 
		            if ( de==1 .and. dsex==1 ) then 
			            oftemp(1,dw0,de,dsex)=exp(psio(9)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                        realpar(9)=oftemp(1,np1,1,1)
                    else if ( de==2 .and. dsex==1 ) then 
			            oftemp(1,dw0,de,dsex)=exp(psio(10) )
                        realpar(10)=oftemp(1,np1,2,1)
                    else if ( de==1 .and. dsex==2 ) then 
			            oftemp(1,dw0,de,dsex)=exp(psio(11)) 
                        realpar(11)=oftemp(1,np1,1,2)
		            else if ( de==2 .and. dsex==2 ) then 
			            oftemp(1,dw0,de,dsex)=exp(psio(12)) 
                        realpar(12)=oftemp(1,np1,2,2)
                    end if 
		            oftemp(2,dw0,de,dsex)=0.0_dp		! 0 since you can't get laid off if you don't have a job! 
		            oftemp(3,dw0,de,dsex)=exp(0.0_dp)		! nothing happens												
	            else  
		            print*, "in fnprof: dw0 > np1 which doesnt' make sense as that's a state variable " 
		            stop
	            end if 
	            oftemp(1:3,dw0,de,dsex)=oftemp(:,dw0,de,dsex)/sum( oftemp(:,dw0,de,dsex)     )
            end do 
        end do 
    end do 
                       realpar(1:2)=oftemp(1:2,1,1,1)  !emp,noed,male
                        realpar(3:4)=oftemp(1:2,1,2,1)  !emp,ed,male
                        realpar(5:6)=oftemp(1:2,1,1,2)  !emp,noed,fem
                        realpar(7:8)=oftemp(1:2,1,2,2)  !emp,ed,fem
                        realpar(9)=oftemp(1,np1,1,1)  !unemp,noed,male
                        realpar(10)=oftemp(1,np1,2,1)  !unemp,ed,male
                        realpar(11)=oftemp(1,np1,1,2)  !unemp,noed,fem
                        realpar(12)=oftemp(1,np1,2,2)  !unemp,ed,fem
    
        realpartemp=realpar
	end subroutine getpars

	subroutine getsteps(par,step)
	real(dp), dimension(:), intent(in) :: par 
	!character(len=15), dimension(:), intent(out) :: name ! names of each parameter
	real(dp), dimension(:), intent(out) :: step 
	integer(i4b) :: i,j
	!name='' 
	!step=0.0_dp
    !name(1)='sig_o'		        ; step(1)=0.5_dp ; if (nepsmove==1) step(1)=0.0_dp 
    !name(2)='uhome(m)'		    ; step(2)=1.0_dp*par(2)
	!name(3)='uhome(f)'		    ; step(3)=1.0_dp*par(3)
	!name(4)='cst(1)'			; step(4)=-1.0_dp*par(4)
	!name(5)='kcst'		        ; step(5)=1.0_dp*par(5)
	!name(6)='divpenalty'		; step(6)=1.0_dp*par(6) ; if (onlysingles) step(6)=0.0_dp
	!name(7)='alphaed(m,1)'		; step(7)=1.0_dp*par(7)
	!name(8)='alphaed(f,1)'		; step(8)=1.0_dp*par(8)
	!name(9)='alphaed(m,2)'		; step(9)=1.0_dp*par(9) 
	!name(10)='alphaed(f,2)'		; step(10)=1.0_dp*par(10) 
	!name(11)='pkid'		        ; step(11)=1.0_dp   
	!name(12)='pmeet'		    ; step(12)=-1.0_dp ; if (onlysingles) step(12)=0.0_dp
	!j=13
    !do i=1,nl
    !    name(j)='uloc'	; step(j)=1.0_dp*par(j)    ; if (i==2) step(j)=0.0_dp ; j=j+1
    !end do 	
    
	!do i=1,nl
    !    name(j)='alf10' ; step(j)=0.5_dp*par(j)  ; j=j+1
	!end do 
	!name(j)='alf11'		; step(j)=0.5_dp*par(j) ; j=j+1
	!name(j)='alf12'     ; step(j)=0.5_dp*par(j) ; j=j+1
    !name(j)='alf13'		; step(j)=0.5_dp*par(j) ; j=j+1
    
    
	!do i=1,nl
    !    name(j)='alf20'	; step(j)=0.5_dp*par(j) ; j=j+1
	!end do
	!name(j)='alf21'		; step(j)=0.5_dp*par(j) ; j=j+1
    !name(j)='alf22'		; step(j)=0.5_dp*par(j) ; j=j+1
    !name(j)='alf23'		; step(j)=0.5_dp*par(j) ; j=j+1
    !name(j)='sig_wge(1)'	; step(j)=0.5_dp ; j=j+1
    !name(j)='sig_wge(2)'	; step(j)=0.5_dp ; j=j+1
    !name(j)='ro'			; step(j)=0.5_dp  ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psio'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psil'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psil'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psil'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psih'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psih'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psih'		; step(j)=0.5_dp ; j=j+1
    !name(j)='psih'		; step(j)=0.5_dp ; j=j+1
    !name(j)='mu_mar(1)'		    ; step(j)=1.0_dp*par(j)    ; if (onlysingles) step(j)=0.0_dp ; j=j+1
    !name(j)='mu_mar(ntypp)'		; step(j)=1.0_dp*par(j)    ; if (onlysingles) step(j)=0.0_dp ; j=j+1
    !name(j)='sig_mar'		; step(j)=0.0_dp    ; if (onlysingles) step(j)=0.0_dp ; j=j+1
    !name(j)='ptypehs(1)'		; step(j)=-1.0_dp    ; j=j+1  
    !name(j)='ptypecol(1)'		; step(j)=-1.0_dp    ; j=j+1 
    !name(j)='alf1t(2)'		    ; step(j)=0.5_dp    ; j=j+1 
    !name(j)='alf2t(2)'		    ; step(j)=0.5_dp    ; j=j+1
    !name(j)='cst(2)'            ; step(j)=1.0_dp*par(j) ; j=j+1    
    !name(j)='ecst'              ; step(j)=-1.0_dp*par(j) ; j=j+1     
    !name(j)='alphakid(m,2)'     ; step(j)=1.0_dp*par(j) ; j=j+1 
    !name(j)='alphakid(f,2)'     ; step(j)=1.0_dp*par(j) ; j=j+1
    !name(j)='scst'              ; step(j)=0.0_dp    !*par(j)     
	!if (j.ne.npars) then 
    !    print*, 'j not equal to npars',j,npars
    !    stop
    !end if 
    step=0.0_dp
	end subroutine getsteps

	subroutine getdistpop
	integer(i4b) :: i,j
	distance=0
	popsize=0.0_dp 
    !ahu 030717: redid this adjacence thing. see map and appendix of the draft. 
	distance(1,2)=1 
	distance(2,1)=1 
	distance(2,5)=1 
	distance(2,3)=1 
	distance(3,2)=1 
	distance(3,4)=1 
	distance(3,6)=1 
	distance(3,5)=1 
	distance(4,3)=1 
	distance(4,8)=1 
	distance(4,7)=1 
	distance(4,6)=1 
	distance(5,8)=1 
	distance(5,3)=1 
	distance(5,6)=1 
	distance(6,5)=1 
	distance(6,7)=1 
	distance(6,4)=1 
	distance(6,3)=1 
	distance(7,6)=1 
	distance(7,8)=1 
	distance(7,4)=1 
	distance(8,7)=1 
	distance(8,9)=1 
	distance(8,4)=1 
	distance(9,8)=1 
    do i=1,nl 
		do j=1,nl 
			if (j==i) then 
				if (distance(j,i)==1) then 		! ahu 061513: for some locaitons, you did have that their distance(i,i) was 1 so corrected this
                    print*, 'distance(i,i) should not be 1 because that is not adjacence'
                    stop
                end if 
			end if 
		end do 
	end do 

	popsize(1)=0.9112_dp	!new england
	popsize(2)=2.695_dp		!middle atlantic 
	popsize(3)=2.9598_dp	!east north central
	popsize(4)=1.2468_dp	!west north central
	popsize(5)=2.5736_dp	!south atlantic
	popsize(6)=1.0089_dp	!east south central
	popsize(7)=1.6121_dp	!west south central
	popsize(8)=0.7661_dp	!mountain
	popsize(9)=2.2757_dp	!pacific
	end subroutine getdistpop

	subroutine getones
		ones%co=-99
		ones%sexr=-99
		ones%hme=-99
		ones%endage=-99
		ones%edr=-99
		ones%expr=-99
        ones%kidr=-99
		ones%hhr=-99
        ones%logwr=-99.0_dp
        ones%wr=-99.0_dp
        ones%wsp=-99.0_dp
        ones%l=-99
        ones%rel=-99
        ones%rellen=-99
        ones%edsp=-99
        ones%expsp=-99
        ones%kidsp=-99
        ones%hhsp=-99
        ones%logwsp=-99.0_dp
		ones%lsp=-99        
		ones%nomiss=-99
        ones%nn=-99
        ones%mm=-99
        ones%r=-99
        ones%typ=-99
        
        ones_init%id=-99
        ones_init%co=-99
        ones_init%sexr=-99
        ones_init%hme=-99
        ones_init%endage=-99
        ones_init%edr=-99
        end subroutine getones	
        
        
end module params
