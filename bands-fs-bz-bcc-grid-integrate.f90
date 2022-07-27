program bccKpointsIntegration
  !*
  !* This program generates the (kx,ky,kz) grid in asending order
  !* for a BCC-IBZ using kmin < kc < kmax and calculates
  !* the area of the sphere of radius kc
  !*
  !* script: fs-bz-bcc-grid-integrate.sh
  !*
  !* input files
  !*
  !* fort.1(=out)=$sym/$case.kband_$band\_$Nk OK
  !* fort.123 OK
  !* fort.244 & fort.245 OK
  !*
  !* output files
  !* fort.12 => ideal FS OK
  !* fort.71 => Ideal FS @ $sym/$case.ideala-ibz\_$knia OK
  !* fort.200 => grad=-1 & BF=0 for FS
  !* fort.201 => grad=-1 & BF=1 for FS
  !* fort.202 => grad=+1 & BF=0 for FS
  !* fort.203 => grad=+1 & BF=1 for FS
  !* fort.8 => kcBand37 OK
  !* fort.246 => FS $case.kflist_$ene with cartesian k-points & data for integration OK
  !* fort.35 => area-vs-k-fs.dat OK
  !*
  !* OLD-output files
  !*
  !* fort.(69,70)   => weight=48   (Ideal,FS)
  !* fort.(71,72)   => all weights (Ideal,FS)
  !* fort.(500,501) => weight=48   Ideal (in->out,out->in)
  !* fort.(502,503) => weight=48   FS    (in->out,out->in)
  !* fort.(400,401) => all weights Ideal (in->out    ,    out->in)
  !*                                 (1,2,3)->(4,5,6) (1,2,3)->(4,5,6)
  !* => Ideal Spherical Fermi Surface k-points are given by
  !* => fort.401@(4,5,6)-all-in  k-points @ideal-bz-bcc-grid-integrate.sh:  kin=$sym/$case.ideala.in
  !* => fort.400@(4,5,6)-all-out k-points @ideal-bz-bcc-grid-integrate.sh: kout=$sym/$case.ideala.out
  !*
  !* fort.(402,403) => all weights FS (in->out,out->in)
  !*
  !* fort.(201,200) => all weights FS (in,out) we chose fort.201 "in" as the Fermi-surface
  !* fort.201       => all weights FS where $10=1 => grad=+1 kpoints or $10=-1 => grad=-1 kpoints
  !* fort.(202,203) => all weights FS (in-out,out->in) vectors for grad=-1
  !* fort.(204,205) => all weights FS (in-out,out->in) vectors for grad=+1
  !*
  !*  fort.12 => all kmin < k < kmax
  !*  fort.29 = (fort.12 = fort.30 = fort.3 = fort.4 not printed)
  !*  fort.803, fort.22 => not printed
  !*  fort.69=fort.70 => weight=48 or inner points
  !*  fort.71 => has GPH points missing
  !*  fort.4 = fort.12 => all kmin < k < kmax
  !*  fort.899 not needed
  !*  fort.61 => gradf(x,y,z)= -1 or 0 & w=48 & w=24 only @ GNP-plane
  !*  fort.60 => gradf(x,y,z)= -1  & w=48 & w=24 only @ GNP-plane

  IMPLICIT NONE
  !REAL(8), ALLOCATABLE, DIMENSION(:) :: kvecin,kvecout
  type ragged_array
     real(8),allocatable :: kvin(:),kvout(:),kvinf(:),kvoutf(:)
     real(8),allocatable :: kvina(:),kvouta(:),kvinfa(:),kvinfap(:),kvoutfa(:),kvoutfap(:)
  end type ragged_array
  type(ragged_array),allocatable:: ikvin(:) ,ikvout(:) ,ikvinf(:) ,ikvoutf(:)
  type(ragged_array),allocatable:: ikvina(:),ikvouta(:),ikvinfa(:),ikvinfap(:),ikvoutfa(:),ikvoutfap(:)
  REAL(8), ALLOCATABLE :: d(:)
  !
  REAL(8), ALLOCATABLE :: energy(:)
  REAL(8), ALLOCATABLE :: rx(:),ry(:),rz(:)
  REAL(8), ALLOCATABLE :: rxx(:),ryy(:),rzz(:)
  REAL(8), ALLOCATABLE :: B(:,:,:),BF(:,:,:)
  REAL(8), ALLOCATABLE :: weight(:,:,:)
  INTEGER, ALLOCATABLE :: tagi(:,:,:)
  REAL(8), ALLOCATABLE :: kpi(:),kpo(:)
  REAL(8), ALLOCATABLE :: rxm(:),rym(:),rzm(:)
  REAL(8), ALLOCATABLE :: rxp(:),ryp(:),rzp(:)
  INTEGER :: N,i,j,k,w,tag
  INTEGER :: Nk,aNk
  INTEGER :: numm,nump,dummy
  INTEGER :: conta,contall,contaf,contafa
  INTEGER :: ini,nf,inia,nfa
  INTEGER :: contin,contout,maxvec,mini
  INTEGER :: contina,contouta
  INTEGER :: continf,contoutf,continfa,contoutfa
  INTEGER :: contap,contam,nm,np,l,m,which
  REAL(8) :: min,gtest
  REAL(8) :: pi,a,g
  REAL(8) :: maxX,maxY,maxZ
  REAL(8) :: km,kc,p
  REAL(8) :: kcBand,band
  REAL(8) :: kmin,kmax
  REAL(8) :: dumy,dx,dy,dz,circun,circuna,circunf,circunfa
  REAL(8) :: gradx,grady,gradz,grad
  REAL(8) :: gradxf,gradyf,gradzf,gradf
  REAL(8) :: g2x,g2y,g2z
  REAL(8) :: g2xf,g2yf,g2zf
  REAL(8) :: aexact,error,errorall
  REAL(8) :: efermi,x,y,z,BB
  REAL(8) :: dum,peso,am,ap
  REAL(8) :: xx,yy,zz,wei,Be
  REAL(8) :: vec,tol
!  REAL(8) :: zm,xm,ym1,ym2
!  REAL(8) :: xp,yp1,yp2,zp
  REAL(8) :: yminp,ymaxp
! REAL(8) :: xminm,xmaxm,yminm,ymaxm,zminm,zmaxm
!  REAL(8) :: z1,z2
!  REAL(8) :: kxyz
  !!                                                      w=48
  !!         wxyz: w=(c->counter,a->area), x=(a->all k,i->inner k), y=(i->below E_F,o->above E_F), z=(m->grad=-1,p->grad=+1)
  REAL(8) :: caom,aaom,caim,aaim
  REAL(8) :: caop,aaop,caip,aaip
  INTEGER, ALLOCATABLE :: ind(:),jnd(:),knd(:),yes(:,:,:)
  INTEGER, ALLOCATABLE :: ii(:),jj(:),kk(:)
  pi=acos(-1.)
  !reads data
  !        1 2  3 4  5 6      7      8    9  10 11  12
  read(*,*)N,Nk,a,kc,p,efermi,kcBand,band,nm,np,tol,which
  !write(*,*)'############################################'
  !write(*,*)N,Nk,a,kc,p,efermi,kcBand,band,nm,np,tol,which
  !write(*,*)'############################################'
  allocate (rxm(nm),rym(nm),rzm(nm))
  allocate (rxp(np),ryp(np),rzp(np))
  aNk=10*Nk
  !allocate (energy(-aNk:aNk),ii(-aNk:aNk),jj(-aNk:aNk),kk(-aNk:aNk),rxx(-aNk:aNk),ryy(-aNk:aNk),rzz(-aNk:aNk))
  allocate (energy(Nk))
  allocate (ii(-aNk:aNk),jj(-aNk:aNk),kk(-aNk:aNk),rxx(-aNk:aNk),ryy(-aNk:aNk),rzz(-aNk:aNk))
  !in and out k-vector
  maxvec=100000
  allocate (d(maxvec),ikvin(maxvec) ,ikvout(maxvec) ,ikvinf(maxvec) ,ikvoutf(maxvec))
  allocate (          ikvina(maxvec),ikvouta(maxvec),ikvinfa(maxvec),ikvinfap(maxvec),ikvoutfa(maxvec),ikvoutfap(maxvec))
  do i=1,maxvec
     allocate (ikvin(i)%kvin(3),  ikvout(i)%kvout(3)  ,ikvinf(i)%kvinf(3)  ,ikvoutf(i)%kvoutf(3))
     allocate (ikvina(i)%kvina(3),ikvouta(i)%kvouta(3),ikvinfa(i)%kvinfa(3),ikvinfap(i)%kvinfap(3),ikvoutfa(i)%kvoutfa(3),ikvoutfap(i)%kvoutfap(3))
  end do
  !
  !Grid lower (ini) and upper (nf) bounds
  ini=-1 !same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  nf=N+1 !same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  !
  !allocate inia=ini-1 and nfa=nf+1 to have enough space
  inia=ini-1 !=-2 same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  nfa=nf+1   !=N+2 same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  allocate (rx(inia:nfa),ry(inia:nfa),rz(inia:nfa))
  allocate (B(inia:nfa,inia:nfa,inia:nfa),BF(inia:nfa,inia:nfa,inia:nfa))
  allocate (weight(inia:nfa,inia:nfa,inia:nfa))
  allocate (yes(inia:nfa,inia:nfa,inia:nfa))
  allocate (ind(inia:nfa),jnd(inia:nfa),knd(inia:nfa))
  allocate (tagi(inia:nfa,inia:nfa,inia:nfa))
  !initialises
  continfa=0
  contoutfa=0
  contin=0
  contout=0
  continf=0
  contoutf=0
  continfa=0
  contoutfa=0
  !
  weight(:,:,:)=0.
  rx(:)=0.
  ry(:)=0.
  rz(:)=0.
  B(:,:,:)=10.
  BF(:,:,:)=10.
  yes(:,:,:)=10
  tagi(:,:,:)=0
  ii(:)=-2
  jj(:)=-2
  kk(:)=-2
  ind(:)=-2
  jnd(:)=-2
  knd(:)=-2
  rxx(:)=0.
  ryy(:)=0.
  rzz(:)=0.
  ! g-values
  g=2*pi/a
  maxX=g
  maxY=g
  maxZ=g
  !k_min and k_max
  kmin=kc-kc*p/100.
  kmax=kc+kc*p/100.
  !counters
  contam=0
  contap=0
  if ( 1 .eq. 1 ) then
     !
     !Generates Data-D
     !
     ! B(i,j,k) & BF(i,j,k) are used for the mask that separates inside (1) from outside (0) of the surface
     !
     do kpi=acos(-1.)
  !reads data
  !        1 2  3 4  5 6      7      8    9  10 11  12
  read(*,*)N,Nk,a,kc,p,efermi,kcBand,band,nm,np,tol,which
  !write(*,*)'############################################'
  !write(*,*)N,Nk,a,kc,p,efermi,kcBand,band,nm,np,tol,which
  !write(*,*)'############################################'
  allocate (rxm(nm),rym(nm),rzm(nm))
  allocate (rxp(np),ryp(np),rzp(np))
  aNk=10*Nk
  !allocate (energy(-aNk:aNk),ii(-aNk:aNk),jj(-aNk:aNk),kk(-aNk:aNk),rxx(-aNk:aNk),ryy(-aNk:aNk),rzz(-aNk:aNk))
  allocate (energy(Nk))
  allocate (ii(-aNk:aNk),jj(-aNk:aNk),kk(-aNk:aNk),rxx(-aNk:aNk),ryy(-aNk:aNk),rzz(-aNk:aNk))
  !in and out k-vector
  maxvec=100000
  allocate (d(maxvec),ikvin(maxvec) ,ikvout(maxvec) ,ikvinf(maxvec) ,ikvoutf(maxvec))
  allocate (          ikvina(maxvec),ikvouta(maxvec),ikvinfa(maxvec),ikvinfap(maxvec),ikvoutfa(maxvec),ikvoutfap(maxvec))
  do i=1,maxvec
     allocate (ikvin(i)%kvin(3),  ikvout(i)%kvout(3)  ,ikvinf(i)%kvinf(3)  ,ikvoutf(i)%kvoutf(3))
     allocate (ikvina(i)%kvina(3),ikvouta(i)%kvouta(3),ikvinfa(i)%kvinfa(3),ikvinfap(i)%kvinfap(3),ikvoutfa(i)%kvoutfa(3),ikvoutfap(i)%kvoutfap(3))
  end do
  !
  !Grid lower (ini) and upper (nf) bounds
  ini=-1 !same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  nf=N+1 !same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  !
  !allocate inia=ini-1 and nfa=nf+1 to have enough space
  inia=ini-1 !=-2 same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  nfa=nf+1   !=N+2 same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
  allocate (rx(inia:nfa),ry(inia:nfa),rz(inia:nfa))
  allocate (B(inia:nfa,inia:nfa,inia:nfa),BF(inia:nfa,inia:nfa,inia:nfa))
  allocate (weight(inia:nfa,inia:nfa,inia:nfa))
  allocate (yes(inia:nfa,inia:nfa,inia:nfa))
  allocate (ind(inia:nfa),jnd(inia:nfa),knd(inia:nfa))
  allocate (tagi(inia:nfa,inia:nfa,inia:nfa))
  !initialises
  continfa=0
  contoutfa=0
  contin=0
  contout=0
  continf=0
  contoutf=0
  continfa=0
  contoutfa=0
                                                                                                                                                                                         101,24         6%

=ini,nf
        rz(k)=real(maxZ*k)/real(N)
        do i=ini,nf
           rx(i)=real(maxX*i)/real(N)
           do j=ini,nf
              ry(j)=real(maxY*j)/real(N)
              !this signals all grid points
              !yes(i,j,k)=0
              ! This complies with BCC
              ! planes GP-GH & GN-GH & GP-GN
              if ( ( (rx(i) - ry(j) ) .le. 0 ) .and. ( ( rz(k)-rx(i) ) .le. 0 ) .and. ( ( rx(i)+ry(j) ) .le. g ) ) then
                 km=sqrt(rx(i)**2+ry(j)**2+rz(k)**2)
                 ! This selects a small region of k-values
                 if ((km.ge.kmin).and.(km.le.kmax)) then
                    !Gets B function
                    ! km < kc inside k-points
                    if ( km .lt. kc ) then
                       B(i,j,k)=1. !inside
                    end if
                    ! km > kc outside k-points
                    if ( km .gt. kc ) then
                       B(i,j,k)=0. !outside
                    end if
                    !
                    !weights
                    !tag=1 => inner points (w=48)
                    !tag=3 => GNH plane    (w=24)
                    !tag=5 => GNP plane    (w=24)
                    !tag=6 => GPH plane    (w=12)
                    !
                    if ((rz(k).gt.0).and.((ry(j)-rx(i)).gt.0).and.(rz(k)-rx(i)).lt.0) then
                       w=48
                       weight(i,j,k)=w
                       tag=1
                       tagi(i,j,k)=tag
                       yes(i,j,k)=1
                       write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                    end if
                    !PLANES-Down
                    !GN-GH-plane including the GN & GH lines: nov-18 OK!
                    if(rz(k).eq.0.) then
                       w=24
                       weight(i,j,k)=w
                       tag=3 !OK
                       tagi(i,j,k)=tag
                       yes(i,j,k)=1
                       write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                    end if
                    !GN-GP-plane including the GP line and excluding the GN line: nov-18-2021 OK!
                    if(((ry(j)-rx(i)).eq.0).and.(rz(k).gt.0)) then
                       w=24
                       weight(i,j,k)=w
                       tag=5 !OK
                       tagi(i,j,k)=tag
                       yes(i,j,k)=1
                       write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                    end if
                    !GP-GH-plane without the GP & GH lines: nov-18-2021 OK!
                    if(((rz(k)-rx(i)).eq.0).and.(rz(k).gt.0).and.((ry(j)-rx(i)).gt.0)) then
                       w=12
                       weight(i,j,k)=w
                       tag=6 !OK
                       yes(i,j,k)=1
                       tagi(i,j,k)=tag
                       write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                    end if
                    ! we don't consider the different weights of the lines, since
                    ! it is only a point that will contribute at most to the gradient
                    ! and the error of not taking the correct (reduced) weight is negligible
                    if ( 1 .eq. 2 ) then
                       !Gamma Point => No points since all spherical surfaces have a finite radii: oct-29 OK!
                       if((rx(i).eq.0.).and.(ry(j).eq.0.).and.(rz(k).eq.0.)) then
                          w=1
                          weight(i,j,k)=w
                          tag=2 !OK
                          write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !PLANES-Down
                       !NP-NH-plane => all spherical surfaces are within the IBZ, thus NO points touch this plane: oct-29 OK!
                       if((rx(i)+ry(j)).eq.g) then
                          w=24
                          weight(i,j,k)=w
                          tag=4
                          !write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !GN-GP-plane: oct-29 OK!
                       !PLANES-Up
                       !LINES-Down
                       !GN-line: included in GN-GP plane. The weight of the line is smaller,
                       !but the number of points is so small that should make a tiny difference the is neglected: oct-29 OK!
                       if( ((rz(k)-rx(i)+ry(j)).eq.0).and.(rz(k).eq.0)) then
                          w=12
                          weight(i,j,k)=w
                          tag=7
                          !write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !NH-line has no points: oct-29 OK!
                       if(((rx(i)+ry(j)).eq.g).and.(rz(k).eq.0)) then
                          w=24
                          weight(i,j,k)=w
                          tag=8
                          !write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !PH-line has no points: oct-29 OK!
                       if(((rx(i)+ry(j)).eq.g).and.(rz(k).eq.rx(i))) then
                          w=24
                          weight(i,j,k)=w
                          tag=9
                          !write(12,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !
                       !mv fort.12 $sym/$out12=$case.kcartesian_$nk
                       !
                       !This signals the valid BCC k-points for the already calculated k-points
                       !in the previous call of new-fermi-surface-via-kf.sh and in the 1st part of this program
                       !that generates the k-points according to
                       !kmin < kc < kmax
                       !yes(i,j,k)=1
                       !
                    end if !( 1 .eq. 2)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-34 for Grad=-1  --DOWN
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if ( band .eq. 34 ) then
                       !THIS REFINEMENT WORKS ONLY FOR EF=-1.78
                       !A DAERING YOUNG SOULD MUST DO THE SAME TRICK
                       !OF BISECTING WITH PLANES AND BOUNDARIES
                       !AS A FUNCTION OF EF
                       !grad=-1
                       !kxyz=sqrt(rx(i)**2+ry(j)**2+rz(k)**2)
                       !if ( (kxyz.ge..04) .and. (kxyz.le..047) .and. (rz(k).ge.0.) ) then
                       if ( rz(k).ge.0. ) then
                          contam=contam+1
                          !fort.21 => OriginalMesh k-points for grad=-1 @Pockets
                          write(21,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !
                    end if !( band .eq. 34 )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-34 for Grad=-1  --UP
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-35 for Grad=-1  --DOWN
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if ( band .eq. 35 ) then
                       !THIS REFINEMENT WORKS ONLY FOR EF=-1.78
                       !A DAERING YOUNG SOULD MUST DO THE SAME TRICK
                       !OF BISECTING WITH PLANES AND BOUNDARIES
                       !AS A FUNCTION OF EF
                       !zx-plane & spheres
                       !kxyz=sqrt(rx(i)**2+ry(j)**2+rz(k)**2)
                       !z1=.01-1.2*(rx(i)-.045)
                       !z2=.045-1.2*(ry(j)-.03)
                       !rep f u ($3<=z2($2)&&$3<=z1($1)&&k($1,$2,$3)<=.063&&k($1,$2,$3)>=.045?$1:1/0):2:3 w p pt 7 ps .3 lc rgb 'red'  t 'p&s'
                       !if ( (rz(k).le.z2) .and. (rz(k).le.z1) .and. (kxyz.le..063) .and. (kxyz.ge..045) .and. (rz(k).ge.0.) ) then
                       if ( rz(k).ge.0. ) then
                          contam=contam+1
                          !fort.21 => OriginalMesh k-points for grad=-1 @Pockets
                          write(21,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !
                    end if !( band .eq. 35 )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-35 for Grad=-1  --UP
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-36 for Grad= +/- 1  --DOWN
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if ( band .eq. 36 ) then
                       !THIS REFINEMENT WORKS ONLY FOR EF=-1.78
                       !A DAERING YOUNG SOULD MUST DO THE SAME TRICK
                       !OF BISECTING WITH PLANES AND BOUNDARIES
                       !AS A FUNCTION OF EF
                       !grad=-1
                       !zm=.027
                       !xm=.045
                       !ym1=.15
                       !ym2=.21
                       !if ( (rz(k).le.zm).and.(rx(i).le.xm).and.(ry(j).ge.ym1).and.(ry(j).le.ym2) ) then
                       if ( (rx(i).ge.0).and.(ry(j).ge.0).and.(rz(k).ge.0) ) then
                          contam=contam+1
                          !fort.21 => OriginalMesh k-points for grad=-1 @Pockets
                          write(21,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !end if
                       !grad=+1
                       !xp=.06
                       !yp1=.0475
                       !yp2=.0725
                       !zp=.06
                       !if ( (ry(j).ge.yp1).and.(ry(j).le.yp2).and.(rx(i).le.xp).and.(rz(k).le.zp).and.(rz(k).ge.0) )then
                       if ( (rx(i).ge.0).and.(ry(j).ge.0).and.(rz(k).ge.0) ) then
                          contap=contap+1
                          !fort.22 => OriginalMesh k-points for grad=+1 @GammaPoint
                          write(22,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !end if
                       !
                    end if !( band .eq. 36 )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-36 for Grad= +/- 1  --UP
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-37 for Grad= +/- 1  --DOWN
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if ( band .eq. 37 ) then
                       !THIS REFINEMENT WORKS ONLY FOR EF=-1.78
                       !A DAERING YOUNG SOULD MUST DO THE SAME TRICK
                       !OF BISECTING WITH PLANES AND BOUNDARIES
                       !AS A FUNCTION OF EF
                       !grad=+1
                       !y-range-grad=+1
                       yminp=0.05
                       ymaxp=0.08
                       !
                       !if( ( ry(j) .ge. yminp ).and.( ry(j) .le. ymaxp) .and. ( rz(k) .ge. 0. ) ) then
                          if( rz(k) .ge. 0. ) then
                          contap=contap+1
                          !fort.22 => OriginalMesh k-points for grad=+1 @baby duck
                          write(22,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !
                       !grad=-1
                       !x-range-grad=-1
                       !xminm=0.00
                       !xmaxm=0.17
                       !y-range-grad=-1
                       !yminm=0.07
                       !ymaxm=0.19
                       !z-range-grad=-1
                       !zminm=0.00
                       !zmaxm=0.08
                       !Bounding planes
                       !z1=.06-(.07/(.16-.08))*(ry(j)-.07)
                       !z2=.12-(.07/(.16-.08))*(ry(j)-.07)
                       !form gnuplot
                       !(  $3>=z1($2)    &&  $3<=z2($2)    &&  $1>=xminm       &&$1<=xmaxm          &&   $2>=yminm        &&  $2<=ymaxm        &&  $3>=zminm        &&  $3<=zmaxm?$1:1/0)
                       !if( (rz(k).gt.z1).and.(rz(k).le.z2).and.(rx(i).ge.xminm).and.(rx(i).le.xmaxm).and.(ry(j).ge.yminm).and.(ry(j).le.ymaxm).and.(rz(k).ge.zminm).and.(rz(k).le.zmaxm) ) then
                       if( rz(k) .ge. 0. ) then
                          contam=contam+1
                          !fort.21 => OriginalMesh k-points for grad=-1 @GammaPoint=>mama duck
                          write(21,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k
                       end if
                       !
                    end if !( band .eq. 37 )
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! selects the k-points around band-37 for Grad= +/- 1  --UP
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 end if ! kmin < km < kmax
              end if !plane GP-GH & GN-GH & GP&GN
           end do !j
        end do !i
     end do !k
     close(22)
     close(21)
     !
     !Generates Data-Up
     !
  end if
  !write(*,*)'stopping @bands-fs-bz-bcc-grid-integrate.f90'
  !stop
  !
  !reads data-----Down
  !From the coarse-grained first FS data generated by
  !fs-bz-bcc-grid-integrate.sh
  !and stored @
  !  symmetries/y2c3.fs-ibz-in-cg-?-r-?-k-?-1-? @ Grad=-1  -> fort.244
  !  symmetries/y2c3.fs-ibz-in-cg-?-r-?-k-?+1-? @ Grad=+1  -> fort.245
  !
  !A denser grid of k-points will be calculated
  !JUST around above coarse-grained FS-k-points.
  !
  !Then, we must run_tiniba.sh with these points,
  !to generate E(k) from which a refined FS would be obtained.
  !
  ! valid for all bands
  open (244, file='fort.244-Ncg-151-N-201')

  do l=1,nm
     !fort.244 -> symmetries/y2c3.fs-ibz-in-?-1-N1 (57) @ Grad=-1
     read(244,*)rxm(l),rym(l),rzm(l),dumy,dumy,dx,dy,dz
  end do
  close(244)
  write(*,*)'fort.244 with',nm,'was read'
! only bands 36 and 37 have two FS
  if ( (band .eq. 36) .or. (band .eq. 37) ) then
    do l=1,np
       fort.245 -> symmetries/y2c3.fs-ibz-in-?+1-N1 (59) @ Grad=+1
       read(245,*)rxp(l),ryp(l),rzp(l),dumy,dumy,dx,dy,dz
    end do
    close(245)
    write(*,*)'fort.245 with',np,'was read'
  end if
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !Read fort.21 (grad=+1) and fort.22 (grad=-1) to calculate the vector
  write(*,*)'        *****'
  write(*,*)'        @bands-fs-bz-bcc-grid-integrate.f90'
  write(*,*)N,'divisions =>',contap,'k-p grad=+1 @CorseGFS',np
  write(*,*)'                      and',contam,'k-p grad=-1 @CorseGFS',nm
  !write(*,*)'stopping @bands-fs-bz-bcc-grid-integrate.f90'
  !stop
  !
  !grad=-1-DOWN
  !valid for all bands
  !
  !New tol value:
  !we use
  !radius=sqrt(dx^2 + dy^2 + dz^2)
  !which gives the radius along the
  !diagonal of a boxel of the Cartesian BCC grid.
  !We take
  !tol=radius
  !which would be equivalent to an small number of boxels around
  !the (xx,yy,zz)-kpoint of the CG-grid. As we increase the refinement of
  !the grid, tol=r would be smaller than the original fixed value of .0075 that
  !works just fine, but has the disadvantage that as we increase the refinement
  !it would include k-points that would be farther away from the CG-FS, thus
  !wasting computational precious time!
  dx=abs(rx(1)-rx(2))
  dy=abs(ry(1)-ry(2))
  dz=abs(rz(1)-rz(2))
  !Comment bellow if you want to externally read the value of tol
  tol=tol*sqrt(dx**2 + dy**2 + dz**2)
  !
  !write(321,*)tol,sqrt(dx**2 + dy**2 + dz**2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! Takes the previous CG as a seed for the next grid                         !!!!!
  !!!! Each k-point of the next grid is within tol-distance of the previous grid !!!!!
  !!!! @@@@@@@@@@@@PARALELIZE@@@@@@@@@@@
  !the rxm,rym,rzm are read from fort.244
  !read(244,*)rxm(l),rym(l),rzm(l),dumy,dumy,dx,dy,dz
  do l=1,contam
     !reads CG grid
     read(21,*)xx,yy,zz,tag,wei,Be,i,j,k
     do m=1,nm
        vec=sqrt( (xx-rxm(m))**2 + (yy-rym(m))**2 + (zz-rzm(m))**2 )
        write(851,*)l,m,vec
        if ( vec .le. tol ) then
           rx(i)=xx
           ry(j)=yy
           rz(k)=zz
           weight(i,j,k)=wei
           B(i,j,k)=Be
           !fort.751 => grad=-1 new k-points around FS@pocket
           write(751,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k,dx,dy,dz
        end if
     end do
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!
  write(*,*)'        @bands-fs-bz-bcc-grid-integrate.f90:end read21&write751=1 for all bands'
  close(751)
  close(851)
  !
  !grad=-1-UP
  !
  !grad=+1-DOWN
  !only bands 36 and 37 have two FS
  !
  if ( (band .eq. 36) .or. (band .eq. 37) ) then
     do l=1,contap
        !reads CG grid
        read(22,*)xx,yy,zz,tag,wei,Be,i,j,k
        do m=1,np
           vec=sqrt( (xx-rxp(m))**2 + (yy-ryp(m))**2 + (zz-rzp(m))**2 )
           write(850,*)l,m,vec
           if ( vec .le. tol ) then
              rx(i)=xx
              ry(j)=yy
              rz(k)=zz
              weight(i,j,k)=wei
              B(i,j,k)=Be
              !fort.750 => grad=+1 new k-points around FS@GammaPoint
              write(750,76)rx(i),ry(j),rz(k),tag,weight(i,j,k),B(i,j,k),i,j,k,dx,dy,dz
           end if
        end do
     end do
     write(*,*)'        @band-36-fs-bz-bcc-integrate.f90:end read-22&write-750 band 36 or 37'
     close(750)
     close(850)
  end if
  !write(*,*)'@bands-fs-bz-bcc-grid-integrate.f90: aquiBOYnas#1'
  !write(*,*)'stopping @bands-fs-bz-bcc-grid-integrate.f90'
  !stop
  !
  !grad=+1-UP
  !
  !
  !reads data from-----Up
  !
  ! 2nd BLOCK-DOWN
  !
  ! reads data from symmetries/case.kband_band_Nk  -----Down
  !
  write(*,*)'        @bands-fs-bz-bcc-grid-integrate.f90: creating refined grid: step 1'
  if ( which .eq. 2 ) then
     BF(:,:,:)=10
     yes(:,:,:)=0
     do i=1,Nk
        !fort.1=out=$sym/$case.kband_$band\_$Nk has all weights as it must = $sym/$case.kcartesian\_$Nk !
        !tag=1 => inner points (w=48)
        !tag=3 => GNH plane    (w=24)
        !tag=5 => GNP plane    (w=24)
        !tag=6 => GPH plane    (w=12)
        !dumy stands for dx dy dz which are calculated below
        !        1 2 3 4   5 6  7     8     9     10   11   12   13
        read(1,*)x,y,z,tag,w,BB,ii(i),jj(i),kk(i),dumy,dumy,dumy,energy(i)
        !fort.1 = fort.751 for x,y,z,tag,w,BB
        !WARNING: BB comes from km < kc < km
        !(x,y,z)
        rxx(ii(i))=x
        ryy(jj(i))=y
        rzz(kk(i))=z
        !(i,j,z=k)
        ind(ii(i))=ii(i)
        jnd(jj(i))=jj(i)
        knd(kk(i))=kk(i)
        !
        weight(ii(i),jj(i),kk(i))=w
        !write(*,*)weight(ii(i),jj(i),kk(i)),ii(i),jj(i),kk(i)
        !valid BCC points determined @ Generates Data (above)
        yes(ii(i),jj(i),kk(i))=1
!!!!!!! FOR THE IDEAL SURFACE-Down !!
        km=sqrt(x**2+y**2+z**2)
        if ( km .lt. kc ) then                !
           B(ii(i),jj(i),kk(i)) = 1. !inside  !
        end if                                !
        if ( km .gt. kc ) then                !
           B(ii(i),jj(i),kk(i)) = 0. !outside !
        end if                                !
!!!!!!! FOR THE IDEAL SURFACE-Up !!
!!!!!!! FOR THE FERMI SURFACE-Down !!
!!!!!!! Grad_i BF could be +1 or -1 for i=x,y,z
        !write(*,*)i,energy(i)
        if (energy(i).lt.efermi) then    !
           BF(ii(i),jj(i),kk(i)) = 1.    !inside
        end if                           !
        if (energy(i).gt.efermi) then    !
           BF(ii(i),jj(i),kk(i)) = 0.    !outside
        end if                           !
!!!!!!! FOR THE FERMI SURFACE-Up !!
        !bellow are needed for debugging pourposes only
        !fort.29 = fort.1=fort.12 as it must!
        !           1         2         3         4   5                         6                     7     8     9     10                     11
        !write(29,29)rx(ii(i)),ry(jj(i)),rz(kk(i)),tag,weight(ii(i),jj(i),kk(i)),B(ii(i),jj(i),kk(i)) ,ii(i),jj(i),kk(i),yes(ii(i),jj(i),kk(i)),energy(i) !Ideal
        !           1         2         3         4   5                         6                     7     8     9     10                     11
        write(30,29)rx(ii(i)),ry(jj(i)),rz(kk(i)),tag,weight(ii(i),jj(i),kk(i)),BF(ii(i),jj(i),kk(i)),ii(i),jj(i),kk(i),yes(ii(i),jj(i),kk(i)),energy(i) !FS
     end do
  end if !( which .eq. 2 )
  !
  ! 2nd BLOCK-UP
  !
  !write(*,*)'@bands-fs-bz-bcc-grid-integrate.f90: aquiBOYnas#2'
  !write(*,*)'stopping @bands-fs-bz-bcc-grid-integrate.f90'
  !stop
  !
  ! 3rd BLOCK-DOWN
  !
  !
  !Integrates Data-Down
  !
  !reads data from symmetries/case.kband_band_Nk  -----Up
  !
  !Takes the gradient using the
  !central-difference scheme
  !Steps of 2 => the gradient goes from a grid point just bellow the surface to a grid point just above the surface
  !              such that the grid point in between is the surface, i.e.
  !gradx=(B(i+1,j,k)-B(i-1,j,k))/2
  !dx=(rx(i+1)-rx(i-1))/2=rx(i+1)-rx(i), dy=(ry(i+1)-ry(i-1))/2=ry(i+1)-ry(i), dz=(rz(i+1)-rz(i-1))/2=rz(i+1)-rz(i)
  !where we divided by 2 as it must be since we skip 2-grid-points, i.e. from i-1 to i+1
  !
  write(*,*)'        @bands-fs-bz-bcc-grid-integrate.f90: creating refined grid: step 2'
  write(*,*)'        @bands-fs-bz-bcc-grid-integrate.f90: which=',which
  if ( which .eq. 2 ) then
     conta=0
     contafa=0
     contall=0
     contin=0
     contina=0
     contout=0
     contouta=0
     circun=0.
     circuna=0.
     !FS
     contaf=0
     continf=0
     contoutf=0
     circunf=0.
     !
     aexact=4.*pi*kc**2 !area of a sphere
     !this choice makes the surface to be as close as possible
     !to the planes that enclose the IBZ
     ini=ini+1 !=0 same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
     nf=nf-1   !=N same as /home/bms/tiniba/ver6.0/utils/grad-method/3D-central-difference-grid-integrator.f90
     do k=ini,nf
        !rz(k)=real(maxZ*k)/real(N)
        do i=ini,nf
           !rx(i)=real(maxX*i)/real(N)
           do j=ini,nf
              !ry(j)=real(maxY*j)/real(N)
              !
              !aquiBOYnas
              !if (weight(i,j,k).gt.0) write(*,*)'w#3rdBlock#1=',weight(i,j,k),i,j,k
              ! This complies with BCC
              !Nov-19 we must check this, it better work since it is the same as that of ideal-bz-bcc-grid-integrate.f90
              ! planes GP-GH & GN-GH & GP-GN
              if ( ( (rx(i) - ry(j) ) .le. 0 ) .and. ( ( rz(k)-rx(i) ) .le. 0 ) .and. ( ( rx(i)+ry(j) ) .le. g ) ) then
                 km=sqrt(rx(i)**2+ry(j)**2+rz(k)**2)
                 if ((km.ge.kmin).and.(km.le.kmax)) then
                    ! fort.1 (above) = fort.8
                    ! fort.803 = fort.804 @ ideal-bz-bcc-grid-integrate.f90
                    !             1     2     3     4           5             6          7
                    !write(803,803)rx(i),ry(j),rz(k),tagi(i,j,k),weight(i,j,k),yes(i,j,k),B(i,j,k)
                    !                 end if
                    !              end if
                    !
                    ! we divide over two since (dx,dy,dz) are the grid spacings
                    ! and we are taking (i+1) & (i-1), => skiping 2 grid points!
                    !dx=(rx(i+1)-rx(i-1))/2.
                    !dy=(ry(j+1)-ry(j-1))/2.
                    !dz=(rz(k+1)-rz(k-1))/2.
                    !Better like this since is equivalent and  there is NO confusion
                    dx=rx(i+1)-rx(i)
                    dy=ry(j+1)-ry(j)
                    dz=rz(k+1)-rz(k)
                    !
                    !the gradient is taken from inside  to outside of the ellipsoid of revolution
                    !since B&BF => B&BF=1 inside & B&BF=0 outside => grad B&BF < 0
                    !
                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ideal-Down => fort.69 WARNING: WE SHOULDN'T NEED THE IDEAL SURFACE FOR THE FS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOT USED-Down
                    if ( 1 .eq. 2 ) then
                       gradx=(B(i+1,j,k)-B(i-1,j,k))
                       grady=(B(i,j+1,k)-B(i,j-1,k))
                       gradz=(B(i,j,k+1)-B(i,j,k-1))
                       if (tagi(i,j,k).eq.5) then
                          !for the GN-GP plane we have to rotate -45^o the gradient=> gradx=grady
                          gradx=B(i+1,j+1,k)-B(i-1,j-1,k)
                          grady=B(i+1,j+1,k)-B(i-1,j-1,k)
                          !checkup
                          !write(3,3)gradx,grady,gradz,weight(i,j,k),i,j,k,rx(i),ry(j),rz(k),tagi(i,j,k),B(i,j,k)
                       end if
                       !checkup: once is working remove below
                       !tag=1 => inner points (w=48)
                       !tag=3 => GNH plane    (w=24)
                       !tag=5 => GNP plane    (w=24)
                       !tag=6 => GPH plane    (w=12)
                       !         1     2     3     4             5 6 7 8     9     10    11          12
                       !write(3,3)gradx,grady,gradz,weight(i,j,k),i,j,k,rx(i),ry(j),rz(k),tagi(i,j,k),B(i,j,k)
                       !         1 2 3   4 5   6          7          8
                       !write(7,8)i,j,k-1,k,k+1,B(i,j,k+1),B(i,j,k-1),gradz
                       !The ideal FS that has spherical points is generated such that Grad B = -1 always as it crosses the surface
                       gtest=-1. ! .eq.
                       if ( (gradx .eq. gtest).or.(grady .eq. gtest).or.(gradz .eq. gtest) ) then !=> option1
                          !finite difference approximation => divide over 2
                          g2x=gradx/2.
                          g2y=grady/2.
                          g2z=gradz/2.
                          grad=sqrt((g2x/dx)**2+(g2y/dy)**2+(g2z/dz)**2)
                          !
                          ! weight=48 => inner points
                          !
                          ! only inner points: this is a very good approximation since
                          !                    most of the points are inside the IBZ
                          !                    and its number scales as ~ kc^2.
                          !                    The points at the perimeter of the surface
                          !                    scale as ~ kc, and their contribution to the
                          !                    total surface should be less.
                          !
                          if ( 1 .eq. 1 ) then
                             if (weight(i,j,k) .eq. 48.) then
                                if (grad .ne. 0 ) then
                                   conta=conta+1
                                   !the integration is given by
                                   circun=circun+weight(i,j,k)*grad*(dx*dy*dz)
                                   if ( B(i,j,k) .eq. 0. ) then
                                      contout=contout+1
                                      !fort.69 => $sym/$out69 => symmetries/$case.ibz-ideal_$kn ($5==0 => out)
                                      !           1     2     3     4             5  6     7     8     9
                                      write(69,96)rx(i),ry(j),rz(k),weight(i,j,k),0.,gradx,grady,gradz,grad
                                      ikvout(contout)%kvout(1)=rx(i)
                                      ikvout(contout)%kvout(2)=ry(j)
                                      ikvout(contout)%kvout(3)=rz(k)
                                   end if
                                   if ( B(i,j,k) .eq. 1. ) then
                                      contin=contin+1
                                      !fort.69 => $sym/$out69 => symmetries/$case.ibz-ideal_$kn ($5==1 => in)
                                      !           1     2     3     4             5  6     7     8     9
                                      write(69,96)rx(i),ry(j),rz(k),weight(i,j,k),1.,gradx,grady,gradz,grad
                                      ikvin(contin)%kvin(1)=rx(i)
                                      ikvin(contin)%kvin(2)=ry(j)
                                      ikvin(contin)%kvin(3)=rz(k)
                                   end if
                                end if ! (grad .ne. 0)
                             end if ! (weight(i,j,k) .eq. 48.)
                          end if ! ( 1 .eq. 1 )
                          !
                          ! all weights => all-k-points-DOWN
                          !
                          if (weight(i,j,k) .gt. 1.) then
                             if (grad .ne. 0 ) then
                                contall=contall+1
                                !the integration is given by
                                circuna=circuna+weight(i,j,k)*grad*dx*dy*dz
                                if ( B(i,j,k) .eq. 0. ) then
                                   contouta=contouta+1
                                   !fort.71 => $sym/$out71 => symmetries/$case.ibz-ideala_$kn ($5==0 => out)
                                   !           1     2     3     4             5        6     7     8     9
                                   write(71,17)rx(i),ry(j),rz(k),weight(i,j,k),B(i,j,k),gradx,grady,gradz,grad
                                   ikvouta(contouta)%kvouta(1)=rx(i)
                                   ikvouta(contouta)%kvouta(2)=ry(j)
                                   ikvouta(contouta)%kvouta(3)=rz(k)
                                end if
                                if ( B(i,j,k) .eq. 1. ) then
                                   contina=contina+1
                                   !fort.71 => $sym/$out71 => symmetries/$case.ibz-ideala_$kn ($5==1 => in)
                                   !           1     2     3     4             5        6     7     8     9
                                   write(71,17)rx(i),ry(j),rz(k),weight(i,j,k),B(i,j,k),gradx,grady,gradz,grad
                                   ikvina(contina)%kvina(1)=rx(i)
                                   ikvina(contina)%kvina(2)=ry(j)
                                   ikvina(contina)%kvina(3)=rz(k)
                                end if
                             end if ! ( grad .ne. 0 )
                          end if ! (weight(i,j,k) .gt. 1.)
                          !
                       end if ! ( (gradx .eq. gtest).or.(grady .eq. gtest).or.(gradz .eq. gtest) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ideal-Up => fort.69 (w=48) & fort.71 (w!=1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOT USED-UP
                    end if !( 1 .eq. 2 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !
                    if ( 1 .eq. 1 ) then
                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FS-Down => fort.70
                       !if (weight(i,j,k).gt.0) write(*,*)'w3rd#2=',weight(i,j,k),i,j,k
                       !
                       gradxf=(BF(i+1,j,k)-BF(i-1,j,k))
                       gradyf=(BF(i,j+1,k)-BF(i,j-1,k))
                       gradzf=(BF(i,j,k+1)-BF(i,j,k-1))
                       if (tagi(i,j,k).eq.5) then
                          !for the GN-GP plane we have to rotate -45^o the gradient=> gradx=grady
                          gradxf=(BF(i+1,j+1,k)-BF(i-1,j-1,k))
                          gradyf=(BF(i+1,j+1,k)-BF(i-1,j-1,k))
                          !checkup
                          !write(3,3)gradxf,gradyf,gradzf,weight(i,j,k),i,j,k,rx(i),ry(j),rz(k),tagi(i,j,k),B(i,j,k)
                       end if
                       !checkup: once is working remove below
                       !tag=1 => inner points (w=48)
                       !tag=3 => GNH plane    (w=24)
                       !tag=5 => GNP plane    (w=24)
                       !tag=6 => GPH plane    (w=12)
                       !          1      2      3      4             5 6 7 8     9     10    11          12
                       !write(4,3)gradxf,gradyf,gradzf,weight(i,j,k),i,j,k,rx(i),ry(j),rz(k),tagi(i,j,k),BF(i,j,k)
                       !         1     2     3     4           5             6
                       write(4,4)rx(i),ry(j),rz(k),tagi(i,j,k),weight(i,j,k),BF(i,j,k)
                       write(41,41)rx(i),ry(j),rz(k),gradxf,gradyf,gradzf
                       ! fort.4 = fort.12 => original FULL GRID!!!
                       !         1 2 3   4 5   6          7          8
                       !write(7,8)i,j,k-1,k,k+1,B(i,j,k+1),B(i,j,k-1),gradz
                       !For the FS Grad B cloud be +1 or -1
                       ! => gtest=+/- 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! gtest=-1.-DOWN !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       gtest=-1. !^- => gradient going from BF(i,j,k)=1 to BF(i,j,k)=0 thus grad=-1.
                       if ( (gradxf .eq. gtest) .or. (gradyf .eq. gtest) .or. (gradzf .eq. gtest) ) then !=> option#1
                          !finite difference approximation => divide over 2
                          g2xf=gradxf/2.
                          g2yf=gradyf/2.
                          g2zf=gradzf/2.
                          gradf=sqrt((g2xf/dx)**2+(g2yf/dy)**2+(g2zf/dz)**2)
                          !
                          ! weight=48 => inner points
                          !
                          ! only inner points: this is a very good approximation since
                          !                    most of the points are inside the IBZ
                          !                    and its number scales as ~ kc^2.
                          !                    The points at the perimeter of the surface
                          !                    scale as ~ kc, and their contribution to the
                          !                    total surface should be less.
                          !
                          !!!!!!!!!!!!!!!!!!!!NOT USED-DOWN
                          if ( 1 .eq. 2 ) then
                             if (weight(i,j,k) .eq. 48.) then
                                if (gradf .ne. 0 ) then
                                   contaf=contaf+1
                                   !the integration is given by
                                   circunf=circunf+weight(i,j,k)*gradf*dx*dy*dz
                                   if ( BF(i,j,k) .eq. 0. ) then
                                      contoutf=contoutf+1
                                      !fort.70 => $sym/$out70 => symmetries/$case.ibz_$kn ($5==0 => out)
                                      !           1     2     3     4             5         6      7      8      9
                                      write(70,96)rx(i),ry(j),rz(k),weight(i,j,k),BF(i,j,k),gradxf,gradyf,gradzf,gradf
                                      ikvoutf(contoutf)%kvoutf(1)=rx(i)
                                      ikvoutf(contoutf)%kvoutf(2)=ry(j)
                                      ikvoutf(contoutf)%kvoutf(3)=rz(k)
                                   end if
                                   if ( BF(i,j,k) .eq. 1. ) then
                                      continf=continf+1
                                      !fort.70 => $sym/$out70 => symmetries/$case.ibz_$kn ($5==0 => out)
                                      !           1     2     3     4             5         6      7      8      9
                                      write(70,96)rx(i),ry(j),rz(k),weight(i,j,k),BF(i,j,k),gradxf,gradyf,gradzf,gradf
                                      ikvinf(continf)%kvinf(1)=rx(i)
                                      ikvinf(continf)%kvinf(2)=ry(j)
                                      ikvinf(continf)%kvinf(3)=rz(k)
                                   end if
                                end if ! (gradf .ne. 0)
                             end if ! (weight(i,j,k) .eq. 48)
                          end if !( 1 .eq. 2  )
                          !!!!!!!!!!!!!!!!!!!!NOT USED-UP
                          !
                          ! all weights => all-k-points & gtest=-1-Down
                          !
                          !write(*,*)'w3rd#3=',weight(i,j,k),i,j,k
                          if (weight(i,j,k) .gt. 0.) then !all
                             if (gradf .ne. 0) then
                                contafa=contafa+1
                                !the integration is given by
                                circunfa=circunfa+weight(i,j,k)*gradf*dx*dy*dz
                                if ( BF(i,j,k) .eq. 0. ) then !out
                                   !OJO#1
                                   if ( (band.eq.36) .or. (band.eq.37) ) then
                                      !A_{all,out}^- => BF(i,j,k)=0 (out) & gtest=-1 (^-) => aom
                                      caom=caom+1
                                      aaom=aaom+weight(i,j,k)*gradf*dx*dy*dz
                                      !out:grad=-1
                                      !gtest=-1. !^- => gradient going from BF(i,j,k)=1 to BF(i,j,k)=0 thus grad=-1.
                                      !            1     2     3     4             5     6  7  8  9
                                      write(200,96)rx(i),ry(j),rz(k),weight(i,j,k),gradf,dx,dy,dz,BF(i,j,k)
                                      ikvoutfa(caom)%kvoutfa(1)=rx(i)
                                      ikvoutfa(caom)%kvoutfa(2)=ry(j)
                                      ikvoutfa(caom)%kvoutfa(3)=rz(k)
                                      !
                                   end if
                                end if
                                if ( BF(i,j,k) .eq. 1. ) then !in
                                   !OJO#2
                                   if ( (band.eq.36) .or. (band.eq.37) ) then
                                      !A_{all,in}^- => BF(i,j,k)=1 (in) & gtest=-1 (^-) => aim
                                      caim=caim+1
                                      aaim=aaim+weight(i,j,k)*gradf*dx*dy*dz
                                      !in:grad=-1
                                      !            1     2     3     4             5     6  7  8  9
                                      write(201,96)rx(i),ry(j),rz(k),weight(i,j,k),gradf,dx,dy,dz,BF(i,j,k)
                                      !WARNING => split for gradf=+/-1??
                                      ikvinfa(caim)%kvinfa(1)=rx(i)
                                      ikvinfa(caim)%kvinfa(2)=ry(j)
                                      ikvinfa(caim)%kvinfa(3)=rz(k)
                                      !
                                   end if
                                end if
                             end if ! (gradf .ne. 0)
                          end if ! (weight(i,j,k) .gt. 0)
                          !
                          ! all weights & gtest=-1-Up
                          !
                       end if ! ( (gradxf .eq. gtest).or.(gradyf .eq. gtest).or.(gradzf .eq. gtest) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! gtest=-1.-UP !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! gtest=1.-DOWN !!!!!!!!
                       gtest=1. !^+ => gradient going from BF(i,j,k)=0 to BF(i,j,k)=1 thus grad=+1.
                       ! band 36 and 37 for Y2C3 has both grad=-1 and grad=+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       if ( (gradxf .eq. gtest).or.(gradyf .eq. gtest).or.(gradzf .eq. gtest) ) then ! => option1
                          !finite difference approximation => divide over 2
                          g2xf=gradxf/2.
                          g2yf=gradyf/2.
                          g2zf=gradzf/2.
                          gradf=sqrt((g2xf/dx)**2+(g2yf/dy)**2+(g2zf/dz)**2)
                          !
                          ! weight=48 => inner points
                          !
                          ! only inner points: this is a very good approximation since
                          !                    most of the points are inside the IBZ
                          !                    and its number scales as ~ kc^2.
                          !                    The points at the perimeter of the surface
                          !                    scale as ~ kc, and their contribution to the
                          !                    total surface should be less.
                          !
                          if ( 1 .eq. 2 ) then
                             if (weight(i,j,k) .eq. 48.) then
                                if ( gradf .ne. 0 ) then
                                   contaf=contaf+1
                                   circunf=circunf+weight(i,j,k)*gradf*dx*dy*dz
                                   if ( BF(i,j,k) .eq. 0. ) then !in
                                      contoutf=contoutf+1
                                      !fort.70 => $sym/$out70 => symmetries/$case.ibz_$kn ($5==0 => out)
                                      !           1     2     3     4             5         6      7      8      9
                                      write(72,96)rx(i),ry(j),rz(k),weight(i,j,k),BF(i,j,k),gradxf,gradyf,gradzf,gradf
                                      ikvoutf(contoutf)%kvoutf(1)=rx(i)
                                      ikvoutf(contoutf)%kvoutf(2)=ry(j)
                                      ikvoutf(contoutf)%kvoutf(3)=rz(k)
                                   end if
                                   if ( BF(i,j,k) .eq. 1. ) then !in
                                      continf=continf+1
                                      !fort.70 => $sym/$out70 => symmetries/$case.ibz_$kn ($5==1 => in)
                                      !           1     2     3     4             5         6      7      8      9
                                      write(72,96)rx(i),ry(j),rz(k),weight(i,j,k),BF(i,j,k),gradxf,gradyf,gradzf,gradf
                                      ikvinf(continf)%kvinf(1)=rx(i)
                                      ikvinf(continf)%kvinf(2)=ry(j)
                                      ikvinf(continf)%kvinf(3)=rz(k)
                                   end if
                                end if !(gradf .ne. 0)
                             end if !(weight(i,j,k) .eq. 48)
                          end if !( 1 .eq. 2 )
                          !
                          ! all weights => all-k-points & gtest=+1-Down
                          !
                          if (weight(i,j,k) .gt. 0.) then !all
                             if ( gradf .ne. 0 ) then
                                contafa=contafa+1
                                !the integration is given by
                                circunfa=circunfa+weight(i,j,k)*gradf*dx*dy*dz
                                if ( BF(i,j,k) .eq. 0. ) then
                                   !A_{all,out}^+ => BF(i,j,k)=0 (out) & gtest=+1 (^+) => caop,aaop
                                   caop=caop+1
                                   aaop=aaop+weight(i,j,k)*gradf*dx*dy*dz
                                   !output-file
                                   !out:grad=+1
                                   !            1     2     3     4             5     6  7  8  9
                                   write(202,96)rx(i),ry(j),rz(k),weight(i,j,k),gradf,dx,dy,dz,BF(i,j,k)
                                   ikvoutfap(caop)%kvoutfap(1)=rx(i)
                                   ikvoutfap(caop)%kvoutfap(2)=ry(j)
                                   ikvoutfap(caop)%kvoutfap(3)=rz(k)
                                end if
                                !
                                if ( BF(i,j,k) .eq. 1. ) then !in
                                   !A_{all,in}^+ => BF(i,j,k)=1 (in) & gtest=+1 (^+) => caip,aaip
                                   caip=caip+1
                                   aaip=aaip+weight(i,j,k)*gradf*dx*dy*dz
                                   !in:grad=+1
                                   !            1     2     3     4             5     6  7  8  9
                                   write(203,96)rx(i),ry(j),rz(k),weight(i,j,k),gradf,dx,dy,dz,BF(i,j,k)
                                   ikvinfap(caip)%kvinfap(1)=rx(i)
                                   ikvinfap(caip)%kvinfap(2)=ry(j)
                                   ikvinfap(caip)%kvinfap(3)=rz(k)
                                end if
                             end if !( gradf .ne. 0 )
                          end if !(weight(i,j,k) .gt. 0.)
                          !
                          ! all weights & gtest=+1-Up
                          !
                       end if ! ( (gradxf .eq. gtest).or.(gradyf .eq. gtest).or.(gradzf .eq. gtest) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! gtest=1.-UP !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if !(1 .eq. 1 )
                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FS-Up => fort.70
                    !
                    !
                    ! Integration-Up
                    !
!!!!!
                    !end if !(yes(i,j,k).eq.1)
                 end if ! if ((km.ge.kmin).and.(km.le.kmax)) then
              end if ! if ( ( (rx(i) - ry(j) ) .le. 0 ) .and. ( ( rz(k)-rx(i) ) .le. 0 ) .and. ( ( rx(i)+ry(j) ) .le. g ) ) then
           end do !j
        end do !i
     end do !k
     !!system call-Down
     !!
     !! We read data from fort.200 and fort.201 & fort.202 and fort.203
     !! to generate the FS k-points just-bellow and just-above,
     !! taking the just-bellow k-points as the FS
     !! We use a system call to data-grad-1+1.sh -r run
     !! to generate the FS k-points.
     !! WARNING: band 36 & 37 have a FS that is naturally splited in two parts
     !! and we work out the integration for each part, just for fun, since
     !! the total contribution is the sum of both parts.
     !! in unit=8 we have the value of kc that we need to run the script
     !! data-grad-1+1.sh -r run
     !!
     if ( 1 .eq. 1 ) then
        open(unit=8,file='fort.8')
        write(8,88)kcBand
        close(unit=8)
        write(*,*)'       *****                                            *****'
        write(*,*)'       ***** Begin system call: data-grad-1+1.sh -r run *****'
        !stop !aquiBOYnas
        call system("data-grad-1+1.sh -r run")
        !CALL SYSTEMQQ("C:\foldername\myprogram.exe argument1 argument2")
        !call execute_command_line ("mkdir example_folder")
        open(unit=123,file='fort.123')
        read(123,*)dummy,numm,nump
        close(unit=123)
        if ( kcBand .ne. 0 ) then
           write(*,*)'       ***** grad=+/- 1',nump,numm,' k-points'
        else
           write(*,*)'       ***** FS= ',numm,' Just above FS= ',nump,' k-points'
        end if
        write(*,*)'       ***** End system call: data-grad-1+1.sh -r run *****************'
        write(*,*)'       *****'
     end if
     !! fort.244 => grad=-1
     !! fort.245 => grad=+1
     !!system call-Up
     !!
     !! New integration over the FS with the
     !! correct k-points
     !! grad=+1-DOWN
     ap=0
     allocate (kpo(nump))
     do i=1,nump
        !          1 2 3 4    5     6  7  8  9
        read(245,*)x,y,z,peso,gradf,dx,dy,dz,dum
        ap=ap+peso*gradf*dx*dy*dz
        write(246,446)x,y,z,peso,gradf,dx,dy,dz,1
        !looking for kmin & kmax so new k-points are generated between these values
        kpo(i)=sqrt(x**2+y**2+z**2)
     end do
     !looking from kmax and kmin
     kmin=minval(kpo,dim=1)
     kmax=maxval(kpo,dim=1)
     write(247,446)kmin,kmax
     !! grad=+1-UP
     !! grad=-1-DOWN
     am=0
     allocate (kpi(numm))
     do i=1,numm
        !          1 2 3 4    5     6  7  8  9
        read(244,*)x,y,z,peso,gradf,dx,dy,dz,dum
        am=am+peso*gradf*dx*dy*dz
        write(246,446)x,y,z,peso,gradf,dx,dy,dz,-1
        !looking for kmin & kmax so new k-points are generated between these values
        kpi(i)=sqrt(x**2+y**2+z**2)
     end do
     !looking from kmax and kmin
     kmin=minval(kpi,dim=1)
     kmax=maxval(kpi,dim=1)
     write(248,446)kmin,kmax
     !! grad=-1-UP
     !!
     if ( 1 .eq. 1 )  then
        error=100*(aexact-circun)/aexact
        if ( 1 .eq. 2) then
           write(*,*)'       ***** # of ideal points*****'
           write(*,73)'    ',conta,'inner IBZ surface points around ',kc, ' +/- ',p,' %'
           write(*,*)contin,'  inner in-points'
           write(*,*)contout,'  inner out-points'
           write(*,73)'    ',contall,'all IBZ surface points around ',kc, ' +/- ',p,' %'
           write(*,*)contina,'    all in-points'
           write(*,*)contouta,'    all out-points'
        end if
        write(*,*)'       ***** all FS points k-points *****'
        if(kcBand.gt.0) then
           write(*,73)'    ',nump+numm,' all FermiSu points around ',kc, ' +/- ',p,' %'
           write(*,*)'    ',numm,'grad=-1',nump,'grad=+1'
        else
           !           write(*,73)'    ',nump+numm,' all FermiSu points around ',kc, ' +/- ',p,' %'
           write(*,*)'    ',numm,'FS',' and ',nump,'Just above FS'
        end if
        write(*,*)'       ***** Area (1/Bohr^2)*****'
        write(*,*)'       ***** inner k-points *********'
        write(*,72)'      IDArea= ',circun,' exact= ',aexact,' %error= ',error
        write(*,*)'       *************************'
        errorall=100*(aexact-circuna)/aexact
        write(*,*)'       ***** all k-points *********'
        write(*,72)'      IDArea= ',circuna,' exact= ',aexact,' %error= ',errorall
        if(kcBand.gt.0) then
           write(*,27)'        FS= grad+1',ap,'+ grad-1',am,'=',ap+am
        else
           write(*,27)'FS Area=',ap,' | aboveFS',am!,'=',ap+am
        end if
        write(*,*)'       **************'
        !           1 2  3       4    5    6      7       8    9    10    11     12    13
        write(35,53)N,Nk,contina,numm,nump,circun,circuna,ap  ,am  ,ap+am,aexact,error,errorall
        !                #idin   #-1  #+1  aIDin  aIDall  FS+1 FS-1 FS    ideal  in    all
     end if
     !
     !Integrates Data-Up
     !
     !
     if ( 1 .eq. 2 ) then
        !Plots for Grad B-Down
        !
        !
        ! ideal-Down
        !
        ! inner-kpoints-Down
        !
        ! 201-in:-1 -> 200-out:-1-Down
        d(:)=1000.
        do i=1,caim
           do j=1,caom
              d(j)=sqrt(sum((ikvoutfa(j)%kvoutfa(:)-ikvinfa(i)%kvinfa(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.202 $sym/out202=$case.ibz-fs-vectors-in-out\_$nk202
           write(206,314)ikvinfa(i)%kvinfa(1),ikvinfa(i)%kvinfa(2),ikvinfa(i)%kvinfa(3),ikvoutfa(mini)%kvoutfa(1),ikvoutfa(mini)%kvoutfa(2),ikvoutfa(mini)%kvoutfa(3)
        end do
        ! 201-in:-1 -> 200-out:-1-Up
        !
        ! 200-out:-1 -> 201-in:-1-Down
        d(:)=1000.
        do i=1,caom
           do j=1,caim
              d(j)=sqrt(sum((ikvinfa(j)%kvinfa(:)-ikvoutfa(i)%kvoutfa(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !aqui boy
           !mv fort.207 $sym/out203=$case.ibz-fs-vectors-out-in\_$nk207
           write(207,314)ikvoutfa(i)%kvoutfa(1),ikvoutfa(i)%kvoutfa(2),ikvoutfa(i)%kvoutfa(3),ikvinfa(mini)%kvinfa(1),ikvinfa(mini)%kvinfa(2),ikvinfa(mini)%kvinfa(3)
        end do
        ! 200-out:-1 -> 201-in:-1-Up
        !
        ! 204-in:+1 -> 204-out:+1-Down
        d(:)=1000.
        do i=1,caip
           do j=1,caop
              d(j)=sqrt(sum((ikvoutfap(j)%kvoutfap(:)-ikvinfap(i)%kvinfap(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.204 $sym/out204=$case.ibz-fs-vectors-in-out+1\_$nk202
           write(204,314)ikvinfap(i)%kvinfap(1),ikvinfap(i)%kvinfap(2),ikvinfap(i)%kvinfap(3),ikvoutfap(mini)%kvoutfap(1),ikvoutfap(mini)%kvoutfap(2),ikvoutfap(mini)%kvoutfap(3)
        end do
        ! 204-in:+1 -> 204-out:+1-Up
        !
        ! 204-out:+1 -> 205-in:+1-Down
        d(:)=1000.
        do i=1,caop
           do j=1,caip
              d(j)=sqrt(sum((ikvinfap(j)%kvinfap(:)-ikvoutfap(i)%kvoutfap(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.204 $sym/out204=$case.ibz-fs-vectors-in-out+1\_$nk202
           write(205,314)ikvoutfap(i)%kvoutfap(1),ikvoutfap(i)%kvoutfap(2),ikvoutfap(i)%kvoutfap(3),ikvinfap(mini)%kvinfap(1),ikvinfap(mini)%kvinfap(2),ikvinfap(mini)%kvinfap(3)
        end do
        ! 204-out:+1 -> 205-in:+1-Up
        ! in -> out-Up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! in -> out-Down
        d(:)=1000.
        do i=1,contin
           do j=1,contout
              d(j)=sqrt(sum((ikvout(j)%kvout(:)**2-ikvin(i)%kvin(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !aqui boy
           !mv fort.500 $sym/out500=$case.ibz-ideal-vectors-in-out\_$nk500
           write(500,314)ikvin(i)%kvin(1),ikvin(i)%kvin(2),ikvin(i)%kvin(3),ikvout(mini)%kvout(1),ikvout(mini)%kvout(2),ikvout(mini)%kvout(3)
        end do
        ! in -> out-Up
        ! out -> in-Down
        d(:)=1000.
        do i=1,contout
           do j=1,contin
              d(j)=sqrt(sum((ikvin(j)%kvin(:)-ikvout(i)%kvout(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.501 $sym/$out3=$case.ibz-vectors-out-in_$kn
           write(501,314)ikvout(i)%kvout(1),ikvout(i)%kvout(2),ikvout(i)%kvout(3),ikvin(mini)%kvin(1),ikvin(mini)%kvin(2),ikvin(mini)%kvin(3)
        end do
        ! out -> in-Up
        !
        ! inner-kpoints-Up
        !
        !
        ! all-kpoints-Down
        !
        ! in -> out-Down
        d(:)=1000.
        do i=1,contina
           do j=1,contouta
              d(j)=sqrt(sum((ikvouta(j)%kvouta(:)-ikvina(i)%kvina(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.400 $sym/out400=$case.ibz-ideala-vectors-in-out\_$nk400
           write(400,314)ikvina(i)%kvina(1),ikvina(i)%kvina(2),ikvina(i)%kvina(3),ikvouta(mini)%kvouta(1),ikvouta(mini)%kvouta(2),ikvouta(mini)%kvouta(3)
        end do
        ! in -> out-Up
        ! out -> in-Down
        d(:)=1000.
        do i=1,contouta
           do j=1,contina
              d(j)=sqrt(sum((ikvina(j)%kvina(:)-ikvouta(i)%kvouta(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.401 $sym/out401=$case.ibz-ideala-vectors-out-in\_$nk401
           write(401,314)ikvouta(i)%kvouta(1),ikvouta(i)%kvouta(2),ikvouta(i)%kvouta(3),ikvina(mini)%kvina(1),ikvina(mini)%kvina(2),ikvina(mini)%kvina(3)
        end do
        ! out -> in-Up
        !
        ! all-kpoints-Up
        !
        ! ideal-Up
        !
        ! Fermi-Down
        !
        ! inner-kpoints-Down
        !
        ! in -> out-Down
        d(:)=1000.
        do i=1,continf
           do j=1,contoutf
              d(j)=sqrt(sum((ikvoutf(j)%kvoutf(:)-ikvinf(i)%kvinf(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.502 $sym/$out3=$case.fs-vectors-in-out_$kn
           write(502,314)ikvinf(i)%kvinf(1),ikvinf(i)%kvinf(2),ikvinf(i)%kvinf(3),ikvoutf(mini)%kvoutf(1),ikvoutf(mini)%kvoutf(2),ikvoutf(mini)%kvoutf(3)
        end do
        ! in -> out-Up
        ! out -> in-Down
        d(:)=1000.
        do i=1,contoutf
           do j=1,continf
              d(j)=sqrt(sum((ikvinf(j)%kvinf(:)-ikvoutf(i)%kvoutf(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.503 $sym/$out3=$case.fs-vectors-out-in_$kn
           write(503,314)ikvoutf(i)%kvoutf(1),ikvoutf(i)%kvoutf(2),ikvoutf(i)%kvoutf(3),ikvinf(mini)%kvinf(1),ikvinf(mini)%kvinf(2),ikvinf(mini)%kvinf(3)
        end do
        ! out -> in-Up
        !
        ! inner-kpoints-Up
        !
        ! all-kpoints-Down
        !
        ! in -> out-Down
        d(:)=1000.
        do i=1,continfa
           do j=1,contoutfa
              d(j)=sqrt(sum((ikvoutfa(j)%kvoutfa(:)-ikvinfa(i)%kvinfa(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.402 $sym/$out402=$case.fsa-vectors-in-out_$kn
           write(402,314)ikvinfa(i)%kvinfa(1),ikvinfa(i)%kvinfa(2),ikvinfa(i)%kvinfa(3),ikvoutfa(mini)%kvoutfa(1),ikvoutfa(mini)%kvoutfa(2),ikvoutfa(mini)%kvoutfa(3)
        end do
        ! in -> out-Up
        ! out -> in-Down
        d(:)=1000.
        do i=1,contoutfa
           do j=1,continfa
              d(j)=sqrt(sum((ikvinfa(j)%kvinfa(:)-ikvoutfa(i)%kvoutfa(:))**2))
           end do
           min=minval(d,dim=1)
           mini=minloc(d,dim=1)
           !mv fort.403 $sym/$out3=$case.fsa-vectors-out-in_$kn
           write(403,314)ikvoutfa(i)%kvoutfa(1),ikvoutfa(i)%kvoutfa(2),ikvoutfa(i)%kvoutfa(3),ikvinfa(mini)%kvinfa(1),ikvinfa(mini)%kvinfa(2),ikvinfa(mini)%kvinfa(3)
        end do
        ! out -> in-Up
        !
        ! all-kpoints-Up
        !
        ! Fermi-Up
        !
        !Plots for Grad B-Up
        !
     end if !( 1 .eq. 2 )
  end if !( which .eq. 2 )
  !
  ! 3rd BLOCK-UP
  !
!!!!!!!!!!!!!!!!!!!!!!!!
2 format(3i6,3f12.6)
3 format(4f6.1,3i5,3f10.5,1i4,1f4.0)
4 format(3f11.5,1i5,3f5.0,1i4,1f4.0)
41 format(6f11.5)
8 format(5i5,3f6.0)
22 format(3f11.6,1f4.0)
29 format(3f11.6,1i4,2f6.1,3i6,1i4,1f10.4)
35 format(6i8,5f15.8,2f10.4)
53 format(5i8,6f15.8,2f10.4)
67 format(4f11.6,4f6.1)
60 format(1f5.0,7f11.6,3f4.0)
72 format(1A16,1f11.5,1A9  ,1f11.5,1A10   ,1f11.2)
27 format(1A18,1f11.5,1A10,1f11.5,1A3,1f11.5)
87 format(1A23,1f6.0,1A13,1f6.0,1A3,1f6.0)
83 format(1A23,2f11.5,1A13,1f11.5,1A3,1f11.5)
73 format(1A6,1i6,1A35,1f6.4,1A5,1f6.2,1A3)
76 format(3f13.8,1i5,2f4.0,3i6,3e13.5)
!67 format(3f13.8,1i5,2f4.0,3i6)
88 format(1f6.4)
96 format(3f11.6,1f5.0,1f15.6,3f12.6,2f5.0)
17 format(3f11.6,5f5.0,1f15.6)
199 format(8i4)
314 format(6f12.8)
446 format(3f11.6,1f7.0,1f15.6,3f11.6,i4)
803 format(3f12.8,1i4,1f6.0,1i4,1f4.0)

end program bccKpointsIntegration


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! chidas juan exec @band 36
! We can change N=85&tol=.005 to get a denser refined mesh around the FS
!
! STEP#0 Coarse-Grained IBZ k-points
!
!                                 (WE CAN PLAY WITH -k 71 as the starting point)
!> new-fermi-surface-via-kf.sh -k "71" -f .128 -w bcc -d 80 -b 36 -v -1.78 -g 0 -p 2
! generates
! eigen_6664_14-nospin & y2c3.klist_6664
! from where the Coarse-Grained FS for band 36 is obtained
! plotting the results (guide.pdf may be followed...we need to document this part with detail)
! we obtain the
!  #grad=-1
!  zm=.027
!  xm=.045
!  ym1=.15
!  ym2=.21
!  #grad=+1
!  xp=.06
!  yp1=.0475
!  yp2=.0725
!  zp=.06
! that restrict the k-points around the
! two FS of band 36
!
! STEP#1 Coarse-Grained FS
!
!The  6664-k obtained @STEP are used by
!                                 71 must be the sames as that used @STEP#0
!> fs-bz-bcc-grid-integrate.sh -N 71 -b 36 -c .128 -p 80 -e -1.78 -k 0.12
! generates:
!  symmetries/y2c3.fs-ibz-in-116-1-N1 (57) @ Grad=-1 -> fort.244
!  symmetries/y2c3.fs-ibz-in-116+1-N1 (59) @ Grad=+1 -> fort.245
! that give the coarse-grained FS
! the info is stored @ .coarse-grained-N
!
! STEP#2 Refinement
!
!band-36-fs-bz-bcc-grid-integrate.sh -N 95 -M 71 -b 36 -c .128 -p 80 -e -1.78 -k 0.12
!
!We use -M 71 same value used @STEPS #0 & #1
!The value used @ -N is used to increase the number of fine-grained k-points around the FS
!
! STEP#3
!
! run refined FS with the set of FS k-points obtained  @SETP#2
!
!FIRST AUTOMATISE ABOVE, THEN WE CONTINUE WITH THE

!3RD BLOCK WHERE WE INTEGRATE OVER THE FS.

!Mighty Fine Folks

!gnuplot> sp 'symmetries/y2c3.fs-ibz-in-116-1-57' u 1:2:3 w p pt 7 ps 1 lc rgb "black" t 'FSCG:grad=-1'
!gnuplot> rep 'symmetries/y2c3.fs-ibz-in-209-1-104' u 1:2:3 w p pt 7 ps 1 lc rgb "blue" t 'FS:grad=-1'
!gnuplot> rep 'symmetries/y2c3.fs-ibz-in-209+1-105' u 1:2:3 w p pt 7 ps 1 lc rgb "red" t 'FS:grad=+1'
!gnuplot> rep 'symmetries/y2c3.fs-ibz-in-116+1-59' u 1:2:3 w p pt 7 ps 1 lc rgb "violet" t 'FSCG:grad=+1'

!WARNING: NEW KPOINTS < ORIGINAL KPOINTS !!!!!!! WRONG !!!!!!!
!gnuplot> sp 'symmetries/y2c3.fs-ibz-in-97+1-51' u 1:2:3 w p pt 7 ps 1 lc rgb "red" t '97+1-51'       NEW KPOINTS
!gnuplot> rep 'symmetries/y2c3.fs-ibz-in-116+1-59' u 1:2:3 w p pt 7 ps 1 lc rgb "blue" t '116+1-59'   ORIGINAL KPOINTS
!gnuplot> rep 'symmetries/y2c3.fs-ibz-in-97-1-46' u 1:2:3 w p pt 7 ps 1 lc rgb "green" t '97-1-46'    NEW KPOINTS
!gnuplot> rep 'symmetries/y2c3.fs-ibz-in-116-1-57' u 1:2:3 w p pt 7 ps 1 lc rgb "black" t '116-1-57'  ORIGINAL KPOINTS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! chidas juan plots
!!$gnuplot> sp 'symmetries/y2c3.fs-ibz-in-116-1-57' u 1:2:3 w p pt 7 ps 1.2 lc rgb 'red' t 'FS:-1'
!!$gnuplot> rep 'fort.751' u 1:2:3 w p pt 7 ps 1 lc rgb 'brown' t 'New:-1'
!!$gnuplot> rep 'symmetries/y2c3.fs-ibz-in-116+1-59' u 1:2:3 w p pt 7 ps 1.2 lc rgb 'blue' t 'FS:+1'
!!$gnuplot> rep 'fort.750' u 1:2:3 w p pt 7 ps 1 lc rgb 'black' t 'New:+1'
!!$gnuplot> rep 'fort.22' u 1:2:3 w p pt 7 ps .6 lc rgb 'violet' t 'all:+1'
!!$gnuplot> rep 'fort.21' u 1:2:3 w p pt 7 ps .6 lc rgb 'green' t 'all:-1'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! doubious-down

!!$#gnuplot>  sp 'symmetries/y2c3.fs-ibz-in-116+1-59' u 1:2:3 w p pt 7 ps 1.2 lc rgb "red" t 'FS:grad=+1'
!!$#gnuplot> rep 'fort.22' u 1:2:3 w p pt 7 ps .8 t 'NEW grad=+1'
!!$#gnuplot> rep 'symmetries/y2c3.fs-ibz-in-116-1-57' u 1:2:3 w p pt 7 ps 1.2 lc rgb "blue" t 'FS:grad=-1'
!!$#gnuplot> rep 'fort.21' u 1:2:3 w p pt 7 ps .8 lc rgb "brown"   t 'NEW grad=-1'
!!$gnuplot> set xlabel 'x'
!!$gnuplot> set ylabel 'y'
!!$gnuplot> set zlabel 'z


!!$gnuplot> sp 'fort.30' u ($6==0?$1:1/0):2:3 w p pt 7 ps .6 lc rgb "red" t 'BF=0'
!!$gnuplot> rep 'fort.30' u ($6==1?$1:1/0):2:3 w p pt 7 ps .6 lc rgb "blue" t 'BF=1'
!!$gnuplot> rep 'y2c3.kflist_189' u 1:2:3 w p pt 7 ps 1

!!$gnuplot> sp 'fort.4' u ($12==1?$8:1/0):9:10 w p pt 7 ps .6 lc rgb "black" t 'BF=1-fort.4'
!!$gnuplot> rep 'fort.4' u ($12==0?$8:1/0):9:10 w p pt 7 ps .6 lc rgb "green" t 'BF=0-fort.4'
!!$gnuplot> rep 'y2c3.kflist_189' u 1:2:3 w p pt 7 ps 1


!!$**************** we have a winner-Down ***********************************************************
!!$data for inner FS k-points
!!$cat fort.201 fort.203 > 13.dat
!!$sort -n -k 1 13.dat > 13.dat-s
!!$sort -u 13.dat-s > 13.dat-su
!!$data for outer FS k-points
!!$cat fort.200 fort.202 > 02.dat
!!$sort -n -k 1 02.dat > 02.dat-s
!!$sort -u 02.dat-s > 02.dat-su
!!$plot to graphically find kc
!!$gnuplot>  k(x,y,z)=sqrt(x**2+y**2+z**2)
!!$gnuplot> p '13.dat-su' u :(k($1,$2,$3)) w p
!!how can we obtain kc automatically?
!!$gnuplot> kc=0.12
!!$gnuplot> rep kc
!!in FS-plots for grad=+/-1
!!$gnuplot> sp '13.dat-su' u (k($1,$2,$3)<kc?$1:1/0):2:3 w p pt 7 ps .8 lc rgb "red" t 'in:+1'
!!$gnuplot> rep '13.dat-su' u (k($1,$2,$3)>kc?$1:1/0):2:3 w p pt 7 ps .8 lc rgb "blue" t 'in:-1'
!!out FS-plots for grad=+/-1
!!$gnuplot> sp '02.dat-su' u (k($1,$2,$3)<kc?$1:1/0):2:3 w p pt 7 ps .8 lc rgb "green" t 'out:+1'
!!$gnuplot> rep '02.dat-su' u (k($1,$2,$3)>kc?$1:1/0):2:3 w p pt 7 ps .8 lc rgb "brown" t 'out:-1'
!!$**************** we have a winner-Up   ***********************************************************



!!$Dec-2-2021
!!$gnuplot> sp 'fort.30' u ($6==1?$1:1/0):2:3 w p pt 7 ps .6
!!$gnuplot> rep 'fort.30' u ($6==0?$1:1/0):2:3 w p pt 7 ps .6
!!$gnuplot> rep 'fort.201' u 1:2:3 w p pt 7 ps .8 lc rgb "red" t 'in:-1'
!!$gnuplot> rep 'fort.203' u 1:2:3 w p pt 7 ps .8 lc rgb "blue" t 'in:+1'

!!$gnuplot> rep 'fort.200' u 1:2:3 w p pt 1 ps .8 t 'out:-1'
!!$gnuplot> rep 'fort.200' u 1:2:3 w p pt 1 ps 1 t 'out:-1'
!!$gnuplot> rep 'fort.202' u 1:2:3 w p pt 7 ps 1 t 'out:+1'




!!$  PLOTS-D
!!$gnuplot> sp 'fort.201' u 6:7:8 w p pt 7 ps 1 t 'Fermi Surface'
!!$gnuplot> sp 'fort.201' u ($10==-1?$6:1/0):7:8 w p pt 7 ps 1 t 'FS grad=-1'
!!$gnuplot> rep 'fort.201' u ($10==1?$6:1/0):7:8 w p pt 7 ps 1 t 'FS grad=+1'
!!$gnuplot> rep 'fort.200' u ($10==-1?$6:1/0):7:8 w p pt 7 ps .8 t 'out grad=-1'
!!$gnuplot> rep 'fort.200' u ($10==1?$6:1/0):7:8 w p pt 7 ps .8 t 'out grad=+1'
!!$  PLOTS-U

!!$ looking for full IBZ
!!$gnuplot> sp 'fort.60' u ($10==-1&&$1!=48?$6:1/0):7:8  w p pt 7 ps .8  t '=> GNP  points @ grad=-1'
!!$gnuplot> rep 'fort.60' u ($10==1&&$1!=48?$6:1/0):7:8  w p pt 7 ps .8  t '=> GNP  points @ grad=+1'
!!$gnuplot> rep 'fort.60' u ($10==1&&$1==48?$6:1/0):7:8  w p pt 7 ps .8  t '=> w=48 points @ grad=+1'
!!$gnuplot> rep 'fort.60' u ($10==-1&&$1==48?$6:1/0):7:8  w p pt 7 ps .8 t '=> w=48 points @ grad=-1'

!!$gnuplot> sp 'fort.61' u ($10==-1&&$1!=48?$6:1/0):7:8  w p pt 7 ps .8 => GNP points    = fort.60
!!$gnuplot> rep 'fort.61' u ($10==-1&&$1==48?$6:1/0):7:8  w p pt 7 ps .8 => w=48 points  = fort.60
!!$gnuplot> rep 'fort.61' u ($10==1&&$1==48?$6:1/0):7:8  w p pt 7 ps .8 => NO POINTS
!!$gnuplot> rep 'fort.61' u ($10==1&&$1!=48?$6:1/0):7:8  w p pt 7 ps .8 => NO POINTS




!!$ attempts so far:
!!$gnuplot> sp 'fort.201' u ($10==-1?$6:1/0):7:8 w p pt 7 ps 1 t '201-in:-1'
!!$gnuplot> rep 'fort.200' u ($10==-1?$6:1/0):7:8 w p pt 7 ps 1 t '200-out:-1'
!!$gnuplot> rep 'fort.202' u 1:2:3:($4-$1):($5-$2):($6-$3) w vectors head filled
!!$gnuplot>  rep 'fort.203' u 1:2:3:($4-$1):($5-$2):($6-$3) w vectors head filled
!!$gnuplot>  rep 'fort.204' u 1:2:3:($4-$1):($5-$2):($6-$3) w vectors head filled
!!$gnuplot>  rep 'fort.205' u 1:2:3:($4-$1):($5-$2):($6-$3) w vectors head filled

!!$ below covers all k-points
!!$gnuplot> sp 'fort.4' u ($4==12?$8:1/0):9:10  w p pt 7 ps .8 t 'GPH-w=12'
!!$gnuplot> rep 'fort.4' u ($4==24?$8:1/0):9:10  w p pt 7 ps .8 t 'GPN&GNH-w=24'
!!$gnuplot> rep 'fort.4' u ($4==48?$8:1/0):9:10  w p pt 7 ps .6 t 'inner-w=48'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$Final files&plots
!!#BF(i,j,k)=1 in k-points => FERMI-SURFACE for Grad = +1 and -1. So far only band 37 has both
!!$gnuplot> rep 'fort.201' u ($10==-1?$6:1/0):7:8 w p pt 7 ps .8 t '200-in:-1'
!!$gnuplot> rep 'fort.201' u ($10==1?$6:1/0):7:8 w p pt 7 ps .8 t '200-in:+1'
!!$
!!#BF(i,j,k)=0 out k-points => above FERMI-SURFACE for Grad = +1 and -1. So far only band 37 has both
!!$gnuplot> rep 'fort.200' u ($10==-1?$6:1/0):7:8 w p pt 7 ps .8 t '200-in:-1'
!!$gnuplot> rep 'fort.200' u ($10==1?$6:1/0):7:8 w p pt 7 ps .8 t '200-in:+1'


!!$supre john-chidas-Down
!!$FS
!!$gnuplot> sp 'fort.60' u ($9==1&&$10==-1?$6:1/0):7:8 w p pt 7 ps .8 t '60-in:-1'
!!$gnuplot> rep 'fort.60' u ($9==1&&$10==+1?$6:1/0):7:8 w p pt 7 ps .8 t '60-in:+1'
!!$gnuplot> rep 'fort.60' u ($9==0&&$10==-1?$6:1/0):7:8 w p pt 7 ps .8 t '60-out:-1'
!!$gnuplot> rep 'fort.60' u ($9==0&&$10==1?$6:1/0):7:8 w p pt 7 ps .8 t '60-out:+1'
!!$supre john-chidas-Up


!!$john-chidas-Down
!!$FS
!!$gnuplot> sp 'fort.30' u ($6==0&&$4==1?$1:1/0):2:3  w p pt 7 ps .8 t 'fs-in'
!!$gnuplot> rep 'fort.30' u ($6==1&&$4==1?$1:1/0):2:3  w p pt 7 ps .8 t 'fs-out'
!!$gnuplot> rep 'fort.30' u ($6==1&&$4==3?$1:1/0):2:3  w p pt 7 ps 1 t 'fs-out'
!!$gnuplot> rep 'fort.30' u ($6==1&&$4==5?$1:1/0):2:3  w p pt 7 ps 1 t 'fs-out'
!!$gnuplot> rep 'fort.30' u ($6==1&&$4==6?$1:1/0):2:3  w p pt 7 ps 1 t 'fs-out'
!!$gnuplot> rep 'fort.30' u ($6==0&&$4==3?$1:1/0):2:3  w p pt 7 ps .8 t 'fs-out'
!!$gnuplot> rep 'fort.30' u ($6==0&&$4==5?$1:1/0):2:3  w p pt 7 ps .8 t 'fs-out'
!!$gnuplot> rep 'fort.30' u ($6==0&&$4==6?$1:1/0):2:3  w p pt 7 ps .8 t 'fs-out'
!!$Ideal
!!$gnuplot> sp 'fort.29' u ($6==0&&$4==1?$1:1/0):2:3  w p pt 7 ps .8 t 'ideal-in'
!!$gnuplot> rep 'fort.29' u ($6==1&&$4==1?$1:1/0):2:3  w p pt 7 ps .8 t 'ideal-out'
!!$gnuplot> rep 'fort.29' u ($6==1&&$4==3?$1:1/0):2:3  w p pt 7 ps 1 t 'ideal-out'
!!$gnuplot> rep 'fort.29' u ($6==1&&$4==5?$1:1/0):2:3  w p pt 7 ps 1 t 'ideal-out'
!!$gnuplot> rep 'fort.29' u ($6==1&&$4==6?$1:1/0):2:3  w p pt 7 ps 1 t 'ideal-out'
!!$gnuplot> rep 'fort.29' u ($6==0&&$4==3?$1:1/0):2:3  w p pt 7 ps .8 t 'ideal-out'
!!$gnuplot> rep 'fort.29' u ($6==0&&$4==5?$1:1/0):2:3  w p pt 7 ps .8 t 'ideal-out'
!!$gnuplot> rep 'fort.29' u ($6==0&&$4==6?$1:1/0):2:3  w p pt 7 ps .8 t 'ideal-out'
!!$Plus
!!$gnuplot> p 'fort.30' u ($6==1?$11:1/0) w p pt 7 ps .7 t 'B=1'
!!$gnuplot> rep 'fort.30' u ($6==0?$11:1/0) w p pt 7 ps .7 t 'B=0'
!!$gnuplot> rep efermi w l lw 4
!!$john-chidas-Up



!gnuplot> pwd
!/home/bms/tiniba/ver6.0/utils/brillouin-zone/ibz/versions-bcc/y2c3
!gnuplot> load 'symmetries/ibz-bcc.g'
!gnuplot> rep 'fort.21' u 1:2:5 w p pt 7 ps .8
!gnuplot> rep '/home/bms/tiniba/ver6.0/utils/grad-method/res/data-753.d' u ($1-$2<=0?$1:1/0):2:3 w p pt 7 ps .9

!!$              !!!!!!!!!!!!!!!!!
!!$              if ( 1 .eq. 2) then
!!$                 if ( ( (rx(i) - ry(j) ) .le. 0 ) .and. ( ( rz(k)-rx(i) ) .le. 0 ) .and. ( ( rx(i)+ry(j) ) .le. g ) ) then
!!$                    km=sqrt(rx(i)**2+ry(j)**2+rz(k)**2)
!!$                    ! This selects a small region of k-values
!!$                    if ((km.ge.kmin).and.(km.le.kmax)) then
!!$                       write(200,199)i,j,k,ind(i),jnd(j),knd(k),yes(i,j,k),yes(ind(i),jnd(j),knd(k))
!!$                       if ( (ind(i).ne.-20) .and. (jnd(j).ne.-20) .and. (knd(k).ne.-20) ) then
!!$                          write(201,199)i,j,k,ind(i),jnd(j),knd(k)
!!$                          if (yes(i,j,k).ne.10) write(202,199)i,j,k,yes(i,j,k)
!!$                       end if
!!$                    end if
!!$                    ! This makes sure that (i,j,k) are equivalent to the original (i,j,k) chosen above
!!$                    ! HOWEVER WE HAVE LESS POINTS, SO SMETHING IS UTTERLY WRONG !!!!!!!!!!!!!!!!!
!!$                    if ( yes(ind(i),jnd(j),knd(k)) .eq. 1 ) then
!!$                       indc=indc+1
!!$                       write(199,199)i,j,k,ind(i),jnd(j),knd(k)
!!$                    end if
!!$                 end if
!!$              end if
!!$              !!!!!!!!!!!!!!!!!


!!$gnuplot> sp 'fort.803' u ($4!=48?$1:1/0):2:3  w p pt 7 ps .6
!!$gnuplot> rep 'fort.3' u ($4==24&&($1==-1||$2==-1)?$8:1/0):9:10  w p pt 7 ps 1
!!$gnuplot> rep 'fort.3' u ($4==12&&($2==-1||$3==-1)?$8:1/0):9:10  w p pt 7 ps 1

!!$gnuplot> sp 'fort.803' u ($4!=48?$1:1/0):2:3  w p pt 7 ps .6
!!$gnuplot> rep 'fort.3' u ($4==24&&($2==-1||$3==-1)?$8:1/0):9:10  w p pt 7 ps 1
!!$gnuplot> rep 'fort.3' u ($4==12&&($2==-1||$3==-1)?$8:1/0):9:10  w p pt 7 ps 1
!!$gnuplot> rep 'fort.3' u ($4==12&&($1==-1||$2==-1||$3==-1)?$8:1/0):9:10  w p pt 7 ps 1
!!$gnuplot> rep 'fort.3' u ($4==24&&($1==-1||$2==-1||$3==-1)?$8:1/0):9:10  w p pt 7 ps 1
!!$gnuplot> rep 'fort.3' u ($4==48&&($1==-1||$2==-1||$3==-1)?$8:1/0):9:10  w p pt 7 ps 1
!!$gnuplot> rep 'fort.3' u ($4==48&&($1==-1||$2==-1||$3==-1)?$8:1/0):9:10  w p pt 7 ps 1

!/Users/bms/chamba/research/spin-current/2d-monochalcogenides/new-results/sns/bz
