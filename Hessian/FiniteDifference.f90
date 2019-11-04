!Compute n-th state energy Hessian by finite difference of analytical gradient
!Usage : ./FiniteDifference.exe n
!Input : intcfl, geom (the geometry to calculate Hessian)
!        geom.all, cartgrd.drt1.state+str(n)+.all (loop information)
!Output: hessian (for Columbus7)
!        freq.out, L.out (for NormalModeLoop.exe)
!        VibrationalFrequency.txt, geom.log (for user)
program main
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization; use Chemistry
    use GeometryTransformation
    implicit none
!Command line input
    character*32::cmd
!Molecule information
    integer::NAtoms,intdim,cartdim
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::ElementNumber,mass,r0,q0
!Loop information
    real*8,allocatable,dimension(:,:)::r,q,cartgrad,intgrad
!Vibration analysis
    real*8,allocatable,dimension(:)::freq
    real*8,allocatable,dimension(:,:)::B,L,Linv,cartmode
!Work variable
    character*32::chartemp; integer::i,j,k; real*8::dbletemp,dbletemp1
    real*8,allocatable,dimension(:,:)::Hessian
    call getarg(1,cmd)!Command line input
!Initialize
    open(unit=99,file='geom',status='old')!Read molecule detail
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r0(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,*)ElementSymbol(i),ElementNumber(i),r0(3*i-2:3*i),mass(i)
            ElementSymbol(i)=trim(adjustl(ElementSymbol(i)))
        end do
        mass=mass*AMUInAU; call StandardizeGeometry(r0,mass,NAtoms,1)
    close(99)
    cartdim=3*NAtoms
    chartemp='Columbus7'; call DefineInternalCoordinate(chartemp,intdim)
    allocate(B(intdim,cartdim)); allocate(q0(intdim))
    call WilsonBMatrixAndInternalCoordinateq(B,q0,r0,intdim,cartdim)
!Read loop information
    allocate(r(cartdim,2*intdim))
    open(unit=99,file='geom.all',status='old')
        do i=1,2*intdim
            do j=1,NAtoms
                read(99,*)chartemp,dbletemp,r(3*j-2:3*j,i),dbletemp1
            end do
        end do
    close(99)
    allocate(cartgrad(cartdim,2*intdim))
    open(unit=99,file='cartgrd.drt1.state'//trim(adjustl(cmd))//'.all',status='old')
        do i=1,2*intdim
            do j=1,NAtoms
                read(99,*)cartgrad(3*j-2:3*j,i)
            end do
        end do
    close(99)
    allocate(q(intdim,2*intdim)); allocate(intgrad(intdim,2*intdim))
    do i=1,2*intdim
        call Cartesian2Internal(r(:,i),cartdim,q(:,i),intdim,1,cartgrad=cartgrad(:,i),intgrad=intgrad(:,i))
    end do
!Calculate Hessian
    allocate(Hessian(intdim,intdim))
    forall(i=1:intdim)
        Hessian(:,i)=(intgrad(:,2*i-1)-intgrad(:,2*i))/(q(i,2*i-1)-q(i,2*i))
    end forall
    open(unit=99,file='hessian',status='replace')!Output for Columbus7
        do i=1,intdim
            do j=1,intdim,8
                do k=j,min(j+7,intdim)
                    write(99,'(F13.6)',advance='no')Hessian(k,i)
                end do
                write(99,*)
            end do
        end do
        write(99,*)Hessian
    close(99)
!Vibration analysis
    !Run Wilson GF method for vibrational frequency and internal normal mode, output vibrational frequency
    allocate(freq(intdim)); allocate(L(intdim,intdim)); allocate(Linv(intdim,intdim))
    call WilsonGFMethod(freq,L,Linv,Hessian,intdim,B,mass,NAtoms)
    open(unit=99,file='freq.out',status='replace')!Output for NormalModeLoop.exe
        write(99,*)freq
    close(99)
    open(unit=99,file='L.out',status='replace')
        write(99,*)L
    close(99)
    open(unit=99,file='VibrationalFrequency.txt',status='replace')!Output for user
        write(99,'(A4,A1,A15)')'Mode',char(9),'Frequency/cm^-1'
        do i=1,intdim; write(99,'(I4,A1,F14.8)')i,char(9),freq(i)/cm_1InAu; end do
    close(99)
    !Convert internal normal mode to Cartesian, output visualization
    allocate(cartmode(cartdim,intdim))
    call InternalMode2CartesianMode(freq,L,InternalDimension,B,cartmode,CartesianDimension)
    chartemp='geom.log'
    call Avogadro_Vibration(NAtoms,ElementSymbol,r0/AInAU,intdim,freq/cm_1InAu,cartmode,FileName=chartemp)
    write(*,'(1x,A77)')'To visualize the molecular structure and vibration, open geom.log in Avogadro'
end program main