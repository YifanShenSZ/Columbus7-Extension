!Check if the real mode loop has gradient component along imaginary mode
!Input: intcfl, geom (the geometry calculated Hessian)
!       VibrationalFrequency.txt, NormalMode.txt (normal mode information)
!       geom.all, cartgrd.drt1.state1.all (real frequency mode loop information)
program main
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization; use Chemistry
    use GeometryTransformation
    implicit none
!Molecule information
    integer::NAtoms,intdim,cartdim
    integer,allocatable,dimension(:)::ElementNumber
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::mass,r0
!Normal mode
    integer::NImagMode
    real*8,allocatable,dimension(:)::freq
    real*8,allocatable,dimension(:,:)::mode,L
!Loop information
    real*8,allocatable,dimension(:,:)::r,cartgrad,intgrad
!Work variable
    character*32::chartemp; real*8::dbletemp; integer::i,j
    real*8,allocatable,dimension(:)::q
!Initialize
    call BetterRandomSeed()
    open(unit=99,file='geom',status='old')!Read molecule detail
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r0(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,*)ElementSymbol(i),dbletemp,r0(3*i-2:3*i),mass(i)
            ElementSymbol(i)=trim(adjustl(ElementSymbol(i)))
            ElementNumber(i)=int(dbletemp)
        end do
        mass=mass*AMUInAU; call StandardizeGeometry(r0,mass,NAtoms,1)
    close(99)
    cartdim=3*NAtoms
    chartemp='Columbus7'; call DefineInternalCoordinate(chartemp,intdim)
!Read vibrational frequency and normal mode
    allocate(freq(intdim))
    open(unit=99,file='VibrationalFrequency.txt',status='old')
        read(99,*)
        do i=1,intdim
            read(99,*)j,freq(i)
        end do
        freq=freq*cm_1InAu
    close(99)
    if(freq(1)<0d0) then!Exist imaginary frequency mode, identify the number of such modes
        do i=1,intdim
            if(freq(i)>0d0) exit
        end do
        NImagMode=i-1
    else
        NImagMode=0
    end if
    allocate(mode(intdim,intdim))
	open(unit=99,file='NormalMode.txt',status='old')
        read(99,*)
		do i=1,intdim
			read(99,'(A9)',advance='no')chartemp
			do j=1,intdim-1
				read(99,'(F18.8,A1)',advance='no')mode(j,i),chartemp
			end do
			read(99,'(F18.8)')mode(intdim,i)
		end do
    close(99)
    allocate(L(intdim,intdim)); L=mode; call My_dgetri(L,intdim)
!Only check real frequency mode loop
    allocate(r(cartdim,2*(intdim-NImagMode)))
    open(unit=99,file='geom.all',status='old')
        do i=1,2*(intdim-NImagMode)
            do j=1,NAtoms
                read(99,*)chartemp,r(3*j-2:3*j,i),dbletemp
            end do
        end do
    close(99)
    allocate(cartgrad(cartdim,2*(intdim-NImagMode)))
    open(unit=99,file='cartgrd.drt1.state1.all',status='old')
        do i=1,2*(intdim-NImagMode)
            do j=1,NAtoms
                read(99,*)cartgrad(3*j-2:3*j,i)
            end do
        end do
    close(99)
    allocate(q(intdim)); allocate(intgrad(intdim,2*(intdim-NImagMode)))
    do i=1,2*(intdim-NImagMode)
        call Cartesian2Internal(r(:,i),cartdim,q,intdim,1,cartgrad=cartgrad(:,i),intgrad=intgrad(:,i))
        intgrad(:,i)=matmul(L,intgrad(:,i))
        do j=1,NImagMode
            if(dAbs(intgrad(j,i))>1d-6) then
                write(*,*)'Point',i,', displacing along real mode',NImagMode+(i+1)/2,', has gradient component along imaginary mode',j
            end if
        end do
    end do
end program main