!Compute ground state Hessian by finite difference of analytical gradient
!Input : intcfl, geom (the geometry to calculate Hessian)
!        geom.all, cartgrd.drt1.state1.all (loop information)
!Output: Hessian, VibrationalFrequency.txt, NormalMode.txt
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
!Loop information
    real*8,allocatable,dimension(:,:)::r,cartgrad,intgrad
!Vibration analysis
    real*8,allocatable,dimension(:)::freq
    real*8,allocatable,dimension(:,:)::B,mode,L
!Work variable
    character*32::chartemp; real*8::dbletemp; integer::i,j
    real*8,allocatable,dimension(:)::q
    real*8,allocatable,dimension(:,:)::Hessian
!Initialize
    call BetterRandomSeed()
    open(unit=99,file='geom',status='old')!Read molecule detail
        write(*,*)'geom format should be A2,I8,3F14.8,F14.8:'
        write(*,*)'    2 space for element symbol'
        write(*,*)'    8 space for integer element number'
        write(*,*)'    3 * 14 space with 8 decimals for Cartesian coordinate in atomic unit'
        write(*,*)'    14 space with 8 decimals for atomic mass in atomic mass unit'
        write(*,*)'Some Columbus version may use different format, please correct it manually'
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r0(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,'(A2,I8,3F14.8,F14.8)')ElementSymbol(i),ElementNumber(i),r0(3*i-2:3*i),mass(i)
            ElementSymbol(i)=trim(adjustl(ElementSymbol(i)))
        end do
        mass=mass*AMUInAU; call StandardizeGeometry(r0,mass,NAtoms,1)
    close(99)
    cartdim=3*NAtoms
    chartemp='Columbus7'; call DefineInternalCoordinate(chartemp,intdim)
!Read loop information
    allocate(r(cartdim,2*intdim))
    open(unit=99,file='geom.all',status='old')
        do i=1,2*intdim
            do j=1,NAtoms
                read(99,'(A10,3F14.8,F14.8)')chartemp,r(3*j-2:3*j,i),dbletemp
            end do
        end do
    close(99)
    allocate(cartgrad(cartdim,2*intdim))
    open(unit=99,file='cartgrd.drt1.state1.all',status='old')
        do i=1,2*intdim
            do j=1,NAtoms
                read(99,*)cartgrad(3*j-2:3*j,i)
            end do
        end do
    close(99)
    allocate(q(intdim)); allocate(intgrad(intdim,2*intdim))
    do i=1,2*intdim
        call Cartesian2Internal(r(:,i),cartdim,q,intdim,1,cartgrad=cartgrad(:,i),intgrad=intgrad(:,i))
    end do
!Calculate Hessian
    allocate(Hessian(intdim,intdim))
    do i=1,intdim
        Hessian(:,i)=intgrad(:,2*i-1)-intgrad(:,2*i)
        if(GeometryTransformation_IntCDef(i).motion(1).type=='stretching') then!Bond length: 0.01A
            Hessian(:,i)=Hessian(:,i)/(0.02d0*AInAU)
        else!Angle: 0.01
            Hessian(:,i)=Hessian(:,i)/0.02d0
        end if
    end do
    open(unit=99,file='Hessian',status='replace')
        write(99,*)Hessian
    close(99)
!Vibration analysis
    allocate(B(intdim,cartdim))
    call WilsonBMatrixAndInternalCoordinateq(B,q,r0,intdim,cartdim)
    allocate(freq(intdim)); allocate(mode(intdim,intdim)); allocate(L(intdim,intdim))
    call WilsonGFMethod(freq,mode,L,Hessian,intdim,B,mass,NAtoms)
    open(unit=99,file='VibrationalFrequency.txt',status='replace')
        write(99,'(A4,A1,A15)')'Mode',char(9),'Frequency/cm^-1'
        do i=1,intdim
            write(99,'(I4,A1,F14.8)')i,char(9),freq(i)/cm_1InAu
        end do
    close(99)
	open(unit=99,file='NormalMode.txt',status='replace')
	    write(99,'(A6,A1)',advance='no')'q\Mode',char(9)
		do i=1,intdim-1
			write(99,'(I6,A1)',advance='no')i,char(9)
		end do
		write(99,'(I6)')intdim
		do i=1,intdim
			write(99,'(I8,A1)',advance='no')i,char(9)
			do j=1,intdim-1
				write(99,'(F18.8,A1)',advance='no')mode(j,i),char(9)
			end do
			write(99,'(F18.8)')mode(intdim,i)
		end do
    close(99)
end program main