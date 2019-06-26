!Generate loop around the geometry calculated Hessian before along each normal mode
!For real frequency modes, will only generate 1 point with approximately 1d-5 Hatree energy higher
!For imaginary frequency modes, will generate 10 points with approximately 1d-5 to 1d-3 Hatree
program main
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization; use Chemistry
    use GeometryTransformation
    implicit none
!Molecule information
    integer::NAtoms,intdim,cartdim
    integer,allocatable,dimension(:)::ElementNumber
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::mass,r0,q0
!Normal mode
    integer::NImagMode
    real*8,allocatable,dimension(:)::freq,Normalq0
    real*8,allocatable,dimension(:,:)::mode,L
!Work variable
    character*32::chartemp
    real*8::dbletemp
    integer::i,j,k
    real*8,allocatable,dimension(:)::q,r
!Initialize
    call BetterRandomSeed()
    open(unit=99,file='geom',status='old')!Read molecule detail
        write(*,*)'geom format should be A2,I8,3F14.8,F14.8:'
        write(*,*)'    2 space for element symbol'
        write(*,*)'    8 space for integer element number'
        write(*,*)'    3 * 14 space with 8 decimals for Cartesian coordinate in atomic unit'
        write(*,*)'    14 space with 8 decimals for atomic mass in atomic mass unit'
        write(*,*)'Some Columbus version may use different format, please correct it manually'
        NAtoms=0
        do
            read(99,*,iostat=i); if(i/=0) exit
            NAtoms=NAtoms+1
        end do
        rewind 99
        allocate(ElementSymbol(NAtoms))
        allocate(ElementNumber(NAtoms))
        allocate(r0(3*NAtoms))
        allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,'(A2,I8,3F14.8,F14.8)')ElementSymbol(i),ElementNumber(i),r0(3*i-2:3*i),mass(i)
            ElementSymbol(i)=trim(adjustl(ElementSymbol(i)))
        end do
        mass=mass*AMUInAU!Convert to atomic unit
        call StandardizeGeometry(r0,mass,NAtoms,1)
    close(99)
    cartdim=3*NAtoms
    chartemp='Columbus7'; call DefineInternalCoordinate(chartemp,intdim)
    allocate(q0(intdim)); q0=InternalCoordinateq(r0,intdim,cartdim)
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
    allocate(Normalq0(intdim)); Normalq0=matmul(mode,q0)
!Generate loop geometries
    allocate(q(intdim)); allocate(r(cartdim))
    open(unit=99,file='geom.imag',status='replace')
        do i=1,NImagMode
            dbletemp=dSqrt(1d-5*2d0/(freq(i)*freq(i)))
            do k=-10,-1
                q=Normalq0; q(i)=q(i)+dbletemp*dble(k); q=matmul(L,q); call WriteGeom()
            end do
            do k=1,10
                q=Normalq0; q(i)=q(i)+dbletemp*dble(k); q=matmul(L,q); call WriteGeom()
            end do
        end do
    close(99)
    open(unit=99,file='geom.real',status='replace')
        do i=NImagMode+1,intdim
            dbletemp=dSqrt(1d-5*2d0/(freq(i)*freq(i)))
            q=Normalq0; q(i)=q(i)+dbletemp; q=matmul(L,q); call WriteGeom()
            q=Normalq0; q(i)=q(i)-dbletemp; q=matmul(L,q); call WriteGeom()
        end do
    close(99)
    contains
    subroutine WriteGeom()
        r=CartesianCoordinater(q,cartdim,intdim,mass=mass,r0=r0)
        do j=1,NAtoms
            write(99,'(A2,I8,3F14.8,F14.8)')ElementSymbol(j),ElementNumber(j),r(3*j-2:3*j),mass(j)/AMUInAU
        end do
    end subroutine WriteGeom
end program main