!Generate loop around the geometry calculated Hessian before along each normal mode
!For imaginary frequency modes, will generate 10 points with approximately 1d-5 to 1d-3 Hatree
!For real frequency modes, will only generate 1 point with approximately 1d-5 Hatree
!Input : intcfl, geom (the geometry calculated Hessian)
!        freq.out, L.out (normal mode information)
!Output: geom.imag, geom.real (loop geometries along imaginary and real frequency modes)
program main
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization; use Chemistry
    use GeometryTransformation
    implicit none
!Molecule information
    integer::NAtoms,intdim,cartdim
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::ElementNumber,mass,r0,q0
!Normal mode
    integer::NImagMode
    real*8,allocatable,dimension(:)::freq
    real*8,allocatable,dimension(:,:)::L
!Work variable
    character*32::chartemp; real*8::dbletemp; integer::i,j,k
    real*8,allocatable,dimension(:)::q,r
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
    allocate(q0(intdim)); q0=InternalCoordinateq(r0,intdim,cartdim)
!Read vibrational frequency and normal mode
    allocate(freq(intdim))
    open(unit=99,file='freq.out',status='old')
        read(99,*)freq
    close(99)
    if(freq(1)<0d0) then!Exist imaginary frequency mode, identify the number of such modes
        do i=1,intdim; if(freq(i)>0d0) exit; end do
        NImagMode=i-1
    else
        NImagMode=0
    end if
    allocate(L(intdim,intdim))
	open(unit=99,file='L.out',status='old')
        read(99,*)L
    close(99)
!Generate loop geometries
    allocate(q(intdim)); allocate(r(cartdim))
    open(unit=99,file='geom.imag',status='replace')
        do i=1,NImagMode
            dbletemp=dSqrt(1d-5*2d0/(freq(i)*freq(i)))
            do k=-10,-1
                q=q0+dbletemp*dble(k)*L(:,i); call WriteGeom()
            end do
            do k=1,10
                q=q0+dbletemp*dble(k)*L(:,i); call WriteGeom()
            end do
        end do
    close(99)
    open(unit=99,file='geom.real',status='replace')
        do i=NImagMode+1,intdim
            dbletemp=dSqrt(1d-5*2d0/(freq(i)*freq(i)))
            q=q0+dbletemp*L(:,i); call WriteGeom()
            q=q0-dbletemp*L(:,i); call WriteGeom()
        end do
    close(99)
    contains
    subroutine WriteGeom()
        r=CartesianCoordinater(q,cartdim,intdim,mass=mass,r0=r0)
        do j=1,NAtoms
            write(99,'(1x,A2,2x,F5.1,4F14.8)')ElementSymbol(j),ElementNumber(j),r(3*j-2:3*j),mass(j)/AMUInAU
        end do
    end subroutine WriteGeom
end program main