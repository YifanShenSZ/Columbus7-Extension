!Generate a linear path from geometry 1 to geometry 2
!Number of steps along the path is controled by NStep
!Input : internal coordinate definition, geom1, geom2
!Output: geom.all
!Dependency: my Fortran-Library, as written in makefile
program main
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization; use Chemistry
    use GeometryTransformation
    implicit none
!Number of steps along the path
    integer::NStep=10
!Molecule information
    integer::NAtoms,intdim,cartdim
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::ElementNumber,mass,r0,q0,r1,q1
!Work variable
    character*32::chartemp; logical::flag; integer::i,j
    real*8::dbletemp,dbletemp1
    real*8,allocatable,dimension(:)::q,r,dq
!Initialize
    open(unit=99,file='geom1',status='old')!Read geometry 1
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r0(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,*)ElementSymbol(i),ElementNumber(i),r0(3*i-2:3*i),mass(i)
        end do
        mass=mass*AMUInAU; call StandardizeGeometry(r0,mass,NAtoms,1)
    close(99)
    cartdim=3*NAtoms
    inquire(file='InternalCoordinateDefinition',exist=flag)
    if(flag) then
        chartemp=''
        call DefineInternalCoordinate(chartemp,intdim)
    else
        chartemp='Columbus7'
        call DefineInternalCoordinate(chartemp,intdim)
    end if
    allocate(q0(intdim)); q0=InternalCoordinateq(r0,intdim,cartdim)
    open(unit=99,file='geom2',status='old')!Read geometry 2
        allocate(r1(3*NAtoms))
        do i=1,NAtoms
            read(99,*)chartemp,dbletemp,r1(3*i-2:3*i),dbletemp1
        end do
        call StandardizeGeometry(r1,mass,NAtoms,1)
    close(99)
    allocate(q1(intdim)); q1=InternalCoordinateq(r1,intdim,cartdim)
!Generate a linear path from geometry 1 to geometry 2
    allocate(q(intdim)); allocate(r(cartdim))
    allocate(dq(intdim)); dq=(q1-q0)/dble(NStep+1)
    open(unit=99,file='geom.all',status='replace')
        do i=1,NStep
            q=q0+dq*dble(i)
            call WriteGeom()
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