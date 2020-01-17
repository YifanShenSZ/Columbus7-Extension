!Generate a linear path from geometry 1 to geometry 2
!Number of steps along the path is controled by NStep
!Input : internal coordinate definition, geom1, geom2, (optional) NStep
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
    real*8,allocatable,dimension(:)::ElementNumber,mass,r1,r2,q1,q2
!Work variable
    character*32::chartemp; logical::flag; integer::i,j
    real*8::dbletemp,dbletemp1
    real*8,allocatable,dimension(:)::q,dq,r,rsave
!Initialize
    inquire(file='InternalCoordinateDefinition',exist=flag)
    if(flag) then
        chartemp=''
        call DefineInternalCoordinate(chartemp,intdim)
    else
        chartemp='Columbus7'
        call DefineInternalCoordinate(chartemp,intdim)
    end if
    open(unit=99,file='geom1',status='old')!Read geometry 1
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r1(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,*)ElementSymbol(i),ElementNumber(i),r1(3*i-2:3*i),mass(i)
        end do
        mass=mass*AMUInAU
    close(99)
    open(unit=99,file='geom2',status='old')!Read geometry 2
        allocate(r2(3*NAtoms))
        do i=1,NAtoms
            read(99,*)chartemp,dbletemp,r2(3*i-2:3*i),dbletemp1
        end do
    close(99)
    inquire(file='NStep',exist=flag)
    if(flag) then
        open(unit=99,file='NStep',status='old'); read(99,*)NStep; close(99)
    end if
    cartdim=3*NAtoms
    allocate(q1(intdim)); q1=InternalCoordinateq(r1,intdim,cartdim)
    allocate(q2(intdim)); q2=InternalCoordinateq(r2,intdim,cartdim)
!Generate a linear path from geometry 1 to geometry 2
    chartemp='assimilate'
    allocate(q(intdim))
    allocate(dq(intdim)); dq=(q2-q1)/dble(NStep+1)
    allocate(r(intdim))
    allocate(rsave(cartdim))
    open(unit=99,file='geom.all',status='replace')
            q=q1+dq
            r=CartesianCoordinater(q,cartdim,intdim,uniquify=chartemp,mass=mass,r0=r1)
            do j=1,NAtoms
                write(99,'(1x,A2,2x,F5.1,4F14.8)')ElementSymbol(j),ElementNumber(j),r(3*j-2:3*j),mass(j)/AMUInAU
            end do
            rsave=r
        do i=2,NStep
            q=q1+dq*dble(i)
            r=CartesianCoordinater(q,cartdim,intdim,uniquify=chartemp,mass=mass,r0=rsave)
            do j=1,NAtoms
                write(99,'(1x,A2,2x,F5.1,4F14.8)')ElementSymbol(j),ElementNumber(j),r(3*j-2:3*j),mass(j)/AMUInAU
            end do
            rsave=r
        end do
    close(99)
end program main