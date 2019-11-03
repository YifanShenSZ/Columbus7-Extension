!Standard input files: intcfl, geom (reference geometry)
!Standard output file: geom.all
!May take more io according to user modification
program main
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization; use Chemistry
    use GeometryTransformation
    implicit none
!Molecule information
    integer::NAtoms,intdim,cartdim
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::ElementNumber,mass,r0,q0
!Work variable
    character*32::chartemp; integer::i,j
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
    allocate(q0(intdim)); q0=InternalCoordinateq(r0,intdim,cartdim)!Reference internal geometry
    allocate(q(intdim)); allocate(r(cartdim))
open(unit=99,file='geom.all',status='replace')!Modification starts here
    do i=-5,5!Example: start from q0, displace 1st internal coordinate
        q=q0
        q(1)=q(1)+dble(i)*0.04d0
        r=CartesianCoordinater(q,cartdim,intdim,mass=mass,r0=r0)
        do j=1,NAtoms
            write(99,'(1x,A2,2x,F5.1,4F14.8)')ElementSymbol(j),ElementNumber(j),r(3*j-2:3*j),mass(j)/AMUInAU
        end do
    end do
close(99)
end program main