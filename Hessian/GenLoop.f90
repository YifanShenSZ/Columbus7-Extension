!Prepare finite difference geometries for computing ground state Hessian
!Input files: intcfl, geom (the geometry to calculate Hessian)
!Output file: geom.all (finite difference geometries)
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
!Work variable
    character*32::chartemp
    integer::i,j
    real*8,allocatable,dimension(:)::q0,q,r
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
        !We only use mass to determine centre of mass, so no need to convert to atomic unit
        call StandardizeGeometry(r0,mass,NAtoms,1)
    close(99)
    cartdim=3*NAtoms
    chartemp='Columbus7'; call DefineInternalCoordinate(chartemp,intdim)
!Generate loop geometries
    allocate(q0(intdim)); q0=InternalCoordinateq(r0,intdim,cartdim)!The geometry to calculate Hessian
    allocate(q(intdim)); allocate(r(cartdim))
    open(unit=99,file='geom.all',status='replace')
        do i=1,intdim!Start from q0, displace +- along each internal coordinate
            if(GeometryTransformation_IntCDef(i).motion(1).type=='stretching') then!Bond length: 0.01A
                q=q0; q(i)=q(i)+0.01d0*AInAU; call WriteGeom()
                q=q0; q(i)=q(i)-0.01d0*AInAU; call WriteGeom()
            else!Angle: 0.01
                q=q0; q(i)=q(i)+0.01d0; call WriteGeom()
                q=q0; q(i)=q(i)-0.01d0; call WriteGeom()
            end if
        end do
    close(99)
    contains
    subroutine WriteGeom()
        r=CartesianCoordinater(q,cartdim,intdim,mass=mass,r0=r0)
        do j=1,NAtoms
            write(99,'(A2,I8,3F14.8,F14.8)')ElementSymbol(j),ElementNumber(j),r(3*j-2:3*j),mass(j)
        end do
    end subroutine WriteGeom
end program main