!Generate a scan along a internal coordinate
!Input : internal coordinate definition, geom
!Output: geom.data
!Dependency: my Fortran-Library, as written in makefile
program main
    use FortranLibrary
    implicit none
!Molecule information
    integer::NAtoms,intdim,cartdim
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::ElementNumber,mass,r0,q0
!Work variable
    character*32::chartemp; logical::flag; integer::i,j
    real*8,allocatable,dimension(:)::q,r
!Initialize
    inquire(file='InternalCoordinateDefinition',exist=flag)
    if(flag) then
        chartemp=''
        call DefineInternalCoordinate(chartemp,intdim)
    else
        chartemp='Columbus7'
        call DefineInternalCoordinate(chartemp,intdim)
    end if
    open(unit=99,file='geom',status='old')!Read reference geometry
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r0(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,*)ElementSymbol(i),ElementNumber(i),r0(3*i-2:3*i),mass(i)
        end do
        mass=mass*AMUInAU
    close(99)
    cartdim=3*NAtoms
    allocate(q0(intdim)); q0=InternalCoordinateq(r0,intdim,cartdim)
!Generate a scan along a internal coordinate
    chartemp='assimilate'
    allocate(q(intdim)); allocate(r(cartdim))
    open(unit=99,file='geom.data',status='replace')
        !Example: start from q0, scan 1st internal coordinate
        do i=-5,5
            q=q0
            q(1)=q(1)+dble(i)*0.01d0
            r=CartesianCoordinater(q,cartdim,intdim,uniquify=chartemp,mass=mass,r0=r0)
            do j=1,NAtoms
                write(99,'(1x,A2,2x,F5.1,4F14.8)')ElementSymbol(j),ElementNumber(j),r(3*j-2:3*j),mass(j)/AMUInAU
            end do
        end do
    close(99)
end program main
