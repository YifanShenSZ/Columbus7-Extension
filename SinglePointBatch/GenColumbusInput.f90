!Standard input files: molecule.xyz, intcfl
!Standard output file: geom.all
!May take more io according to user modification
program main
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization; use Chemistry
    use GeometryTransformation
    implicit none
    type MoleculeDetails!Store the details of an example molecule
        integer::NAtoms
        character*2,allocatable,dimension(:)::ElementSymbol
        real*8,allocatable,dimension(:,:)::RefConfig!Short for REFerence CONFIGuration
        real*8,allocatable,dimension(:)::mass
    end type MoleculeDetails
!Molecule information
    type(MoleculeDetails)::MoleculeDetail
    integer::intdim,cartdim
    real*8,allocatable,dimension(:)::rstandard
!Work variable
    character*32::chartemp
    integer::i,j
    real*8,allocatable,dimension(:)::q0,q,r
!Initialize
    call BetterRandomSeed()
    open(unit=99,file='molecule.xyz',status='old')!Read molecule detail
        read(99,*)MoleculeDetail.NAtoms
            allocate(MoleculeDetail.ElementSymbol(MoleculeDetail.NAtoms))
            allocate(MoleculeDetail.RefConfig(3,MoleculeDetail.NAtoms))
            allocate(MoleculeDetail.mass(MoleculeDetail.NAtoms))
        read(99,*)
        do i=1,MoleculeDetail.NAtoms
            read(99,'(A2,3F20.15)')MoleculeDetail.ElementSymbol(i),MoleculeDetail.RefConfig(:,i)
            MoleculeDetail.ElementSymbol(i)=trim(adjustl(MoleculeDetail.ElementSymbol(i)))
        end do
        read(99,*)
        do i=1,MoleculeDetail.NAtoms
            read(99,*)MoleculeDetail.mass(i)
        end do
    close(99)
    cartdim=3*MoleculeDetail.NAtoms
    chartemp='Columbus7'; call DefineInternalCoordinate(chartemp,intdim)
    allocate(rstandard(cartdim)); rstandard=reshape(MoleculeDetail.RefConfig,[cartdim])
    call StandardizeGeometry(rstandard,MoleculeDetail.mass,MoleculeDetail.NAtoms,1)
!Modification starts here
    open(unit=99,file='q0',status='old')!Reference internal geometry
        allocate(q0(intdim))
        read(99,*)q0
    close(99)
    allocate(q(intdim)); allocate(r(cartdim))
    open(unit=99,file='geom.all',status='replace')
        do i=-4,-1!Start from q0, displace 1st internal coordinate
            q=q0
            q(1)=q(1)+dble(i)*0.04d0
            r=CartesianCoordinater(q,cartdim,intdim,mass=MoleculeDetail.mass,r0=rstandard)
            do j=1,MoleculeDetail.NAtoms
                write(99,'(A2,I8,3F14.8,F14.8)')MoleculeDetail.ElementSymbol(j),Symbol2Number(MoleculeDetail.ElementSymbol(j)),r(3*j-2:3*j),MoleculeDetail.mass(j)
            end do
        end do
    close(99)
end program main