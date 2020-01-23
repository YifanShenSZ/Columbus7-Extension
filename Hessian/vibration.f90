!Analyze vibration based on the internal coordinate Hessian
!Input : internal coordinate definition, geom, hessian
!Output: geom.log
!Dependency: my Fortran-Library, as written in makefile
program main
    use FortranLibrary
    implicit none
!Molecule information
    integer::NAtoms,intdim,cartdim
    character*2,allocatable,dimension(:)::ElementSymbol
    real*8,allocatable,dimension(:)::ElementNumber,mass,r,q
!Wilson GF method information
    real*8,allocatable,dimension(:)::freq
    real*8,allocatable,dimension(:,:)::B,L,Linv,H,cartmode
!Work variable
    character*32::chartemp; logical::flag; integer::i,j,jstart,jstop,jj
!Initialize
    inquire(file='InternalCoordinateDefinition',exist=flag)
    if(flag) then
        chartemp=''
        call DefineInternalCoordinate(chartemp,intdim)
    else
        chartemp='Columbus7'
        call DefineInternalCoordinate(chartemp,intdim)
    end if
    open(unit=99,file='geom',status='old')
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,*)ElementSymbol(i),ElementNumber(i),r(3*i-2:3*i),mass(i)
        end do
        mass=mass*AMUInAU
    close(99)
    open(unit=99,file='hessian',status='old')
        allocate(H(intdim,intdim))
        do i=1,intdim
            do j=1,intdim,8
                jstart=j
                jstop=j+8; if(jstop>intdim) jstop=intdim
                do jj=jstart,jstop-1
                    read(99,'(F13.6)',advance='no')H(i,jj)
                end do
                read(99,'(F13.6)')H(i,jstop)
            end do
        end do
    close(99)
!Do the job
    !The internal coordinate and vibration routines of Columbus use weird unit:
    !    energy in 10^-18 J, length in A (to be continued)
    H=H/4.35974417d0! 1 Hatree = 4.35974417 * 10^-18 J
    do i=1,intdim
        if(GeometryTransformation_IntCDef(i).motion(1).type=='stretching') then
            H(:,i)=H(:,i)/AInAU
            H(i,:)=H(i,:)/AInAU
        end if
    end do
    !Get B matrix
    cartdim=3*NAtoms
    allocate(q(intdim)); allocate(B(intdim,cartdim))
    call WilsonBMatrixAndInternalCoordinateq(B,q,r,intdim,cartdim)
    !Calculate internal coordinate vibration
    allocate(freq(intdim)); allocate(L(intdim,intdim)); allocate(Linv(intdim,intdim))
    call WilsonGFMethod(freq,L,Linv,H,intdim,B,mass,NAtoms)
    !Transform to Cartesian coordinate
    allocate(cartmode(cartdim,intdim))
    call InternalMode2CartesianMode(freq,L,intdim,B,cartmode,cartdim)
    !Output for visualization
    chartemp='geom.log'
    call Avogadro_Vibration(NAtoms,ElementSymbol,r/AInAU,intdim,freq/cm_1InAu,cartmode,FileName=chartemp)
end program main