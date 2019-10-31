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
    real*8,allocatable,dimension(:,:)::B,L,Linv,cartmode
!Work variable
    character*32::chartemp; integer::i,j; real*8::dbletemp
    real*8,allocatable,dimension(:)::q
    real*8,allocatable,dimension(:,:)::Hessian
!Initialize
    open(unit=99,file='geom',status='old')!Read molecule detail
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do; rewind 99
        allocate(ElementSymbol(NAtoms)); allocate(ElementNumber(NAtoms))
        allocate(r0(3*NAtoms)); allocate(mass(NAtoms))
        do i=1,NAtoms
            read(99,*)ElementSymbol(i),dbletemp,r0(3*i-2:3*i),mass(i)
            ElementSymbol(i)=trim(adjustl(ElementSymbol(i)))
            ElementNumber(i)=int(dbletemp)
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
                read(99,*)chartemp,r(3*j-2:3*j,i),dbletemp
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
    !Construct Wilson B matrix
    allocate(B(intdim,cartdim))
    call WilsonBMatrixAndInternalCoordinateq(B,q,r0,intdim,cartdim)
    !Run Wilson GF method for vibrational frequency and internal normal mode, output vibrational frequency
    allocate(freq(intdim)); allocate(L(intdim,intdim)); allocate(Linv(intdim,intdim))
    call WilsonGFMethod(freq,L,Linv,Hessian,intdim,B,mass,NAtoms)
    open(unit=99,file='VibrationalFrequency.txt',status='replace')
        write(99,'(A4,A1,A15)')'Mode',char(9),'Frequency/cm^-1'
        do i=1,intdim; write(99,'(I4,A1,F14.8)')i,char(9),freq(i)/cm_1InAu; end do
    close(99)
    open(unit=99,file=)
    !Convert internal normal mode to Cartesian, output visualization
    allocate(cartmode(cartdim,intdim))
    call InternalMode2CartesianMode(freq,L,InternalDimension,B,cartmode,CartesianDimension)
    chartemp='geom.log'
    call Avogadro_Vibration(NAtoms,ElementSymbol,r0/AInAU,intdim,freq/cm_1InAu,cartmode,FileName=chartemp)
    write(*,'(1x,A77)')'To visualize the molecular structure and vibration, open geom.log in Avogadro'
end program main