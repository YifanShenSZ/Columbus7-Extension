!This is to reformat collected Columbus7 data. Take 1 argument $1 = NState
!Please note that this program is only able to deal with up to 9 states
program main
    implicit none
    integer::NState,NAtoms,NPoint
    integer::i,istate,jstate,ip,ia
    character*1::chartemp,chartemp2
    character*128::cmd
    real*8::dbtemp
    real*8,allocatable,dimension(:,:)::energy
    real*8,allocatable,dimension(:,:,:,:)::grad
    real*8,allocatable,dimension(:,:,:,:,:)::nac
    call getarg(1,cmd); read(cmd,*)NState
    open(unit=99,file='energy.temp',status='old')
        NPoint=0; do; read(99,*,iostat=i); if(i/=0) exit; NPoint=NPoint+1; end do; rewind 99
        allocate(energy(NState,NPoint))
        do ip=1,NPoint
            read(99,*)energy(:,ip),dbtemp
        end do
    close(99)
    open(unit=99,file='cartgrd.drt1.state1.temp',status='old')
        NAtoms=0; do; read(99,*,iostat=i); if(i/=0) exit; NAtoms=NAtoms+1; end do
        NAtoms=NAtoms/NPoint
    close(99)
    allocate(grad(3,NAtoms,NState,NPoint))
    allocate( nac(3,NAtoms,NState,NState,NPoint))
    do istate=1,NState
        write(chartemp,'(I1)')istate
        open(unit=99,file='cartgrd.drt1.state'//chartemp//'.temp',status='old')
            do ip=1,NPoint
                do ia=1,NAtoms
                    read(99,*)grad(:,ia,istate,ip)
                end do
            end do
        close(99)
        do jstate=istate+1,NState
            write(chartemp2,'(I1)')jstate
            open(unit=99,file='cartgrd.nad.drt1.state'//chartemp//'.drt1.state'//chartemp2//'.temp',status='old')
                do ip=1,NPoint
                    do ia=1,NAtoms
                        read(99,*)nac(:,ia,istate,jstate,ip)
                    end do
                    nac(:,:,istate,jstate,ip)=nac(:,:,istate,jstate,ip)*(energy(jstate,ip)-energy(istate,ip))
                end do
            close(99)
        end do
    end do
    open(unit=99,file='energy.all',status='replace')
        do ip=1,NPoint
            write(99,*)energy(:,ip)
        end do
    close(99)
    do istate=1,NState
        write(chartemp,'(I1)')istate
        open(unit=99,file='cartgrd.drt1.state'//chartemp//'.all',status='replace')
            do ip=1,NPoint
                do ia=1,NAtoms
                    write(99,*)grad(:,ia,istate,ip)
                end do
            end do
        close(99)
        do jstate=istate+1,NState
            write(chartemp2,'(I1)')jstate
            open(unit=99,file='cartgrd.nad.drt1.state'//chartemp//'.drt1.state'//chartemp2//'.all',status='replace')
                do ip=1,NPoint
                    do ia=1,NAtoms
                        write(99,*)nac(:,ia,istate,jstate,ip)
                    end do
                end do
            close(99)
        end do
    end do
end program main