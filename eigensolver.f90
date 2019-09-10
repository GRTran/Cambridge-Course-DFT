program eigensolver
  !-------------------------------------------------!
  ! This program is designed to illustrate the use  !
  ! of numerical methods and optimised software     !
  ! libraries to solve large Hermitian eigenvalue   !
  ! problems where only the lowest few eigenstates  !
  ! are required. Such problems arise in many       !
  ! electronic structure calculations, where we     !
  ! solve a large Schroedinger equation for just    !
  ! enough states to accommodate the number of      !
  ! electrons.                                      !
  !                                                 !
  ! In this program we will use the special case    !
  ! where the Hamiltonian is real and symmetric.    !
  ! This is a good approximation for large systems. !
  !-------------------------------------------------!
  ! Written by Phil Hasnip (University of York)     !
  ! with contributions from Dave Quigley            !
  ! (University of Warwick)                         !
  !-------------------------------------------------!
  ! Version 0.8, last modified 3rd Sept. 2019       !
  !-------------------------------------------------!
  implicit none

  ! Define what we mean by double-precision
  integer, parameter                              :: dp=selected_real_kind(15,300)

  real(kind=dp),    dimension(:),   allocatable   :: H_kinetic  ! the kinetic energy operator
  real(kind=dp),    dimension(:),   allocatable   :: H_local    ! the local potential operator
  integer                                         :: num_wavevectors ! no. wavevectors
  integer                                         :: num_pw     ! no. plane-waves = 2*wavevectors+1
  integer                                         :: num_states ! no. eigenstates req'd

  complex(kind=dp), dimension(:,:), allocatable   :: trial_wvfn,gradient,search_direction
  real(kind=dp),    dimension(:),   allocatable   :: eigenvalue
  integer                                         :: iter,max_iter
  real(kind=dp)                                   :: energy,prev_energy,energy_tol
  integer                                         :: status

  real(kind=dp), dimension(:),   allocatable      :: products

  integer                                         :: i,j,k,np,nb

  real(kind=dp), dimension(:),   allocatable      :: full_eigenvalues

  ! Declare timing variables. NB intrinsic call is single-precision
  real(kind=kind(0.0))                            :: init_cpu_time
  real(kind=kind(0.0))                            :: curr_cpu_time

  ! Determine whether there should be a preconditioner
  integer :: precon
  character(len=20) :: num1char
  character(len=20) :: num2char

  ! BLAS
  complex(kind=dp),external :: zdotc

  ! No. nonzero wavevectors "G" in our wavefunction expansion

  ! first argument is the number of wavevectors
  call GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values
  ! second argument is the preconditioner presence
  call GET_COMMAND_ARGUMENT(2,num2char)

  read(num1char,*) num_wavevectors                    !then, convert them to REALs
  read(num2char,*) precon
  !num_wavevectors = 200
  !precon = 0

  ! No. plane-waves in our wavefunction expansion. One plane-wave has
  ! wavevector 0, and for all the others there are plane-waves at +/- G
  num_pw = 2*num_wavevectors+1

  ! No. eigenstates to compute
  num_states = 1

  ! Catch any nonsensical combinations of parameters
  if(num_states>=num_pw) stop 'Error, num_states must be less than num_pw'

  ! Desired tolerance on the eigenvalue sum when using an iterative search.
  ! The iterative search will stop when the change in the eigenvalue sum
  ! per iteration is less than this tolerance.
  energy_tol = 1.0e-10_dp

  ! Initialise the random number generator
  call init_random_seed()

  ! Now we need to allocate space for the Hamiltonian terms, and call
  ! the initialisation subroutine
  allocate(H_kinetic(num_pw),stat=status)
  if(status/=0) stop 'Error allocating RAM to H_kinetic'

  allocate(H_local(num_pw),stat=status)
  if(status/=0) stop 'Error allocating RAM to H_local'

  write(*,*) 'Initialising Hamiltonian...'

  ! Initialise the Hamiltonian
  call init_H(num_pw,H_kinetic,H_local)

  ! Perform an exact diagonalisation for comparison
  call cpu_time(init_cpu_time)

  allocate(full_eigenvalues(num_pw),stat=status)
  if(status/=0) stop 'Error allocating RAM to full_eigenvalues'

  write(*,*) 'Starting full diagonalisation...'
  full_eigenvalues = 0.0_dp

  ! You need to complete this subroutine:
  !call exact_diagonalisation(num_pw,H_kinetic,H_local,full_eigenvalues)
  call cpu_time(curr_cpu_time)

  !write(*,*) 'Full diagonalisation took ',curr_cpu_time-init_cpu_time,' secs',full_eigenvalues(1)

  ! Allocate memory for iterative eigenvector search

  allocate(trial_wvfn(num_pw,num_states),gradient(num_pw,num_states),search_direction(num_pw,num_states),stat=status)
  if(status/=0) stop 'Error allocating RAM to trial_wvfn, gradient and search_direction'

  allocate(eigenvalue(num_states),stat=status)
  if(status/=0) stop 'Error allocating RAM to eigenvalues'

  allocate(products(num_states),stat=status)
  if(status/=0) stop 'Error allocating RAM to products'

  write(*,*) 'Starting iterative search for eigenvalues'
  write(*,*) ' '
  write(*,*) '+-----------+----------------+-----------------+'
  write(*,*) '| iteration |     energy     |  energy change  |'
  write(*,*) '+-----------+----------------+-----------------+'

  call cpu_time(init_cpu_time)
  
  ! open the file depending on whether preconditioning has been selected
  if (precon == 1) then
  	open(1560,file='preconditioned_data.dat',access='append')
  else
  	open(1560,file='normal_data.dat',access='append')
  endif

  ! We start from a random guess
  call randomise_state(num_pw,num_states,trial_wvfn)

  ! Our states should be orthogonal and normalised
  call orthonormalise(num_pw,num_states,trial_wvfn)

  ! Apply the H to this state
  call apply_H(num_pw,num_states,trial_wvfn,H_kinetic,H_local,gradient)

  ! Compute the eigenvalues
  do nb=1,num_states
    eigenvalue(nb) = zdotc(num_pw,trial_wvfn(1,nb),1,gradient(1,nb),1)
  end do

  energy = sum(eigenvalue)
  write(*,'(a,g15.8,a,t49,a)') ' |  Initial  | ',energy,'|','|'

  ! In case of problems, we cap the total number of iterations
  max_iter = 40000

  main_loop: do iter=1,max_iter

    prev_energy = energy

    ! The constrained gradient is H.wvfn - (wvfn+.H.wvfn)*wvfn
    ! -- i.e. it is orthogonal to wvfn

    ! You need to complete this subroutine:
    call orthogonalise(num_pw,num_states,gradient,trial_wvfn)

    ! The steepest descent direction is minus the gradient
    call zcopy(num_pw*num_states,gradient,1,search_direction,1)
    call zdscal(num_pw*num_states,-1.0_dp,search_direction,1)

    ! Any modifications to the search direction go here

    ! modify the search direction by implementing a preconditioner
    if (precon == 1) then
    	call precondition(num_pw,num_states,search_direction,trial_wvfn,H_kinetic)
    	! make the search direction orthogonal to the trial wavefunction
    	call orthogonalise(num_pw,num_states,search_direction,trial_wvfn)
    endif

    ! Search along this direction for the best approx. eigenvector
    call line_search(num_pw,num_states,trial_wvfn,H_kinetic,H_local,search_direction,gradient,eigenvalue,energy)

     ! Check convergence
    if(abs(prev_energy-energy)<energy_tol) then
      write(*,*) '+-----------+----------------+-----------------+'
      write(*,*) 'Eigenvalues converged'
      exit
    end if

    energy = sum(eigenvalue)

    write(*,'(a,i5,1x,a,g15.8,a,g15.8,a)') ' |    ',iter,' | ',energy,'| ',prev_energy-energy,' |'

  end do main_loop

  call cpu_time(curr_cpu_time)

  write(*,*) 'Iterative search took ',curr_cpu_time-init_cpu_time,' secs'
  
  ! here we write the data for the preconditioned and non preconditioned solution to file with append
  write(1560,*) curr_cpu_time-init_cpu_time, full_eigenvalues(1), eigenvalue(1), iter

  ! Finally summarise the results
  write(*,'(a,t17,a,t38,a)') 'Eigenvalue','Iterative','Exact'
  do nb=1,num_states
    write(*,'(i5,5x,2g20.10)') nb,eigenvalue(nb),full_eigenvalues(nb)
  end do
  
  ! close the output writer file
  close(1560)
  call output_results(num_pw,num_states,H_local,trial_wvfn)

  ! Deallocate memory
  deallocate(full_eigenvalues,stat=status)
  if(status/=0) stop 'Error deallocating RAM from full_eigenvalues'

  deallocate(search_direction,stat=status)
  if(status/=0) stop 'Error deallocating RAM from search_direction'

  deallocate(gradient,stat=status)
  if(status/=0) stop 'Error deallocating RAM from gradient'

  deallocate(trial_wvfn,stat=status)
  if(status/=0) stop 'Error deallocating RAM from trial_wvfn'

end program eigensolver

  ! -- YOU WILL NEED TO FINISH THESE SUBROUTINES --

  subroutine exact_diagonalisation(num_pw,H_kinetic,H_local,full_eigenvalues)
    !-------------------------------------------------!
    ! This subroutine takes a compact representation  !
    ! of the matrix H, constructs the full H, and     !
    ! diagonalises to get the whole eigenspectrum.    !
    !-------------------------------------------------!
    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                    intent(in)  :: num_pw
    real(kind=dp),  dimension(num_pw),          intent(in)  :: H_kinetic
    real(kind=dp),  dimension(num_pw),          intent(in)  :: H_local
    real(kind=dp),  dimension(num_pw),          intent(out) :: full_eigenvalues

    complex(kind=dp), dimension(:,:), allocatable :: full_H
    real(kind=dp),    dimension(:),   allocatable :: lapack_real_work
    complex(kind=dp), dimension(:),   allocatable :: lapack_cmplx_work
    integer                                       :: status,np
    integer                                       :: info

    ! Delete this line once you've coded this subroutine
    ! stop 'Subroutine exact_diagonalisation has not been written yet'

    ! First we allocate and construct the full H
    allocate(full_H(num_pw,num_pw),stat=status)
    if(status/=0) stop 'Error allocating RAM to full_H in exact_diagonalisation'

    call construct_full_H(num_pw,H_kinetic,H_local,full_H)

    ! Use LAPACK to get eigenvalues and eigenvectors
    ! NB H is Hermitian

    ! allocating the lapack work
    allocate(lapack_real_work(3*num_pw-2))
    ! allocate(lapack_real_work(6633))
    allocate(lapack_cmplx_work(2*num_pw-1))
    write(*,*) 'using zheev'
    call zheev('V','U',num_pw,full_H,num_pw,full_eigenvalues,lapack_cmplx_work, &
                size(lapack_cmplx_work),lapack_real_work,info)
    if(info/=0) then
      write(*,*) info
      stop 'error in zheev'
    else
      write(*,*) 'optimal size of work :: ', lapack_cmplx_work(1)
    end if

! character 	JOBZ,
! character 	UPLO,
! integer 	N,
! complex*16, dimension( lda, * ) 	A,
! integer 	LDA,
! double precision, dimension( * ) 	W,
! complex*16, dimension( * ) 	WORK,
! integer 	LWORK,
! double precision, dimension( * ) 	RWORK,
! integer 	INFO
! )





    ! Deallocate memory
    deallocate(full_H,stat=status)
    if(status/=0) stop 'Error deallocating RAM from full_H in exact_diagonalisation'

    return
  end subroutine exact_diagonalisation

  subroutine orthogonalise(num_pw,num_states,state,ref_state)
    !-------------------------------------------------!
    ! This subroutine takes a set of states and       !
    ! orthogonalises them to a set of reference       !
    ! states.                                         !
    !-------------------------------------------------!
    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                        intent(in)    :: num_pw
    integer,                                        intent(in)    :: num_states
    complex(kind=dp), dimension(num_pw,num_states), intent(inout) :: state
    complex(kind=dp), dimension(num_pw,num_states), intent(in)    :: ref_state

    complex(kind=dp)                             :: overlap
    integer                                      :: nb1,nb2,np
    integer                                      :: status

    ! Delete this line once you've coded this subroutine
    ! stop 'Subroutine orthogonalise has not been written yet'

    do nb2=1,num_states
      do nb1=1,num_states
        overlap = cmplx(0.0_dp,0.0_dp,dp)

        ! Compute overlap of state nb2 with ref_state nb1
        overlap = dot_product(state(:,nb1),conjg(ref_state(:,nb1)))

        ! Remove overlapping parts of ref_state nb1 from state nb2
        state(:,nb1) = state(:,nb1) - overlap*ref_state(:,nb1)

      end do
    end do

    return

  end subroutine orthogonalise

  subroutine precondition(num_pw,num_states,search_direction,trial_wvfn,H_kinetic)
    !-------------------------------------------------!
    ! This subroutine takes a search direction and    !
    ! applies a simple kinetic energy-based           !
    ! preconditioner to improve the conditioning of   !
    ! the eigenvalue search.                          !
    !-------------------------------------------------!
    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                        intent(in)    :: num_pw
    integer,                                        intent(in)    :: num_states
    complex(kind=dp), dimension(num_pw,num_states), intent(inout) :: search_direction
    complex(kind=dp), dimension(num_pw,num_states), intent(in)    :: trial_wvfn
    real(kind=dp),    dimension(num_pw),            intent(in)    :: H_kinetic

    integer                                      :: np,nb
    real(kind=dp)                                :: kinetic_eigenvalue,x,temp

    ! Delete this line once you've coded this subroutine
    !stop 'Subroutine precondition has not been written yet'

    kinetic_eigenvalue = 0.0_dp
    do nb=1,num_states
       ! Compute kinetic energy eigenvalue for state nb
       ! Our best guess at this is currently the Rayleigh quotient
       ! wvfn+.H_kinetic.wvfn. H_kinetic is diagonal and stored as a
       ! vector so the ith element of H_kinetic.wvfn is H_kinetic(i)*wvfn(i)
       ! Hence....
	
     
       kinetic_eigenvalue = sum((conjg(trial_wvfn(:,nb))*H_kinetic(:)*trial_wvfn(:,nb))) &
	/ sum((trial_wvfn(:,nb)*conjg(trial_wvfn(:,nb))))


	!print*, kinetic_eigenvalue
 	!kinetic_eigenvalue = H_kinetic(num_pw)
       do np=1, num_pw
          ! Compute and apply the preconditioning.
          
	  x = H_kinetic(np) / (2._dp*kinetic_eigenvalue)
          !print*, x, H_kinetic(np), kinetic_eigenvalue

	  temp = (8._dp + 4._dp*x + 2._dp*x**2 + x**3) / (8._dp + 4._dp*x + 2._dp*x**2 + x**3 + x**4)

	  search_direction(np,nb) = search_direction(np,nb)* temp

	  ! preconditioning is applied to variable search_direction
          
          !print*, temp
       end do
       
       

    end do
	!stop
    return

  end subroutine precondition

  subroutine diagonalise(num_pw,num_states,state,H_state,eigenvalues,rotation)
    !-------------------------------------------------!
    ! This subroutine takes a set of states and       !
    ! H acting on those states, and transforms the    !
    ! states to diagonalise <state|H|state>.          !
    !-------------------------------------------------!
    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                            intent(in)    :: num_pw
    integer,                                            intent(in)    :: num_states
    complex(kind=dp), dimension(num_pw,num_states),     intent(inout) :: state
    complex(kind=dp), dimension(num_pw,num_states),     intent(inout) :: H_state
    real(kind=dp), dimension(num_states),               intent(inout) :: eigenvalues
    complex(kind=dp), dimension(num_states,num_states), intent(out)   :: rotation

    integer                                      :: nb1,nb2
    integer                                      :: optimal_size,status

    real(kind=dp),    dimension(:),  allocatable :: lapack_real_work
    complex(kind=dp), dimension(:),  allocatable :: lapack_cmplx_work

    complex(kind=dp), external :: zdotc

    ! Delete this line once you've coded this subroutine
    stop 'Subroutine diagonalise has not been written yet'

    ! Compute the subspace H matrix and store in rotation array
    do nb2=1,num_states
      do nb1=1,num_states
        rotation(nb1,nb2) = zdotc(num_pw,state(:,nb1),1,H_state(:,nb2),1)
      end do
    end do

    ! Compute diagonalisation
    !
    ! determine optimal size of work arrays

    ! allocate memory for work arrays

    ! zero the work array

    ! Use LAPACK to diagonalise the Hamiltonian in this subspace
    ! NB H is Hermitian

    ! Finally apply the diagonalising rotation

    call transform(num_pw,num_states,state,rotation)
    call transform(num_pw,num_states,H_state,rotation)

    return

  end subroutine diagonalise

  ! -- THE FOLLOWING SUBROUTINES ARE ALREADY WRITTEN --
  !       (you may wish to optimise them though)

  subroutine orthonormalise(num_pw,num_states,state)
    !-------------------------------------------------!
    ! This subroutine takes a set of states and       !
    ! orthonormalises them.                           !
    !-------------------------------------------------!
    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                        intent(in)    :: num_pw
    integer,                                        intent(in)    :: num_states
    complex(kind=dp), dimension(num_pw,num_states), intent(inout) :: state

    complex(kind=dp), dimension(:,:), allocatable   :: overlap
    complex(kind=dp)                                :: tmp_overlap
    integer                                      :: nb2,nb1,np
    integer                                      :: status

    complex(kind=dp), external :: zdotc

    allocate(overlap(num_states,num_states),stat=status)
    if(status/=0) stop 'Error allocating RAM to overlap in orthonormalise'

    ! Compute the overlap matrix
    do nb2=1,num_states
      do nb1=1,num_states
        overlap(nb1,nb2) = zdotc(num_pw,state(:,nb1),1,state(:,nb2),1)
      end do
    end do

    ! Compute orthogonalising transformation

    ! First compute Cholesky (U.U^H) factorisation of the overlap matrix
    call zpotrf('U',num_states,overlap,num_states,status)
    if(status/=0) stop 'zpotrf failed in orthonormalise'

    ! invert this upper triangular matrix
    call ztrtri('U','N',num_states,overlap,num_states,status)
    if(status/=0) stop 'ztrtri failed in orthonormalise'

    ! Set lower triangle to zero
    do nb2 = 1,num_states
      do nb1 = nb2+1,num_states
        overlap(nb1,nb2)=cmplx(0.0_dp,0.0_dp,dp)
      end do
    end do

    ! overlap array now contains the (upper triangular) orthonormalising transformation
    call transform(num_pw,num_states,state,overlap)

    deallocate(overlap,stat=status)
    if(status/=0) stop 'Error deallocating RAM from overlap in orthonormalise'

    return

  end subroutine orthonormalise

  subroutine transform(num_pw,num_states,state,transformation)
    !-------------------------------------------------!
    ! This subroutine takes a set of states and       !
    ! orthonormalises them.                           !
    !-------------------------------------------------!
    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                            intent(in)    :: num_pw
    integer,                                            intent(in)    :: num_states
    complex(kind=dp), dimension(num_pw,num_states),     intent(inout) :: state
    complex(kind=dp), dimension(num_states,num_states), intent(inout) :: transformation

    integer                                        :: nb2,nb1,np
    integer                                        :: status
    complex(kind=dp), dimension(:,:), allocatable  :: new_state

    allocate(new_state(size(state,1),size(state,2)),stat=status)
    if(status/=0) stop 'Error allocating RAM to new_state in transform'

    new_state = cmplx(0.0_dp,0.0_dp,dp)
    do nb2=1,num_states
      do nb1=1,num_states
        do np=1,num_pw
          new_state(np,nb1) = new_state(np,nb1) + state(np,nb2)*transformation(nb2,nb1)
        end do
      end do
    end do

    state = new_state

    deallocate(new_state,stat=status)
    if(status/=0) stop 'Error deallocating RAM from new_state in transform'

    return

  end subroutine transform

  subroutine line_search(num_pw,num_states,approx_state,H_kinetic,H_local,direction, &
                       & gradient,eigenvalue,energy)
    !-------------------------------------------------!
    ! This subroutine takes an approximate eigenstate !
    ! and searches along a direction to find an       !
    ! improved approximation.                         !
    !-------------------------------------------------!
    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                        intent(in)    :: num_pw
    integer,                                        intent(in)    :: num_states
    complex(kind=dp), dimension(num_pw,num_states), intent(inout) :: approx_state
    real(kind=dp), dimension(num_pw),               intent(in)    :: H_kinetic
    real(kind=dp), dimension(num_pw),               intent(in)    :: H_local
    complex(kind=dp), dimension(num_pw,num_states), intent(inout) :: direction
    complex(kind=dp), dimension(num_pw,num_states), intent(inout) :: gradient
    real(kind=dp),                                  intent(inout) :: energy
    real(kind=dp), dimension(num_states),           intent(inout) :: eigenvalue

    real(kind=dp)                                 :: init_energy
    real(kind=dp)                                 :: tmp_energy
!    real(kind=dp), save                           :: trial_step=0.00004_dp
    real(kind=dp), save                           :: trial_step=0.4_dp
    real(kind=dp)                                 :: step
    real(kind=dp)                                 :: opt_step
    complex(kind=dp), dimension(:,:), allocatable :: tmp_state
    real(kind=dp)                                 :: d2E_dstep2
    real(kind=dp)                                 :: best_step
    real(kind=dp)                                 :: best_energy
    real(kind=dp)                                 :: denergy_dstep
    real(kind=dp)                                 :: mean_norm
    integer                                       :: loop
    integer                                       :: nb
    integer                                       :: status

    ! BLAS external functions
    real(kind=dp),    external :: dznrm2
    complex(kind=dp), external :: zdotc

    ! To try to keep a convenient step length, we reduce the size of the search direction
    mean_norm = 0.0_dp
    do nb=1,size(approx_state,2)
      mean_norm = mean_norm + dznrm2(num_pw,direction,1)
    end do
    mean_norm = mean_norm/real(size(approx_state,2),dp)

    call zdscal(num_pw*num_states,1.0_dp/mean_norm,direction,1)

    ! The rate-of-change of the energy is just 2*Re(conjg(direction).gradient)
    denergy_dstep = 0.0_dp
    do nb=1,size(approx_state,2)
      denergy_dstep = denergy_dstep + 2*real(zdotc(num_pw,direction(:,nb),1,gradient(:,nb),1),dp)
    end do

    allocate(tmp_state(size(approx_state,1),size(approx_state,2)),stat=status)
    if(status/=0) stop 'Error allocating RAM to tmp_state in line_search'

    best_step   = 0.0_dp
    best_energy = energy

    ! First take a trial step in the direction
    step = trial_step

    ! We find a trial step that lowers the energy:
    do loop=1,10

      tmp_state = approx_state + step*direction

      call orthonormalise(num_pw,num_states,tmp_state)

      ! Apply the H to this state
      call apply_H(num_pw,num_states,tmp_state,H_kinetic,H_local,gradient)

      ! Compute the new energy estimate
      tmp_energy = 0.0_dp
      do nb=1,num_states
        tmp_energy = tmp_energy + real(zdotc(num_pw,tmp_state(:,nb),1,gradient(:,nb),1),dp)
      end do

      if(tmp_energy<energy) then
        exit
      else
        d2E_dstep2 = (tmp_energy - energy - step*denergy_dstep )/(step**2)
        if(d2E_dstep2<0.0_dp) then
          if(tmp_energy<energy) then
            exit
          else
            step = step/4.0_dp
          end if
        else
          step  = -denergy_dstep/(2*d2E_dstep2)
        end if
      end if

    end do

    if(tmp_energy<best_energy) then
      best_step   = step
      best_energy = tmp_energy
    end if

    ! We now have the initial eigenvalue, the initial gradient, and a trial step
    ! -- we fit a parabola, and jump to the estimated minimum position
    ! Set default step and energy
    d2E_dstep2 = (tmp_energy - energy - step*denergy_dstep )/(step**2)


    if(d2E_dstep2<0.0_dp) then
      ! Parabolic fit gives a maximum, so no good
      write(*,'(a)') '** Warning, parabolic stationary point is a maximum **'

      if(tmp_energy<energy) then
        opt_step = step
      else
        opt_step = 0.1_dp*step
      end if
    else
      opt_step  = -denergy_dstep/(2*d2E_dstep2)
    end if


!    e = e0 + de*x + c*x**2
! => c = (e - e0 - de*x)/x**2
! => min. at -de/(2c)
!
!    de/dx = de + 2*c*x

    call zaxpy(num_pw*num_states,cmplx(opt_step,0.0_dp,dp),direction,1,approx_state,1)

    call orthonormalise(num_pw,num_states,approx_state)

    ! Apply the H to this state
    call apply_H(num_pw,num_states,approx_state,H_kinetic,H_local,gradient)

    ! Compute the new energy estimate
    energy = 0.0_dp
    do nb=1,size(approx_state,2)
      eigenvalue(nb) = real(zdotc(num_pw,approx_state(:,nb),1,gradient(:,nb),1),dp)
      energy = energy + eigenvalue(nb)
    end do

    ! This ought to be the best, but check...
    if(energy>best_energy) then
!      if(best_step>0.0_dp) then
      if(abs(best_step-epsilon(1.0_dp))>0.0_dp) then

        call zaxpy(num_pw*num_states,cmplx(best_step,0.0_dp,dp),direction,1,gradient,1)

        call orthonormalise(num_pw,num_states,approx_state)

        ! Apply the H to this state
        call apply_H(num_pw,num_states,approx_state,H_kinetic,H_local,gradient)

        ! Compute the new energy estimate
        energy = 0.0_dp
        do nb=1,size(approx_state,2)
          eigenvalue(nb) = real(zdotc(num_pw,approx_state(:,nb),1,gradient(:,nb),1),dp)
          energy = energy + eigenvalue(nb)
        end do

      else
        write(*,*) 'Oh dear:',best_step
        stop 'Problem with line search'
      end if
    end if

!    write(*,'(3f25.15,a)') opt_step,step,energy,' <- test2'

    ! We'll use this step as the basis of our trial step next time
!    trial_step = 2*opt_step

    deallocate(tmp_state,stat=status)
    if(status/=0) stop 'Error deallocating RAM from tmp_state in line_search'

    return

  end subroutine line_search

  subroutine output_results(num_pw,num_states,H_local,wvfn)
    !-------------------------------------------------!
    ! This subroutine writes the potential and the    !
    ! wavefunction to two files (pot.dat and wvfn.dat !
    ! respectively).                                  !
    !-------------------------------------------------!
    use fftw3

    implicit none

    ! Define what we mean by double-precision
    integer, parameter                           :: dp=selected_real_kind(15,300)

    integer,                                        intent(in)    :: num_pw
    integer,                                        intent(in)    :: num_states
    real(kind=dp), dimension(num_pw),               intent(in)    :: H_local
    complex(kind=dp), dimension(num_pw,num_states), intent(in)    :: wvfn

    integer                                                       :: np,nb,status
    character(len=255)                                            :: filename

    complex(kind=dp), dimension(:),                allocatable    :: realspace_wvfn
    complex(kind=dp), dimension(:),                allocatable    :: tmp_wvfn

    type(C_PTR) :: plan                       ! defined in fftw3

    ! First write the local potential
    open(unit=10,file='pot.dat',form='formatted')

    do np=1,num_pw
      write(10,*) real(np-1,kind=dp)/real(num_pw,kind=dp),H_local(np)
    end do

    close(10)

    ! Now FFT and write the eigenstates
    allocate(realspace_wvfn(num_pw),stat=status)
    if(status/=0) stop 'Error allocating realspace in output_results'

    allocate(tmp_wvfn(num_pw),stat=status)
    if(status/=0) stop 'Error allocating tmp_wvfn in output_results'

    ! First create a plan. We want a forward transform and we want to estimate
    ! (rather than measure) the most efficient way to perform the transform.
    plan = fftw_plan_dft_1d(num_pw,tmp_wvfn,realspace_wvfn,FFTW_FORWARD,FFTW_ESTIMATE)

    do nb=1,num_states
      write(filename,*) nb
      filename='wvfn_'//trim(adjustl(filename))//'.dat'
      open(unit=10,file=filename,form='formatted')

      tmp_wvfn = wvfn(:,nb)

      ! Compute the FFT of in using the current plan, store in out.
      call fftw_execute_dft(plan,tmp_wvfn,realspace_wvfn)

      do np=1,num_pw
        write(10,*) real(np-1,kind=dp)/real(num_pw,kind=dp),real(realspace_wvfn(np),dp)
      end do

      close(10)
    end do

    call fftw_destroy_plan(plan)
    deallocate(realspace_wvfn,stat=status)
    if(status/=0) stop 'Error deallocating realspace in output_results'

    return

  end subroutine output_results
