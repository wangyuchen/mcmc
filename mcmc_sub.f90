subroutine mcmc(l,g0,g1,g,p,result)
  implicit none
  integer, intent(in)::l, g0, g1, g
  real, dimension(0:2, l), intent(in) :: p

  integer(kind=8), dimension(l, 0:2),intent(out):: result

  interface
    function prior(total)
      real(kind = 8) :: prior
      integer, dimension(0:2), intent(in) :: total
    end function prior

    subroutine substitute(Iold, Iprop, l, sub_flag, pos)
      integer, intent(in) :: l
      integer, intent(in) :: Iold(:)
      logical, intent(in) :: sub_flag

      integer, dimension(l), intent(out) :: Iprop
      integer, dimension(2), intent(out) :: pos
    end subroutine substitute
  end interface



  ! problem-length vector
  integer, dimension(l) :: Iold, Iprop



  ! iteration-length vector
  real(kind=8), dimension(0:g) :: logprob
  integer, dimension(l, 0:2) :: sum_counts

  real(kind=8) p1, p2, mh, delta1, delta2
  real(kind=8), dimension(l) :: lk2
  logical is_sub
  integer, dimension(0:2) :: total_old, total_prop
  integer, dimension(2) :: pos
  integer :: acpt = 0

  ! misc
  integer :: i
  real :: t1, t2, t3, t4, t5, t_s1 = 0, t_s2 = 0, t_s3 = 0

  call cpu_time(t1)
  ! call init_random_seed()

  ! set initial Iold
  ! g0, g1, and l-g0-g1 is the initial count for 1, 2, 0
  Iold(1:g0) = 1
  Iold((g0+1):(g0+g1)) = 2
  Iold((g0+g1+1):l) = 0
  total_old = (/l-g0-g1, g0, g1/)
  where (Iold == 0)
    lk2 = p(0, :)
  else where (Iold == 1)
    lk2 = p(1, :)
  else where (Iold == 2)
    lk2 = p(2, :)
  end where
  logprob(0) = sum(lk2) + prior(total_old)

  sum_counts=0
  do i = 1, g
    call cpu_time(t2)

    call random_number(rand)
    is_sub = rand < .5
    call substitute(Iold, Iprop, l, is_sub, pos)

    ! if (modulo(i, l*30) == 0) then
    !   print *, total
    !   print *, i
    ! end if

    call cpu_time(t3)
    t_s1 = t_s1 + t3 - t2


    total_prop = total_old

    if (is_sub) then
      total_prop(Iold(pos(1))) = total_prop(Iold(pos(1))) - 1
      total_prop(Iprop(pos(1))) = total_prop(Iprop(pos(1))) + 1
      p1 = 1. / total_old(Iold(pos(1)))
      p2 = 1. / total_prop(Iprop(pos(1)))
    else
      p1 = 1.0
      p2 = 1.0
    end if

    if (is_sub) then
      delta1 = p(Iprop(pos(1)), pos(1)) - p(Iold(pos(1)), pos(1))
      delta2 = prior(total_prop) - prior(total_old)
    else
      delta1 = p(Iprop(pos(1)), pos(1)) - p(Iold(pos(1)), pos(1))
      delta2 = p(Iprop(pos(2)), pos(2)) - p(Iold(pos(2)), pos(2))
    end if
    mh = delta1 + delta2 + log(p2) - log(p1)

    call cpu_time(t4)
    t_s2 = t_s2 + t4 - t3

    ! if accepted, Iold = Iprop
    ! if not, Iold doesn't change
    call random_number(rand)
    if (rand < exp(mh)) then
      total_old = total_prop
      Iold = Iprop
      logprob(i) = logprob(i-1) + delta1 + delta2
    else
      logprob(i) = logprob(i-1)
    end if

    if (i > g/3) then
      if (rand < exp(mh)) then
        acpt = acpt + 1
      end if
      where (Iold == 0)
        sum_counts(:, 0) = sum_counts(:, 0) + 1
      elsewhere (Iold == 1)
        sum_counts(:, 1) = sum_counts(:, 1) + 1
      elsewhere(Iold == 2)
        sum_counts(:, 2) = sum_counts(:, 2) + 1
      end where
    end if

    call cpu_time(t5)
    t_s3 = t_s3 + t5 - t4
  end do
  result= sum_counts
  !open(10, file='result_Iprop.txt')
  !write(10, *), maxloc(sum_counts, dim=2) - 1
  !close(10)

  !open(12, file='sumcounts.txt')
  !write(12, *), sum_counts
  !close(12)

  !open(13, file='maxprob.txt')
  !write(13, *), logprob
  !close(13)



  !print '(f8.4, a)', real(acpt)/l/2, "%"

  !call cpu_time(t5)
  !print *, t_s1, t_s2, t_s3
  !print *, 'total time taken', (t5-t1)

end subroutine mcmc





subroutine substitute(Iold, Iprop, l, sub_flag, pos)
  !---------------------------------------------------------------------
  !
  !  this subroutine randomly substitute/switch one element in Iold.
  !  input: Iold, l,
  !    1. sub_flag (indicator of substitute/switch)
  !  output: p1, p2
  !    1. Iprop: different from Iold in one or two positions
  !    2. pos: if substitute, pos(2) = 0
  !
  !---------------------------------------------------------------------
  implicit none
  integer, intent(in) :: l
  integer, intent(in) :: Iold(:)
  logical, intent(in) :: sub_flag

  integer, dimension(l), intent(out) :: Iprop
  integer, dimension(2), intent(out) :: pos

  real :: rand1, rand2
  integer :: tmk1, tmkn

  ! tmk1 is sampled from (0, 1, 2)
  call random_number(rand1)
  tmk1 = floor(3 * rand1)

  call random_number(rand2)
  ! tmkn is sampled from the other two
  if (rand2 < .5 .and. tmk1 /= 0) then
    tmkn = 0
  else
    if (rand2 > .5 .and. tmk1 /= 2) then
      tmkn = 2
    else
      tmkn = 1
    end if
  end if

  ! find Iold(pos(1)) = tmk1
  call random_number(rand1)
  pos(1) = ceiling(l * rand1)
  do while (Iold(pos(1)) /= tmk1)
    call random_number(rand1)
    pos(1) = ceiling(l * rand1)
  end do

  Iprop = Iold

  ! substitute/switch
  if (sub_flag) then
    pos(2) = 0
    Iprop(pos(1)) = tmkn
  else
    ! find Iold(pos(2)) = tmkn
    call random_number(rand2)
    pos(2) = ceiling(l * rand2)
    do while (Iold(pos(2)) /= tmkn)
      call random_number(rand2)
      pos(2) = ceiling(l * rand2)
    end do

    Iprop(pos(1)) = Iold(pos(2))
    Iprop(pos(2)) = Iold(pos(1))
  end if
end subroutine substitute


subroutine init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine init_random_seed


function prior(total)
  implicit none
  ! interface
  !   function gammln(xx)
  !     REAL(kind=8) gammln
  !     integer j
  !     DOUBLE PRECISION ser, stp, tmp, x,y, cof(6),xx
  !   end function gammln
  ! end interface

  real(kind = 8) :: prior
  integer, dimension(0:2), intent(in) :: total

  real(kind = 8), parameter, dimension(1:3) :: k = (/1, 1, 1/)
  real(kind = 8) :: l1, l2
  integer :: i0, i1, i2

  i0 = total(0)
  i1 = total(1)
  i2 = total(2)

  l1 = dlgama(k(1)+i0) + dlgama(k(2)+i1) + &
       dlgama(k(3)+i2) - dlgama(sum(k)+i0+i1+i2)
  l2 = dlgama(sum(k)) - dlgama(k(1)) - dlgama(k(2)) - dlgama(k(3))

  prior = l1 + l2
end function prior
