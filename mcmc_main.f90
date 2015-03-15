program main
  integer, parameter:: l=2000, g0=100, g1=200, g=l*1000
  real, dimension(0:2, l)::p
  integer(kind=8), dimension(l, 0:2):: result
  integer, dimension(l):: Iprop
  external mcmc

  open(1, file='p_matrix.txt', status='old', action='read')
  read(1,*) p
  print*,p(0:2, 1)
  print*,p(0:2, 2)
  print*, g

  call mcmc(l,g0, g1, g, p, result)
  Iprop=maxloc(result, dim=2) - 1

  open(10, file='sumcounts3.txt')
  write(10, *), result
  close(10)

  open(11, file='result_Iprop3.txt')
  write(11, *), Iprop
  close(11)
end program main
