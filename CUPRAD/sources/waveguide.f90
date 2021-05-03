PROGRAM waveguide
  IMPLICIT NONE

  INTEGER(4) i_x_max, i_z_max, i_x,i_z
  REAL(8), ALLOCATABLE     :: xx(:),zz(:),Indice(:,:)

  i_x_max = 6
  i_z_max = 3

  ALLOCATE(xx(i_x_max),zz(i_z_max),Indice(i_x_max,i_z_max))

  xx(1) = 0.D0
  xx(2) = 10.D0
  xx(3) = 10.01D0
  xx(4) = 20.D0
  xx(5) = 20.01D0
  xx(6) = 30.D0

  zz(1) = 0.D0
  zz(2) = 1000.D0
  zz(3) = 2000.D0

  Indice = 1.5134D0
  Indice(1:2,:) = 1.5147D0
  Indice(3:4,3) = 1.5147D0

  open(12,file='waveguide.inp',action="write",form="unformatted")
  write(12)i_x_max, i_z_max
  write(12)(xx(i_x),i_x=1,i_x_max)
  do i_z = 1, i_z_max
     write(12) zz(i_z)
     write(12) (Indice(i_x,i_z),i_x=1,i_x_max)
  enddo
  close(12)

END PROGRAM waveguide
