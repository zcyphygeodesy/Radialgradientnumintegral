!  Radialgradientnumintegral.f90 
!
!  FUNCTIONS:
!  Radialgradientnumintegral - Entry point of console application.
!
!****************************************************************************
      program Radialgradientnumintegral
      implicit none
	character*800::calcpntfl,dwmhgrdfl,gravgrdfl
	integer knd
	real*8::dr
!---------------------------------------------------------------------
      !输入计算点文件名，程序默认计算点在等位边界面上
      !Input the calculation point file on the equipotential boundary surface.
      write(calcpntfl,*)'calcpnt.txt'
      !输入等位边界面大地高格网文件名
      !Input the ellipsoidal height grid file of the equipotential surface.
      write(dwmhgrdfl,*)'landgeoidhgt.dat'
      !输入等位面上残差场元格网文件名。
      !Input residual gravity field element grid file on the equipotential surface.
      !程序要求等位面大地高格网及其面上残差场元格网有相同的格网规格
      !The same grid specifications required for the ellipsoidal height grid of the equipotential
      !surface and residual field element grid on the surface.
      write(gravgrdfl,*)'resGMlgeoid541_1800.ksi'
      !输入积分半径(m)
      dr=120.d3!Integral radius (m)
      write(*, *)"    Begin compulation......"
      call Radialgradientpnt(calcpntfl,dwmhgrdfl,gravgrdfl,dr)
      pause
      end
