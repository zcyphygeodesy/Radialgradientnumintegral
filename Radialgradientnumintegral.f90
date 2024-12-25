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
      !���������ļ���������Ĭ�ϼ�����ڵ�λ�߽�����
      !Input the calculation point file on the equipotential boundary surface.
      write(calcpntfl,*)'calcpnt.txt'
      !�����λ�߽����ظ߸����ļ���
      !Input the ellipsoidal height grid file of the equipotential surface.
      write(dwmhgrdfl,*)'landgeoidhgt.dat'
      !�����λ���ϲвԪ�����ļ�����
      !Input residual gravity field element grid file on the equipotential surface.
      !����Ҫ���λ���ظ߸����������ϲвԪ��������ͬ�ĸ������
      !The same grid specifications required for the ellipsoidal height grid of the equipotential
      !surface and residual field element grid on the surface.
      write(gravgrdfl,*)'resGMlgeoid541_1800.ksi'
      !������ְ뾶(m)
      dr=120.d3!Integral radius (m)
      write(*, *)"    Begin compulation......"
      call Radialgradientpnt(calcpntfl,dwmhgrdfl,gravgrdfl,dr)
      pause
      end
