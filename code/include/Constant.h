#ifndef CONSTANT
#define CONSTANT

/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ���춯��ѧ�����ʵ����
* ��ϵ��: ��� 1522620129@qq.com
* �ļ���: Constant.h
* ���ݼ�������������ļ������ĳ��� 
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2022-09-15    ���     ������ļ�����ȫ�˲����õ��ĳ���
****************************************************************************/

//��ѧ����	
const double DPI = 3.1415926535897932384626433832795; //Բ����
const double D2PI = 6.283185307179586476925286766559;
const double D2R = 0.017453292519943295769236907684886;
const double R2D = 57.295779513082320876798154814105;
const double S2R = 4.848136811095359935899141E-6;
const double EPSILON = 1.0e-14;

//������س���
const double JD2Second = 86400.0;       //��תs
const double muEarth = 398600.4415e9;	//������������(m^3/s^2)
const double Re = 6378137.0;			//����뾶(m)
const double J2 = 1.08262668e-3;	    //J2�㶯
const double g0 = 9.80665;              //�������ٶ�(m/s^2)

const double OmegaEarth = 7.292115146251017e-005;//7.292158553e-5;//rad/s
const double EarthFlattening = 1.0 / 298.257223563;  //���ʣ�
const double EarthPoleRadius = 6356752.3142;//m
const double EarthRefRadius = 6378136.3;//m

//����������س���
const double muSun = 1.32712440018e20;//m^3/s^2
const double muMoon = 4.902801076e12; //m^3/s^2
const double AU = 1.49597870700e11;//m
const double Rm = 1738200.0;//m

const double LightSpeed = 2.99792458e8;//m/s
const double LuminositySun = 3.823e26;  //̫�����䣿

//ʱ��ת�����
const double JD_MJD = 2400000.5;
const double TT_TAI = 32.184;
const double MJDJ2000 = 51544.5;
const double JC2JD = 36525.0;


#endif

