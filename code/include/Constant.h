#ifndef CONSTANT
#define CONSTANT

/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院航天动力学与控制实验室
* 联系人: 武迪 1522620129@qq.com
* 文件名: Constant.h
* 内容简述：整理各项文件依赖的常数 
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2022-09-15    武迪     整理该文件，补全了部分用到的常数
****************************************************************************/

//数学常数	
const double DPI = 3.1415926535897932384626433832795; //圆周率
const double D2PI = 6.283185307179586476925286766559;
const double D2R = 0.017453292519943295769236907684886;
const double R2D = 57.295779513082320876798154814105;
const double S2R = 4.848136811095359935899141E-6;
const double EPSILON = 1.0e-14;

//地球相关常数
const double JD2Second = 86400.0;       //天转s
const double muEarth = 398600.4415e9;	//地球引力常数(m^3/s^2)
const double Re = 6378137.0;			//地球半径(m)
const double J2 = 1.08262668e-3;	    //J2摄动
const double g0 = 9.80665;              //重力加速度(m/s^2)

const double OmegaEarth = 7.292115146251017e-005;//7.292158553e-5;//rad/s
const double EarthFlattening = 1.0 / 298.257223563;  //扁率？
const double EarthPoleRadius = 6356752.3142;//m
const double EarthRefRadius = 6378136.3;//m

//其他天体相关常数
const double muSun = 1.32712440018e20;//m^3/s^2
const double muMoon = 4.902801076e12; //m^3/s^2
const double AU = 1.49597870700e11;//m
const double Rm = 1738200.0;//m

const double LightSpeed = 2.99792458e8;//m/s
const double LuminositySun = 3.823e26;  //太阳辐射？

//时间转换相关
const double JD_MJD = 2400000.5;
const double TT_TAI = 32.184;
const double MJDJ2000 = 51544.5;
const double JC2JD = 36525.0;


#endif

