#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_
#include <math.h>       /* isnan, sqrt */
#include "Constant.h"
#include "OrbitMath.h"
//#include "MinpackSolver.h"
/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院航天动力学与控制实验室
* 联系人: 武迪 1522620129@qq.com
* 文件名: OrbitFun.h
* 内容简述：整理轨道根数转换和解析递推函数
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2022-09-15    武迪     整理该文件，添加了平瞬根转换和递推
****************************************************************************/

/********************************经典轨道根数三种角度关系********************************/
// E2f 根据偏近点角和偏心率求真近点角
double E2f(int& flag, double E, double e);

// E2M 根据偏近点角和偏心率求平近点角
double E2M(int& flag, double E, double e);

// f2E 根据真近点角和偏心率求偏近点角
double f2E(int& flag, double f, double e);

// M2E 根据平近点角和偏心率求偏近点角
double M2E(int& flag, double M, double e, int MaxIter = 100, double epsilon = 1.0e-14);

// f0dt2ft 根据初始真近点角和演化时间求最终真近点角
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter = 100, double epsilon = 1.0e-14);

// f0ft2dt 根据初始真近点角和最终真近点角求演化时间
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu);

/*******************************经典轨道根数、直角坐标、改进春分点轨道根数的转换*******************************/
// coe2rv 根据经典轨道根数求惯性直角坐标系下的位置和速度分量
void coe2rv(int& flag, double* rv, const double* coe, double mu);

// rv2coe 根据地心惯性直角坐标系下的位置和速度分量求经典轨道根数
void rv2coe(int& flag, double* coe, const double* RV, double mu);

/***********************************已知初值和时间，求末端状态***************************************/
//根据初始时刻状态coe0求末端时刻dt的状态coe1，按二体推进。若计算成功,flag返回1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu);

//根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu);

#endif