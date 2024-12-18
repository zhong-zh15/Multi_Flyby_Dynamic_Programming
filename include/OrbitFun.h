#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_
#include <math.h>       /* isnan, sqrt */
#include "Constant.h"
#include "OrbitMath.h"
//#include "MinpackSolver.h"
/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ���춯��ѧ�����ʵ����
* ��ϵ��: ��� 1522620129@qq.com
* �ļ���: OrbitFun.h
* ���ݼ���������������ת���ͽ������ƺ���
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2022-09-15    ���     ������ļ��������ƽ˲��ת���͵���
****************************************************************************/

/********************************�������������ֽǶȹ�ϵ********************************/
// E2f ����ƫ����Ǻ�ƫ������������
double E2f(int& flag, double E, double e);

// E2M ����ƫ����Ǻ�ƫ������ƽ�����
double E2M(int& flag, double E, double e);

// f2E ���������Ǻ�ƫ������ƫ�����
double f2E(int& flag, double f, double e);

// M2E ����ƽ����Ǻ�ƫ������ƫ�����
double M2E(int& flag, double M, double e, int MaxIter = 100, double epsilon = 1.0e-14);

// f0dt2ft ���ݳ�ʼ�����Ǻ��ݻ�ʱ��������������
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter = 100, double epsilon = 1.0e-14);

// f0ft2dt ���ݳ�ʼ�����Ǻ��������������ݻ�ʱ��
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu);

/*******************************������������ֱ�����ꡢ�Ľ����ֵ���������ת��*******************************/
// coe2rv ���ݾ��������������ֱ������ϵ�µ�λ�ú��ٶȷ���
void coe2rv(int& flag, double* rv, const double* coe, double mu);

// rv2coe ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ����󾭵�������
void rv2coe(int& flag, double* coe, const double* RV, double mu);

/***********************************��֪��ֵ��ʱ�䣬��ĩ��״̬***************************************/
//���ݳ�ʼʱ��״̬coe0��ĩ��ʱ��dt��״̬coe1���������ƽ���������ɹ�,flag����1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu);

//���ݳ�ʼʱ��t0��״̬rv0��ĩ��ʱ��t1��״̬rv1���������ƽ���������ɹ�,flag����1
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu);

#endif