#ifndef _LAMBERT_H_
#define _LAMBERT_H_

/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ���춯��ѧ�����ʵ����
* ��ϵ��: ��� 1522620129@qq.com
* �ļ���: Lambert.h
* ���ݼ�����
*           ����ʦ����Izzo���lambert�����򣨽�Ϊͨ�ã�
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2022-09-15    ���     �����ļ�����Ӳ���ע��
****************************************************************************/

/***********************************Izzo���Lambert����**********************************/
//Subfunction that evaluates the time of flight as a function of x
double x2tof(double x, double s, double c, int lw, int N);

void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2,
	double tf, const double* unith, int& flag, double mu, int way = 0, int N = 0, int branch = 0,
	int Maxiter = 200, double tol = 1.0e-11);

//��Ȧlambert�������
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double GM);

/***********************************Lambert�����**********************************/
//������lambert_solver
//������lambert�������
class lambert_solver
{
public:
	lambert_solver(const double* rv1, const double* rv2, const double tof, const double mu, const int way = 0);
	void solve_single(const int N = 0, const int branch = 0);
	int solve_multi();
	void get_dv1(double* dv1);
	void get_dv2(double* dv2);
	double get_Mdv1();
	double get_Mdv2();


	double m_rv1[6], m_rv2[6], m_tof, m_mu;
	int m_way;
	int m_flag;
	double m_dv1[3], m_dv2[3], m_Mdv1, m_Mdv2;
	int m_N, m_branch;

private:

};

#endif
