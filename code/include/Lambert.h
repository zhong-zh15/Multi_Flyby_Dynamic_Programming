#ifndef _LAMBERT_H_
#define _LAMBERT_H_

/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院航天动力学与控制实验室
* 联系人: 武迪 1522620129@qq.com
* 文件名: Lambert.h
* 内容简述：
*           蒋老师翻译Izzo组的lambert求解程序（较为通用）
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2022-09-15    武迪     创建文件并添加部分注释
****************************************************************************/

/***********************************Izzo组的Lambert程序**********************************/
//Subfunction that evaluates the time of flight as a function of x
double x2tof(double x, double s, double c, int lw, int N);

void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2,
	double tf, const double* unith, int& flag, double mu, int way = 0, int N = 0, int branch = 0,
	int Maxiter = 200, double tol = 1.0e-11);

//多圈lambert遍历求解
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double GM);

/***********************************Lambert求解器**********************************/
//类名：lambert_solver
//描述：lambert求解器；
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
