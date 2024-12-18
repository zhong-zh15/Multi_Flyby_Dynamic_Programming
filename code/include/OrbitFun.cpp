#include "OrbitFun.h"

/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院航天动力学与控制实验室
* 联系人: 武迪 1522620129@qq.com
* 文件名: OrbitFun.cpp
* 内容简述：整理轨道根数转换和解析递推函数
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2022-09-15    武迪     整理该文件，添加了平瞬根转换和递推
****************************************************************************/

//经典轨道根数 coe[5]:
// 	coe[0]半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
//                      其值取为近星距.
// 	coe[1]偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
// 	coe[2]轨道倾角i（弧度）:范围0<=i<=180度.
//	coe[3]升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,可将其值取为0
//	coe[4]近拱点幅角omega（弧度）:当偏心率为0时没有意义,可将其值置为0
//	coe[5]真近点角f（弧度）:当偏心率为0时没有意义,可将其值取为omega+f,即纬度幅角

/*********************************经典轨道根数三种角度关系***************************************************************/
// E2f 根据偏近点角和偏心率求真近点角
//输入参数E, e:
//   E:偏近点角,单位：弧度,对于双曲轨道,指H,有r=a(ecoshH-1)
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//输出参数f:
// 	真近点角（弧度）
double E2f(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }

	double f = 0.0;
	if (e >= 0.0 && e < 1.0)//圆和椭圆轨道
	{
		double E0 = fmod(E, D2PI);
		if (E0 > DPI)
			E0 -= D2PI;
		if (E0 < -DPI)
			E0 += D2PI;
		f = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(0.5 * E0));
		f = f + E - E0;
	}
	else if (e > 1.0)//双曲线轨道
		f = 2.0 * atan(sqrt((e + 1.0) / (e - 1.0)) * tanh(0.5 * E));
	else //抛物线轨道
	{
		f = E; //抛物线轨道没有定义偏近点角.在此将真近点角的值置为E
	}
	flag = 1;
	return f;
}

// E2M 根据偏近点角和偏心率求平近点角
//输入参数E, e:
//   E:偏近点角,单位：弧度,对于双曲轨道,指H
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//输出参数M:
// 	平近点角(弧度),对于双曲轨道,指N
double E2M(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }
	double M = 0.0;
	if (e >= 0.0 && e < 1.0)//圆和椭圆轨道
	{
		double E0 = fmod(E, D2PI);
		M = E0 - e * sin(E0);
		M = M + E - E0;
	}
	else if (e > 1.0)//双曲线轨道
		M = e * sinh(E) - E;
	else //抛物线轨道
	{
		M = E;
		//		cout<<"抛物线轨道没有定义偏近点角.在此将平近点角值置为E."<<endl;
	}
	flag = 1;
	return M;
}

// f2E 根据真近点角和偏心率求偏近点角
//输入参数f, e:
//   f:真近点角,单位：弧度
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//输出参数E:
// 	偏近点角（弧度）,对于双曲轨道,指H,有r=a(ecoshH-1)
double f2E(int& flag, double f, double e)
{
	if (e < 0.0) { flag = 0; return f; }

	double E = 0.0;
	if (e >= 0.0 && e < 1.0)//圆和椭圆轨道
	{
		double f0 = fmod(f, D2PI);
		if (f0 > DPI)
			f0 -= D2PI;
		if (f0 < -DPI)
			f0 += D2PI;
		E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(0.5 * f0));
		E += f - f0;
	}
	else if (e > 1.0)//双曲线轨道
	{
		if (f > DPI - acos(1.0 / e) || f < -DPI + acos(1.0 / e))
		{
			//			cout<<"不可能达到的双曲轨道."<<endl;
			flag = 0;
			return f;
		}
		else
			E = 2.0 * atanh(sqrt((e - 1.0) / (1.0 + e)) * tan(0.5 * f));
	}
	else//if(abs(e-1.0)<epsilon)抛物线轨道
	{
		E = f;
		//		cout<<"抛物线轨道没有定义偏近点角.在此将其值置为f."<<endl;
	}
	flag = 1;
	return E;
}

// M2E 根据平近点角和偏心率求偏近点角
//输入参数M, e, MaxIter, epsilon:
//   M:平近点角,单位：弧度,对于双曲轨道,指N
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   MaxIter:最大迭代次数,默认为60
//   epsilon:迭代绝对误差,默认为1e-14;
//输出参数E:
// 	偏近点角（弧度），对双曲轨道，指H
double M2E(int& flag, double M, double e, int MaxIter, double epsilon)
{
	if (epsilon <= 0.0 || MaxIter < 1 || e < 0.0) { flag = 0; return M; }

	//迭代方法见《Solar System Dynamics》Chapter2,Carl D.Murray and Stanley F.Dermott著
	double E = 0.0, Minus = 0.0, DeMinus = 0.0, DeDeMinus = 0.0, DeDeDeMinus = 0.0, Delta1 = 0.0, Delta2 = 0.0, Delta3 = 0.0;
	int N = 0;
	if (e >= 0.0 && e < 1.0)//圆和椭圆轨道
	{
		double RM = fmod(M, D2PI);
		if (RM < 0.0)
			RM += D2PI;
		double sinRM = sin(RM);
		E = RM + 0.85 * e * Sign(sinRM);
		N = 0;
		Delta3 = 1.0;
		while (fabs(Delta3) >= epsilon && N < MaxIter)
		{
			Minus = E - e * sin(E) - RM;
			DeMinus = 1.0 - e * cos(E);
			DeDeMinus = e * sin(E);
			DeDeDeMinus = e * cos(E);
			Delta1 = -Minus / DeMinus;
			Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
			Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
			E = E + Delta3;
			N = N + 1;
		}
		E = E + M - RM;
	}
	else if (e > 1.0)//双曲线轨道
	{
		E = asinh(M / e);
		Delta3 = 1.0;
		N = 0;
		while (fabs(Delta3) >= epsilon && N < MaxIter)
		{
			Minus = e * sinh(E) - E - M;
			DeMinus = e * cosh(E) - 1.0;
			DeDeMinus = e * sinh(E);
			DeDeDeMinus = e * cosh(E);
			Delta1 = -Minus / DeMinus;
			Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
			Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
			E = E + Delta3;
			N = N + 1;
		}
	}
	else //(abs(e-1.0)<epsilon)抛物线轨道
	{
		E = M;
		//		cout<<"抛物线轨道没有定义偏近点角.在此将其值置为M."<<endl;
	}
	if (((e >= 0.0 && e < 1.0) || (e > 1.0)) && fabs(Delta3) >= 5.0 * epsilon && N >= MaxIter)
	{
		//		cout<<"迭代不收敛,请降低精度epsilon或增加迭代次数限制."<<endl;
		flag = 0;
		return M;
	}
	flag = 1;
	return E;
}

// f0dt2ft 根据初始真近点角和演化时间求最终真近点角
//输入参数f0,t,a,e,mu,MaxIter,epsilon:
//   f0:初始真近点角,单位：弧度
//   t:轨道深化时间,单位：秒
//   a:轨道半长轴，单位：米。对于抛物线轨道，指近星距，即p/2
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   MaxIter:最大迭代次数
//   epsilon:迭代绝对误差,默认为1e-14;
//输出参数ft:
// 	最终真近点角(弧度)
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter, double epsilon)
{
	if (mu <= 0.0 || MaxIter < 1 || a <= 0.0 || e < 0.0) { flag = 0; return f0; }

	double ft = 0.0;
	if ((e >= 0.0 && e < 1.0) || (e > 1.0))//圆,椭圆,双曲轨道
	{
		double E = f2E(flag, f0, e);
		if (flag == 0) return f0;
		double M = E2M(flag, E, e);
		if (flag == 0) return f0;
		M += sqrt(mu / (a * a * a)) * dt;
		E = M2E(flag, M, e, MaxIter, epsilon);
		if (flag == 0) return f0;
		ft = E2f(flag, E, e);
		if (flag == 0) return f0;
	}
	else //(abs(e-1.0)<epsilon)抛物线轨道
	{
		if ((f0 < -DPI) || (f0 > DPI))
		{
			//			cout<<"对于抛物线轨道，初始真近点角应在-180至180度之间."<<endl;
			flag = 0; return f0;
		}
		else if (f0 > DPI || f0 < -DPI)
			ft = f0;
		else
		{
			double B = 0.75 * sqrt(2.0 * mu / (a * a * a)) * dt + 0.5 * tan(0.5 * f0) * ((tan(0.5 * f0)) * (tan(0.5 * f0)) + 3.0);
			double B1B = B + sqrt(1.0 + B * B);
			double tanv = 0.0;
			if (fabs(dt) < D2PI * sqrt((a * a * a) / mu) / 1000.0)//推进时间为小量的情况
			{
				double A = pow(B1B, 2.0 / 3.0);
				tanv = 2.0 * A * B / (1.0 + (1.0 + A) * A);
			}
			else//不是小量的情况
			{
				double temp = pow(B1B, 1.0 / 3.0);
				tanv = temp - 1.0 / temp;
			}
			ft = 2.0 * atan(tanv);
		}
	}
	flag = 1;
	return ft;
}

/*********************************经典轨道根数、直角坐标、改进春分点轨道根数的转换**************************************/
void coe2rv(int& flag, double* rv, const double* coe, double mu)
{
	flag = 0;
	if (mu <= 0.0 || coe[0] <= 0.0 || coe[1] < 0.0 || coe[2]<0.0 || coe[2]>DPI)
		return;
	if ((coe[1] * cos(coe[5])) < -1.0)
	{
		//		cout<<"不可能达到的双曲轨道."<<endl;		
		return;
	}

	double p = coe[0] * fabs(1.0 - coe[1] * coe[1]);//半通径
	if (coe[1] == 1.0)//如果是抛物线轨道,区别对待.
		p = 2.0 * coe[0];

	double sini, cosi, sinO, cosO, sino, coso;
	sini = sin(coe[2]);
	cosi = cos(coe[2]);
	sinO = sin(coe[3]);
	cosO = cos(coe[3]);
	sino = sin(coe[4]);
	coso = cos(coe[4]);

	//轨道平面法向单位矢量,即角动量单位矢量
	double HVector[3] = { sini * sinO, -sini * cosO, cosi };

	//偏心率单位矢量,或叫Laplace矢量
	double PVector[3] = { cosO * coso - sinO * sino * cosi, sinO * coso + cosO * sino * cosi, sino * sini };

	//半通径方向单位矢量,PVector,QVector,HVector构成右手坐标系
	// QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
	double QVector[3];
	V_Cross(QVector, HVector, PVector);

	double r = 0.0;
	if ((coe[1] * cos(coe[5])) + 1.0 <= 0.0)
	{
		//		cout<<"抛物或双曲轨道到达无穷远处."<<endl;
		r = 1.0e308;
	}
	else
		r = p / (1.0 + coe[1] * cos(coe[5]));

	for (int i = 0; i < 3; i++)
	{
		rv[i] = r * (cos(coe[5]) * PVector[i] + sin(coe[5]) * QVector[i]);
		rv[3 + i] = sqrt(mu / p) * (-sin(coe[5]) * PVector[i] + (cos(coe[5]) + coe[1]) * QVector[i]);
	}
	flag = 1;
	return;
}

void rv2coe(int& flag, double* coe, const double* RV, double mu)
{
	int i;
	flag = 0;
	if (mu <= 0.0)
		return;

	double R[3] = { RV[0], RV[1], RV[2] };
	double V[3] = { RV[3], RV[4], RV[5] };
	double radius = V_Norm2(R, 3);//距离
	double velocity = V_Norm2(V, 3);//速度
	if (radius <= 0.0 || velocity <= 0.0)
		return;
	double unitR[3];
	for (i = 0; i < 3; i++) unitR[i] = R[i] / radius;//径向单位矢量    
	double unitV[3];
	for (i = 0; i < 3; i++) unitV[i] = V[i] / velocity;//切向单位矢量
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h = radius * velocity * V_Norm2(hvector, 3);//角动量值
	if (h <= 0.0)
		return;
	double unith[3];
	for (i = 0; i < 3; i++) unith[i] = hvector[i] / V_Norm2(hvector, 3);//轨道面法向单位矢量
	//偏心率矢量
	double evector[3];
	V_Cross(evector, unitV, unith);
	for (i = 0; i < 3; i++) evector[i] = (velocity * h / mu) * evector[i] - unitR[i];
	coe[1] = V_Norm2(evector, 3);//偏心率
	double p = h * h / mu;
	if (coe[1] == 1.0)
		coe[0] = 0.5 * p;//抛物线轨道的近星距
	else
		coe[0] = p / (fabs(1.0 - coe[1] * coe[1]));//半长轴
	bool judge = (coe[1] > 0.0);
	double unite[3] = { 0.0 };
	if (judge)
		for (i = 0; i < 3; i++) unite[i] = evector[i] / coe[1];//偏心率单位矢量
	coe[2] = acos(unith[2]);//轨道倾角

	double unitN[3] = { -unith[1], unith[0], 0.0 };//节线矢量,未归一化

	double temp[3];

	if (V_Norm2(unitN, 3) == 0.0)
	{
		coe[3] = 0.0;//升交点赤经
//		cout<<"轨道倾角接近0或180度,升交点赤经容易奇异.在此将其置为零."<<endl;
		if (!judge)
		{
			coe[4] = 0.0;//近星点幅角
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;        
			coe[5] = atan2(unitR[1] * unith[2], unitR[0]);//真近点角
		}
		else
		{
			V_Cross(temp, unite, unitR);
			coe[4] = atan2(unite[1] * unith[2], unite[0]); //近星点幅角       
			coe[5] = atan2(V_Dot(unith, temp, 3), V_Dot(unite, unitR, 3));
		}
	}
	else
	{
		V_Cross(temp, unitN, unitR);
		coe[3] = atan2(unith[0], -unith[1]);
		coe[5] = atan2(V_Dot(unith, temp, 3), V_Dot(unitN, unitR, 3));
		if (!judge)
		{
			coe[4] = 0.0;
			//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;
		}
		else
		{
			V_Cross(temp, unitN, unite);
			coe[4] = atan2(V_Dot(unith, temp, 3), V_Dot(unite, unitN, 3));
			coe[5] = coe[5] - coe[4];
		}
	}
	//转换到[0,2pi)中
	coe[3] = fmod(coe[3], D2PI);
	if (coe[3] < 0.0)
		coe[3] += D2PI;
	coe[4] = fmod(coe[4], D2PI);
	if (coe[4] < 0.0)
		coe[4] += D2PI;
	coe[5] = fmod(coe[5], D2PI);
	if (coe[1] >= 1.0)
	{
		if (coe[5] > DPI - acos(1.0 / coe[1]))
			coe[5] -= D2PI;
		else if (coe[5] < -DPI + acos(1.0 / coe[1]))
			coe[5] += D2PI;
	}
	flag = 1;
	return;
}


/***********************************已知初值和时间，求末端状态***************************************/
//根据初始时刻状态coe0求末端时刻dt的状态coe1，按二体推进。若计算成功,flag返回1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu)
{
	V_Copy(coe1, coe0, 6);
	coe1[5] = f0dt2ft(flag, coe1[5], dt, coe1[0], coe1[1], mu);
}

//根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu)
{
	double coe[6];
	rv2coe(flag, coe, rv0, mu);
	if (flag == 0) return;
	coe[5] = f0dt2ft(flag, coe[5], dt, coe[0], coe[1], mu);
	if (flag == 0) return;
//	coe[5] = fmod(coe[5], D2PI);
	coe2rv(flag, rv1, coe, mu);
}
