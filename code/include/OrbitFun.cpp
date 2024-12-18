#include "OrbitFun.h"

/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ���춯��ѧ�����ʵ����
* ��ϵ��: ��� 1522620129@qq.com
* �ļ���: OrbitFun.cpp
* ���ݼ���������������ת���ͽ������ƺ���
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2022-09-15    ���     ������ļ��������ƽ˲��ת���͵���
****************************************************************************/

//���������� coe[5]:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,�ɽ���ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵȡΪomega+f,��γ�ȷ���

/*********************************�������������ֽǶȹ�ϵ***************************************************************/
// E2f ����ƫ����Ǻ�ƫ������������
//�������E, e:
//   E:ƫ�����,��λ������,����˫�����,ָH,��r=a(ecoshH-1)
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//�������f:
// 	�����ǣ����ȣ�
double E2f(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }

	double f = 0.0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
	{
		double E0 = fmod(E, D2PI);
		if (E0 > DPI)
			E0 -= D2PI;
		if (E0 < -DPI)
			E0 += D2PI;
		f = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(0.5 * E0));
		f = f + E - E0;
	}
	else if (e > 1.0)//˫���߹��
		f = 2.0 * atan(sqrt((e + 1.0) / (e - 1.0)) * tanh(0.5 * E));
	else //�����߹��
	{
		f = E; //�����߹��û�ж���ƫ�����.�ڴ˽������ǵ�ֵ��ΪE
	}
	flag = 1;
	return f;
}

// E2M ����ƫ����Ǻ�ƫ������ƽ�����
//�������E, e:
//   E:ƫ�����,��λ������,����˫�����,ָH
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//�������M:
// 	ƽ�����(����),����˫�����,ָN
double E2M(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }
	double M = 0.0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
	{
		double E0 = fmod(E, D2PI);
		M = E0 - e * sin(E0);
		M = M + E - E0;
	}
	else if (e > 1.0)//˫���߹��
		M = e * sinh(E) - E;
	else //�����߹��
	{
		M = E;
		//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽�ƽ�����ֵ��ΪE."<<endl;
	}
	flag = 1;
	return M;
}

// f2E ���������Ǻ�ƫ������ƫ�����
//�������f, e:
//   f:������,��λ������
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//�������E:
// 	ƫ����ǣ����ȣ�,����˫�����,ָH,��r=a(ecoshH-1)
double f2E(int& flag, double f, double e)
{
	if (e < 0.0) { flag = 0; return f; }

	double E = 0.0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
	{
		double f0 = fmod(f, D2PI);
		if (f0 > DPI)
			f0 -= D2PI;
		if (f0 < -DPI)
			f0 += D2PI;
		E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(0.5 * f0));
		E += f - f0;
	}
	else if (e > 1.0)//˫���߹��
	{
		if (f > DPI - acos(1.0 / e) || f < -DPI + acos(1.0 / e))
		{
			//			cout<<"�����ܴﵽ��˫�����."<<endl;
			flag = 0;
			return f;
		}
		else
			E = 2.0 * atanh(sqrt((e - 1.0) / (1.0 + e)) * tan(0.5 * f));
	}
	else//if(abs(e-1.0)<epsilon)�����߹��
	{
		E = f;
		//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽���ֵ��Ϊf."<<endl;
	}
	flag = 1;
	return E;
}

// M2E ����ƽ����Ǻ�ƫ������ƫ�����
//�������M, e, MaxIter, epsilon:
//   M:ƽ�����,��λ������,����˫�����,ָN
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   MaxIter:����������,Ĭ��Ϊ60
//   epsilon:�����������,Ĭ��Ϊ1e-14;
//�������E:
// 	ƫ����ǣ����ȣ�����˫�������ָH
double M2E(int& flag, double M, double e, int MaxIter, double epsilon)
{
	if (epsilon <= 0.0 || MaxIter < 1 || e < 0.0) { flag = 0; return M; }

	//������������Solar System Dynamics��Chapter2,Carl D.Murray and Stanley F.Dermott��
	double E = 0.0, Minus = 0.0, DeMinus = 0.0, DeDeMinus = 0.0, DeDeDeMinus = 0.0, Delta1 = 0.0, Delta2 = 0.0, Delta3 = 0.0;
	int N = 0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
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
	else if (e > 1.0)//˫���߹��
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
	else //(abs(e-1.0)<epsilon)�����߹��
	{
		E = M;
		//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽���ֵ��ΪM."<<endl;
	}
	if (((e >= 0.0 && e < 1.0) || (e > 1.0)) && fabs(Delta3) >= 5.0 * epsilon && N >= MaxIter)
	{
		//		cout<<"����������,�뽵�;���epsilon�����ӵ�����������."<<endl;
		flag = 0;
		return M;
	}
	flag = 1;
	return E;
}

// f0dt2ft ���ݳ�ʼ�����Ǻ��ݻ�ʱ��������������
//�������f0,t,a,e,mu,MaxIter,epsilon:
//   f0:��ʼ������,��λ������
//   t:����ʱ��,��λ����
//   a:����볤�ᣬ��λ���ס����������߹����ָ���Ǿ࣬��p/2
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   MaxIter:����������
//   epsilon:�����������,Ĭ��Ϊ1e-14;
//�������ft:
// 	����������(����)
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter, double epsilon)
{
	if (mu <= 0.0 || MaxIter < 1 || a <= 0.0 || e < 0.0) { flag = 0; return f0; }

	double ft = 0.0;
	if ((e >= 0.0 && e < 1.0) || (e > 1.0))//Բ,��Բ,˫�����
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
	else //(abs(e-1.0)<epsilon)�����߹��
	{
		if ((f0 < -DPI) || (f0 > DPI))
		{
			//			cout<<"���������߹������ʼ������Ӧ��-180��180��֮��."<<endl;
			flag = 0; return f0;
		}
		else if (f0 > DPI || f0 < -DPI)
			ft = f0;
		else
		{
			double B = 0.75 * sqrt(2.0 * mu / (a * a * a)) * dt + 0.5 * tan(0.5 * f0) * ((tan(0.5 * f0)) * (tan(0.5 * f0)) + 3.0);
			double B1B = B + sqrt(1.0 + B * B);
			double tanv = 0.0;
			if (fabs(dt) < D2PI * sqrt((a * a * a) / mu) / 1000.0)//�ƽ�ʱ��ΪС�������
			{
				double A = pow(B1B, 2.0 / 3.0);
				tanv = 2.0 * A * B / (1.0 + (1.0 + A) * A);
			}
			else//����С�������
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

/*********************************������������ֱ�����ꡢ�Ľ����ֵ���������ת��**************************************/
void coe2rv(int& flag, double* rv, const double* coe, double mu)
{
	flag = 0;
	if (mu <= 0.0 || coe[0] <= 0.0 || coe[1] < 0.0 || coe[2]<0.0 || coe[2]>DPI)
		return;
	if ((coe[1] * cos(coe[5])) < -1.0)
	{
		//		cout<<"�����ܴﵽ��˫�����."<<endl;		
		return;
	}

	double p = coe[0] * fabs(1.0 - coe[1] * coe[1]);//��ͨ��
	if (coe[1] == 1.0)//����������߹��,����Դ�.
		p = 2.0 * coe[0];

	double sini, cosi, sinO, cosO, sino, coso;
	sini = sin(coe[2]);
	cosi = cos(coe[2]);
	sinO = sin(coe[3]);
	cosO = cos(coe[3]);
	sino = sin(coe[4]);
	coso = cos(coe[4]);

	//���ƽ�淨��λʸ��,���Ƕ�����λʸ��
	double HVector[3] = { sini * sinO, -sini * cosO, cosi };

	//ƫ���ʵ�λʸ��,���Laplaceʸ��
	double PVector[3] = { cosO * coso - sinO * sino * cosi, sinO * coso + cosO * sino * cosi, sino * sini };

	//��ͨ������λʸ��,PVector,QVector,HVector������������ϵ
	// QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
	double QVector[3];
	V_Cross(QVector, HVector, PVector);

	double r = 0.0;
	if ((coe[1] * cos(coe[5])) + 1.0 <= 0.0)
	{
		//		cout<<"�����˫�������������Զ��."<<endl;
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
	double radius = V_Norm2(R, 3);//����
	double velocity = V_Norm2(V, 3);//�ٶ�
	if (radius <= 0.0 || velocity <= 0.0)
		return;
	double unitR[3];
	for (i = 0; i < 3; i++) unitR[i] = R[i] / radius;//����λʸ��    
	double unitV[3];
	for (i = 0; i < 3; i++) unitV[i] = V[i] / velocity;//����λʸ��
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h = radius * velocity * V_Norm2(hvector, 3);//�Ƕ���ֵ
	if (h <= 0.0)
		return;
	double unith[3];
	for (i = 0; i < 3; i++) unith[i] = hvector[i] / V_Norm2(hvector, 3);//����淨��λʸ��
	//ƫ����ʸ��
	double evector[3];
	V_Cross(evector, unitV, unith);
	for (i = 0; i < 3; i++) evector[i] = (velocity * h / mu) * evector[i] - unitR[i];
	coe[1] = V_Norm2(evector, 3);//ƫ����
	double p = h * h / mu;
	if (coe[1] == 1.0)
		coe[0] = 0.5 * p;//�����߹���Ľ��Ǿ�
	else
		coe[0] = p / (fabs(1.0 - coe[1] * coe[1]));//�볤��
	bool judge = (coe[1] > 0.0);
	double unite[3] = { 0.0 };
	if (judge)
		for (i = 0; i < 3; i++) unite[i] = evector[i] / coe[1];//ƫ���ʵ�λʸ��
	coe[2] = acos(unith[2]);//������

	double unitN[3] = { -unith[1], unith[0], 0.0 };//����ʸ��,δ��һ��

	double temp[3];

	if (V_Norm2(unitN, 3) == 0.0)
	{
		coe[3] = 0.0;//������ྭ
//		cout<<"�����ǽӽ�0��180��,������ྭ��������.�ڴ˽�����Ϊ��."<<endl;
		if (!judge)
		{
			coe[4] = 0.0;//���ǵ����
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;        
			coe[5] = atan2(unitR[1] * unith[2], unitR[0]);//������
		}
		else
		{
			V_Cross(temp, unite, unitR);
			coe[4] = atan2(unite[1] * unith[2], unite[0]); //���ǵ����       
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
			//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;
		}
		else
		{
			V_Cross(temp, unitN, unite);
			coe[4] = atan2(V_Dot(unith, temp, 3), V_Dot(unite, unitN, 3));
			coe[5] = coe[5] - coe[4];
		}
	}
	//ת����[0,2pi)��
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


/***********************************��֪��ֵ��ʱ�䣬��ĩ��״̬***************************************/
//���ݳ�ʼʱ��״̬coe0��ĩ��ʱ��dt��״̬coe1���������ƽ���������ɹ�,flag����1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu)
{
	V_Copy(coe1, coe0, 6);
	coe1[5] = f0dt2ft(flag, coe1[5], dt, coe1[0], coe1[1], mu);
}

//���ݳ�ʼʱ��t0��״̬rv0��ĩ��ʱ��t1��״̬rv1���������ƽ���������ɹ�,flag����1
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
