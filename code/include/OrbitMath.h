#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_
#include <math.h>
#include <assert.h>
#include <vector>

/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ���춯��ѧ�����ʵ����
* ��ϵ��: ��� 1522620129@qq.com
* �ļ���: OrbitMath.h
* ���ݼ��������������������;����������ѧ����
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2022-09-14    ���     ����������ļ�
* 01b       2020-09-15    ���     ���ע�ͣ�������ࡢ������������Ӳ��־������㺯��
****************************************************************************/

//*******************************��������***********************************
//�ͷ�����B���ڴ�
template<class T> inline void V_Dele(T* B)
{
	delete[] B;
	B = NULL;
}

//������A��ֵ����B��������ά��ΪN
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=A[I_];
}

//������ֵ���θ�����ά����B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
}

//������ֵ���θ�����ά����B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
	B[3]=vx;
	B[4]=vy;
	B[5]=vz;
}

//��Nά�������෴����B[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=-A[I_];
}

//��Nά����A��B�ĸ���Ԫ����ӵõ�����C[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]+B[I_];	
}

//��Nά����B�е�ÿ��Ԫ�ؾ�������A�õ�����C[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A+B[I_];	
}

//��Nά����A��B�ĸ���Ԫ������õ�����C[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B[I_];	
}

//��Nά����B�е�ÿ��Ԫ�ؾ���ȥ��A�õ�����C[i]=B[i]-A
template<class T> inline void V_Minus(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=B[I_]-A;
}

//��Nά����A��B�ĸ���Ԫ����˵õ�����C[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]*B[I_];	
}

//��Nά����B�е�ÿ��Ԫ�ؾ�������A�õ�����C[i]=B[i]*A
template<class T> inline void V_Multi(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A*B[I_];
}

//��Nά����A��B�ĸ���Ԫ������õ�����C[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B[I_];	
}

//��Nά����A�е�ÿ��Ԫ�ؾ�������B�õ�����C[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B;
}

//��Nά����A��B���ڻ�
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) result+=A[I_]*B[I_];
	return result;
}

//����ά����A��B�����C=AXB,������V_Cross(A,A,B)��V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];
}

//��Nά����B�и���Ԫ�صľ���ֵ�õ�����A[i]=|B[i]|
template<class T> inline void V_Absol(T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			A[I_]=B[I_];
		else
			A[I_]=-B[I_];
	}
}

//��Nά����B��1-����
template<class T> inline T V_Norm1(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			result+=B[I_];
		else
			result-=B[I_];
	}
	return result;
}

//��Nά����B��2-����
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result=V_Dot(B,B,N);
	return sqrt(result);
}

//��Nά����B������-����
template<class T> inline T V_NormInf(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) 
	{
		if(B[I_]>=0) {if(B[I_]>result) result=B[I_];}
		else {if(-B[I_]>result) result=-B[I_];}
	}
	return result;
}

//��Nά����B�ĵ�λ����C
template<class T> inline void V_Vers(T* C, const T* B, int N)
{
	double a = V_Norm2(B, N);
	V_Divid(C, B, a, N);
}

//��Nά����B�е����Ԫ��
template<class T> inline T V_Max(const T* B, int N)
{
	T result = B[0];
	for (int I_ = 0; I_ < N; I_++) if (B[I_] > result) result = B[I_];
	return result;
}

//��Nά����B�е����Ԫ�ؼ���λ��index
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]>maximal) {index=I_;maximal=B[I_];}
	return maximal;
}

//��Nά����B�е���СԪ��
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]<result) result=B[I_];
	return result;
}

//��Nά����B�е���СԪ�ؼ���λ��index
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]<minimal) {index=I_;minimal=B[I_];}
	return minimal;
}

//Nά����B��Aÿ��Ԫ�ض����ʱ������
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for (int I_ = 0; I_ < N; I_++)	if (B[I_] != A[I_])return false;
	return true;
}

//��Nά����p�е�Ԫ�ذ��մ�С�������У�����p
template <class T> void PaiXu(T* p, int N)
{
	T tmp = 0;
	for (int i = 0; i < N; i++)
	{
		bool flag = false;
		for (int j = 0; j < N - 1 - i; j++)
		{
			if (p[j] > p[j + 1])
			{
				tmp = p[j];
				p[j] = p[j + 1];
				p[j + 1] = tmp;

				flag = true;
			}
		}
		if (!flag)
		{
			break;
		}
	}
}

//���ֲ��ң����ش�С������������S��С����e��λ�ã���һ��λ�ô���e��
template <typename T> inline int binSearch(const std::vector<T>& S, T const& e, int lo, int hi)
{
	while (lo < hi)
	{
		int mi = (lo + hi) >> 1;
		e < S[mi] ? hi = mi : lo = mi + 1;
	}
	return --lo;
}

//���ֲ��ң����ش�С������������S��С����e��λ�ã���һ��λ�ô���e��
template <typename T> inline int binSearch(T* S, T const& e, int lo, int hi)
{
	while (lo < hi)
	{
		int mi = (lo + hi) >> 1;
		e < S[mi] ? hi = mi : lo = mi + 1;
	}
	return --lo;
}
//*******************************�����������*********************************

//*******************************��������*************************************
//��������ֵA�ķ��ţ�������1��-1��0��
template<class T> inline int Sign(const T & A)
{
	if(A>0)
		return 1;
	else if(A<0)
		return -1;
	else
		return 0;
}

//����������x��y�Ľϴ�ֵ
template <class T> inline T Max (T x, T y) 
{ 
	return (x>y)?x:y;
}

//����������x��y�Ľ�Сֵ
template <class T> inline T Min (T x, T y) 
{
	return (x<y)?x:y;
}

//���ؽ�y�ķ��Ÿ���x��ֵ
template <class T> inline T CopySign (T x, T y)
{
	return (y < 0) ? ((x < 0) ? x : -x) : ((x > 0) ? x : -x);
}

//ȡ�����x������,����double��.
inline double ANINT(double x)
{
	double left = fmod(x, 1.0);
	if (fabs(left) < 0.5)
		return x - left;
	else if (left >= 0.5)
		return x - left + 1;
	else
		return x - left - 1;
}

//ȡ�����x������,����int��.
inline int NINT(double x)
{
	double left = ANINT(x);
	return (int)left;
}

//����(x,y)��Ӧ�����,[0,2*pi)
inline double NiceAngle(double x, double y)
{
	double temp=atan2(y, x);
	if(temp<0.0)
		temp+= 6.283185307179586476925286766559;
	return temp;
}

//���Ƕ�ֵalphaת����[0,2*PI)��
inline double NiceAngle(double alpha)
{
	double temp=fmod(alpha, 6.283185307179586476925286766559);
	if(temp<0.0)
		temp+= 6.283185307179586476925286766559;
	return temp;
}
//*******************************�����������************************************

//*******************************���������**************************************


//��N��M�о���A��ֵ����B
template<class T> inline void M_Copy(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[I_][J_];
}
template<class T> inline void M_Copy(T* B, T* A, int N, int M)
{
	for (int I_ = 0; I_ < N * M; I_++) B[I_] = A[I_];
}

//���Ÿ�ֵ���θ���3��3�еľ���B
template<class T> inline void M_Copy(T** B, T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
{
	B[0][0]=a11;
	B[0][1]=a12;
	B[0][2]=a13;
	B[1][0]=a21;
	B[1][1]=a22;
	B[1][2]=a23;
	B[2][0]=a31;
	B[2][1]=a32;
	B[2][2]=a33;
}
template<class T> inline void M_Copy(T* B, T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
{
	B[0] = a11;
	B[1] = a12;
	B[2] = a13;
	B[3] = a21;
	B[4] = a22;
	B[5] = a23;
	B[6] = a31;
	B[7] = a32;
	B[8] = a33;
}

//������A��ֵ�������N��M�о���B�����д洢��A��ǰM����ΪB�ĵ�һ��
template<class T> inline void M_Copy(T** B, T* A, int N, int M)
{
	int i=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[i++];
}


//����ת��B[i][j]=A[j][i]
template<class T> inline void M_Tranpose(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[J_][I_];
}
template<class T> inline void M_Tranpose(T* B, T* A, int N, int M)
{
	for (int I_ = 0; I_ < N; I_++) for (int J_ = 0; J_ < M; J_++) B[I_ * M + J_] = A[J_ * N + I_];
}

//C[i][j]=A[i][j]+B[i][j]
template<class T> inline void M_Add(T** C, T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A[I_][J_]+B[I_][J_];	
}
template<class T> inline void M_Add(T* C, T* A, T* B, int N, int M)
{
	for (int I_ = 0; I_ < N * M; I_++) C[I_] = A[I_] + B[I_];
}

//C[i][j]=B[i][j]+A
template<class T> inline void M_Add(T** C, T** B, T A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A+B[I_][J_];	
}
template<class T> inline void M_Add(T* C, T* B, T A, int N, int M)
{
	for (int I_ = 0; I_ < N * M; I_++)  C[I_] = A + B[I_];
}



//��N��K�о������K��M�о���B�õ�N��M�о���C[i][j]=A[i][k]*B[k][j]
template<class T> inline void M_Multi(T** C, T** A, T** B, int N, int M, int K)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		C[I_][J_]=0;
		for(int i=0;i<K;i++){ C[I_][J_]+=A[I_][i]*B[i][J_];	}
	}
}
template<class T> inline void M_Multi(T* C, T* A, T* B, int N, int M, int K)
{
	for (int I_ = 0; I_ < N; I_++) for (int J_ = 0; J_ < M; J_++)
	{
		int tmp = I_ * M + J_;
		C[tmp] = 0;
		for (int i = 0; i < K; i++) { C[tmp] += A[I_ * K + i] * B[i * M + J_]; }
	}
}

//��N��M�о���A����M������B�õ�N������C[i]=A[i][j]*B[j]
template<class T> inline void M_Multi(T* C, T** A, T* B, int N, int M)
{
	for(int I_=0;I_<N;I_++)
	{
		C[I_]=0;
		for(int J_=0;J_<M;J_++)	C[I_]+=A[I_][J_]*B[J_];
	}
}
template<class T> inline void M_Multi(T* C, T* A, T* B, int N, int M)
{
	for (int I_ = 0; I_ < N; I_++)
	{
		int tmp = I_ * M;
		C[I_] = 0;
		for (int J_ = 0; J_ < M; J_++)	C[I_] += A[tmp + J_] * B[J_];
	}
}

//C[i][j]=A[i][j]*B
template<class T> inline void M_Multi(T** C, T** A, T B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=B*A[I_][J_];
}
template<class T> inline void M_Multi(T* C, T* A, T B, int N, int M)
{
	for (int I_ = 0; I_ < N * M; I_++)  C[I_] = B * A[I_];
}

#endif