#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_
#include <math.h>
#include <assert.h>
#include <vector>

/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院航天动力学与控制实验室
* 联系人: 武迪 1522620129@qq.com
* 文件名: OrbitMath.h
* 内容简述：整理针对数、数组和矩阵操作的数学函数
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2022-09-14    武迪     创建整理该文件
* 01b       2020-09-15    武迪     添加注释，整理分类、输入和输出，添加部分矩阵运算函数
****************************************************************************/

//*******************************向量运算***********************************
//释放向量B的内存
template<class T> inline void V_Dele(T* B)
{
	delete[] B;
	B = NULL;
}

//将向量A的值赋给B，向量的维度为N
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=A[I_];
}

//将三个值依次赋给三维向量B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
}

//将六个值依次赋给六维向量B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
	B[3]=vx;
	B[4]=vy;
	B[5]=vz;
}

//求N维向量的相反向量B[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=-A[I_];
}

//将N维向量A和B的各个元素相加得到向量C[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]+B[I_];	
}

//将N维向量B中的每个元素均加上数A得到向量C[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A+B[I_];	
}

//将N维向量A和B的各个元素相减得到向量C[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B[I_];	
}

//将N维向量B中的每个元素均减去数A得到向量C[i]=B[i]-A
template<class T> inline void V_Minus(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=B[I_]-A;
}

//将N维向量A和B的各个元素相乘得到向量C[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]*B[I_];	
}

//将N维向量B中的每个元素均乘以数A得到向量C[i]=B[i]*A
template<class T> inline void V_Multi(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A*B[I_];
}

//将N维向量A和B的各个元素相除得到向量C[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B[I_];	
}

//将N维向量A中的每个元素均除以数B得到向量C[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B;
}

//求N维向量A和B的内积
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) result+=A[I_]*B[I_];
	return result;
}

//求三维向量A和B的外积C=AXB,不能用V_Cross(A,A,B)或V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];
}

//求N维向量B中各个元素的绝对值得到向量A[i]=|B[i]|
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

//求N维向量B的1-范数
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

//求N维向量B的2-范数
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result=V_Dot(B,B,N);
	return sqrt(result);
}

//求N维向量B的无穷-范数
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

//求N维向量B的单位向量C
template<class T> inline void V_Vers(T* C, const T* B, int N)
{
	double a = V_Norm2(B, N);
	V_Divid(C, B, a, N);
}

//求N维向量B中的最大元素
template<class T> inline T V_Max(const T* B, int N)
{
	T result = B[0];
	for (int I_ = 0; I_ < N; I_++) if (B[I_] > result) result = B[I_];
	return result;
}

//求N维向量B中的最大元素及其位置index
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]>maximal) {index=I_;maximal=B[I_];}
	return maximal;
}

//求N维向量B中的最小元素
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]<result) result=B[I_];
	return result;
}

//求N维向量B中的最小元素及其位置index
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]<minimal) {index=I_;minimal=B[I_];}
	return minimal;
}

//N维向量B与A每个元素都相等时返回真
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for (int I_ = 0; I_ < N; I_++)	if (B[I_] != A[I_])return false;
	return true;
}

//将N维向量p中的元素按照从小到大排列，返回p
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

//二分查找，返回从小到大排列数组S中小于数e的位置，下一个位置大于e。
template <typename T> inline int binSearch(const std::vector<T>& S, T const& e, int lo, int hi)
{
	while (lo < hi)
	{
		int mi = (lo + hi) >> 1;
		e < S[mi] ? hi = mi : lo = mi + 1;
	}
	return --lo;
}

//二分查找，返回从小到大排列数组S中小于数e的位置，下一个位置大于e。
template <typename T> inline int binSearch(T* S, T const& e, int lo, int hi)
{
	while (lo < hi)
	{
		int mi = (lo + hi) >> 1;
		e < S[mi] ? hi = mi : lo = mi + 1;
	}
	return --lo;
}
//*******************************向量运算结束*********************************

//*******************************数的运算*************************************
//返回输入值A的符号（整数：1，-1，0）
template<class T> inline int Sign(const T & A)
{
	if(A>0)
		return 1;
	else if(A<0)
		return -1;
	else
		return 0;
}

//返回两个数x和y的较大值
template <class T> inline T Max (T x, T y) 
{ 
	return (x>y)?x:y;
}

//返回两个数x和y的较小值
template <class T> inline T Min (T x, T y) 
{
	return (x<y)?x:y;
}

//返回将y的符号赋予x的值
template <class T> inline T CopySign (T x, T y)
{
	return (y < 0) ? ((x < 0) ? x : -x) : ((x > 0) ? x : -x);
}

//取最近接x的整数,返回double型.
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

//取最近接x的整数,返回int型.
inline int NINT(double x)
{
	double left = ANINT(x);
	return (int)left;
}

//坐标(x,y)对应的相角,[0,2*pi)
inline double NiceAngle(double x, double y)
{
	double temp=atan2(y, x);
	if(temp<0.0)
		temp+= 6.283185307179586476925286766559;
	return temp;
}

//将角度值alpha转换到[0,2*PI)中
inline double NiceAngle(double alpha)
{
	double temp=fmod(alpha, 6.283185307179586476925286766559);
	if(temp<0.0)
		temp+= 6.283185307179586476925286766559;
	return temp;
}
//*******************************数的运算结束************************************

//*******************************矩阵的运算**************************************


//将N行M列矩阵A的值赋给B
template<class T> inline void M_Copy(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[I_][J_];
}
template<class T> inline void M_Copy(T* B, T* A, int N, int M)
{
	for (int I_ = 0; I_ < N * M; I_++) B[I_] = A[I_];
}

//将九个值依次赋给3行3列的矩阵B
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

//将数组A的值逐个赋给N行M列矩阵B，按行存储，A的前M个数为B的第一行
template<class T> inline void M_Copy(T** B, T* A, int N, int M)
{
	int i=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[i++];
}


//矩阵转置B[i][j]=A[j][i]
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



//将N行K列矩阵乘以K行M列矩阵B得到N行M列矩阵C[i][j]=A[i][k]*B[k][j]
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

//将N行M列矩阵A乘以M行向量B得到N行向量C[i]=A[i][j]*B[j]
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