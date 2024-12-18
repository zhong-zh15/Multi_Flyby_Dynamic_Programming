#ifndef GTOC11_problem_h
#define GTOC11_problem_h
#include <vector>
//#include <math.h>
#include <string>
#include <cmath>
#include <functional>

#include "GTOC4_Problem.h"
#include "OrbitFun.h"
extern double data_GTOC11[83453][7];

////////////////////////////////////////////////////////////////
/*                the constant of the GTOC11 problem           */
////////////////////////////////////////////////////////////////

//GTOC11������ʵ�ʲ�����δ���й�һ��
constexpr double GTOC11_mu = 1.32712440018e11 * 1.0e9;														//̫����������(m^3/s^2)
constexpr double GTOC11_AU = 1.49597870691e8 * 1.0e3;														//m
constexpr double GTOC11_Gamma_ATD = 1.0e-4;																	//m/s^2
constexpr double GTOC11_alpha = 6.0e-9;																		//s^-1
constexpr double GTOC11_Day2Second = 86400.0;																//��תs
constexpr double GTOC11_Year2Day = 365.25;																	//��ת��
constexpr double GTOC11_MJD_Init = 95739.0;																	//��ʼʱ��MJD
constexpr double GTOC11_MJD_End = 103044.0;																	//����ʱ��MJD

//��һ����λ����������������һ��Ϊ1������ת�볤��Ϊ1������ת����Ϊ2pi�����й�һ��
const double GTOC11_LUnit = GTOC11_AU;																	//���ȹ�һ�����÷�����һ�����ȣ������٣�= ʵ�ʳ��ȣ�m�� / LUnit
const double GTOC11_TUnit = sqrt(GTOC11_LUnit * GTOC11_LUnit * GTOC11_LUnit / GTOC11_mu);				//ʱ���һ�����÷�����һ��ʱ�䣨�����٣�= ʵ��ʱ�䣨s�� / TUnit
const double GTOC11_MuUnit = GTOC11_LUnit / GTOC11_TUnit * GTOC11_LUnit / GTOC11_TUnit * GTOC11_LUnit;	//����������һ�����÷�����һ�������������ɿ�Ϊ1
const double GTOC11_VUnit = GTOC11_LUnit / GTOC11_TUnit;												//�ٶȹ�һ�����÷�����һ���ٶȣ�m/s��= ʵ���ٶȣ�m/s�� / VUnit
const double GTOC11_AUnit = GTOC11_VUnit / GTOC11_TUnit;												//���ٶȹ�һ�����÷�����һ�����ٶȣ�m/s��= ʵ�ʼ��ٶȣ�m/s�� / AUnit

//�����������
constexpr int GTOC11_Num_Mother = 10;                                                                       //ĸ������
constexpr int GTOC11_Num_Station = 12;                                                                      //����վ����
constexpr int GTOC11_Num_Asteroid = 83453;                                                                  //С���Ǹ���
constexpr double GTOC11_Phase_Station = 30.0 * 0.017453292519943295769236907684886;                         //����վ��λ��

//��һ������
constexpr double GTOC11_Gms_Unit = 1.0;                                                                     //��һ��̫����������
const double GTOC11_ATD_Unit = GTOC11_Gamma_ATD / GTOC11_AUnit;			     							//��һ��С���Ǽ��ٶ�
const double GTOC11_Ini_Unit = GTOC11_MJD_Init * GTOC11_Day2Second / GTOC11_TUnit;						//��һ����ʼʱ��
const double GTOC11_End_Unit = GTOC11_MJD_End * GTOC11_Day2Second / GTOC11_TUnit;						//��һ������ʱ��

//��һ��Լ��
constexpr int GTOC11_Num_Impulse = 4;                                                                       //����������
const double GTOC11_ATD_Active = 30.0 * GTOC11_Day2Second / GTOC11_TUnit;							    //��һ������ʱ
const double GTOC11_Dyson_Dur = 90.0 * GTOC11_Day2Second / GTOC11_TUnit;							    //��һ������
const double GTOC11_Esc_Vel = 6.0e3 / GTOC11_VUnit;                                                     //��һ����v����
const double GTOC11_Rel_Dis = 1.0e3 / GTOC11_LUnit;                                                     //��һ��Խ��Ծ���
const double GTOC11_Rel_Vel = 2.0e3 / GTOC11_VUnit;                                                     //��һ��Խ����ٶ�
constexpr double GTOC11_Dyson_min = 0.65;                                                                   //��һ��С����뾶
constexpr double GTOC11_Mot_min = 0.4;                                                                      //��һĸ��̫����С����
const double GTOC11_Err_pos = 10.0e3 / GTOC11_LUnit;                                                    //��һλ�����
const double GTOC11_Err_vel = 0.01 / GTOC11_VUnit;                                                      //��һ�ٶ����
constexpr double GTOC11_Err_mas = 1.0;                                                                      //����������һ����


constexpr double step_t_GTOC11 = 32.0 * 86400.0; // 0.5 * 86400.0; //0.2
constexpr int num_t_GTOC11 = 31;// 101; //301 //7 801

void read_data_GTOC11();
void read_data_GTOC11(std::string filename);

//���tʱ��С����id��rv��tΪ��һ����λ
void getrv_asteroid_GTOC11(double* rv, int id, double t);

//���tʱ��С����id��coe��tΪ��һ����λ��coeΪ��һ������������
void getcoe_asteroid_GTOC11(double* coe, int id, double t);

void getrv_Earth_GTOC11(double* rv, double t);

//���-2��Ӧ���������ӦС����
//t����ĵ�λ��s����0��ʼ�����������ʼʱ�� 95739.0 (MJD)
//coe �����ĵ�λ��m
void getcoe_GTOC11(double* coe, int id, double t);

//���-2��Ӧ���������ӦС����
//t����ĵ�λ��s����0��ʼ�����������ʼʱ�� 95739.0 (MJD)
//rv �����ĵ�λ��m,m/s
void getrv_GTOC11(double* rv, int id, double t);




void GTOC11_sequence_process();
// ���ڴ洢 GTOC11 ���������׶����к����ʱ������
// ����ת��Ϊs
//void GTOC11_sequence(std::vector<std::vector<int>>& Sequence_GTOC11_multi_temp, std::vector<std::vector<double>> & T_GTOC11_multi_temp, std::vector<std::vector<double>>& Maximum_T_GTOC11_temp);
const std::vector<std::vector<int>>& GTOC11_sequence_get_sequence();
const std::vector<std::vector<double>>& GTOC11_sequence_get_T();
const std::vector<std::vector<double>>& GTOC11_sequence_get_Maximum_T();




inline int compute_v_nominal_GTOC11(int ID_departue, int ID_arrival, double time_departue, double time_arrival, double* rv_departue, double* rv_arrival)
{

	getrv_GTOC11(rv_departue, ID_departue, time_departue);
	getrv_GTOC11(rv_arrival, ID_arrival, time_arrival);

	lambert_solver lambertsolver(rv_departue, rv_arrival, time_arrival - time_departue, GTOC11_mu);
	lambertsolver.solve_multi();
	//lambertsolver.solve_single();
	double dv_total_pulse = lambertsolver.get_Mdv1() + lambertsolver.get_Mdv2();
	double dv0[3], dvt[3];
	lambertsolver.get_dv1(dv0); lambertsolver.get_dv2(dvt);
	V_Add(rv_departue + 3, rv_departue + 3, dv0, 3);
	V_Minus(rv_arrival + 3, rv_arrival + 3, dvt, 3);

	return lambertsolver.m_flag;
}




/****************************************************************************
* ������   : rendezvous2flyby()
* ��  ��   : �ӽ���������Ϊ��Խ����
* �� ��    : rv��С���ǵ�λ���ٶȣ�6ά
*            last_dv_end������С��ǰ����һ���ƶ������٣�ԭʼdv��3ά
*            next_dv_start������ȥ��һ��С���Ǽ��٣��ƶ���ԭʼdv��3ά
* ȫ�ֱ��� : ��
* ��    �� : last_dv_end���������Խ�����ġ�����С��ǰ����һ���ƶ������٣�dv��3ά
*            next_dv_start���������Խ�����ġ�����ȥ��һ��С���Ǽ��٣��ƶ���dv��3ά
****************************************************************************/
void rendezvous2flyby(const double* rv, double* last_dv_end, double* next_dv_start);

// ��ת����ϵ����F1F2Cת��XOYƽ�棬F1F2��X�᷽�򣬲�����ת�ƾ���
// r1r2r3�ֱ�Ϊ���������꣬TMΪת�ƾ���(ƽ���ԭ)
void calc_transF2X(double* r1, double* r2, double* r3, double** TM);
// �ж�F1F2�߶��Ƿ��е���Բ��(��)
bool verifyFinC(const double* rC, double r, double c);
// ��볤��a�Ĵη��̲�ȡ�ڶ����a(a����>a����>c)
double calc_eqnA(double m, double n, double r, double c);
// ����Բ��Բ�����Ĵη���x����
double calc_eqnX(double m, double n, double r, double c, double a);

#endif

