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

//GTOC11参数，实际参数，未进行归一化
constexpr double GTOC11_mu = 1.32712440018e11 * 1.0e9;														//太阳引力常数(m^3/s^2)
constexpr double GTOC11_AU = 1.49597870691e8 * 1.0e3;														//m
constexpr double GTOC11_Gamma_ATD = 1.0e-4;																	//m/s^2
constexpr double GTOC11_alpha = 6.0e-9;																		//s^-1
constexpr double GTOC11_Day2Second = 86400.0;																//天转s
constexpr double GTOC11_Year2Day = 365.25;																	//年转天
constexpr double GTOC11_MJD_Init = 95739.0;																	//初始时刻MJD
constexpr double GTOC11_MJD_End = 103044.0;																	//结束时刻MJD

//归一化单位，按照引力常数归一化为1，地球公转半长轴为1，地球公转周期为2pi，进行归一化
const double GTOC11_LUnit = GTOC11_AU;																	//长度归一化，用法：归一化长度（无量纲）= 实际长度（m） / LUnit
const double GTOC11_TUnit = sqrt(GTOC11_LUnit * GTOC11_LUnit * GTOC11_LUnit / GTOC11_mu);				//时间归一化，用法：归一化时间（无量纲）= 实际时间（s） / TUnit
const double GTOC11_MuUnit = GTOC11_LUnit / GTOC11_TUnit * GTOC11_LUnit / GTOC11_TUnit * GTOC11_LUnit;	//引力常数归一化，用法：归一化的引力常数可看为1
const double GTOC11_VUnit = GTOC11_LUnit / GTOC11_TUnit;												//速度归一化，用法：归一化速度（m/s）= 实际速度（m/s） / VUnit
const double GTOC11_AUnit = GTOC11_VUnit / GTOC11_TUnit;												//加速度归一化，用法：归一化加速度（m/s）= 实际加速度（m/s） / AUnit

//任务参数常数
constexpr int GTOC11_Num_Mother = 10;                                                                       //母船数量
constexpr int GTOC11_Num_Station = 12;                                                                      //发电站数量
constexpr int GTOC11_Num_Asteroid = 83453;                                                                  //小行星个数
constexpr double GTOC11_Phase_Station = 30.0 * 0.017453292519943295769236907684886;                         //发电站相位差

//归一化常数
constexpr double GTOC11_Gms_Unit = 1.0;                                                                     //归一化太阳引力常数
const double GTOC11_ATD_Unit = GTOC11_Gamma_ATD / GTOC11_AUnit;			     							//归一化小行星加速度
const double GTOC11_Ini_Unit = GTOC11_MJD_Init * GTOC11_Day2Second / GTOC11_TUnit;						//归一化初始时刻
const double GTOC11_End_Unit = GTOC11_MJD_End * GTOC11_Day2Second / GTOC11_TUnit;						//归一化结束时刻

//归一化约束
constexpr int GTOC11_Num_Impulse = 4;                                                                       //最大脉冲个数
const double GTOC11_ATD_Active = 30.0 * GTOC11_Day2Second / GTOC11_TUnit;							    //归一激活用时
const double GTOC11_Dyson_Dur = 90.0 * GTOC11_Day2Second / GTOC11_TUnit;							    //归一建造间隔
const double GTOC11_Esc_Vel = 6.0e3 / GTOC11_VUnit;                                                     //归一出发v无穷
const double GTOC11_Rel_Dis = 1.0e3 / GTOC11_LUnit;                                                     //归一飞越相对距离
const double GTOC11_Rel_Vel = 2.0e3 / GTOC11_VUnit;                                                     //归一飞越相对速度
constexpr double GTOC11_Dyson_min = 0.65;                                                                   //归一最小建造半径
constexpr double GTOC11_Mot_min = 0.4;                                                                      //归一母船太阳最小距离
const double GTOC11_Err_pos = 10.0e3 / GTOC11_LUnit;                                                    //归一位置误差
const double GTOC11_Err_vel = 0.01 / GTOC11_VUnit;                                                      //归一速度误差
constexpr double GTOC11_Err_mas = 1.0;                                                                      //质量误差（待归一化）


constexpr double step_t_GTOC11 = 32.0 * 86400.0; // 0.5 * 86400.0; //0.2
constexpr int num_t_GTOC11 = 31;// 101; //301 //7 801

void read_data_GTOC11();
void read_data_GTOC11(std::string filename);

//获得t时刻小行星id的rv，t为归一化单位
void getrv_asteroid_GTOC11(double* rv, int id, double t);

//获得t时刻小行星id的coe，t为归一化单位，coe为归一化经典轨道根数
void getcoe_asteroid_GTOC11(double* coe, int id, double t);

void getrv_Earth_GTOC11(double* rv, double t);

//序号-2对应地球，其余对应小行星
//t传入的单位是s，从0开始，代表任务初始时刻 95739.0 (MJD)
//coe 传出的单位是m
void getcoe_GTOC11(double* coe, int id, double t);

//序号-2对应地球，其余对应小行星
//t传入的单位是s，从0开始，代表任务初始时刻 95739.0 (MJD)
//rv 传出的单位是m,m/s
void getrv_GTOC11(double* rv, int id, double t);




void GTOC11_sequence_process();
// 用于存储 GTOC11 问题的任务阶段序列和最大时间限制
// 最后均转化为s
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
* 函数名   : rendezvous2flyby()
* 功  能   : 从交会条件变为飞越条件
* 输 入    : rv：小行星的位置速度，6维
*            last_dv_end：到达小行前的上一次制动（加速）原始dv，3维
*            next_dv_start：出发去下一颗小行星加速（制动）原始dv，3维
* 全局变量 : 无
* 输    出 : last_dv_end：【满足飞越条件的】到达小行前的上一次制动（加速）dv，3维
*            next_dv_start：【满足飞越条件的】出发去下一颗小行星加速（制动）dv，3维
****************************************************************************/
void rendezvous2flyby(const double* rv, double* last_dv_end, double* next_dv_start);

// 旋转坐标系，将F1F2C转到XOY平面，F1F2沿X轴方向，并计算转移矩阵
// r1r2r3分别为三个点坐标，TM为转移矩阵(平面→原)
void calc_transF2X(double* r1, double* r2, double* r3, double** TM);
// 判断F1F2线段是否有点在圆内(上)
bool verifyFinC(const double* rC, double r, double c);
// 解半长轴a四次方程并取第二大的a(a外切>a内切>c)
double calc_eqnA(double m, double n, double r, double c);
// 解椭圆与圆外切四次方程x坐标
double calc_eqnX(double m, double n, double r, double c, double a);

#endif

