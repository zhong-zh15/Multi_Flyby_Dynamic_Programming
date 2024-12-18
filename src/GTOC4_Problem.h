#ifndef GTOC4_Problem_h
#define GTOC4_Problem_h

#include <string>
#include <vector>
#include "Constant.h"
#include "Lambert.h"
#include "OrbitMath.h"

/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: GTCC4_problem.h
* Description: the constant of the GTOC4 problem
*              the dynamics model, the constraints, the nominal trajectory and the user-defined functions
* Log:
*Version      Date        Author           Description
* 01        2023-03-02    Zhong Zhang       Create
* 01        2024-12-11    Zhong Zhang       Combine the different files to one file
****************************************************************************/

////////////////////////////////////////////////////////////////
/*                the constant of the GTOC4 problem           */
////////////////////////////////////////////////////////////////

//Unit: AU for distance, second for time, rad for angle
const double Initial_MJD_GTOC4 = 57023.0;		//the initial MJD of the GTOC4 problem, Day
const double End_MJD_GTOC4 = 61041.0;			//the end MJD of the GTOC4 problem, Day
const double g0_GTOC4 = 9.80665;				//the gravity of the Earth, m/s^2
const double AU_GTOC4 = 1.49597870691e8;		//the distance of the Earth and the Sun, km
const double Year_GTOC4 = 365.25;               //the year of the GTOC4 problem, Day

const double Mu_GTOC4 = 1.32712440018e11 * 1.0e9; //the gravitational parameter of the Sun, m^3/s^2
const double Earth_coe_GTOC4[6] = { 0.999988049532578 * AU_GTOC4 * 1000.0, 1.671681163160e-2, 0.8854353079654e-3 * D2R, 175.40647696473 * D2R, 287.61577546182 * D2R, 4.4635844709062136 };   //257.60683707535 * D2R
const double Earth_coe_MJD_GTOC4 = 54000.0;       //the initial MJD of the GTOC4 problem, Day

const double step_t = 32 * 86400.0; // 0.5 * 86400.0; //0.2
const int num_t = 31;// 101; //301 //7

extern double asteroid_data[1436][6]; //the data of the asteroids
extern  std::string asteroid_name[1436]; //the name of the asteroids

////////////////////////////////////////////////////////////////
/*                             Read data                      */
////////////////////////////////////////////////////////////////
void read_asteroid_data_GTOC4();
void getcoe_asteroid(double* coe, int id, double t);
void getrv_asteroid(double* rv, int id, double t);
void getcoe_Earth(double* coe, int id, double t);
void getrv_Earth(double* rv, int id, double t);
void getcoe(double* coe, int id, double t);
void getrv(double* rv, int id, double t);

using namespace std;

void GTOC4_sequence_process();

inline int compute_v_nominal(int ID_departue, int ID_arrival, double time_departue, double time_arrival, double* rv_departue, double* rv_arrival)
{

    getrv(rv_departue, ID_departue, time_departue);
    getrv(rv_arrival, ID_arrival, time_arrival);

    lambert_solver lambertsolver(rv_departue, rv_arrival, time_arrival - time_departue, Mu_GTOC4);
    lambertsolver.solve_multi();
    //lambertsolver.solve_single();
    double dv_total_pulse = lambertsolver.get_Mdv1() + lambertsolver.get_Mdv2();
    double dv0[3], dvt[3];
    lambertsolver.get_dv1(dv0); lambertsolver.get_dv2(dvt);
    V_Add(rv_departue + 3, rv_departue + 3, dv0, 3);
    V_Minus(rv_arrival + 3, rv_arrival + 3, dvt, 3);

    return lambertsolver.m_flag;
}

const std::vector<int>& GTOC4_get_result_sequence();

const std::vector<double>& GTOC4_get_result_T();

#endif