#include "GTOC4_Problem.h"

#include "OrbitMath.h"
#include <fstream>
#include <iostream>
#include "Constant.h"
#include "GTOC4_Problem.h"

#include <algorithm>

#include "OrbitFun.h"

//Global variables
double asteroid_data[1436][6]; //the data of the asteroids
std::string asteroid_name[1436]; //the name of the asteroids

std::vector<int> Sequence_GTOC4;
std::vector<double> T_GTOC4;

/****************************************************************************
* Function     : read_asteroid_data_GTOC4
* Description  : read the data of the asteroids
****************************************************************************/
void read_asteroid_data_GTOC4()
{
	std::ifstream fin;
	fin.open("../input/GTOC4_data.txt");

	std::string str;
	int flag;
	char ch[1000];
	fin.getline(ch, 1000);
	fin.getline(ch, 1000);
	//fin >> str >> str >> str >> str >> str >> str >> str >> str;
	//0:a(AU); 1:e; 2:i(rad); 3:RAAN(rad); 4:omega(rad); 5:M(rad)
	for (int i = 0; i < 1436; i++)
	{
		fin >> asteroid_name[i];
		double time_temp;
		fin >> time_temp >> asteroid_data[i][0]
			>> asteroid_data[i][1] >> asteroid_data[i][2]
			>> asteroid_data[i][3] >> asteroid_data[i][4]
			>> asteroid_data[i][5];
		asteroid_data[i][0] *= AU_GTOC4 * 1000.0; //m
		asteroid_data[i][2] *= D2R;
		asteroid_data[i][3] *= D2R;
		asteroid_data[i][4] *= D2R;
		asteroid_data[i][5] *= D2R;

		asteroid_data[i][5] = M2E(flag, asteroid_data[i][5], asteroid_data[i][1]);
		asteroid_data[i][5] = E2f(flag, asteroid_data[i][5], asteroid_data[i][1]);

		if (fabs(time_temp - 54800.0) > 1.0e-5)
		{
			coe02coef(flag, asteroid_data[i], asteroid_data[i], (54800.0 - time_temp) * 86400.0, Mu_GTOC4);
			std::cout << " Time is not initial " << std::endl;
		}
	}



	//Earth_coe_GTOC4[5] = M2E(flag, Earth_coe_GTOC4[5], Earth_coe_GTOC4[1]);
	//Earth_coe_GTOC4[5] = E2f(flag, Earth_coe_GTOC4[5], Earth_coe_GTOC4[1]);

	fin.close();
}

/****************************************************************************
* Function     : getcoe_asteroid
* Description  : get the coe of the asteroid
*                input:
*					coe: the coe of the asteroid (AU and rad)
*					id: Asteroid ID (from 0)
*                   t: time (Second from initial date)
*                ouput:
*					coe: the coe of the asteroid (AU and rad)
****************************************************************************/
void getcoe_asteroid(double* coe, int id, double t)
{
	if (id < 0 || id >= 1436)
		std::cout << "ERROR::id is wrong!!! (getcoe_asteroid)" << std::endl;

	int flag;
	for (int i = 0; i < 6; i++)
		coe[i] = asteroid_data[id][i];
	coe[5] = f0dt2ft(flag, coe[5], t + (Initial_MJD_GTOC4 - 54800.0) * 86400.0, coe[0], coe[1], Mu_GTOC4);
}
/****************************************************************************
* Function     : getrv_asteroid
* Description  : get the rv of the asteroid
*                input:
*					rv: the rv of the asteroid (AU and AU/s)
*					id: Asteroid ID (from 0)
*                   t: time (Second from initial date)
*                ouput:
*					rv: the rv of the asteroid (AU and AU/s)
****************************************************************************/
void getrv_asteroid(double* rv, int id, double t)
{
	if (id < 0 || id >= 1436)
		std::cout << "ERROR::id is wrong!!! (getrv_asteroid)" << std::endl;

	double coe[6];
	int flag;
	getcoe_asteroid(coe, id, t);
	coe2rv(flag, rv, coe, Mu_GTOC4);
}

/****************************************************************************
* Function     : getcoe_Earth
* Description  : get the coe of the Earth
*                input:
*					coe: the coe of the Earth (AU and rad)
*					id: Earth ID (-3)
*                   t: time (Second from initial date)
*                ouput:
*					coe: the coe of the asteroid (AU and rad)
****************************************************************************/
void getcoe_Earth(double* coe, int id, double t)
{
	if (id != -3)
		std::cout << "ERROR::id is wrong!!! (getcoe_Earth)" << std::endl;

	int flag;
	for (int i = 0; i < 6; i++)
		coe[i] = Earth_coe_GTOC4[i];
	coe[5] = f0dt2ft(flag, coe[5], t + (Initial_MJD_GTOC4 - Earth_coe_MJD_GTOC4) * 86400.0, coe[0], coe[1], Mu_GTOC4);
}
/****************************************************************************
* Function     : getrv_asteroid
* Description  : get the rv of the Earth
*                input:
*					rv: the rv of the Earth (AU and AU/s)
*					id: Earth ID ( -3)
*                   t: time (Second from initial date)
*                ouput:
*					rv: the rv of the Earth (AU and AU/s)
****************************************************************************/
void getrv_Earth(double* rv, int id, double t)
{
	if (id != -3)
		std::cout << "ERROR::id is wrong!!! (getrv_Earth)" << std::endl;

	double coe[6];
	int flag;
	getcoe_Earth(coe, id, t);
	coe2rv(flag, rv, coe, Mu_GTOC4);
}


void getcoe(double* coe, int id, double t)
{
	if (id == -3)
		getcoe_Earth(coe, id, t);
	else
	{
		getcoe_asteroid(coe, id, t);
	}
}

void getrv(double* rv, int id, double t)
{
	if (id == -3)
		getrv_Earth(rv, id, t);
	else
	{
		getrv_asteroid(rv, id, t);
	}
}



void chimpionship_result(vector<int>& Sequence_gtoc4, vector<double>& T_sequence_gtoc4)
{
	Sequence_gtoc4 = { 0,1058,1126,1329,1182,938,1078,1290,1001,602,1377,1254,955,1091,1023,753,152,959,1338,1366,801,1357,1246,
		1330,667,927,736,1084,235,816,988,133,834,223,940,1294,710,1239,1380,249,1119,476,652,971,744,1101,556,1386,1174,367 };


	T_sequence_gtoc4 = { 58676.4,58801.45,58973.95,59119.45,59326.95,59419.55,59586.25,59746.75,59901.25,60078.75,60233.95,60422.75,60549.75,
		60665.25,60897.85,61054.35,61110.35,61221.35,61377.65,61486.75,61646.75,61806.65,61921.85,62057.15,62151.65,58731.65,58866.95,
		59084.55,59221.75,59366.25,59520.05,59685.15,59863.25,59955.85,60157.05,60334.75,60507.25,60601.15,60796.75,60938.15,61082.15,61201.65,
		61285.55,61441.75,61562.75,61727.45,61867.75,61979.55,62105.85,62261.25,
	};

	sort(T_sequence_gtoc4.begin(), T_sequence_gtoc4.end());

	for (int i = 0; i < Sequence_gtoc4.size(); i++)
	{
		Sequence_gtoc4[i] = Sequence_gtoc4[i] - 1;
	}
	Sequence_gtoc4[0] = -3;

	for (int i = 0; i < T_sequence_gtoc4.size(); i++)
	{
		T_sequence_gtoc4[i] = (T_sequence_gtoc4[i] - Initial_MJD_GTOC4) * 86400.0;
	}

	//vector<int> Sequence_gtoc4_filter={9,15,30,46};
	//for (auto it = Sequence_gtoc4_filter.rbegin(); it != Sequence_gtoc4_filter.rend(); ++it) {
	//	if (*it < Sequence_gtoc4.size()) {
	//		Sequence_gtoc4.erase(Sequence_gtoc4.begin() + *it);
	//		T_sequence_gtoc4.erase(T_sequence_gtoc4.begin() + *it);
	//	}
	//}


}


//Jena GTOC4 low-thrust result
void jena_result(vector<int>& Sequence_gtoc4, vector<double>& T_sequence_gtoc4)
{
	Sequence_gtoc4 = { -3,612,443,931,475,642,871,927,1312,1352,318,824,633,1377,857,1367,509,821,1019,
686,1253,1435,768,1407,125,461,653,958,1336,1279,307,1362,458,631,379,333,1415,1100,321,964,1419,1156,914,120,781,691,892,
92,1323,1420,800};

	T_sequence_gtoc4 = {
	57938,
	58014.076,58097.728,58119.077,58170.188,58257.831,58355.073,58448.887,58592.5,58657.82,58693.581,58822.871,58844.361,58924.084,
58988.288,59095.175,59152.196,59203.169,59262.899,59359.805,59436.238,59527.581,59563.955,59661.652,59750.498,59804.801,59923.251,
59964.86,60075.508,60114.119,60187.866,60252.911,60275.778,60362.533,60422.974,60474.506,60526.301,60582.821,60657.606,60731.297,60806.039,
60846.344,60909.121,60984.884,61048.13,61112.945,61201.947,61285.019,61308.591,61380.251,61554.717
	};

	for (int i = 0; i < T_sequence_gtoc4.size(); i++)
	{
		T_sequence_gtoc4[i] = (T_sequence_gtoc4[i] - Initial_MJD_GTOC4) * 86400.0;
	}

}


void GTOC4_sequence_process()
{
	vector<int> Sequence_gtoc4_temp;
	vector<double> T_sequence_gtoc4_temp;
	chimpionship_result(Sequence_gtoc4_temp, T_sequence_gtoc4_temp);

	//jena_result(Sequence_gtoc4, T_sequence_gtoc4);
	//auto [epochs, ids] = myresult("DP_result.txt");
	//T_sequence_gtoc4 = epochs;
	//Sequence_gtoc4 = ids;

	Sequence_GTOC4 = Sequence_gtoc4_temp;
	T_GTOC4 = T_sequence_gtoc4_temp;
}


const std::vector<int>& GTOC4_get_result_sequence() {
	return Sequence_GTOC4;
}

const std::vector<double>& GTOC4_get_result_T() {
	return T_GTOC4;
}

