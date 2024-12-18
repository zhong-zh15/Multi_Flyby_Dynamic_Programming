#ifndef DYNAMIC_PROGRAMMING_FLYBY_2IMPULSE_H
#define DYNAMIC_PROGRAMMING_FLYBY_2IMPULSE_H

#include <complex.h>
#include <vector>

using namespace std;

class state_flyby_lambert
{
public:
    double t;
    double t_last;
    double rv[6];
    double total_dv;
    state_flyby_lambert* pre_state;

    // Constructors
    state_flyby_lambert() : t(-1.0e16), t_last(-1.0e16), total_dv(0.0), pre_state(nullptr)
    {
        std::fill(rv, rv + 6, 1.0e16);
    }

    state_flyby_lambert(const state_flyby_lambert& other)
        : t(other.t), t_last(other.t_last), total_dv(other.total_dv), pre_state(other.pre_state)
    {
        std::memcpy(rv, other.rv, 6 * sizeof(double));
    }

    // Assignment Operator
    state_flyby_lambert& operator=(const state_flyby_lambert& state)
    {
        if (this != &state) // Prevent self-assignment
        {
            t = state.t;
            t_last = state.t_last;
            std::memcpy(rv, state.rv, 6 * sizeof(double));
            total_dv = state.total_dv;
            pre_state = state.pre_state;
        }
        return *this;
    }
};

class output_unit_flyby_lambert
{
public:
    double t;
    double rv[6];
    double dv;
    double total_dv;
    int id;

    // Constructors
    output_unit_flyby_lambert() : t(-1.0e16), dv(0.0), total_dv(0.0), id(-1)
    {
        std::fill(rv, rv + 6, 1.0e16);
    }

    output_unit_flyby_lambert(double t_a, double* rv_a, double dv_a, double total_dv_a, int id_a)
        : t(t_a), dv(dv_a), total_dv(total_dv_a), id(id_a)
    {
        std::memcpy(rv, rv_a, 6 * sizeof(double));
    }

    // Assignment Operator
    output_unit_flyby_lambert& operator=(const output_unit_flyby_lambert& unit)
    {
        if (this != &unit) // Prevent self-assignment
        {
            t = unit.t;
            std::memcpy(rv, unit.rv, 6 * sizeof(double));
            dv = unit.dv;
            total_dv = unit.total_dv;
            id = unit.id;
        }
        return *this;
    }
};
int DP_flyby_Lambert(
    const vector<int>& sequence,
	//任务阶段序列
	double& total_dv,
	//返回总dv
	vector<output_unit_flyby_lambert>& result //返回结果，该结构体在头文件中定义
	, double max_dv_threshold = 1.0e10
);

//为了通用性考虑，方便后续扩展，只考虑二脉冲Lambert解情况，状态量只需要当前时刻和前一个时刻，不并行
int DP_flyby_Lambert_iter(
	const vector<int>& sequence,
	//任务阶段序列
	double& total_dv,
	//返回总dv
	vector<output_unit_flyby_lambert>& result,
    bool if_display = false,
	string output_file = "",
	//返回结果，该结构体在头文件中定义
	double max_dv_threshold = 1.0e9,
	double step_t = 10.0 * 86400.0, double t_range = 100.0 * 86400.0, bool if_first = false
);

#endif // !DYNAMIC_PROGRAMMING_FLYBY_2IMPULSE_H


