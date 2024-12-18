#ifndef DP_flyby_2impulse_user_function_h
#define DP_flyby_2impulse_user_function_h
#include "DP_flyby_2impulse.h"
#include "GTOC4_Problem.h"
#include <vector>
#include <ostream>

using namespace std;
extern int problem_case;

// 用户自定义函数
// 返回可能的t值，由问题的特性决定，用户根据自己的问题定义
vector<double> possible_t_values(const vector<int>& sequence, //任务阶段序列
    int stage, std::vector<state_flyby_lambert>& last_layer_state //需要的参数，供扩展用
);
// 用户自定义函数
// 计算当前时刻的rv,dv,以及总dv，由问题的特性决定，用户根据自己的问题定义
void compute_rv_dv(vector<int> sequence,
    state_flyby_lambert& last_state,
    int id_departue, int id_now,
    double state_t,
    double* rv_arrival,
    double& dv,
    double& total_dv
);

// 用户自定义函数
// 检查 state_t 与 last_state_t 是否能够形成一个满足要求的解，可以剪枝一些冗余情况
// 也可以完全不做限制
bool check_t(double last_state_t, double state_t, int ID_departue, int ID_now);

// 用户自定义函数
// 输出最终结果，可以根据自己的问题定义输出内容
void print_optimal_path(const std::vector<output_unit_flyby_lambert>& optimal_path, 
std::ostream& out_stream);


#endif

