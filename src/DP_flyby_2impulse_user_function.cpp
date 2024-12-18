#include "DP_flyby_2impulse_user_function.h"

#include <iomanip>
#include <iostream>
#include <numeric>    //accumulate

#include "GTOC11_problem.h"
#include "OrbitFun.h"
#include "OrbitMath.h"


int problem_case = 4;  // 4: GTOC4, 11: GTOC11

// 用户自定义函数
// 返回可能的t值，由问题的特性决定，用户根据自己的问题定义
vector<double> possible_t_values(const vector<int>& sequence, //任务阶段序列
    int stage, std::vector<state_flyby_lambert>& last_layer_state //需要的参数，供扩展用
) {
    if (problem_case == 4)
    {
        const vector<double>& T_sequence = GTOC4_get_result_T();
        double t_min = T_sequence.front();
        double t_max = T_sequence.back();
		double t_mid = (t_max + t_min) / 2.0;
        t_min = t_mid - 5.0 * Year_GTOC4 * 86400.0;
		t_max = t_mid + 5.0 * Year_GTOC4 * 86400.0;

		double t_step = (t_max - t_min) / (T_sequence.size() - 1);
		double t_center = t_min + (stage +0.5) * t_step;

        //double t_center = T_sequence[stage]; //测试是否原序列计算结果与以往结果的值时一样的，验证正确性，已验证
        //cout << "t_center difference (day): " << fabs(t_center - T_sequence[stage])/86400.0 << endl;
        std::vector<double> possible_t_values;


        if (last_layer_state.size() > 0)
        {
            t_min = last_layer_state[0].t;
			t_max = min(t_max,last_layer_state.back().t + Year_GTOC4 * 86400.0);
        }


        int num_points_left = num_t / 2;
        int num_points_right = num_t - num_points_left;

        for (int i = -num_points_left; i < num_points_right; ++i) {
            double x = t_center + i * step_t;// std::min(step_t, arrival_time_feasible[1] - arrival_time_feasible[0]);
            if (t_min <= x && x <= t_max) {
                possible_t_values.push_back(x);
            }
        }

        return possible_t_values;
	}
	else if (problem_case == 11)
	{

        const auto& sequence_all = GTOC11_sequence_get_sequence();
        const auto& T_sequence_all = GTOC11_sequence_get_T();
        const auto& T_MAX_sequence_all = GTOC11_sequence_get_Maximum_T();


        //判断sequence在sequence_all的第几个里，给出index
		int index = -1;
		for (int i = 0; i < sequence_all.size(); i++)
		{
			if (sequence_all[i] == sequence)
			{
				index = i;
				break;
			}
		}
		const auto & T_sequence = T_sequence_all[index];
        const auto& T_MAX_sequence = T_MAX_sequence_all[index];


		//double t_center = T_sequence[stage]; //测试是否原序列计算结果与以往结果的值时一样的，验证正确性，已验证


		//t_center 应该时从0到最后时刻均分的
        double t_step_center = (T_MAX_sequence.back() - 0.0) / (sequence_all[index].size() - 1);
        double t_center =   (stage+0.5) * t_step_center;

		//cout << "t_center difference (day): " << (t_center - T_sequence[stage])/86400.0 << endl;

		double t_min = 0.0; 
        double t_max = T_MAX_sequence[stage];
		if (last_layer_state.size() > 0)
		{
			t_min = last_layer_state[0].t;
            t_max = min(t_max, last_layer_state.back().t + GTOC11_Year2Day * 2.*86400.0);
		}

		std::vector<double> possible_t_values;
		int num_points_left = num_t_GTOC11 / 2;
		int num_points_right = num_t_GTOC11 - num_points_left;
		for (int i = -num_points_left; i < num_points_right; ++i) {
			double x = t_center + i * step_t_GTOC11;
			if (t_min <= x && x <= t_max) {
				possible_t_values.push_back(x);
			}
		}

		return possible_t_values;
	}

    return  std::vector<double> {};
}

// 用户自定义函数
// 计算当前时刻的rv,dv,以及总dv，由问题的特性决定，用户根据自己的问题定义
void compute_rv_dv(vector<int> sequence,
    state_flyby_lambert& last_state,
    int id_departue, int id_now,
    double state_t,
    double* rv_arrival,
    double& dv,
    double& total_dv
) {
    if (problem_case == 4)
    {
        double rv_departue[6];
        compute_v_nominal(id_departue, id_now, last_state.t, state_t, rv_departue, rv_arrival);
        double dv_vector[3];

        V_Minus(dv_vector, rv_departue + 3, last_state.rv + 3, 3);
        dv = V_Norm2(dv_vector, 3);

        //GTOC4: first give C3 maximum 4.0 km/s
        if (sequence.front() == id_departue)
        {
            getrv(last_state.rv, id_departue, last_state.t);
        	V_Minus(dv_vector, rv_departue + 3, last_state.rv + 3, 3);
            dv = V_Norm2(dv_vector, 3);
            dv = dv - 4000.0 < 0.0 ? 0.0 : dv - 4000.0;
        }

        //GTOC4: last target is rendezvous
        if (sequence.back() == id_now)
        {
            double rv_arrival_target[6];
            getrv(rv_arrival_target, id_now, state_t);
            double dv_vector_rendezvous[3];
            V_Minus(dv_vector_rendezvous, rv_arrival_target + 3, rv_arrival + 3, 3);
            dv += V_Norm2(dv_vector_rendezvous, 3);
        }

        total_dv = last_state.total_dv + dv;
        
    }
	else if (problem_case == 11)
    {
        double rv_departue[6];
        int flag = compute_v_nominal_GTOC11(id_departue, id_now, last_state.t, state_t, rv_departue, rv_arrival);
        double dv_vector[3];

        double coe[6];
		rv2coe(flag, coe, rv_departue, GTOC11_mu);

        if (flag != 1)
        {
            dv = 1.e10;
            total_dv = 1.0e10;
            return;
        }

        double rv_body[6];
		getrv_GTOC11(rv_body, id_departue, last_state.t);

		//last_state.rv 本质上就是上一个状态求完Lambert后的rv_arrival
        double last_dv_end[3];
		double next_dv_start[3];
		V_Minus(last_dv_end, rv_body + 3, last_state.rv + 3, 3);
		V_Minus(next_dv_start, rv_departue + 3, rv_body + 3, 3);

        rendezvous2flyby(rv_body, last_dv_end, next_dv_start); //飞越2.0km/s 由该函数控制

        double dv_test = V_Norm2(last_dv_end, 3);

        dv = V_Norm2(last_dv_end, 3) + V_Norm2(next_dv_start, 3);

        //GTOC11: first give C3 maximum 6.0 km/s
        if (sequence.front() == id_departue)
        {
            getrv_GTOC11(last_state.rv, id_departue, last_state.t);
            V_Minus(dv_vector, rv_departue + 3, last_state.rv + 3, 3);
            dv = V_Norm2(dv_vector, 3);

            dv = dv - 6000.0 < 0.0 ? 0.0 : dv - 6000.0;
        }

		//GTOC11: last target flyby speed cannot exceed 2.0 km/s
        if (sequence.back() == id_now)
        {
            double rv_arrival_target[6];
            getrv_GTOC11(rv_arrival_target, id_now, state_t);
            double dv_vector_rendezvous[3];
            V_Minus(dv_vector_rendezvous, rv_arrival_target + 3, rv_arrival + 3, 3);

            double dv_rendezvous = V_Norm2(dv_vector_rendezvous, 3);
            dv_rendezvous  = dv_rendezvous -2000.0 < 0.0 ? 0.0 : dv_rendezvous-2000.0;

            dv += dv_rendezvous;
        }

        total_dv = last_state.total_dv + dv;
    }


}


// 用户自定义函数
// 检查 state_t 与 last_state_t 是否能够形成一个满足要求的解，可以剪枝一些冗余情况
// 也可以完全不做限制
bool check_t(double last_state_t, double state_t, int ID_departue, int ID_now) {

    if (problem_case == 4)
    {
        if (last_state_t > -1.e6)
        {
            //如果state_t 与 last_state_t 相差小于1天，返回false
            if (state_t - last_state_t < 1.0 * 86400.) {
                return false;
            }
            if ((state_t - last_state_t) > 365.25 * 86400.0 * 1.0) {
                //如果state_t 与 last_state_t 相差大于1年，返回false
                return false;
            }
        }

        const vector<double>& T_sequence_gtoc4 = GTOC4_get_result_T();
        const vector<double>& T_sequence = GTOC4_get_result_T();
        double t_min = T_sequence.front();
        double t_max = T_sequence.back();
        double t_mid = (t_max + t_min) / 2.0;
        t_min = t_mid - 5.0 * GTOC11_Year2Day * 86400.0;
        t_max = t_mid + 5.0 * GTOC11_Year2Day * 86400.0;

        if (state_t  > t_max) {
            return false;
        }
        if (state_t < t_min) {
            return false;
        }

        //double rv_departure[6];
        //double rv_arrival[6];
        //double coe[6];

        //compute_v_nominal(ID_departue, ID_now, last_state_t, state_t, rv_departure, rv_arrival);
        //int flag;
        //rv2coe(flag, coe, rv_departure, Mu_GTOC4);

        ////如果rv_departure的coe不满足要求(半长轴大于1.8 AU 或 偏心率大于0.3)，返回false
        //if (0.5 * AU_GTOC4 * 1000.0 > coe[0] ||
        //    coe[0] > 1.6 * AU_GTOC4 * 1000.0 ||
        //    coe[1] > 0.3)
        //{
        //    return false;
        //}

        return true;
    }
    else if (problem_case == 11)
    {
        if (state_t <0.0)
        {
            return false;
        }
        if (last_state_t > -1.e6)
        {
            //如果state_t 与 last_state_t 相差小于2天，返回false
            if (state_t - last_state_t < 1.0 * 86400.) {
                return false;
            }
            if ((state_t - last_state_t) > 365.25 * 86400.0 * 2.0) {
                //如果state_t 与 last_state_t 相差大于2年，返回false
                return false;
            }
        }
        const auto& sequence_all = GTOC11_sequence_get_sequence();
        //const auto& T_sequence_all = GTOC11_sequence_get_T();
        const auto& T_MAX_sequence_all = GTOC11_sequence_get_Maximum_T();

        //判断ID_now在sequence_all的第几个里（查找，可能在某条序列的中间），给出index
        int index_i = -1;
        int index_j = -1;
		for (int i = 0; i < sequence_all.size(); i++) //判断ID_now是否在sequence_all[i]中，如果在，返回i
		{
			for (int j = 0; j < sequence_all[i].size(); j++)
			{
				if (sequence_all[i][j] == ID_now)
				{
					index_i = i;
					index_j = j;
					break;
				}
			}
		}
        //const auto& T_sequence = T_sequence_all[index_i];
        const auto& T_MAX_sequence = T_MAX_sequence_all[index_i];

		if (state_t > T_MAX_sequence[index_j])
		{
			return false;
		}

        return true;
    }

	return false;
}


// 用户自定义函数
// 输出最终结果，可以根据自己的问题定义输出内容
void print_optimal_path(const std::vector<output_unit_flyby_lambert>& optimal_path, std::ostream& out_stream) {

    if (problem_case == 4)
    {
        // 输出 dv
        double total_dv = std::accumulate(optimal_path.begin(), optimal_path.end(), 0.0,
            [](double sum, const output_unit_flyby_lambert& unit) { return sum + unit.dv; });
        out_stream << "Total DV: " << std::fixed << std::setprecision(3)
            << total_dv << "   Mass: "<< 1500.0 * exp(-total_dv/g0_GTOC4/3000.0)<< std::endl;

        // 输出标题
        out_stream << std::setw(6) << "Stage" << std::setw(19) << "Epoch"
            << std::setw(19) << "X(m)" << std::setw(19) << "Y(m)" << std::setw(19) << "Z(m)"
            << std::setw(13) << "VX(m/s)" << std::setw(13) << "VY(m/s)" << std::setw(13) << "VZ(m/s)"
            << std::setw(11) << "DV(m/s)" << std::setw(11) << "TOTAL_DV" << std::setw(11) << "MASS(kg)"
            << std::setw(6) << "ID" << std::endl;

        // 输出 optimal_path 中的数据
        int stage = 0;
        for (const auto& unit : optimal_path) {
            out_stream << std::setw(6) << stage++;
            out_stream << std::setw(19) << std::fixed << std::setprecision(2) << unit.t; // Epoch
            out_stream << std::setw(19) << std::fixed << std::setprecision(1) << unit.rv[0]; // X
            out_stream << std::setw(19) << std::fixed << std::setprecision(1) << unit.rv[1]; // Y
            out_stream << std::setw(19) << std::fixed << std::setprecision(1) << unit.rv[2]; // Z
            out_stream << std::setw(13) << std::fixed << std::setprecision(3) << unit.rv[3]; // VX
            out_stream << std::setw(13) << std::fixed << std::setprecision(3) << unit.rv[4]; // VY
            out_stream << std::setw(13) << std::fixed << std::setprecision(3) << unit.rv[5]; // VZ
            out_stream << std::setw(11) << std::fixed << std::setprecision(2) << unit.dv; // DV
            out_stream << std::setw(11) << std::fixed << std::setprecision(2) << unit.total_dv; // total_DV
            out_stream << std::setw(11) << std::fixed << std::setprecision(2) << 1500.0 * exp(-unit.total_dv / g0_GTOC4 / 3000.0); // total_DV
            out_stream << std::setw(6) << std::noshowpoint << unit.id << std::endl; // ID
        }
    }
    else if (problem_case == 11)
    {
        // 输出 dv
        out_stream << "Total DV: " << std::fixed << std::setprecision(3)
            << std::accumulate(optimal_path.begin(), optimal_path.end(), 0.0,
                [](double sum, const output_unit_flyby_lambert& unit) { return sum + unit.dv; }) << std::endl;

        // 输出标题
        out_stream << std::setw(6) << "Stage" << std::setw(19) << "Epoch"
            << std::setw(19) << "X(m)" << std::setw(19) << "Y(m)" << std::setw(19) << "Z(m)"
            << std::setw(13) << "VX(m/s)" << std::setw(13) << "VY(m/s)" << std::setw(13) << "VZ(m/s)"
            << std::setw(11) << "DV(m/s)" << std::setw(11) << "TOTAL_DV"
            << std::setw(6) << "ID" << std::endl;

        // 输出 optimal_path 中的数据
        int stage = 0;
        for (const auto& unit : optimal_path) {
            out_stream << std::setw(6) << stage++;
            out_stream << std::setw(19) << std::fixed << std::setprecision(2) << unit.t; // Epoch
            out_stream << std::setw(19) << std::fixed << std::setprecision(1) << unit.rv[0]; // X
            out_stream << std::setw(19) << std::fixed << std::setprecision(1) << unit.rv[1]; // Y
            out_stream << std::setw(19) << std::fixed << std::setprecision(1) << unit.rv[2]; // Z
            out_stream << std::setw(13) << std::fixed << std::setprecision(3) << unit.rv[3]; // VX
            out_stream << std::setw(13) << std::fixed << std::setprecision(3) << unit.rv[4]; // VY
            out_stream << std::setw(13) << std::fixed << std::setprecision(3) << unit.rv[5]; // VZ
            out_stream << std::setw(11) << std::fixed << std::setprecision(2) << unit.dv; // DV
            out_stream << std::setw(11) << std::fixed << std::setprecision(2) << unit.total_dv; // total_DV
            out_stream << std::setw(6) << std::noshowpoint << unit.id << std::endl; // ID
        }
    }
}

