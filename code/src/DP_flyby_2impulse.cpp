#include "DP_flyby_2impulse.h"
#include "DP_flyby_2impulse_user_function.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>

#include "OrbitFun.h"

using namespace std;

//初始化状态矩阵
std::vector<std::vector<state_flyby_lambert>> dp_state_all_stages_flyby_lambert;


//为了通用性考虑，方便后续扩展，只考虑二脉冲Lambert解情况，状态量只需要当前时刻和前一个时刻，不并行
int DP_flyby_Lambert(
	const vector<int> &sequence,
	//任务阶段序列
	double & total_dv,
	//返回总dv
	vector<output_unit_flyby_lambert> &result //返回结果，该结构体在头文件中定义
	, double max_dv_threshold
){
    //double max_dv_threshold = 1.0e10;

    int Sequence_N = static_cast<int>(sequence.size());
    dp_state_all_stages_flyby_lambert.clear();
    dp_state_all_stages_flyby_lambert.resize(Sequence_N);

    // Initialize stage 0
    {
        int stage = 0;
        auto t_list = possible_t_values(sequence, stage, dp_state_all_stages_flyby_lambert[stage]);
        for(int i = 0; i < t_list.size(); i++)
        {
	        state_flyby_lambert temp;
			temp.t = t_list[i];
            //getrv(temp.rv, sequence[stage], temp.t);
			dp_state_all_stages_flyby_lambert[stage].push_back(temp);
        }
    }

    // Compute the solution for subproblems from stage 1 to stage N-1
    for (int stage = 1; stage < Sequence_N; ++stage) {
        auto start_t = std::chrono::steady_clock::now();
        std::cout << " stage:" << stage << endl; 

        auto t_values = possible_t_values(sequence, stage, dp_state_all_stages_flyby_lambert[stage-1]);
        // Compute last_stage_t
        auto& last_stage_states = dp_state_all_stages_flyby_lambert[stage - 1];

        int ID_departue = sequence[stage - 1];
        int ID_now = sequence[stage];

        // 计算每个状态，先对t进行遍历，再对上一阶段的状态进行遍历
		// 计划的结构是按照同一个t排序，然后前一阶段的t排序，这样可以减少计算量
        // 计算复杂度为O(n^3) 其中n为每层可能的时间数量（即一层的状态数为n^2）
        for (int i = 0; i < t_values.size(); i++)
        {
            double state_t = t_values[i];

            //接下来对上一阶段的t遍历，上一阶段同样的t存在多个解，其中选择dv最小的解，作为DP的递归排序
            //因此要分出同样的上一阶段t，然后对dv进行排序取最好的
            //last_stage_states 已经按照 last_state_t 排序
            double last_state_t_previous = -1.0e16;
			std::vector<state_flyby_lambert*> last_stage_states_same_t; //只放地址，不拷贝，减少内存占用

            for (int j = 0; j <= last_stage_states.size(); j++)
            {
                if (j == last_stage_states.size() || fabs(last_stage_states[j].t - last_state_t_previous) > 1.0e-5)
                {
					if (!last_stage_states_same_t.empty()) //对应初始时刻为空的情况，不需要处理
                    {
                        // 对于同一组 last_state_t，处理 t_values，找到 total_dv 最小的状态
                        state_flyby_lambert min_state;
                        double min_total_dv = max_dv_threshold;
                        min_state.total_dv = min_total_dv;

                        for (const auto last_state: last_stage_states_same_t)
                        {
							if (last_state->total_dv < min_total_dv) //剪枝，不影响最优性
                            {
                                double rv_temp[6];
                                double dv;
                                double total_dv_temp;
                                // 计算当前时刻的 rv, dv, 以及 total_dv
                                compute_rv_dv(sequence, (*last_state), ID_departue, ID_now, state_t, rv_temp, dv, total_dv_temp);

                                // 如果 total_dv 更小，更新最小状态
                                if (total_dv_temp < min_total_dv)
                                {
                                    min_total_dv = total_dv_temp;
                                    min_state.t = state_t;
                                    min_state.t_last = last_state->t;
                                    memcpy(min_state.rv, rv_temp, 6 * sizeof(double));
                                    min_state.total_dv = total_dv_temp;
                                    min_state.pre_state = last_state;
                                }
                            }
                        }

						if (min_state.total_dv < max_dv_threshold) //未小于预先设定的阈值，剪枝，不影响最优性
						{
                            dp_state_all_stages_flyby_lambert[stage].push_back(min_state);
						}
                        
                        last_stage_states_same_t.clear();
                    }
                    if (j == last_stage_states.size())
                        break;

                    last_state_t_previous = last_stage_states[j].t;
                }
                
                // 检查 state_t 与 last_state_t 是否能够形成一个满足要求的解，可以剪枝一些冗余情况
                if (check_t(last_stage_states[j].t, state_t, ID_departue, ID_now))
                {
                    last_stage_states_same_t.push_back(&last_stage_states[j]);
                }
            }


        }
        if (dp_state_all_stages_flyby_lambert[stage].size() == 0) //means no result exist!
        {
            break;
        }
        cout << " stage:" << stage << "  End!" << "   state size:" << dp_state_all_stages_flyby_lambert[stage].size() << endl;
        auto end_t = std::chrono::steady_clock::now();
        double duration_t = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count() / 1000.0;
        std::cout << "Time is " << duration_t << "s ; which is " << duration_t / 3600.0 << " hours" << std::endl;
    }

    // Find the last stage with non-empty states
    int last_non_empty_stage = Sequence_N-1;
    while (last_non_empty_stage >= 0 && dp_state_all_stages_flyby_lambert[last_non_empty_stage].empty()) {
        --last_non_empty_stage;
    }

    // If all stages are empty, no solution exists
    if (last_non_empty_stage < 0) {
        cout << "No solution exists!" << endl;
        return -1; // Return an empty path with a negative result
    }
    else if (last_non_empty_stage != Sequence_N - 1) {
        cout << " Solution Not Complete!" << endl;
		return 0; // Return an incomplete path with a zero result
    }

    // Extract the optimal path
    vector<output_unit_flyby_lambert> optimal_path(last_non_empty_stage + 1);

	// find the lowest total_dv in the last non-empty stage, which is dp_state_all_stages_flyby_lambert[last_non_empty_stage]
	 auto min_iter = min_element(
         dp_state_all_stages_flyby_lambert[last_non_empty_stage].begin(), 
         dp_state_all_stages_flyby_lambert[last_non_empty_stage].end(),
		[](const state_flyby_lambert& a, const state_flyby_lambert& b) 
         {return a.total_dv < b.total_dv; });

	 state_flyby_lambert min_state = *min_iter;

    // Find the state with the maximum metric in the last non-empty stage
    int max_index = 0;
	total_dv = min_state.total_dv;
    cout << "Total dv  " << min_state.total_dv << endl;
    ////Todo: delete,test
    //cout << "mass:  " << 1500.0*exp(-min_state.total_dv/3000.0/g0_GTOC4) << endl;


    // Backtrack the optimal solution using the saved indices
    for (int stage = last_non_empty_stage; stage >= 0; --stage) {
        //double rv_target[6];
        int id_target = sequence[stage];
        double t = min_state.t;
        //getrv(rv_target, id_target, t);
        if (stage != 0) {
            optimal_path[stage] = output_unit_flyby_lambert(t, min_state.rv, min_state.total_dv - min_state.pre_state->total_dv, min_state.total_dv, id_target);
            min_state = *min_state.pre_state;
		}
		else {
			optimal_path[stage] = output_unit_flyby_lambert(t, min_state.rv, min_state.total_dv, min_state.total_dv, id_target);
		}
    }

    result = optimal_path;

    print_optimal_path(result, cout);

    std::ofstream file_output("DP_result.txt", std::ios::app);
    print_optimal_path(result, file_output);
    file_output.close();
    
    return 1;
}

int DP_flyby_Lambert_iter(const vector<int>& sequence, double& total_dv, vector<output_unit_flyby_lambert>& result,
                          bool if_display, string output_file, double max_dv_threshold, double step_t, double t_range, bool if_first)
{
	int Sequence_N = static_cast<int>(sequence.size());
	dp_state_all_stages_flyby_lambert.clear();
	dp_state_all_stages_flyby_lambert.resize(Sequence_N);

	// Initialize stage 0
	{
		int stage = 0;
		vector<double> t_values;
		if (if_first)
		{
			t_values = possible_t_values(sequence, stage, dp_state_all_stages_flyby_lambert[stage]);
		}
		else
		{
			double t_center = result[stage].t;
			std::vector<double> possible_t_values;
			int num_points_left = static_cast<int>(t_range / step_t);
			int num_points_right = num_points_left;
			for (int i = -num_points_left; i <= num_points_right; ++i) {
				double x = t_center + i * step_t;
				if(check_t(-1e9, x, sequence[0], sequence[1]))
				t_values.push_back(x);
			}
		}

		for (int i = 0; i < t_values.size(); i++)
		{
			state_flyby_lambert temp;
			temp.t = t_values[i];
			dp_state_all_stages_flyby_lambert[stage].push_back(temp);
		}
	}

	// Compute the solution for subproblems from stage 1 to stage N-1
	for (int stage = 1; stage < Sequence_N; ++stage) {
		auto start_t = std::chrono::steady_clock::now();
		if(if_display) std::cout << " stage:" << stage << endl;

		vector<double> t_values;
		if (if_first)
		{
			t_values = possible_t_values(sequence, stage, dp_state_all_stages_flyby_lambert[stage - 1]);
		}
		else
		{
			double t_center = result[stage].t;
			std::vector<double> possible_t_values;
			int num_points_left = static_cast<int>(t_range / step_t);
			int num_points_right = num_points_left;
			for (int i = -num_points_left; i <= num_points_right; ++i) {
				double x = t_center + i * step_t;
				t_values.push_back(x);
			}
		}
		// Compute last_stage_t
		auto& last_stage_states = dp_state_all_stages_flyby_lambert[stage - 1];

		int ID_departue = sequence[stage - 1];
		int ID_now = sequence[stage];

		// 计算每个状态，先对t进行遍历，再对上一阶段的状态进行遍历
		// 计划的结构是按照同一个t排序，然后前一阶段的t排序，这样可以减少计算量
		// 计算复杂度为O(n^3) 其中n为每层可能的时间数量（即一层的状态数为n^2）
		for (int i = 0; i < t_values.size(); i++)
		{
			double state_t = t_values[i];

			//接下来对上一阶段的t遍历，上一阶段同样的t存在多个解，其中选择dv最小的解，作为DP的递归排序
			//因此要分出同样的上一阶段t，然后对dv进行排序取最好的
			//last_stage_states 已经按照 last_state_t 排序
			double last_state_t_previous = -1.0e16;
			std::vector<state_flyby_lambert*> last_stage_states_same_t; //只放地址，不拷贝，减少内存占用

			for (int j = 0; j <= last_stage_states.size(); j++)
			{
				if (j == last_stage_states.size() || fabs(last_stage_states[j].t - last_state_t_previous) > 1.0e-5)
				{
					if (!last_stage_states_same_t.empty()) //对应初始时刻为空的情况，不需要处理
					{
						// 对于同一组 last_state_t，处理 t_values，找到 total_dv 最小的状态
						state_flyby_lambert min_state;
						double min_total_dv = max_dv_threshold;
						min_state.total_dv = min_total_dv;

						for (const auto last_state : last_stage_states_same_t)
						{
							if (last_state->total_dv < min_total_dv) //剪枝，不影响最优性
							{
								double rv_temp[6];
								double dv;
								double total_dv_temp;
								// 计算当前时刻的 rv, dv, 以及 total_dv
								compute_rv_dv(sequence, (*last_state), ID_departue, ID_now, state_t, rv_temp, dv, total_dv_temp);

								// 如果 total_dv 更小，更新最小状态
								if (total_dv_temp < min_total_dv)
								{
									min_total_dv = total_dv_temp;
									min_state.t = state_t;
									min_state.t_last = last_state->t;
									memcpy(min_state.rv, rv_temp, 6 * sizeof(double));
									min_state.total_dv = total_dv_temp;
									min_state.pre_state = last_state;
								}
							}
						}

						if (min_state.total_dv < max_dv_threshold) //未小于预先设定的阈值，剪枝，不影响最优性
						{
							dp_state_all_stages_flyby_lambert[stage].push_back(min_state);
						}

						last_stage_states_same_t.clear();
					}
					if (j == last_stage_states.size())
						break;

					last_state_t_previous = last_stage_states[j].t;
				}

				// 检查 state_t 与 last_state_t 是否能够形成一个满足要求的解，可以剪枝一些冗余情况
				if (check_t(last_stage_states[j].t, state_t, ID_departue, ID_now))
				{
					last_stage_states_same_t.push_back(&last_stage_states[j]);
				}
			}


		}
		if (dp_state_all_stages_flyby_lambert[stage].size() == 0) //means no result exist!
		{
			break;
		}
		if (if_display) cout << " stage:" << stage << "  End!" << "   state size:" << dp_state_all_stages_flyby_lambert[stage].size() << endl;
		auto end_t = std::chrono::steady_clock::now();
		double duration_t = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count() / 1000.0;
		if (if_display) std::cout << "Time is " << duration_t << "s ; which is " << duration_t / 3600.0 << " hours" << std::endl;
	}

	// Find the last stage with non-empty states
	int last_non_empty_stage = Sequence_N - 1;
	while (last_non_empty_stage >= 0 && dp_state_all_stages_flyby_lambert[last_non_empty_stage].empty()) {
		--last_non_empty_stage;
	}

	// If all stages are empty, no solution exists
	if (last_non_empty_stage < 0) {
		cout << "No solution exists!" << endl;
		return -1; // Return an empty path with a negative result
	}
	else if (last_non_empty_stage != Sequence_N - 1) {
		cout << " Solution Not Complete!" << endl;
		return 0; // Return an incomplete path with a zero result
	}

	// Extract the optimal path
	vector<output_unit_flyby_lambert> optimal_path(last_non_empty_stage + 1);

	// find the lowest total_dv in the last non-empty stage, which is dp_state_all_stages_flyby_lambert[last_non_empty_stage]
	auto min_iter = min_element(
		dp_state_all_stages_flyby_lambert[last_non_empty_stage].begin(),
		dp_state_all_stages_flyby_lambert[last_non_empty_stage].end(),
		[](const state_flyby_lambert& a, const state_flyby_lambert& b)
		{return a.total_dv < b.total_dv; });

	state_flyby_lambert min_state = *min_iter;

	// Find the state with the maximum metric in the last non-empty stage
	int max_index = 0;
	total_dv = min_state.total_dv;
	if (if_display) cout << "Total dv  " << min_state.total_dv << endl;

	// Backtrack the optimal solution using the saved indices
	for (int stage = last_non_empty_stage; stage >= 0; --stage) {
		int id_target = sequence[stage];
		double t = min_state.t;

		if (stage != 0) {
			optimal_path[stage] = output_unit_flyby_lambert(t, min_state.rv, min_state.total_dv - min_state.pre_state->total_dv, min_state.total_dv, id_target);
			min_state = *min_state.pre_state;
		}
		else {
			optimal_path[stage] = output_unit_flyby_lambert(t, min_state.rv, min_state.total_dv, min_state.total_dv, id_target);
		}
	}

	result = optimal_path;

	if (step_t > 100.0) //100 经测试不会对结果造成1cm/s 的影响
	{
		DP_flyby_Lambert_iter(
			sequence,
			total_dv,
			result,
			false,
			"",
			total_dv+1.0e-3,
			step_t / 2.0, t_range / 2.0, false);
	}

	//if (if_display) print_optimal_path(result, cout);

	// 当output_file 不是nullptr时，输出到文件
	if (output_file != "")
	{
		print_optimal_path(result, cout);
		std::ofstream file_output(output_file);
		print_optimal_path(result, file_output);
		file_output.close();
	}

	return 1;
}





