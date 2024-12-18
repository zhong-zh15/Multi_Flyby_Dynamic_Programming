#include "GTOC4_Problem.h"
#include "DP_flyby_2impulse.h"
#include "GTOC11_problem.h"
#include <chrono>
#include <fstream>
#include <iostream>

#include "DP_flyby_2impulse_user_function.h"

int main()
{
	read_data_GTOC11();
	GTOC11_sequence_process();

    read_asteroid_data_GTOC4();
    GTOC4_sequence_process();

	//{
 //       double last_dv_end[3]= { -1018.0583418161496, -1287.8309532006251, 1145.6853084687846 };
 //       double next_dv_start[3]= { 957.93728541547353, 1313.2104821549692, -1166.7200050782531 };
 //       double rv_body[6];

 //       rendezvous2flyby(rv_body, last_dv_end, next_dv_start); //飞越2.0km/s 由该函数控制
	//}


    const vector<int>& Sequence_gtoc4= GTOC4_get_result_sequence();
    const vector<double>& T_sequence_gtoc4 = GTOC4_get_result_T();

    problem_case = 4;
	{
        std::cout << "Computing GTOC4: "  << "  Starting!!!!!" << std::endl;
        auto start_t = std::chrono::steady_clock::now();
        double dv;
        double dv_threhold = 1.0e8;
        vector<output_unit_flyby_lambert> result;
        //DP_flyby_Lambert(Sequence_gtoc4, dv, result, dv_threhold);
        //std::string output_file_name = "DP_result_GTOC4.txt";
        std::string output_file_name = "";
        DP_flyby_Lambert_iter(Sequence_gtoc4, dv, result, false, output_file_name, dv_threhold, step_t, step_t *12.0, true);
        auto end_t = std::chrono::steady_clock::now();
        double duration_t = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count() / 1000.0;
        std::cout << "Computing GTOC4:: " << "  END!!!!!" << "   Time is " << duration_t << "s ; which is " << duration_t / 60.0 << " mins;  which is " << duration_t / 3600.0 << " hours" << std::endl;
	}


    const auto& sequence_gtoc11_all = GTOC11_sequence_get_sequence();
    const auto& T_sequence_all = GTOC11_sequence_get_T();
    const auto& T_MAX_gtoc11_sequence_all = GTOC11_sequence_get_Maximum_T();

    std::ofstream file_output("JUST_T_RESULT.txt");
    problem_case = 11;
    for(int i = 0; i < sequence_gtoc11_all.size(); i++)
    {
        std::cout << "Computing Mission: " << i << ".  Starting!!!!!" << std::endl;
        auto start_t = std::chrono::steady_clock::now();
        double dv;
		double dv_threhold = 1.0e8;
        vector<output_unit_flyby_lambert> result;
        //DP_flyby_Lambert(sequence_gtoc11_all[i], dv, result, dv_threhold);
        //std::string output_file_name = "DP_result_" + std::to_string(i) + ".txt";
        std::string output_file_name = "";
        DP_flyby_Lambert_iter(sequence_gtoc11_all[i], dv, result, false, output_file_name, dv_threhold, step_t_GTOC11, step_t_GTOC11*6.0, true);
        auto end_t = std::chrono::steady_clock::now();
        double duration_t = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count() / 1000.0;
        std::cout << "Computing Mission: " << i << ".  END!!!!!" << "   Time is " << duration_t << "s ; which is " << duration_t / 60.0 << " mins;  which is " << duration_t / 3600.0 << " hours" << std::endl;
        double t_last = 0.0;
        for(int transfer_id = 0; transfer_id < result.size() ; transfer_id++)
        {
            file_output << std::fixed << std::setprecision(18) << (result[transfer_id].t - t_last) / (T_MAX_gtoc11_sequence_all[i][transfer_id] - t_last) << "  ";
			t_last = result[transfer_id].t;
        }
        file_output << std::endl;
    }

    file_output.close();

	system("pause");

	return 0;
}
