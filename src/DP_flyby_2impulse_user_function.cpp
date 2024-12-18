#include "DP_flyby_2impulse_user_function.h"

#include <iomanip>
#include <iostream>
#include <numeric>    //accumulate

#include "GTOC11_problem.h"
#include "OrbitFun.h"
#include "OrbitMath.h"


int problem_case = 4;  // 4: GTOC4, 11: GTOC11

// �û��Զ��庯��
// ���ؿ��ܵ�tֵ������������Ծ������û������Լ������ⶨ��
vector<double> possible_t_values(const vector<int>& sequence, //����׶�����
    int stage, std::vector<state_flyby_lambert>& last_layer_state //��Ҫ�Ĳ���������չ��
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

        //double t_center = T_sequence[stage]; //�����Ƿ�ԭ���м����������������ֵʱһ���ģ���֤��ȷ�ԣ�����֤
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


        //�ж�sequence��sequence_all�ĵڼ��������index
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


		//double t_center = T_sequence[stage]; //�����Ƿ�ԭ���м����������������ֵʱһ���ģ���֤��ȷ�ԣ�����֤


		//t_center Ӧ��ʱ��0�����ʱ�̾��ֵ�
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

// �û��Զ��庯��
// ���㵱ǰʱ�̵�rv,dv,�Լ���dv������������Ծ������û������Լ������ⶨ��
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

		//last_state.rv �����Ͼ�����һ��״̬����Lambert���rv_arrival
        double last_dv_end[3];
		double next_dv_start[3];
		V_Minus(last_dv_end, rv_body + 3, last_state.rv + 3, 3);
		V_Minus(next_dv_start, rv_departue + 3, rv_body + 3, 3);

        rendezvous2flyby(rv_body, last_dv_end, next_dv_start); //��Խ2.0km/s �ɸú�������

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


// �û��Զ��庯��
// ��� state_t �� last_state_t �Ƿ��ܹ��γ�һ������Ҫ��Ľ⣬���Լ�֦һЩ�������
// Ҳ������ȫ��������
bool check_t(double last_state_t, double state_t, int ID_departue, int ID_now) {

    if (problem_case == 4)
    {
        if (last_state_t > -1.e6)
        {
            //���state_t �� last_state_t ���С��1�죬����false
            if (state_t - last_state_t < 1.0 * 86400.) {
                return false;
            }
            if ((state_t - last_state_t) > 365.25 * 86400.0 * 1.0) {
                //���state_t �� last_state_t ������1�꣬����false
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

        ////���rv_departure��coe������Ҫ��(�볤�����1.8 AU �� ƫ���ʴ���0.3)������false
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
            //���state_t �� last_state_t ���С��2�죬����false
            if (state_t - last_state_t < 1.0 * 86400.) {
                return false;
            }
            if ((state_t - last_state_t) > 365.25 * 86400.0 * 2.0) {
                //���state_t �� last_state_t ������2�꣬����false
                return false;
            }
        }
        const auto& sequence_all = GTOC11_sequence_get_sequence();
        //const auto& T_sequence_all = GTOC11_sequence_get_T();
        const auto& T_MAX_sequence_all = GTOC11_sequence_get_Maximum_T();

        //�ж�ID_now��sequence_all�ĵڼ�������ң�������ĳ�����е��м䣩������index
        int index_i = -1;
        int index_j = -1;
		for (int i = 0; i < sequence_all.size(); i++) //�ж�ID_now�Ƿ���sequence_all[i]�У�����ڣ�����i
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


// �û��Զ��庯��
// ������ս�������Ը����Լ������ⶨ���������
void print_optimal_path(const std::vector<output_unit_flyby_lambert>& optimal_path, std::ostream& out_stream) {

    if (problem_case == 4)
    {
        // ��� dv
        double total_dv = std::accumulate(optimal_path.begin(), optimal_path.end(), 0.0,
            [](double sum, const output_unit_flyby_lambert& unit) { return sum + unit.dv; });
        out_stream << "Total DV: " << std::fixed << std::setprecision(3)
            << total_dv << "   Mass: "<< 1500.0 * exp(-total_dv/g0_GTOC4/3000.0)<< std::endl;

        // �������
        out_stream << std::setw(6) << "Stage" << std::setw(19) << "Epoch"
            << std::setw(19) << "X(m)" << std::setw(19) << "Y(m)" << std::setw(19) << "Z(m)"
            << std::setw(13) << "VX(m/s)" << std::setw(13) << "VY(m/s)" << std::setw(13) << "VZ(m/s)"
            << std::setw(11) << "DV(m/s)" << std::setw(11) << "TOTAL_DV" << std::setw(11) << "MASS(kg)"
            << std::setw(6) << "ID" << std::endl;

        // ��� optimal_path �е�����
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
        // ��� dv
        out_stream << "Total DV: " << std::fixed << std::setprecision(3)
            << std::accumulate(optimal_path.begin(), optimal_path.end(), 0.0,
                [](double sum, const output_unit_flyby_lambert& unit) { return sum + unit.dv; }) << std::endl;

        // �������
        out_stream << std::setw(6) << "Stage" << std::setw(19) << "Epoch"
            << std::setw(19) << "X(m)" << std::setw(19) << "Y(m)" << std::setw(19) << "Z(m)"
            << std::setw(13) << "VX(m/s)" << std::setw(13) << "VY(m/s)" << std::setw(13) << "VZ(m/s)"
            << std::setw(11) << "DV(m/s)" << std::setw(11) << "TOTAL_DV"
            << std::setw(6) << "ID" << std::endl;

        // ��� optimal_path �е�����
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

