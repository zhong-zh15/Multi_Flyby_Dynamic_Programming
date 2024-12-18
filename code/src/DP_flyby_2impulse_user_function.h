#ifndef DP_flyby_2impulse_user_function_h
#define DP_flyby_2impulse_user_function_h
#include "DP_flyby_2impulse.h"
#include "GTOC4_Problem.h"
#include <vector>
#include <ostream>

using namespace std;
extern int problem_case;

// �û��Զ��庯��
// ���ؿ��ܵ�tֵ������������Ծ������û������Լ������ⶨ��
vector<double> possible_t_values(const vector<int>& sequence, //����׶�����
    int stage, std::vector<state_flyby_lambert>& last_layer_state //��Ҫ�Ĳ���������չ��
);
// �û��Զ��庯��
// ���㵱ǰʱ�̵�rv,dv,�Լ���dv������������Ծ������û������Լ������ⶨ��
void compute_rv_dv(vector<int> sequence,
    state_flyby_lambert& last_state,
    int id_departue, int id_now,
    double state_t,
    double* rv_arrival,
    double& dv,
    double& total_dv
);

// �û��Զ��庯��
// ��� state_t �� last_state_t �Ƿ��ܹ��γ�һ������Ҫ��Ľ⣬���Լ�֦һЩ�������
// Ҳ������ȫ��������
bool check_t(double last_state_t, double state_t, int ID_departue, int ID_now);

// �û��Զ��庯��
// ������ս�������Ը����Լ������ⶨ���������
void print_optimal_path(const std::vector<output_unit_flyby_lambert>& optimal_path, 
std::ostream& out_stream);


#endif

