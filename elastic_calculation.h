//弾性計算メイン
#include "Rigidbody.h"
;//セミコロンを入れないとコンパイルエラーになる・・・
void calc_elastic(vector<mpselastic> &PART, elastic &ELAST, int t, double **F);

//内部ひずみ・応力から加速度を計算する
//void calc_accel_for_2D(vector<mpselastic> &PART, elastic &ELAST);
void calc_accel_for_3D(vector<mpselastic> &PART, elastic &ELAST);

//体積ひずみ・圧力から加速度を計算する／接触による影響も同時に計算する
void calc_pressure_and_contact(vector<mpselastic> &PART, elastic &ELAST);
void calc_contact(vector<mpselastic> &PART, elastic &ELAST);

//接触判定
void contact_judge(vector<mpselastic> &PART, elastic &ELAST);

//加速度の修正
void modify_acceleration(vector<mpselastic> &PART, elastic &ELAST, double **F);

//速度・位置の計算
void calc_pre_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, double **F);
void calc_post_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST,int t);
void calc_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, double **F);

//剛体の位置修正
void calc_rigid(vector<mpselastic> &PART,vector<Rigidbody> &terminal);

//角速度の計算
void calc_angular_velocity(elastic &ELAST, vector<mpselastic> &PART);
void calc_angular_velocity_using_vectorSTL(elastic &ELAST, vector<mpselastic> &PART);

//クォータニオンの計算
void calc_quaternion(elastic &ELAST, vector<mpselastic> &PART);
void calc_quaternion_using_vectorSTL(elastic &ELAST, vector<mpselastic> &PART);
bool check_q_norm(const double *q, const double EPSILON); //クォータニオンの差分をチェックする
bool check_q_norm(const vector<double> &q, const double EPSILON); //クォータニオンの差分をチェックする

//クォータニオン用の回転行列
inline double rotx(double *r, double *r0, double *q);
inline double roty(double *r, double *r0, double *q);
inline double rotz(double *r, double *r0, double *q);
inline double rotn(double *q);

//並列計算用にqにポインタを使わない
inline double rotx(const double *r, const double *r0, const vector<double> &q);
inline double roty(const double *r, const double *r0, const vector<double> &q);
inline double rotz(const double *r, const double *r0, const vector<double> &q);
inline double rotn(const vector<double> &q);

//回転行列の計算
void calc_r0_ij(vector<mpselastic> &PART);
void rotate_r0(const double *rInit, const double *q, double *result);
void quaternion(double &R, double axis[3], double angle);

//密度の更新(flag_modify_densityで制御する)
void calc_modified_density(vector<mpselastic> &PART, elastic &ELAST);

//エネルギーの計算
void calc_hamiltonian(elastic &ELAST, int t);

//剛体の計算
void rigid_calc(vector<mpselastic> &PART, elastic &ELAST);

//非線形計算
void calc_nonlinear_accel_for_3D_test(vector<mpselastic> &PART, elastic &ELAST);
void calc_nonlinear_elastic_constants(vector<mpselastic> &PART, double *E_ij, double *E_ji);
void calc_nonlinear_elastic(vector<mpselastic> &PART, elastic &ELAST, int t, double **F);//I/F
void calc_nonlinear_accel_for_3D(vector<mpselastic> &PART, elastic &ELAST);
void calc_nonlinear_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, vector<vector<double>> &residual_acceleration);
void calc_residual_acceleration(vector<mpselastic> &PART, elastic &ELAST, vector<vector<double>> &res_accel, int t, double **F);

//クォータニオン用のヤコビアン
inline double rotx_x(double *r, double *r0, double *q){return (2*(r[1]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])-r[2]*(r0[0]*q[1]-r0[1]*q[0]-r0[2]*q[3])));}
inline double rotx_y(double *r, double *r0, double *q){return (2*(r[1]*(-r0[0]*q[3]+r0[1]*q[2]-r0[2]*q[1])-r[2]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])));}
inline double rotx_z(double *r, double *r0, double *q){return (2*(r[1]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])-r[2]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])));}
inline double rotx_s(double *r, double *r0, double *q){return (2*(r[1]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])-r[2]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])));}

inline double roty_x(double *r, double *r0, double *q){return (2*(r[2]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])-r[0]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])));}
inline double roty_y(double *r, double *r0, double *q){return (2*(r[2]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])-r[0]*(-r0[0]*q[3]+r0[1]*q[2]-r0[2]*q[1])));}
inline double roty_z(double *r, double *r0, double *q){return (2*(r[2]*(-r0[0]*q[2]-r0[1]*q[3]+r0[2]*q[0])-r[0]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])));}
inline double roty_s(double *r, double *r0, double *q){return (2*(r[2]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])-r[0]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])));}

inline double rotz_x(double *r, double *r0, double *q){return (2*(r[0]*(r0[0]*q[1]-r0[1]*q[0]-r0[2]*q[3])-r[1]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])));}
inline double rotz_y(double *r, double *r0, double *q){return (2*(r[0]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])-r[1]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])));}
inline double rotz_z(double *r, double *r0, double *q){return (2*(r[0]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])-r[1]*(-r0[0]*q[2]-r0[1]*q[3]+r0[2]*q[0])));}
inline double rotz_s(double *r, double *r0, double *q){return (2*(r[0]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])-r[1]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])));}

inline double rotn_x(double *q){return 2*q[0];}
inline double rotn_y(double *q){return 2*q[1];}
inline double rotn_z(double *q){return 2*q[2];}
inline double rotn_s(double *q){return 2*q[3];}

//並列計算用
inline double rotx_x(const double *r, const double *r0, const vector<double> &q){return (2*(r[1]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])-r[2]*(r0[0]*q[1]-r0[1]*q[0]-r0[2]*q[3])));}
inline double rotx_y(const double *r, const double *r0, const vector<double> &q){return (2*(r[1]*(-r0[0]*q[3]+r0[1]*q[2]-r0[2]*q[1])-r[2]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])));}
inline double rotx_z(const double *r, const double *r0, const vector<double> &q){return (2*(r[1]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])-r[2]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])));}
inline double rotx_s(const double *r, const double *r0, const vector<double> &q){return (2*(r[1]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])-r[2]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])));}

inline double roty_x(const double *r, const double *r0, const vector<double> &q){return (2*(r[2]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])-r[0]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])));}
inline double roty_y(const double *r, const double *r0, const vector<double> &q){return (2*(r[2]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])-r[0]*(-r0[0]*q[3]+r0[1]*q[2]-r0[2]*q[1])));}
inline double roty_z(const double *r, const double *r0, const vector<double> &q){return (2*(r[2]*(-r0[0]*q[2]-r0[1]*q[3]+r0[2]*q[0])-r[0]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])));}
inline double roty_s(const double *r, const double *r0, const vector<double> &q){return (2*(r[2]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])-r[0]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])));}

inline double rotz_x(const double *r, const double *r0, const vector<double> &q){return (2*(r[0]*(r0[0]*q[1]-r0[1]*q[0]-r0[2]*q[3])-r[1]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])));}
inline double rotz_y(const double *r, const double *r0, const vector<double> &q){return (2*(r[0]*(r0[0]*q[0]+r0[1]*q[1]+r0[2]*q[2])-r[1]*(-r0[0]*q[1]+r0[1]*q[0]+r0[2]*q[3])));}
inline double rotz_z(const double *r, const double *r0, const vector<double> &q){return (2*(r[0]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])-r[1]*(-r0[0]*q[2]-r0[1]*q[3]+r0[2]*q[0])));}
inline double rotz_s(const double *r, const double *r0, const vector<double> &q){return (2*(r[0]*(r0[0]*q[2]+r0[1]*q[3]-r0[2]*q[0])-r[1]*(r0[0]*q[3]-r0[1]*q[2]+r0[2]*q[1])));}

inline double rotn_x(const vector<double> &q){return 2*q[0];}
inline double rotn_y(const vector<double> &q){return 2*q[1];}
inline double rotn_z(const vector<double> &q){return 2*q[2];}
inline double rotn_s(const vector<double> &q){return 2*q[3];}

//その他
void check_elastic_config(elastic &ELAST);
void check_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, const int i, double **F);
void vector_product(double *a, double *b, double *result); //ベクトル積

/****************************************************************/
/*********************以下旧バージョンの関数*********************/
/****************************************************************/

//加速度の計算
void calc_accel_for_3D_ver_1(vector<mpselastic> &PART, elastic &ELAST);
void calc_accel_for_3D_ver_2(vector<mpselastic> &PART, elastic &ELAST);
void calc_accel_for_3D_ver_3(vector<mpselastic> &PART, elastic &ELAST);
void stress_3D_original(vector<mpselastic> &PART, elastic &ELAST, double n0);
void stress_3D_test(vector<mpselastic> &PART, elastic &ELAST, double n0);

void calc_pressure_and_contact_ver_3(vector<mpselastic> &PART, elastic &ELAST);

void velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, double n0, int *INDEX, int **MESH, int t, double **F, double &kinetic_v, double &potential);

//角速度と角度の更新
void calc_angular_velocity_and_angle(vector<mpselastic> &PART, elastic &ELAST);

//二次元シンプレクティックスキーム用のトルク・角度・角速度計算
void calc_torque_for_2D(vector<mpselastic> &PART, elastic &ELAST);

//その他
void output_stress_for_GNU(vector<mpselastic> &PART, elastic &ELAST, int t);
void numerical_jacobian(double **J, double *r, double *r0, double *q);