#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

double Normal(double exp, double var)
{
	static int sw = 0;
	static double n1;
	double n2;

	double r1 = ((double)rand()) / RAND_MAX;
	double r2 = ((double)rand()) / RAND_MAX;
	if(sw == 0){
		n1 = sqrt(-2.0 * log(r1)) * cos(2.0 * M_PI * r2);
		n2 = sqrt(-2.0 * log(r1)) * sin(2.0 * M_PI * r2);
		sw = 1;
		return sqrt(var) * n2 + exp;
	}
	else{
		sw = 0;
		return sqrt(var) * n1 + exp;
	}
}

#define STV_NUM	3				// 状態数
#define NOI_NUM	2				// ノイズ源数
#define MES_NUM	1				// 観測数
#define Fs		100.0			// サンプリング周期
#define Ts		(1.0/Fs)		// サンプリング間隔
#define T		2.0			// 信号周期
#define OMEGA	(2.0*M_PI/T)	// 信号角速度
#define BIAS	1.0
#define LOOP	1500
#define SIG_DELAY	100
#define SIG_TIME	400
#define SIG_WAIT	400

class Sensor
{

};

class KalmanFilter
{
protected:
	// name for debug
	const int stv_num;
	const int inp_num;
	const int mes_num;
	std::string name;
	// status space model
	Matrix<double> F_m;
	Matrix<double> x_v;
	Matrix<double> G_m;
	Matrix<double> u_v;
	Matrix<double> Q_m;
	// measurement
	Matrix<double> z_v;
	Matrix<double> H_m;
	Matrix<double> R_m;
	// Kalman Filter
	Matrix<double> P_m;
	Matrix<double> K_m;

	virtual void init() {}
	virtual void updateModel() {}

public:
	KalmanFilter(int stv_num, int inp_num, int mes_num, std::string name="")
		:stv_num(stv_num), inp_num(inp_num), mes_num(mes_num), name(name)
		// status space model
		,F_m(stv_num,stv_num,name+":F_m")
		,x_v(stv_num,1,name+":x_v")
		,G_m(stv_num,inp_num,name+":G_m")
		,u_v(inp_num,1,name+":u_v")
		,Q_m(stv_num,stv_num,name+":Q_m")
		// measurement
		,z_v(mes_num,1,name+":z_v")
		,H_m(mes_num,stv_num,name+":H_m")
		,R_m(mes_num,mes_num,name+":R_m")
		// Kalman Filter
		,P_m(stv_num,stv_num,name+":P_m")
		,K_m(stv_num,mes_num,name+":K_m")
		{
			// init
			init();
		}
	virtual ~KalmanFilter() {
	}
	void update(double ts) {
		Matrix<double> PHI_m(F_m.Row(),F_m.Col(),name+":I + dtxF_m");
		Matrix<double> PSI_m(G_m.Row(),G_m.Col(),name+":dtxG_m");
		// update status model for linearize
		updateModel();
		PHI_m = F_m.ide() + F_m.mul(ts);
		PSI_m = G_m.mul(ts);
		// transition
		x_v = PHI_m * x_v + PSI_m * u_v;
		P_m = PHI_m * P_m * PHI_m.tra() + Q_m;
		// gain
		K_m = P_m * H_m.tra() * (H_m * P_m * H_m.tra() + R_m).inv();
		// least squares
		x_v = x_v + K_m * (z_v - H_m * x_v);
		// update error co-variance
		P_m = ( Matrix<double>::ide(K_m.Row()) - K_m * H_m ) * P_m;
	}
};

double x[STV_NUM] = {
	0.0,	// accelerometer noise
	0.0,	// gyro bias(k)
	0.0,	// gyro bias diff
};
double f[STV_NUM*STV_NUM] = {
	0.0,	0.0,	0.0,	// accelerometer noise
	0.0,	1.0,	0.0,	// gyro bias(k)
	0.0,	0.0,	1.0,	// gyro bias diff
};
double b[STV_NUM*NOI_NUM] = {
	1.0,	0.0,	// accelerometer
	0.0,	1.0,	// gyro bias(k)
	0.0,	0.0,	// gyro bias diff
};
double h[MES_NUM*STV_NUM] = {
	1.0,	-Ts,	-1.0,	// (noise of accelerometer) - total_time * (differential of gyro bias)
};
// var
double var_e[STV_NUM] = {
	1.0,	// about accelerometer
	1.0,	// about gyro bias(k)
	1.0,	// about gyro bias diff
};
double var_v[NOI_NUM] = {
	0.1,	// on accelerometer
	0.00001,	// on gyro bias
};
double var_w[MES_NUM] = {
	0.0,	// on measure
};

double getTheta(int i)
{
	double theta=0.0;
	double noise=0.0;
	if(i < SIG_DELAY){
		i = SIG_DELAY;
	}
	i -= SIG_DELAY;
	if(SIG_TIME < i){
		i = SIG_TIME;
	}
	theta = 0.25 * (OMEGA*Ts*i - sin(OMEGA*Ts*i));
//	noise = Normal(0.0, 0.1);
	return theta + noise;
}

double getOmega(int i)
{
	double omega=0.0;
	double noise=0.0;
	if(i < SIG_DELAY){
		i = SIG_DELAY;
	}
	i -= SIG_DELAY;
	if(SIG_TIME < i){
		i = SIG_TIME;
	}
	omega = 0.25*OMEGA*(1.0 - cos(OMEGA*Ts*i));
//	noise = Normal(0.0, 0.1);
	return omega + noise;
}

int main(int argc, char* argv[])
{
	// status equation
	Matrix<double> X(STV_NUM, 1, "X");
	Matrix<double> F(STV_NUM, STV_NUM, "F");
	Matrix<double> B(STV_NUM, NOI_NUM, "B");
	Matrix<double> V(NOI_NUM, 1, "V");
	Matrix<double> Q(NOI_NUM, NOI_NUM, "Q");
	// measure
	Matrix<double> Z(MES_NUM, 1, "Z");
	Matrix<double> H(MES_NUM, STV_NUM, "H");
	Matrix<double> R(MES_NUM, MES_NUM, "R");
	// Kalman
	Matrix<double> P(STV_NUM, STV_NUM, "P");
	Matrix<double> K(STV_NUM, MES_NUM, "K");
	// sensor
	double theta_a, theta_w;
	// integral
	double sig_t, sig_o;

	// init
	X.Set(x);
	H.Set(h);
	F.Set(f);
	B.Set(b);
	// var
	Q.Dia(var_v);
	R.Dia(var_w);
	P.Dia(var_e);
	// integral
	theta_a = theta_w = 0.0;
	sig_t = sig_o = 0.0;

	unsigned int loop = LOOP;
	if(1 < argc) loop = atoi(argv[1]);

	for(int i = 0; i < loop; i++){
		// true parameter
		double theta = getTheta(i);
		double omega = getOmega(i);
		// sensor white noise
		for(int n = 0; n < NOI_NUM; n++){
			V[n][0] = Normal(0.0, var_v[n]);
		}
		// sensor
		double sen_a = theta + V[0][0];
		double sen_w = omega + BIAS + V[1][0];
		// measure
		theta_a = sen_a;
		theta_w += Ts*sen_w;
		Z[0][0] = theta_a - theta_w;
		// P previous
		P = F * P * F.tra() + B * Q * B.tra();
		// gain
		K = P * H.tra() * (H * P * H.tra() + R).inv();
		K.print();
		// update
		X = F * X + K * (Z - H * F * X);
		// P next
		P = ( Matrix<double>::ide(K.Row()) - K * H) * P;
		// feed back
		theta_a -= X[0][0];
		theta_w -= Ts*X[1][0];
		// output
		printf("%lf %lf %lf %lf %lf %lf %lf %lf %s",
				theta,			// theta
				omega,			// omega
				sen_a,			// accelerometer
				sen_w,			// gyro
				theta_a,		// from accl
				theta_w,		// from gyro
				X[1][0],			// true tyro bias
				X[2][0],		// estimate gyro bias
				"\n"
			  );
		// status update
	}

	return 0;
}

