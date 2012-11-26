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

#define STV_NUM	2				// 状態数
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

double s[STV_NUM] = {
	0.0,	// accelerometer noise
	BIAS,	// gyro bias
};
double x[STV_NUM] = {
	0.0,	// accelerometer noise
	0.0,	// gyro bias
};
double f[STV_NUM*STV_NUM] = {
	0.0,	0.0,	// accelerometer noise
	0.0,	1.0,	// gyro bias
};
double b[STV_NUM*NOI_NUM] = {
	1.0,	0.0,	// accelerometer noise
	0.0,	1.0,	// gyro bias drift
};
double h[MES_NUM*STV_NUM] = {
	1.0,	1.0,	// 
};
// var
double var_e[STV_NUM] = {
	1.0,	// about accelerometer
	1.0,	// about gyro bias
};
double var_v[NOI_NUM] = {
	0.1,		// on accelerometer
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
		double sen_w = omega + BIAS + V[1][0] + Normal(0.0, 0.01);
		// measure
		theta_a = sen_a;
		theta_w += Ts * sen_w;
		Z[0][0] = theta_a - theta_w;
		// P previous
		P = F * P * F.tra() + B * Q * B.tra();
		// gain
		H[0][1] = - (i+1) * Ts;
		K = P * H.tra() * (H * P * H.tra() + R).inv();
		// update
		X = F * X + K * (Z - H * F * X);
		// P next
		P = ( (K * H).ide() - K * H) * P;
		// feed back
		theta_a -= X[0][0];
		double theta_ww = theta_w - (i+1)*Ts*X[1][0];
//		theta_w -= Ts*(X[1][0] + X[2][0]);
		// output
		printf("%lf %lf %lf %lf %lf %lf %lf %lf %s",
				theta,			// theta
				omega,			// omega
				sen_a,			// accelerometer
				sen_w,			// gyro
				theta_a,		// from accl
				theta_ww,		// from gyro
				BIAS,			// true tyro bias
				X[1][0],		// estimate gyro bias
				"\n"
			  );
		// status update
	}
	P.print();

	return 0;
}

