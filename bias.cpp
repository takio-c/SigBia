#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

float Uniform( void ){
	static int x=10;
	int a=1103515245,b=12345,c=2147483647;
	x = (a*x + b)&c;

	return (float)(((double)x+1.0) / ((double)c+2.0));
}

float Normal( void ){
	return sqrt( -log(Uniform())*2.0 ) * sin( 2.0*M_PI*Uniform() );
}

#define VAR(s)	sqrt(var_##s)
#define STV_NUM	3				// 状態数
#define SEN_NUM	2				// 状態数
#define MES_NUM	1				// 観測数
#define Fs		100.0			// サンプリング周期
#define Ts		(1.0/Fs)		// サンプリング間隔
#define T		10.0			// 信号周期
#define OMEGA	(2.0*M_PI/T)	// 信号角速度
#define BIAS	20.0

double s[STV_NUM] = {
	0.0,	// low freq noise
	BIAS,	// bias
	0.0,	// high freq noise
};
double x[STV_NUM] = {
	0.0,	// signal
	0.0,	// differencial
	0.0,	// bias
};
double f[STV_NUM*STV_NUM] = {
	0.0,	0.0,	0.0,	// for signal
	0.0,	1.0,	0.0,	// for differencial
	0.0,	0.0,	0.0,	// for bias
};
double b[STV_NUM*SEN_NUM] = {
	1.0,	0.0,	// for sensor1 noise
	0.0,	0.0,	// for sensor1 bias
	0.0,	1.0,	// for sensor1 noise
};
double h[MES_NUM*STV_NUM] = {
	1.0,	1.0,	-1.0,	// 
};
// var
double var_e[STV_NUM*STV_NUM] = {
	1000.0,	0.0,	0.0,	// 
	0.0,	1000.0,	0.0,	// 
	0.0,	0.0,	1000.0,	// 
};
double var_v[SEN_NUM] = {
	60.,	// on sensor
	60.,	// on sensor
};
double var_w[MES_NUM] = {
	0.00001,	// on measure
};

double getSig(int i)
{
	double sig;
	double noise;
	sig = 40*sin(0.005*i) + 10;
	noise = sqrt(10)*Normal();
	return sig + noise;
}

int main(int argc, char* argv[])
{
	// status equation
	Matrix<double> S(STV_NUM, 1, "S");
	Matrix<double> X(STV_NUM, 1, "X");
	Matrix<double> F(STV_NUM, STV_NUM, "F");
	Matrix<double> B(STV_NUM, SEN_NUM, "B");
	Matrix<double> V(SEN_NUM, 1, "V");
	Matrix<double> Q(SEN_NUM, SEN_NUM, "Q");
	// measure
	Matrix<double> Z(MES_NUM, 1, "Z");
	Matrix<double> H(MES_NUM, STV_NUM, "H");
	Matrix<double> W(MES_NUM, 1, "W");
	Matrix<double> R(MES_NUM, MES_NUM, "R");
	// Kalman
	Matrix<double> P(STV_NUM, STV_NUM, "P");
	Matrix<double> K(STV_NUM, MES_NUM, "K");

	// init
	S.Set(s);
	X.Set(x);
	H.Set(h);
	F.Set(f);
	B.Set(b);
	// var
	Q.Dia(var_v);
	R.Dia(var_w);
	P.Set(var_e);

	unsigned int loop = 1000;
	if(1 < argc) loop = atoi(argv[1]);

	for(int i = 0; i < loop; i++){
		double sig = getSig(i);
		for(int n = 0; n < SEN_NUM; n++){
			V[n][0] = VAR(v[n])*Normal();
		}
		for(int n = 0; n < MES_NUM; n++){
			W[n][0] = VAR(w[n])*Normal();
		}
		// make error on sensor
		S = F * S + B * V;
		// set sensor
		double s1 = sig + S[0][0] + S[1][0];
		double s2 = sig + S[2][0];
		Z = H * S + W;
		// P previous
		P = F * P * F.tra() + B * Q * B.tra();
		// gain
		K = P * H.tra() * (H * P * H.tra() + R).inv();
		// update
		X = F * X + K * (Z - H * F * X);
		// P next
		P = ( (K * H).ide() - K * H) * P;
		// set estimate
		double x1 = s1 - X[0][0] - X[1][0];
		double x2 = s2 - X[2][0];
		// output
		printf("%lf %lf %lf %lf %lf %lf %lf%s",
				sig,			// signal
				s1,				// sensor1
				x1,				// estimate1
				s2,				// sensor2
				x2,				// estimate2
				S[1][0],		// sensor1 bias
				X[1][0],		// estimate bias
				"\n"
			  );
		// status update
	}
	P.print();

	return 0;
}

