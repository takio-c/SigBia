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
#define MES_NUM	1				// 観測数
#define Fs		100.0			// サンプリング周期
#define Ts		(1.0/Fs)		// サンプリング間隔
#define T		10.0			// 信号周期
#define OMEGA	(2.0*M_PI/T)	// 信号角速度

double sig[STV_NUM] = {
	0.0,	// signal
	0.0,	// differencial
	9.0,	// bias
};
double init[STV_NUM] = {
	0.0,	// signal
	0.0,	// differencial
	0.0,	// bias
};
double phi[STV_NUM*STV_NUM] = {
	1.0,	Ts,		0.0,	// for signal
	0.0,	1.0,	0.0,	// for differencial
	0.0,	0.0,	1.0,	// for bias
};
double h[MES_NUM*STV_NUM] = {
	1.0,	0.0,	1.0,	// 
};
// var
double var_e[STV_NUM*STV_NUM] = {
	0.0,	0.0,	0.0,	// 
	0.0,	0.0,	0.0,	// 
	0.0,	0.0,	1.0,	// 
};
double var_w[STV_NUM] = {
	0.0,	// on signal
	0.0,	// on differencial
	0.0,	// on bias
};
double var_v[MES_NUM] = {
	0.1,	// on sensor
};

int main(int argc, char* argv[])
{
	// true
	Matrix<double> S(STV_NUM, 1);
	// estimate
	Matrix<double> X(STV_NUM, 1);
	// sensor
	Matrix<double> H(MES_NUM, STV_NUM);
	Matrix<double> Z(MES_NUM, 1);
	Matrix<double> V(MES_NUM, 1);
	Matrix<double> R(MES_NUM, MES_NUM);
	// status vector
	Matrix<double> PHI(STV_NUM, STV_NUM);
	Matrix<double> W(STV_NUM, 1);
	Matrix<double> Q(STV_NUM, STV_NUM);
	// Kalman
	Matrix<double> P(STV_NUM, STV_NUM);
	Matrix<double> K(STV_NUM, MES_NUM);

	// init
	S.Set(sig);
	X.Set(init);
	H.Set(h);
	PHI.Set(phi);
	// var
	R.Dia(var_v);
	Q.Dia(var_w);
	P.Set(var_e);

	unsigned int loop = 10000;
	if(1 < argc) loop = atoi(argv[1]);

	for(int i = 0; i < loop; i++){
		for(int n = 0; n < STV_NUM; n++){
			W[n][0] = VAR(w[n])*Normal();
		}
		for(int n = 0; n < MES_NUM; n++){
			V[n][0] = VAR(v[n])*Normal();
		}
		// set sensor
		S = PHI * S + W;
		Z = H * S + V;
		// P previous
		P = PHI * P * PHI.tra() + Q;
		// gain
		K = P * H.tra() * (H * P * H.tra() + R).inv();
		// update
		X = PHI * X + K * (Z - H * PHI * X);
		// P next
		P = ( (K * H).ide() - K * H) * P;
		// output
		printf("%lf %lf %lf %lf %lf %s",
				Z[0][0],		// sensor
				S[0][0],		// true signal
				X[0][0],		// estimate signal
				P[0][0],		// error variance for signal
				P[2][2],		// error variance for bias
				"\n"
			  );
		// status update
#if 0
		if(2500 < i && i <= 5000){
			S[1][0] = OMEGA*cos(OMEGA*i*Ts);
			X[1][0] = OMEGA*cos(OMEGA*i*Ts);
		}
#endif
	}

	return 0;
}

