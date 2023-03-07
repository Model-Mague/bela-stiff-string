/*
 ____  _____ _        _
| __ )| ____| |      / \
|  _ \|  _| | |     / _ \
| |_) | |___| |___ / ___ \
|____/|_____|_____/_/   \_\

In render() you'll see a nested for loop structure. You'll see this in all Bela projects.
The first for loop cycles through 'audioFrames', the second through 'audioChannels' (in this case left 0 and right 1).
------------------------------------
Finite Difference Time Domain Stiff String, adapted to Real-Time C++ from Stefan Bilbao's MATLAB example.
extracted from: https://www2.ph.ed.ac.uk/~sbilbao/matlabpage.html

*/
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <vector>
#include <array>
#include "Eigen/Sparse"
#include "Eigen/LU"

#define DESKTOP_BUILD

#ifndef DESKTOP_BUILD
#include <Bela.h>
#endif

float inharmonicity = 0.001f;                                 // inharmonicity parameter (>0)
float f0 = 100.f;                                  // fundamental(Hz)
float TF = 2.f;                                    // duration of simulation(s)
float ctr = 0.1f;
float wid = 0.05f;                   			 // center location/width of excitation
float u0 = 1.f;
float v0 = 0.f;                           		 // maximum initial displacement/velocity
float rp[2] = { 0.3f, 0.7f };                        // positions of readout(0-1)
float loss[2][2] = { {100.f, 10.f}, {1000.f, 8.f} };       // loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
float theta = 1.0f;
float SR;

float NF;
float k;
float K;
float gumo;

float h;
float N;
float mu;
float lambda;
Eigen::Vector2f rp_int;
Eigen::Vector2f rp_frac;

float zeta1;
float zeta2;
float sig0;
float sig1;

Eigen::MatrixXf A;
Eigen::MatrixXf B;
Eigen::MatrixXf C;
Eigen::VectorXf u2;
Eigen::VectorXf u1;
Eigen::MatrixXf u;
Eigen::MatrixXf out;


void export_matrix(Eigen::MatrixXf mx, const std::string& path)
{
	std::ofstream out(path);
	out << mx;
}


Eigen::MatrixXf toeplitz(Eigen::VectorXf r, Eigen::MatrixXf zeros_1d_matrix)
{
	const auto size = r.size() + zeros_1d_matrix.cols();
	Eigen::MatrixXf out(size, size);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == 0)
			{
				if (j < r.size())
					out(i, j) = r(j);
				else
					out(i, j) = 0;//zeros_1d_matrix(0, j - r.size());
			}
			else if (j == 0)
			{
				if (i < r.size())
					out(i, j) = r(i);
				else
					out(i, j) = 0;//zeros_1d_matrix(0, i - r.size());
			}
			else out(i, j) = out(i - 1, j - 1);
		}
	}
	return out;
}

// Zeros(int n, m) -> ArrayXXf a3 = ArrayXXf::Zero(n, m);
Eigen::MatrixXf zeros(int n, int m)
{
	return Eigen::MatrixXf::Zero(n, m);
}

void max_set(Eigen::VectorXf& inout, float threshold)
{
	// This modifies the inout array
	for (int i = 0; i < inout.size(); i++)
	{
		if (inout(i) < threshold)
		{
			inout(i) = threshold;
		}
	}
}

void sign_set(Eigen::VectorXf& inout)
{
	for (int i = 0; i < inout.size(); i++)
	{
		if (inout(i) < 0.0f) inout(i) = -1.f;
		else if (inout(i) > 0.0f) inout(i) = 1.f;
	}
}

Eigen::VectorXf gen_consecutive_vec(int from, int to)
{
	int count = to - from;
	Eigen::VectorXf result(count);
	for (int i = 0; i < count; i++) result(i) = (float)(from + i);
	return result;
}


#ifdef DESKTOP_BUILD
bool setup()
#else
bool setup(BelaContext* context, void* userData)
#endif
{
	Eigen::MatrixXd m(2, 2);
	m(0, 0) = 3;
	m(1, 0) = 2.5;
	m(0, 1) = -1;
	m(1, 1) = m(1, 0) + m(0, 1);

#ifdef DESKTOP_BUILD
	SR = 44100;
#else
	SR = context->audioSampleRate;
#endif
	NF = floor(SR * TF);
	k = 1 / SR;                                  // time step
	gumo = 2 * f0;
	K = sqrt(inharmonicity) * (gumo / (float)M_PI);      // set parameters

	// stability conditions

	h = sqrt((powf(gumo, 2) * powf(k, 2) + sqrt(powf(gumo, 4) * powf(k, 4) + 16 * powf(K, 2) * powf(k, 2) * (2 * theta - 1))) / (2 * (2 * theta - 1)));
	N = floor(1 / h);
	h = 1 / N;
	mu = K * k / powf(h, 2);
	lambda = gumo * k / h;


	// readout interpolation parameters
	for (int i = 0; i < 2; i++)
	{
		rp_int(i) = 1 + floor(N * rp[i]);  // rounded grid index for readout
	}

	for (int i = 0; i < 2; i++)
	{
		rp_frac(i) = 1 + rp[i] / h - rp_int(i);  // fractional part of readout location
	}
	// set scheme loss 
	float gammaSq2 = powf(gumo, 2);
	float gammaSq4 = powf(gumo, 4);
	float Ksq2 = powf(K, 2);
	float piLoss00 = powf(2 * (float)M_PI * loss[0][0], 2);
	float piLoss10 = powf(2 * (float)M_PI * loss[1][0], 2);

	zeta1 = (-gammaSq2 + sqrt(gammaSq4 + 4 * Ksq2 * piLoss00)) / (2 * Ksq2);
	zeta2 = (-gammaSq2 + sqrt(gammaSq4 + 4 * Ksq2 * piLoss10)) / (2 * Ksq2);
	sig0 = 6 * (float)log(10) * (-zeta2 / loss[0][1] + zeta1 / loss[1][1]) / (zeta1 - zeta2);
	sig1 = 6 * (float)log(10) * (1 / loss[0][1] - 1 / loss[1][1]) / (zeta1 - zeta2);


	// create update matrices
	// zeros(1, N-3) 
	// [0]
	// or 
	// [0 0]
	// [0 0]
	// Commonly-used zero mx

	// First vector 2; zeros vector N - 3
	// Total dimension 2 + N - 3 = N - 1
	// Matrix element count (N-1) * (N-1)

	Eigen::MatrixXf zerosNMinus3 = zeros(1, (int)N - 3);
	Eigen::MatrixXf zerosNMinus4 = zeros(1, (int)N - 4);


	// M = sparse(toeplitz([theta (1-theta)/2 zeros(1,N-3)]));
	Eigen::VectorXf mComponent(2);
	mComponent(0) = theta;
	mComponent(1) = (1 - theta) / 2;

	Eigen::MatrixXf M = toeplitz(mComponent, zerosNMinus3);

	Eigen::VectorXf aComponent(2);
	aComponent(0) = sig1 * k / (powf(h, 2)) + sig0 * k / 2;
	aComponent(1) = -sig1 * k / (2 * powf(h, 2));
	A = toeplitz(aComponent, zerosNMinus3);

	// A = M + A;
	A = M + A;

	// C = M+sparse(toeplitz([-sig1*k/(h^2)-sig0*k/2 sig1*k/(2*h^2) zeros(1,N-3)]));
	Eigen::VectorXf cComponent(2);
	cComponent(0) = -sig1 * k / (powf(h, 2)) - sig0 * k / 2;
	cComponent(1) = sig1 * k / (2 * powf(h, 2));
	C = toeplitz(cComponent, zerosNMinus3);
	// C = M + C
	C = M + C;


	// float B = 2*M+sparse(toeplitz([-2*lambda^2-6*mu^2 lambda^2+4*mu^2 -mu^2 zeros(1,N-4)]));
	Eigen::VectorXf bComponent(3);
	bComponent(0) = -2 * powf(lambda, 2) - 6 * powf(mu, 2);
	bComponent(1) = powf(lambda, 2) + 4 * powf(mu, 2);
	bComponent(2) = -powf(mu, 2);

	B = 2 * M + toeplitz(bComponent, zerosNMinus4);


	// Create raised cosine
	// xax = [1:N-1]'*h;
	Eigen::VectorXf xax = gen_consecutive_vec(1, (int)N); // maybe needs tanspose here
	Eigen::VectorXf xaxMinusWid = gen_consecutive_vec(1, (int)N);
	Eigen::VectorXf xaxPlusWid = gen_consecutive_vec(1, (int)N);
	Eigen::VectorXf rc = gen_consecutive_vec(1, (int)N);

	for (int i = 0; i < (int)N - 1; i++) {
		xax(i) = xax(i) * h;
		// Copy values into soon-to-be-used matrices
		xaxMinusWid(i) = xax(i) - ctr - wid / 2;
		xaxPlusWid(i) = xax(i) - ctr + wid / 2;

		rc(i) = 1 + (float)cos(2 * M_PI * (xax(i) - ctr) / wid);
	}

	int n_as_int = (int)roundf(N);
	// ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
	Eigen::VectorXf ind(n_as_int - 1);
	// Multiply element-wise (xax-ctr-wid/2).*(xax-ctr+wid/2)
	for (int i = 0; i < (int)N - 1; i++)	ind(i) = -1 * xaxMinusWid(i) * xaxPlusWid(i);

	max_set(ind, 0);
	sign_set(ind);

	// rc = 0.5*ind.*(rc);
	for (int i = 0; i < (int)N - 1; i++) rc(i) = 0.5f * ind(i) * rc(i);

	/*
		u2 = u0*rc;
		u1 = (u0+k*v0)*rc;
		u = zeros(N+1,1);
		out = zeros(NF,2);
	*/
	u2 = Eigen::VectorXf(n_as_int - 1);
	u2 = u0 * rc;

	u1 = Eigen::VectorXf(n_as_int - 1);
	u1 = (u0 + k * v0) * rc;

	u = zeros((int)N + 1, 1);
	out = zeros((int)NF, 2);

	//export_matrix(u1, "u1.txt");
	//export_matrix(u2, "u2.txt");
	//export_matrix(xax, "xax.txt"); // maybe needs tanspose here
	//export_matrix(xaxMinusWid, "xaxMinusWid.txt");
	//export_matrix(xaxPlusWid, "xaxPlusWid.txt");
	//export_matrix(rc, "rc.txt");
	//export_matrix(ind, "ind.txt");

	return true;
}


#ifdef DESKTOP_BUILD
int main(int argc, char** argv)
#else
void render(BelaContext* context, void* userData)
#endif
{
#ifdef DESKTOP_BUILD
	setup();
#endif

	int first_index = (int)roundf(rp_int(0));
	int second_index = (int)roundf(rp_int(1));

	Eigen::Vector2f urpint;
	Eigen::Vector2f urpintplus;
	Eigen::Vector2f rp_frac_minus1 = Eigen::Vector2f::Ones(2) - rp_frac;

#ifdef DESKTOP_BUILD
	for (unsigned int n = 0; n < 20; n++)
	{
#else
	for (unsigned int n = 0; n < context->audioFrames; n++)
	{
#endif
		Eigen::MatrixXf Bmember = B * u1 - C * u2;
		Eigen::VectorXf u = A.lu().solve(Bmember);
		urpint(0) = u(first_index);
		urpint(1) = u(second_index);
		urpintplus(0) = u(first_index + 1);
		urpintplus(1) = u(second_index + 1);
		// (1-rp_frac).*u(rp_int)'+rp_frac.*u(rp_int+1)'
		Eigen::Vector2f out = (rp_frac_minus1).cwiseProduct(urpint) + rp_frac.cwiseProduct(urpintplus); // this is two numbers
		//export_matrix(rp_frac_minus1, "rp_frac_minus1.txt"); // OK
		//export_matrix(urpint, "urpint.txt");  // Precision issue
		//export_matrix(rp_frac, "rp_frac.txt"); // OK
		//export_matrix(urpintplus, "urpintplus.txt"); // Precision issue

#ifndef DESKTOP_BUILD
		for (unsigned int channel = 0; context->audioOutChannels; channel++)
		{
			audioWrite(context, n, channel, out(channel % 2));
		}
#endif
		// Update
		u2 = u1;
		u1 = u;
	}
}

#ifndef DESKTOP_BUILD
void cleanup(BelaContext* context, void* userData)
{

}
#endif
