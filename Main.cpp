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
#include <vector>
#include <iostream>
#include <vector>
#include <array>
#include "Eigen/Sparse"
#include "Eigen/LU"

float inharmonicity = 0.001;                                 // inharmonicity parameter (>0)
float f0 = 100;                                  // fundamental(Hz)
float TF = 2;                                    // duration of simulation(s)
float ctr = 0.1;
float wid = 0.05;                   			 // center location/width of excitation
float u0 = 1;
float v0 = 0;                           		 // maximum initial displacement/velocity
float rp[2] = { 0.3, 0.7 };                        // positions of readout(0-1)
float loss[2][2] = { {100, 10}, {1000, 8} };       // loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
float theta = 1.0;
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




Eigen::MatrixXf toeplitz(Eigen::VectorXf r, Eigen::MatrixXf zeros_1d_matrix)
{
	// @HERE I really don't think the size calculation here is right
	int size = r.size() + zeros_1d_matrix.cols();
	Eigen::MatrixXf out(size, size);
	int n_as_int = (int)roundf(N);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == 0)
			{
				if (j < n_as_int)
					out(i, j) = r(j);
				else
					out(i, j) = zeros_1d_matrix(0, j - n_as_int);
			}
			else if (j == 0)
			{
				if (i < n_as_int)
					out(i, j) = r(i);
				else
					out(i, j) = zeros_1d_matrix(0, i - n_as_int);
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

void max_set(Eigen::VectorXf inout, float threshold)
{
	// This modifies the inout array
	for (int i = 0; i < inout.size(); i++)
	{
		if (inout(i) < threshold) inout(i) = threshold;
	}
}

void sign_set(Eigen::VectorXf inout)
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

bool setup()
{
	Eigen::MatrixXd m(2, 2);
	m(0, 0) = 3;
	m(1, 0) = 2.5;
	m(0, 1) = -1;
	m(1, 1) = m(1, 0) + m(0, 1);

	SR = 44100;//context->audioSampleRate;
	NF = floor(SR * TF);
	k = 1 / SR;                                  // time step
	gumo = 2 * f0;
	K = sqrt(inharmonicity) * (gumo / M_PI);      // set parameters

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
	// set scheme loss parameters
	zeta1 = (-powf(gumo, 2) + sqrt(powf(gumo, 4) + 4 * powf(K, 2) * (2 * M_PI * powf(loss[0][0], 2)))) / (2 * powf(K, 2));
	zeta2 = (-powf(gumo, 2) + sqrt(powf(gumo, 4) + 4 * powf(K, 2) * (2 * M_PI * powf(loss[1][0], 2)))) / (2 * powf(K, 2));
	sig0 = 6 * log10(10) * (-zeta2 / loss[0][1] + zeta1 / loss[1][1]) / (zeta1 - zeta2);
	sig1 = 6 * log10(10) * (1 / loss[0][1] - 1 / loss[1][1]) / (zeta1 - zeta2);


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

	Eigen::MatrixXf zerosNMinus3 = zeros(1, N - 3);
	Eigen::MatrixXf zerosNMinus4 = zeros(1, N - 4);


	// M = sparse(toeplitz([theta (1-theta)/2 zeros(1,N-3)]));
	Eigen::VectorXf mComponent(2);
	mComponent(0) = theta;
	mComponent(1) = (1 - theta) / 2;
	Eigen::MatrixXf M = toeplitz(mComponent, zerosNMinus3); // @HERE Does this really return a matrix of elementCountMinus elements?  Function returns 25 here
	//float* M = (float*)malloc(elementCountMinus * sizeof(float));


	// A = M+sparse(toeplitz([sig1*k/(h^2)+sig0*k/2 -sig1*k/(2*h^2) zeros(1,N-3)]));
	Eigen::VectorXf aComponent(2);
	aComponent(0) = sig1 * k / (powf(h, 2)) + sig0 * k / 2;
	aComponent(1) = -sig1 * k / (2 * powf(h, 2));
	A = toeplitz(aComponent, zerosNMinus3);

	// A = M + A;
	A = M + A;
	//for (int i = 0; i < elementCountMinus; i++) A[i] += M[i];

	// C = M+sparse(toeplitz([-sig1*k/(h^2)-sig0*k/2 sig1*k/(2*h^2) zeros(1,N-3)]));
	Eigen::VectorXf cComponent(2);
	cComponent(0) = -sig1 * k / (powf(h, 2)) - sig0 * k / 2;
	cComponent(1) = sig1 * k / (2 * powf(h, 2));
	C = toeplitz(cComponent, zerosNMinus3);

	// C = M + C
	C = M + C;
	//for (int i = 0; i < elementCountMinus; i++) C[i] += M[i];

	// float B = 2*M+sparse(toeplitz([-2*lambda^2-6*mu^2 lambda^2+4*mu^2 -mu^2 zeros(1,N-4)]));
	Eigen::VectorXf bComponent(3);
	bComponent(0) = -2 * powf(lambda, 2) - 6 * powf(mu, 2);
	bComponent(1) = powf(lambda, 2) + 4 * powf(mu, 2);
	bComponent(2) = -powf(mu, 2);
	B = toeplitz(bComponent, zerosNMinus4);


	// Create raised cosine
	// xax = [1:N-1]'*h;
	Eigen::VectorXf xax = gen_consecutive_vec(1, N - 1); // maybe needs tanspose here
	Eigen::VectorXf xaxMinusWid = gen_consecutive_vec(1, N - 1);
	Eigen::VectorXf xaxPlusWid = gen_consecutive_vec(1, N - 1);
	Eigen::VectorXf rc = gen_consecutive_vec(1, N - 1);

	for (int i = 0; i < N - 2; i++) {
		xax(i) = xax(i) * h;
		// Copy values into soon-to-be-used matrices
		xaxMinusWid(i) = xax(i) - ctr - wid / 2;
		xaxPlusWid(i) = xax(i) - ctr + wid / 2;
		rc(i) = 1 + cos(2 * M_PI * (xax(i) - ctr) / wid);
	}

	int n_as_int = (int)roundf(N);
	// ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
	Eigen::VectorXf ind(n_as_int - 2);
	// Multiply element-wise (xax-ctr-wid/2).*(xax-ctr+wid/2)
	for (int i = 0; i < N - 2; i++)	ind(i) = -1 * xaxMinusWid(i) * xaxPlusWid(i);

	max_set(ind, 0);
	sign_set(ind);

	// rc = 0.5*ind.*(rc);
	for (int i = 0; i < N - 2; i++) rc(i) = 0.5 + ind(i) * rc(i);

	/*
		u2 = u0*rc;
		u1 = (u0+k*v0)*rc;
		u = zeros(N+1,1);
		out = zeros(NF,2);
	*/
	u2 = Eigen::VectorXf(n_as_int - 2);
	// @HERE Should be u2 *= rc; but check how to do this in Eigen
	u2 = u0 * rc;
	//for (int i = 0; i < N - 2; i++) u2(i) = rc(i) * u0;

	u1 = Eigen::VectorXf(n_as_int - 2);
	u1 = (u0 + k * v0) * rc;
	// @HERE Should be u1 = (u0 + k * v0) * rc; but check how to do this in Eigen
	//for (int i = 0; i < N - 2; i++) u1[i] = (u0 + k * v0) * rc[i]; // @HERE was bug, u1[i]

	u = zeros(N + 1, 1);
	out = zeros(NF, 2);

	std::cout << "SETUP DONE" << std::endl;

	return true;
}



int main(int argc, char** argv)
{
	// B
	// u1
	// C
	// u2
	// A
	// rp_int
	// u
	// rp_frac
	//std::cout << "B: " << B << std::endl;
	//std::cout << "u1 " << u1 << std::endl;
	//std::cout << "C " << C << std::endl;
	//std::cout << "u2 " << u2 << std::endl;
	//std::cout << "A " << A << std::endl;
	//std::cout << "rp_int " << rp_int << std::endl;
	//std::cout << "u " << u << std::endl;
	//std::cout << "rp_frac " << rp_frac << std::endl;
	//Eigen::VectorXf result(-1); // crash
	setup();
	for (unsigned int n = 0; n < 10; n++)
	{
		Eigen::MatrixXf Bmember = B * u1 - C * u2;
		Eigen::VectorXf u = A.lu().solve(Bmember);

		std::cout << "After call to solve" << std::endl << u << std::endl;

		// u(rp_int)'

		// New vector
		// First elements is u((int)roundf(rp_int(0))) second element is u((int)roundf(rp_int(1)))
		int first_index = (int)roundf(rp_int(0));
		int second_index = (int)roundf(rp_int(1));
		Eigen::Vector2f urpint;
		urpint(0) = u(first_index);
		urpint(1) = u(second_index);
		std::cout << "urpint: " << urpint << std::endl;

		auto urpintT = urpint; // tried removing T
		rp_int.array() += 1;

		Eigen::Vector2f urpintplus;
		urpintplus(0) = u(first_index);
		urpintplus(1) = u(second_index);
		std::cout << "After urpintplus" << std::endl << urpintplus;

		auto urpintplusT = urpintplus; // tried removing T
		// out(n, :) = (1 - rp_frac).*u(rp_int)'+rp_frac.*u(rp_int+1)';
		// lookatme = (1-rp_frac).*u(rp_int)'+rp_frac.*u(rp_int+1)';  RETURNS Vector2f

		// a([1, 3]) will select elements 1 and 3 from the array

		// (1-rp_frac)  .*  u(rp_int)'  +  rp_frac  .*  u(rp_int+1)';

		//                                    1         - rp_frac      .*      u(rp_int)'  rp_frac          .*      urpintplusT
		Eigen::Vector2f out = (Eigen::Vector2f::Ones(2) - rp_frac).cwiseProduct(urpintT) + rp_frac.cwiseProduct(urpintplusT); // this is two numbers
		//for (unsigned int channel = 0; context->audioOutChannels; channel++)
		//{
		//	audioWrite(context, n, channel, out(channel % 2));
		//}
		//Eigen::VectorXf result(-1); // crash

		// Update
		u2 = u1;
		u1 = u;
	}
}
