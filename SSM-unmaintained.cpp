//~ Solar System Model
//~ By Behzad Zahedi
//~ behzadxahedi@gmail.com

/*

 The following functions are parts of Vallado ASTIOD code (http://celestrak.org/software/vallado-sw.php) and are slightly adapted to be used in this project:
 - lambertB
 - seebattin
 - kbattin
 - lambhodograph

 To use this code you'd need the following C++ Libraries:
  - boost
  - eigen3

 */


#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen/Dense>
#include <boost/multiprecision/eigen.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <string>
#include <random>
//#include <boost/multiprecision/gmp.hpp>

using namespace std;
using namespace Eigen;
using namespace boost::multiprecision;

//typedef boost::multiprecision::mpf_float_50 gmpf;
typedef boost::multiprecision::cpp_dec_float_50 gmpf;
typedef Matrix<gmpf, Dynamic, Dynamic> MatrixXz;
typedef Matrix<gmpf, 3, 3> Matrix3z;
typedef Matrix<gmpf, 3, 1> Matrix31z;
typedef Matrix<gmpf, 6, 1> Matrix61z;

class Sun {
private:
	const gmpf MU = 132712000000;
	const gmpf R = 696000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
};

class Mercury {
private:
	const gmpf MU = 22030;
	const gmpf R = 2440;
	const gmpf SOI = 112000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Venus {
private:
	const gmpf MU = 324900;
	const gmpf R = 6052;
	const gmpf SOI = 616000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Earth {
private:
	const gmpf MU = 398600.4418;
	const gmpf R = 6378;
	const gmpf SOI = 925000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Moon {
private:
	const gmpf MU = 4903;
	const gmpf R = 1737;
	const gmpf SOI = 66100;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Mars {
private:
	const gmpf MU = 42828;
	const gmpf R = 3396;
	const gmpf SOI = 577000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Jupiter {
private:
	const gmpf MU = 126686000;
	const gmpf R = 71490;
	const gmpf SOI = 48200000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Saturn {
private:
	const gmpf MU = 37931000;
	const gmpf R = 60270;
	const gmpf SOI = 54800000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Uranus {
private:
	const gmpf MU = 5794000;
	const gmpf R = 25560;
	const gmpf SOI = 51800000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Neptune {
private:
	const gmpf MU = 6835100;
	const gmpf R = 24760;
	const gmpf SOI = 86600000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Pluto {
private:
	const gmpf MU = 830;
	const gmpf R = 1195;
	const gmpf SOI = 3080000;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
	gmpf soi() { return SOI; };
};

class Dione {
private:
	const gmpf MU = 73.116;
	const gmpf R = 562.5;
public:
	gmpf mu() { return MU; };
	gmpf r() { return R; };
};

class Uti {
private:
	const gmpf PI = 3.141592653589793;
	const gmpf DTR = PI / 180;
	const gmpf RTD = 180 / PI;
	const gmpf AU2KM = 149597870.7;
	const gmpf SEC2DAY = 1.0 / (60.0 * 60.0 * 24.0);
public:
	gmpf pi() { return PI; };
	gmpf dtr() { return DTR; };
	gmpf rtd() { return RTD; };
	gmpf au2km() { return AU2KM; };
	gmpf km2au(gmpf x) { return x * (1 / AU2KM); };
	gmpf sec2day() { return SEC2DAY; };
	void DispMat(MatrixXz x, string s) {
		//cout << s << ":" << endl << setprecision(numeric_limits<mpf_float_50>::max_digits10) << x << endl << endl;
		cout << s << ":" << endl << setprecision(numeric_limits<cpp_dec_float_50>::max_digits10) << x << endl << endl;
	};
	void DispNum(gmpf x, string s) {
		//cout << s << ":" << endl << setprecision(numeric_limits<mpf_float_50>::max_digits10) << x << endl << endl;
		cout << s << ":" << endl << setprecision(numeric_limits<cpp_dec_float_50>::max_digits10) << x << endl << endl;
	};
	gmpf nrm31(Matrix31z M) {
		gmpf out = sqrt(pow(M(0, 0), 2) + pow(M(1, 0), 2) + pow(M(2, 0), 2));
		return out;
	};
	Matrix3z xRot(gmpf d) {
		Matrix3z out;
		out << 1, 0, 0,
			0, cos(d), sin(d),
			0, -sin(d), cos(d);
		return out;
	};
	Matrix3z yRot(gmpf d) {
		Matrix3z out;
		out << cos(d), 0, -sin(d),
			0, 1, 0,
			sin(d), 0, cos(d);
		return out;
	};
	Matrix3z zRot(gmpf d) {
		Matrix3z out;
		out << cos(d), sin(d), 0,
			-sin(d), cos(d), 0,
			0, 0, 1;
		return out;
	};
	Matrix61z matError(Matrix61z mj, Matrix61z mn) {
		Matrix61z out;
		out(0, 0) = (mj(0, 0) - mn(0, 0)) * (100 / mj(0, 0));
		out(1, 0) = (mj(1, 0) - mn(1, 0)) * (100 / mj(1, 0));
		out(2, 0) = (mj(2, 0) - mn(2, 0)) * (100 / mj(2, 0));
		out(3, 0) = (mj(3, 0) - mn(3, 0)) * (100 / mj(3, 0));
		out(4, 0) = (mj(4, 0) - mn(4, 0)) * (100 / mj(4, 0));
		out(5, 0) = (mj(5, 0) - mn(5, 0)) * (100 / mj(5, 0));
		return out.cwiseAbs();
	};
	Matrix61z mat36(Matrix31z a, Matrix31z b) {
		Matrix61z out;
		out(0, 0) = a(0, 0);
		out(1, 0) = a(1, 0);
		out(2, 0) = a(2, 0);
		out(3, 0) = b(0, 0);
		out(4, 0) = b(1, 0);
		out(5, 0) = b(2, 0);
		return out;
	};
	gmpf it360(gmpf angle) {
		angle = angle * rtd();
		if (angle >= 360.0) {
			angle = angle - floor(angle / 360) * 360;
		}
		else if (angle < 0.0) {
			angle = angle - (ceil(angle / 360) - 1) * 360;
		}
		return angle * dtr();
	};
	gmpf sgn(gmpf x) {
		if (x < 0.0) {
			return -1.0;
		}
		else {
			return 1.0;
		}
	};
	Matrix61z rCOE(Matrix61z a) {
		a(1, 0) = a(1, 0) * rtd();
		a(3, 0) = a(3, 0) * rtd();
		a(4, 0) = a(4, 0) * rtd();
		a(5, 0) = a(5, 0) * rtd();
		return a;
	};
	gmpf angle3d(Matrix31z V1, Matrix31z V2) {
		gmpf ct, mag1, mag2, angle;
		mag1 = nrm31(V1);
		mag2 = nrm31(V2);
		ct = (V1 / mag1).dot(V2 / mag2);
		angle = acos(ct);
		DispNum(angle * rtd(), "Angle");	//in radians
		return angle;
	};
};

class Orbital {
public:
	gmpf ae2h(gmpf a, gmpf e, gmpf mu) {
		gmpf h;
		if (e > 0.0 && e < 1.0) {
			h = sqrt(a * mu * (1.0 - (e * e)));
		}
		if (e > 1.0) {
			h = sqrt(a * mu * ((e * e) - 1.0));
		}
		return h;
	};
	Matrix3z T_P2I(gmpf RAAN, gmpf inc, gmpf AOP) {
		Uti u;
		Matrix3z T_I2P;
		T_I2P = (u.zRot(AOP) * u.xRot(inc)) * u.zRot(RAAN);
		return T_I2P.transpose();
	};
	Matrix61z SV2COE(Matrix61z sv, gmpf mu) {
		Uti u;
		Matrix31z R, V, H, N, tmpZ, E;
		R << sv(0, 0), sv(1, 0), sv(2, 0);
		V << sv(3, 0), sv(4, 0), sv(5, 0);
		tmpZ << 0, 0, 1;
		gmpf r = u.nrm31(R);
		gmpf Vr = R.dot(V) / r;
		H = R.cross(V);
		gmpf h = u.nrm31(H);
		gmpf inc = acos(H(2, 0) / h);
		N = tmpZ.cross(H);
		gmpf n = u.nrm31(N);
		gmpf RAAN;
		if (N(1, 0) >= 0) {
			RAAN = acos(N(0, 0) / n);
		}
		else {
			RAAN = (2 * u.pi()) - acos(N(0, 0) / n);
		}
		E = (1 / mu) * ((V.cross(H)) - (mu * (R / r)));
		gmpf e = u.nrm31(E);
		gmpf AOP;
		if (E(2, 0) >= 0) {
			AOP = acos((N / n).dot(E / e));
		}
		else {
			AOP = (2 * u.pi()) - acos((N / n).dot(E / e));
		}
		gmpf f;
		if (Vr >= 0) {
			f = acos((E / e).dot(R / r));
		}
		else {
			f = (2 * u.pi()) - acos((E / e).dot(R / r));
		}
		gmpf a = ((h * h) / mu) * (1 / (1 - (e * e)));
		//~ ADD a Test Thingy HERE!
		Matrix61z out;
		out << a, inc, e, RAAN, AOP, f;
		return out;
	};
	Matrix61z COE2SV(const Matrix61z coe, gmpf mu) {
		Uti u;
		const gmpf a = coe(0, 0);
		const gmpf inc = coe(1, 0);
		const gmpf e = coe(2, 0);
		const gmpf RAAN = coe(3, 0);
		const gmpf AOP = coe(4, 0);
		const gmpf f = coe(5, 0);
		const gmpf h = ae2h(a, e, mu);
		Matrix31z RP(3, 1);
		RP << cos(f), sin(f), 0;
		RP = ((h * h) / mu) * (1 / (1 + (e * cos(f)))) * RP;
		Matrix31z VP(3, 1);
		VP << -sin(f), (e + cos(f)), 0;
		VP = (mu / h) * VP;
		Matrix3z TP2I = T_P2I(RAAN, inc, AOP);
		Matrix31z RI, VI;
		RI = TP2I * RP;
		VI = TP2I * VP;
		Matrix61z out = u.mat36(RI, VI);
		return out;
	};
	gmpf Me2E(gmpf M, gmpf e) {
		Uti u;
		gmpf E, ratio, FE, FPE;
		if (M <= u.pi()) {
			E = M + (e / 2);
		}
		else if (M > u.pi()) {
			E = M - (e / 2);
		}
		ratio = 1;
		while (abs(ratio) > 1e-20) {
			FE = E - (e * sin(E)) - M;
			FPE = 1 - (e * cos(E));
			ratio = FE / FPE;
			E = E - ratio;
		}
		return E;
	};
	Matrix61z Eph(int Planet, gmpf JD) {
		Uti u;
		Sun sun;
		gmpf a, adot, e, edot, i, idot, RAAN, RAANdot, LOP, LOPdot, ML, MLdot, T0;
		switch (Planet) {
			//Mercury DE405
		case 1:
			a = 0.38709927;
			adot = 0.00000037;
			e = 0.20563593;
			edot = 0.00001906;
			i = 7.00497902;
			idot = -0.00590158;
			RAAN = 48.33076593;
			RAANdot = -0.12534081;
			LOP = 77.45779628;
			LOPdot = 0.16047689;
			ML = 252.25032350;
			MLdot = 149472.67411175;
			break;
			//Venus DE405
		case 2:
			a = 0.72333566;
			adot = 0.00000390;
			e = 0.00677672;
			edot = -0.00004107;
			i = 3.39467605;
			idot = -0.00078890;
			RAAN = 76.67984255;
			RAANdot = -0.27769418;
			LOP = 131.60246718;
			LOPdot = 0.00268329;
			ML = 181.97909950;
			MLdot = 58517.81538729;
			break;
			//EMB DE405
		case 3:
			a = 1.00000261;
			adot = 0.00000562;
			e = 0.01671123;
			edot = -0.00004392;
			i = -0.00001531;
			idot = -0.01294668;
			RAAN = -5.11260389;
			RAANdot = -0.24123856;
			LOP = 102.93768193;
			LOPdot = 0.32327364;
			ML = 100.46457166;
			MLdot = 35999.37244981;
			break;
			//Mars DE405
		case 4:
			a = 1.52371034;
			adot = 0.00001847;
			e = 0.09339410;
			edot = 0.00007882;
			i = 1.84969142;
			idot = -0.00813131;
			RAAN = 49.55953891;
			RAANdot = -0.29257343;
			LOP = -23.94362959;
			LOPdot = 0.44441088;
			ML = -4.55343205;
			MLdot = 19140.30268499;
			break;
			//Jupiter DE405
		case 5:
			a = 5.20288700;
			adot = -0.00011607;
			e = 0.04838624;
			edot = -0.00013253;
			i = 1.30439695;
			idot = -0.00183714;
			RAAN = 100.47390909;
			RAANdot = 0.20469106;
			LOP = 14.72847983;
			LOPdot = 0.21252668;
			ML = 34.39644051;
			MLdot = 3034.74612775;
			break;
			//Saturn DE405
		case 6:
			a = 9.53667594;
			adot = -0.00125060;
			e = 0.05386179;
			edot = -0.00050991;
			i = 2.48599187;
			idot = 0.00193609;
			RAAN = 113.66242448;
			RAANdot = -0.28867794;
			LOP = 92.59887831;
			LOPdot = -0.41897216;
			ML = 49.95424423;
			MLdot = 1222.49362201;
			break;
			//Uranus DE405
		case 7:
			a = 19.18916464;
			adot = -0.00196176;
			e = 0.04725744;
			edot = -0.00004397;
			i = 0.77263783;
			idot = -0.00242939;
			RAAN = 74.01692503;
			RAANdot = 0.04240589;
			LOP = 170.95427630;
			LOPdot = 0.40805281;
			ML = 313.23810451;
			MLdot = 428.48202785;
			break;
			//Neptune DE405
		case 8:
			a = 30.06992276;
			adot = 0.00026291;
			e = 0.00859048;
			edot = 0.00005105;
			i = 1.77004347;
			idot = 0.00035372;
			RAAN = 131.78422574;
			RAANdot = -0.00508664;
			LOP = 44.96476227;
			LOPdot = -0.32241464;
			ML = -55.12002969;
			MLdot = 218.45945325;
			break;
			//Pluto DE405
		case 9:
			a = 39.48211675;
			adot = -0.00031596;
			e = 0.24882730;
			edot = 0.00005170;
			i = 17.14001206;
			idot = 0.00004818;
			RAAN = 110.30393684;
			RAANdot = -0.01183482;
			LOP = 224.06891629;
			LOPdot = -0.04062942;
			ML = 238.92903833;
			MLdot = 145.20780515;
			break;
		default:
			cout << "Invalid Planet!" << endl << endl;
		}
		Matrix61z Ele0, Rates, Ele, SV, COE;
		Ele0(0, 0) = a * u.au2km();
		Ele0(1, 0) = e;
		Ele0(2, 0) = i;
		Ele0(3, 0) = RAAN;
		Ele0(4, 0) = LOP;
		Ele0(5, 0) = ML;
		Rates(0, 0) = adot * u.au2km();
		Rates(1, 0) = edot;
		Rates(2, 0) = idot;
		Rates(3, 0) = RAANdot;
		Rates(4, 0) = LOPdot;
		Rates(5, 0) = MLdot;

		T0 = (JD - 2451545.0) / (36525.0);
		Ele = Ele0 + (Rates * T0);
		Ele(2, 0) = u.it360(Ele(2, 0) * u.dtr());
		Ele(3, 0) = u.it360(Ele(3, 0) * u.dtr());
		Ele(4, 0) = u.it360(Ele(4, 0) * u.dtr());
		Ele(5, 0) = u.it360(Ele(5, 0) * u.dtr());

		//~ Matrix61z zzz = Ele;
		//~ zzz(2,0) *= u.rtd();
		//~ zzz(3,0) *= u.rtd();
		//~ zzz(4,0) *= u.rtd();
		//~ zzz(5,0) *= u.rtd();
		//~ u.DispMat(zzz, "Planet Elements");

		gmpf AOP, MA, E, TA;
		AOP = u.it360(Ele(4, 0) - Ele(3, 0));
		MA = u.it360(Ele(5, 0) - Ele(4, 0));
		E = Me2E(MA, Ele(1, 0));
		TA = u.it360(2 * atan(sqrt((1 + Ele(1, 0)) / (1 - Ele(1, 0))) * tan(E / 2)));
		COE << Ele(0, 0), Ele(2, 0), Ele(1, 0), Ele(3, 0), AOP, TA;
		SV = COE2SV(COE, sun.mu());
		return SV;
	};
	gmpf JD(Matrix61z Time) {
		gmpf J0, UT, JD, T0;
		J0 = 367.0 * Time(0, 0) - floor(7 * (Time(0.0) + floor((Time(1, 0) + 9.0) / 12.0)) / 4.0) + floor(275.0 * Time(1, 0) / 9.0) + Time(2, 0) + 1721013.5;
		UT = Time(3, 0) + (Time(4, 0) / 60.0) + (Time(5, 0) / 3600.0);
		JD = J0 + (UT / 24.0);
		return JD;
	};
};

class Trajectory {
public:
	Matrix61z Govern(Matrix61z M, gmpf mu) {
		Uti u;
		Matrix61z out(6, 1);
		gmpf rxDot = M(3, 0);
		gmpf ryDot = M(4, 0);
		gmpf rzDot = M(5, 0);
		gmpf r = sqrt(pow(M(0, 0), 2) + pow(M(1, 0), 2) + pow(M(2, 0), 2));
		gmpf rxDDot = -(mu / pow(r, 3)) * M(0, 0);
		gmpf ryDDot = -(mu / pow(r, 3)) * M(1, 0);
		gmpf rzDDot = -(mu / pow(r, 3)) * M(2, 0);
		out << rxDot, ryDot, rzDot, rxDDot, ryDDot, rzDDot;
		return out;
	};
	Matrix61z RK4(Matrix61z in, gmpf t0, gmpf t, gmpf stp, gmpf mu, string SVName, string COEName) {
		Uti u;
		Orbital o;
		int n = (int)((t - t0) / stp);
		Matrix61z k1, k2, k3, k4, out, COETemp;
		ofstream SVFile;
		ofstream COEFile;
		SVFile.precision(25);
		COEFile.precision(25);
		SVFile.open(SVName);
		COEFile.open(COEName);
		SVFile << t0 << " " << in(0, 0) << " " << in(1, 0) << " " << in(2, 0) << " " << in(3, 0) << " " << in(4, 0) << " " << in(5, 0) << endl;
		COETemp = o.SV2COE(in, mu);
		COEFile << t0 << " " << COETemp(0, 0) << " " << COETemp(1, 0) * u.rtd() << " " << COETemp(2, 0) << " " << COETemp(3, 0) * u.rtd() << " " << COETemp(4, 0) * u.rtd() << " " << COETemp(5, 0) * u.rtd() << endl;
		for (int i = 1.0; i <= n; i++) {
			t0 = t0 + stp;
			k1 = stp * Govern(in, mu);
			k2 = stp * Govern(in + 0.5 * k1, mu);
			k3 = stp * Govern(in + 0.5 * k2, mu);
			k4 = stp * Govern(in + k3, mu);
			out = in + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
			SVFile << t0 << " " << out(0, 0) << " " << out(1, 0) << " " << out(2, 0) << " " << out(3, 0) << " " << out(4, 0) << " " << out(5, 0) << endl;
			COETemp = o.SV2COE(out, mu);
			COEFile << t0 << " " << COETemp(0, 0) << " " << COETemp(1, 0) * u.rtd() << " " << COETemp(2, 0) << " " << COETemp(3, 0) * u.rtd() << " " << COETemp(4, 0) * u.rtd() << " " << COETemp(5, 0) * u.rtd() << endl;
			in = out;
		}
		SVFile.close();
		COEFile.close();
		return out;
	};
	Matrix61z Governp(Matrix61z M, gmpf mu, gmpf JD) {
		Uti u;
		Orbital o;
		Mercury mercury;
		Venus venus;
		Earth earth;
		Mars mars;
		Jupiter jupiter;
		Saturn saturn;
		Uranus uranus;
		Neptune neptune;
		Pluto pluto;

		Matrix61z out, SV1, SV2, SV3, SV4, SV5, SV6, SV7, SV8, SV9;
		Matrix31z R, SRPa, R1, R2, R3, R4, R5, R6, R7, R8, R9, TPPT, TPP1, TPP2, TPP3, TPP4, TPP5, TPP6, TPP7, TPP8, TPP9;

		R(0) = M(0);
		R(1) = M(1);
		R(2) = M(2);

		SRPa = SRP(R);

		SV1 = o.Eph(1, JD);
		R1(0) = SV1(0);
		R1(1) = SV1(1);
		R1(1) = SV1(2);
		TPP1 = TPP(R, R1, mercury.mu());
		//~ u.DispNum(u.nrm31(TPP1),"1p");

		SV2 = o.Eph(2, JD);
		R2(0) = SV2(0);
		R2(1) = SV2(1);
		R2(1) = SV2(2);
		TPP2 = TPP(R, R2, venus.mu());
		//~ u.DispNum(u.nrm31(TPP2),"2p");

		SV3 = o.Eph(3, JD);
		R3(0) = SV3(0);
		R3(1) = SV3(1);
		R3(1) = SV3(2);
		TPP3 = TPP(R, R3, earth.mu());
		//~ u.DispNum(u.nrm31(TPP3),"3p");

		SV4 = o.Eph(4, JD);
		R4(0) = SV4(0);
		R4(1) = SV4(1);
		R4(1) = SV4(2);
		TPP4 = TPP(R, R4, mars.mu());
		//~ u.DispNum(u.nrm31(TPP4),"4p");

		SV5 = o.Eph(5, JD);
		R5(0) = SV5(0);
		R5(1) = SV5(1);
		R5(1) = SV5(2);
		TPP5 = TPP(R, R5, jupiter.mu());
		//~ u.DispNum(u.nrm31(TPP5),"5p");

		SV6 = o.Eph(6, JD);
		R6(0) = SV6(0);
		R6(1) = SV6(1);
		R6(1) = SV6(2);
		TPP6 = TPP(R, R6, saturn.mu());
		//~ u.DispNum(u.nrm31(TPP6),"6p");

		SV7 = o.Eph(7, JD);
		R7(0) = SV7(0);
		R7(1) = SV7(1);
		R7(1) = SV7(2);
		TPP7 = TPP(R, R7, uranus.mu());
		//~ u.DispNum(u.nrm31(TPP7),"7p");

		SV8 = o.Eph(8, JD);
		R8(0) = SV8(0);
		R8(1) = SV8(1);
		R8(1) = SV8(2);
		TPP8 = TPP(R, R8, neptune.mu());
		//~ u.DispNum(u.nrm31(TPP8),"8p");

		SV9 = o.Eph(9, JD);
		R9(0) = SV9(0);
		R9(1) = SV9(1);
		R9(1) = SV9(2);
		TPP9 = TPP(R, R9, pluto.mu());
		//~ u.DispNum(u.nrm31(TPP9),"9p");

	//~ TPPT << 0, 0, 0;
	//~ TPPT = TPP1;
	//~ TPPT = TPP1 + TPP2;
	//~ TPPT = TPP1 + TPP2 + TPP3;
	//~ TPPT = TPP1 + TPP2 + TPP3 + TPP4;
	//~ TPPT = TPP1 + TPP2 + TPP3 + TPP4 + TPP5;
	//~ TPPT = TPP1 + TPP2 + TPP3 + TPP4 + TPP5 + TPP6;
	//~ TPPT = TPP1 + TPP2 + TPP3 + TPP4 + TPP5 + TPP6 + TPP7;
	//~ TPPT = TPP1 + TPP2 + TPP3 + TPP4 + TPP5 + TPP6 + TPP7 + TPP8;
	//~ TPPT = TPP1 + TPP2 + TPP3 + TPP4 + TPP5 + TPP6 + TPP7 + TPP8 + TPP9;
		TPPT = TPP1 + TPP2 + TPP3 + TPP4 + TPP5 + TPP6 + TPP7 + TPP8 + TPP9 + SRPa;

		gmpf rxDot = M(3, 0);
		gmpf ryDot = M(4, 0);
		gmpf rzDot = M(5, 0);
		gmpf r = sqrt(pow(M(0, 0), 2) + pow(M(1, 0), 2) + pow(M(2, 0), 2));
		gmpf rxDDot = -(mu / pow(r, 3)) * M(0, 0) + TPPT(0);
		gmpf ryDDot = -(mu / pow(r, 3)) * M(1, 0) + TPPT(1);
		gmpf rzDDot = -(mu / pow(r, 3)) * M(2, 0) + TPPT(2);
		out << rxDot, ryDot, rzDot, rxDDot, ryDDot, rzDDot;
		return out;
	};
	Matrix31z SRP(Matrix31z Rh) {			// Solar Radiation Pressure P.695 Cutris
		Uti u;															// R heliocentric in kms
		gmpf nu = 1.0; 													//always in sunlight
		gmpf S0 = 63.15e12; 												//Solar Constant W/Km2
		gmpf c = 2.988e5; 												//speed of light  Km/s
		//~ gmpf Rsc = 0.0079; 													//s/c radius km
		//~ gmpf Asc = u.pi()*Rsc*Rsc;											// area of sphere km2 
		//~ gmpf m = 100.0; 													//s/c mass kg
		gmpf Apm = 1.0e-6;
		gmpf Cr = 1.3;										// 1 for black body, 2 for 100% reflection
		Matrix31z p = nu * (S0 * pow(696.0 / u.nrm31(Rh), 2.0)) * Cr * (1.0 / c) * Apm * (Rh / u.nrm31(Rh));
		return p;
	};
	Matrix31z TPP(Matrix31z R, Matrix31z Rp, gmpf MUp) {		// Third planet perurbation P.713 Cutris
		Uti u;															// R heliocentric in kms, Rp Planet position
		Matrix31z Rrel = Rp - R; 										// R relative
		Matrix31z p = MUp * ((Rrel) / pow(u.nrm31(Rrel), 3) - (Rp) / pow(u.nrm31(Rp), 3));
		return p;
	};
	Matrix61z RK4p(Matrix61z in, gmpf t0, gmpf t, gmpf stp, gmpf mu, gmpf JD) {
		Uti u;
		Orbital o;
		int n = (int)((t - t0) / stp);
		//~ cout << "Number of Iterations: " << n << endl;
		//~ cout << "This could take a while depending on specified step size." << endl << endl;
		cout << "Two Body Function With Perturbations Running.(" << n << " iterations)" << endl << endl;
		Matrix61z k1, k2, k3, k4, out, COETemp;
		ofstream SVpFile;
		ofstream COEpFile;
		SVpFile.precision(25);
		COEpFile.precision(25);
		SVpFile.open("SVp.dat");
		COEpFile.open("COEp.dat");
		SVpFile << t0 << " " << in(0, 0) << " " << in(1, 0) << " " << in(2, 0) << " " << in(3, 0) << " " << in(4, 0) << " " << in(5, 0) << endl;
		COETemp = o.SV2COE(in, mu);
		COEpFile << t0 << " " << COETemp(0, 0) << " " << COETemp(1, 0) * u.rtd() << " " << COETemp(2, 0) << " " << COETemp(3, 0) * u.rtd() << " " << COETemp(4, 0) * u.rtd() << " " << COETemp(5, 0) * u.rtd() << endl;
		for (int i = 1.0; i <= n; i++) {
			//~ cout << "Progress: " << (i*100)/n << "%" << endl;
			t0 = t0 + stp;
			JD = JD + stp * u.sec2day();
			k1 = stp * Governp(in, mu, JD);
			k2 = stp * Governp(in + 0.5 * k1, mu, JD);
			k3 = stp * Governp(in + 0.5 * k2, mu, JD);
			k4 = stp * Governp(in + k3, mu, JD);
			out = in + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
			SVpFile << t0 << " " << out(0, 0) << " " << out(1, 0) << " " << out(2, 0) << " " << out(3, 0) << " " << out(4, 0) << " " << out(5, 0) << endl;
			COETemp = o.SV2COE(out, mu);
			COEpFile << t0 << " " << COETemp(0, 0) << " " << COETemp(1, 0) * u.rtd() << " " << COETemp(2, 0) << " " << COETemp(3, 0) * u.rtd() << " " << COETemp(4, 0) * u.rtd() << " " << COETemp(5, 0) * u.rtd() << endl;
			in = out;
		}
		SVpFile.close();
		COEpFile.close();
		//~ cout << endl;
		return out;
	};
};

class Lambert {
public:
	Matrix61z lambertG(Matrix31z R1, Matrix31z R2, gmpf t, gmpf mu) {
		//Gauss Method
		Uti u;
		gmpf r1 = u.nrm31(R1);
		gmpf r2 = u.nrm31(R2);
		gmpf cosDF = R1.dot(R2) / (r1 * r2);
		gmpf sinDF = 1 * sqrt(1 - (cosDF * cosDF));
		gmpf L = (r1 + r2) / (4 * sqrt(r1 * r2) * cos(acos(cosDF) / 2)) - 0.5;
		gmpf m = (mu * t * t) / pow(2 * sqrt(r1 * r2) * cos(acos(cosDF) / 2), 3);
		gmpf y = 1;
		int n = 0;
		gmpf x1, x2;
		//~ NEEDS PROPER STOPPER HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		while (n <= 50) {
			n++;
			x1 = m / y / y - L;
			x2 = (4.0 / 3.0) * (1 + (6.0 * x1 / 5.0) + (48.0 * x1 * x1 / 35.0) + (480.0 * x1 * x1 * x1 / 315.0));
			y = 1 + x2 * (L + x1);
		}
		gmpf p = (r1 * r2 * (1 - cosDF)) / (r1 + r2 - (2 * sqrt(r1 * r2) * cos(acos(cosDF) / 2) * (1 - (2 * x1))));
		gmpf f = 1 - (r2 / p) * (1 - cosDF);
		gmpf g = (r1 * r2 * sinDF) / sqrt(mu * p);
		gmpf gdot = 1 - (r1 / p) * (1 - cosDF);
		Matrix31z V1 = (R2 - f * R1) / g;
		Matrix31z V2 = (gdot * R2 - R1) / g;
		Matrix61z out = u.mat36(V1, V2);
		return out;
	};
	Matrix61z lambhodograph(Matrix31z r1, Matrix31z v1, Matrix31z r2, gmpf p, gmpf ecc, gmpf dnu, gmpf dtsec, gmpf mu) {
		Uti u;

		Matrix31z nvec, rcrv, rcrr, v1t, v2t;
		gmpf eps, magr1, magr2, a, b, x1, x2, y2a, y2b, ptx, pi, twopi;
		pi = u.pi();
		twopi = 2 * u.pi();
		int i;

		eps = 1.0e-8;  // -14

		magr1 = u.nrm31(r1);
		magr2 = u.nrm31(r2);

		a = mu * (1.0 / magr1 - 1.0 / p);  // not the semi - major axis
		b = pow(mu * ecc / p, 2) - a * a;
		if (b <= 0.0)
			x1 = 0.0;
		else
			x1 = -sqrt(b);

		// 180 deg, and multiple 180 deg transfers
		if (fabs(sin(dnu)) < eps)
		{
			rcrv = r1.cross(v1);
			for (i = 0; i < 3; i++)
				nvec[i] = rcrv[i] / u.nrm31(rcrv);
			if (ecc < 1.0)
			{
				ptx = twopi * sqrt(p * p * p / pow(mu * (1.0 - ecc * ecc), 3));
				if (fmod(dtsec, ptx) > ptx * 0.5)
					x1 = -x1;
			}
		}
		else
		{
			y2a = mu / p - x1 * sin(dnu) + a * cos(dnu);
			y2b = mu / p + x1 * sin(dnu) + a * cos(dnu);
			if (fabs(mu / magr2 - y2b) < fabs(mu / magr2 - y2a))
				x1 = -x1;

			// depending on the cross product, this will be normal or in plane,
			// or could even be a fan
			rcrr = r1.cross(r2);
			for (i = 0; i < 3; i++)
				nvec[i] = rcrr[i] / u.nrm31(rcrr); // if this is r1, v1, the transfer is coplanar!
			if (fmod(dnu, twopi) > pi)
			{
				for (i = 0; i < 3; i++)
					nvec[i] = -nvec[i];
			}
		}

		rcrv = nvec.cross(r1);
		rcrr = nvec.cross(r2);
		x2 = x1 * cos(dnu) + a * sin(dnu);
		for (i = 0; i < 3; i++)
		{
			v1t[i] = (sqrt(mu * p) / magr1) * ((x1 / mu) * r1[i] + rcrv[i]) / magr1;
			v2t[i] = (sqrt(mu * p) / magr2) * ((x2 / mu) * r2[i] + rcrr[i]) / magr2;
		}
		return u.mat36(v1t, v2t);
	};
	static gmpf kbattin(gmpf v) {
		gmpf d[21] = {
			1.0 / 3.0, 4.0 / 27.0,
			8.0 / 27.0, 2.0 / 9.0,
			22.0 / 81.0, 208.0 / 891.0,
			340.0 / 1287.0, 418.0 / 1755.0,
			598.0 / 2295.0, 700.0 / 2907.0,
			928.0 / 3591.0, 1054.0 / 4347.0,
			1330.0 / 5175.0, 1480.0 / 6075.0,
			1804.0 / 7047.0, 1978.0 / 8091.0,
			2350.0 / 9207.0, 2548.0 / 10395.0,
			2968.0 / 11655.0, 3190.0 / 12987.0,
			3658.0 / 14391.0
		};
		gmpf del, delold, term, termold, sum1;
		int i;

		/* ---- process forwards ---- */
		sum1 = d[0];
		delold = 1.0;
		termold = d[0];
		i = 1;
		while ((i <= 20) && (fabs(termold) > 0.00000001))
		{
			del = 1.0 / (1.0 - d[i] * v * delold);
			term = termold * (del - 1.0);
			sum1 = sum1 + term;
			i++;
			delold = del;
			termold = term;
		}
		//return sum1;

		int ktr = 20;
		gmpf sum2 = 0.0;
		gmpf term2 = 1.0 + d[ktr] * v;
		for (i = 1; i <= ktr - 1; i++)
		{
			sum2 = d[ktr - i] * v / term2;
			term2 = 1.0 + sum2;
		}

		return (d[0] / term2);
	};
	static gmpf seebattin(gmpf v2) {
		gmpf c[21]{
			0.2,
			9.0 / 35.0, 16.0 / 63.0,
			25.0 / 99.0, 36.0 / 143.0,
			49.0 / 195.0, 64.0 / 255.0,
			81.0 / 323.0, 100.0 / 399.0,
			121.0 / 483.0, 144.0 / 575.0,
			169.0 / 675.0, 196.0 / 783.0,
			225.0 / 899.0, 256.0 / 1023.0,
			289.0 / 1155.0, 324.0 / 1295.0,
			361.0 / 1443.0, 400.0 / 1599.0,
			441.0 / 1763.0, 484.0 / 1935.0
		};
		// first term is diff, indices are offset too
		gmpf c1[20]{
			9.0 / 7.0, 16.0 / 63.0,
			25.0 / 99.0, 36.0 / 143.0,
			49.0 / 195.0, 64.0 / 255.0,
			81.0 / 323.0, 100.0 / 399.0,
			121.0 / 483.0, 144.0 / 575.0,
			169.0 / 675.0, 196.0 / 783.0,
			225.0 / 899.0, 256.0 / 1023.0,
			289.0 / 1155.0, 324.0 / 1295.0,
			361.0 / 1443.0, 400.0 / 1599.0,
			441.0 / 1763.0, 484.0 / 1935.0
		};

		gmpf term, termold, del, delold, sum1, eta, sqrtopv;
		int i;

		sqrtopv = sqrt(1.0 + v2);
		eta = v2 / pow(1.0 + sqrtopv, 2);

		/* ---- process forwards ---- */
		delold = 1.0;
		termold = c[0];  // * eta
		sum1 = termold;
		i = 1;
		while ((i <= 20) && (fabs(termold) > 0.000001)) {
			del = 1.0 / (1.0 + c[i] * eta * delold);
			term = termold * (del - 1.0);
			sum1 = sum1 + term;
			i++;
			delold = del;
			termold = term;
		}

		//   return ((1.0 / (8.0 * (1.0 + sqrtopv))) * (3.0 + sum1 / (1.0 + eta * sum1)));
		gmpf seebatt = 1.0 / ((1.0 / (8.0 * (1.0 + sqrtopv))) * (3.0 + sum1 / (1.0 + eta * sum1)));

		int ktr = 19;
		gmpf sum2 = 0.0;
		gmpf term2 = 1.0 + c1[ktr] * eta;
		for (i = 0; i <= ktr - 1; i++)
		{
			sum2 = c1[ktr - i] * eta / term2;
			term2 = 1.0 + sum2;
		}

		return (8.0 * (1.0 + sqrtopv) /
			(3.0 +
				(1.0 /
					(5.0 + eta + ((9.0 / 7.0) * eta / term2)))));
	};
	Matrix61z lambertB(Matrix31z r1, Matrix31z r2, Matrix31z v1, char dm, char df, int nrev, gmpf dtsec, const gmpf mu) {
		//Batin Method
		Uti uu;
		gmpf pi = uu.pi();
		const gmpf small = 0.000001;

		Matrix31z rcrossr, v1dvl, v2dvl, v2;
		int i, loops;

		Matrix31z v1t, v2t;
		Matrix61z v1tv2t;

		gmpf   u, b, x, xn, y, L, m, cosdeltanu, sindeltanu, dnu, a,
			ror, h1, h2, tempx, eps, denom, chord, k2, s,
			p, ecc, f, A, y1, bigt;
		gmpf magr1, magr2, magrcrossr, lam, temp, temp1, temp2;

		//~ error = 0;
		magr1 = uu.nrm31(r1);
		magr2 = uu.nrm31(r2);

		cosdeltanu = r1.dot(r2) / (magr1 * magr2);
		// make sure it's not more than 1.0
		if (abs(cosdeltanu) > 1.0) {
			cosdeltanu = 1.0 * uu.sgn(cosdeltanu);
		}
		rcrossr = r1.cross(r2);
		magrcrossr = uu.nrm31(rcrossr);
		if (dm == 's')
			sindeltanu = magrcrossr / (magr1 * magr2);
		else
			sindeltanu = -magrcrossr / (magr1 * magr2);

		dnu = atan2(sindeltanu, cosdeltanu);
		// the angle needs to be positive to work for the long way
		if (dnu < 0.0)
			dnu = 2.0 * pi + dnu;

		// these are the same
		chord = sqrt(magr1 * magr1 + magr2 * magr2 - 2.0 * magr1 * magr2 * cosdeltanu);
		//chord = mag(r2 - r1);

		s = (magr1 + magr2 + chord) * 0.5;
		ror = magr2 / magr1;
		eps = ror - 1.0;

		lam = 1.0 / s * sqrt(magr1 * magr2) * cos(dnu * 0.5);
		L = pow((1.0 - lam) / (1.0 + lam), 2);
		m = 8.0 * mu * dtsec * dtsec / (s * s * s * pow(1.0 + lam, 6));

		// initial guess
		if (nrev > 0)
			xn = 1.0 + 4.0 * L;
		else
			xn = L;   //l    // 0.0 for par and hyp, l for ell

		// alt approach for high energy(long way, retro multi - rev) case
		if ((df == 'r') && (nrev > 0))
		{
			xn = 1e-20;  // be sure to reset this here!!
			x = 10.0;  // starting value
			loops = 1;
			while ((abs(xn - x) >= small) && (loops <= 20))
			{
				x = xn;
				temp = 1.0 / (2.0 * (L - x * x));
				temp1 = sqrt(x);
				temp2 = (nrev * pi * 0.5 + atan(temp1)) / temp1;
				h1 = temp * (L + x) * (1.0 + 2.0 * x + L);
				h2 = temp * m * temp1 * ((L - x * x) * temp2 - (L + x));

				b = 0.25 * 27.0 * h2 / (pow(temp1 * (1.0 + h1), 3));
				if (b < -1.0) // reset the initial condition
					f = 2.0 * cos(1.0 / 3.0 * acos(sqrt(b + 1.0)));
				else
				{
					A = pow(sqrt(b) + sqrt(b + 1.0), (1.0 / 3.0));
					f = A + 1.0 / A;
				}

				y = 2.0 / 3.0 * temp1 * (1.0 + h1) * (sqrt(b + 1.0) / f + 1.0);
				xn = 0.5 * ((m / (y * y) - (1.0 + L)) - sqrt(pow(m / (y * y) - (1.0 + L), 2) - 4.0 * L));
				//~ HERE! fprintf(outfile," %3i yh %11.6f x %11.6f h1 %11.6f h2 %11.6f b %11.6f f %11.7f \n", loops, y, x, h1, h2, b, f);
				loops = loops + 1;
			}  // while
			x = xn;
			a = s * pow(1.0 + lam, 2) * (1.0 + x) * (L + x) / (8.0 * x);
			p = (2.0 * magr1 * magr2 * (1.0 + x) * pow(sin(dnu * 0.5), 2)) / (s * pow(1.0 + lam, 2) * (L + x));  // thompson
			ecc = sqrt(1.0 - p / a);
			v1tv2t = lambhodograph(r1, v1, r2, p, ecc, dnu, dtsec, mu);

			v1t << v1tv2t[0], v1tv2t[1], v1tv2t[2];
			v2t << v1tv2t[3], v1tv2t[4], v1tv2t[5];

			//~ HERE! fprintf(outfile,"high v1t %16.8f %16.8f %16.8f %16.8f\n", v1t, uu.nrm31(v1t));
		}
		else
		{
			// standard processing
			// note that the dr nrev = 0 case is not represented
			loops = 1;
			y1 = 0.0;
			x = 10.0;  // starting value
			while ((abs(xn - x) >= small) && (loops <= 30))
			{
				if (nrev > 0)
				{
					x = xn;
					temp = 1.0 / ((1.0 + 2.0 * x + L) * (4.0 * x));
					temp1 = (nrev * pi * 0.5 + atan(sqrt(x))) / sqrt(x);
					h1 = temp * pow(L + x, 2) * (3.0 * pow(1.0 + x, 2) * temp1 - (3.0 + 5.0 * x));
					h2 = temp * m * ((x * x - x * (1.0 + L) - 3.0 * L) * temp1 + (3.0 * L + x));
				}
				else
				{
					x = xn;
					tempx = seebattin(x);
					denom = 1.0 / ((1.0 + 2.0 * x + L) * (4.0 * x + tempx * (3.0 + x)));
					h1 = pow(L + x, 2) * (1.0 + 3.0 * x + tempx) * denom;
					h2 = m * (x - L + tempx) * denom;
				}

				// ---------------------- - evaluate cubic------------------
				b = 0.25 * 27.0 * h2 / (pow(1.0 + h1, 3));

				u = 0.5 * b / (1.0 + sqrt(1.0 + b));
				k2 = kbattin(u);
				y = ((1.0 + h1) / 3.0) * (2.0 + sqrt(1.0 + b) / (1.0 + 2.0 * u * k2 * k2));
				xn = sqrt(pow((1.0 - L) * 0.5, 2) + m / (y * y)) - (1.0 + L) * 0.5;

				y1 = sqrt(m / ((L + x) * (1.0 + x)));
				loops = loops + 1;
				//        fprintf(1, ' %3i yb %11.6f x %11.6f k2 %11.6f b %11.6f u %11.6f y1 %11.7f \n', loops, y, x, k2, b, u, y1);
			}  // while

		}

		if (loops < 30)
		{
			// blair approach use y from solution
			//       lam = 1.0 / s * sqrt(magr1*magr2) * cos(dnu*0.5);
			//       m = 8.0*mu*dtsec*dtsec / (s ^ 3 * (1.0 + lam) ^ 6);
			//       L = ((1.0 - lam) / (1.0 + lam)) ^ 2;
			//a = s*(1.0 + lam) ^ 2 * (1.0 + x)*(lam + x) / (8.0*x);
			// p = (2.0*magr1*magr2*(1.0 + x)*sin(dnu*0.5) ^ 2) ^ 2 / (s*(1 + lam) ^ 2 * (lam + x));  % loechler, not right ?
			p = (2.0 * magr1 * magr2 * y * y * pow(1.0 + x, 2) * pow(sin(dnu * 0.5), 2)) / (m * s * pow(1.0 + lam, 2));  // thompson
			ecc = sqrt((eps * eps + 4.0 * magr2 / magr1 * pow(sin(dnu * 0.5), 2) * pow((L - x) / (L + x), 2)) / (eps * eps + 4.0 * magr2 / magr1 * pow(sin(dnu * 0.5), 2)));
			v1tv2t = lambhodograph(r1, v1, r2, p, ecc, dnu, dtsec, mu);

			v1t << v1tv2t[0], v1tv2t[1], v1tv2t[2];
			v2t << v1tv2t[3], v1tv2t[4], v1tv2t[5];
			//            fprintf(1, 'oldb v1t %16.8f %16.8f %16.8f %16.8f\n', v1dv, mag(v1dv));

			// Battin solution to orbital parameters(and velocities)
			// thompson 2011, loechler 1988
			if (dnu > pi)
				lam = -sqrt((s - chord) / s);
			else
				lam = sqrt((s - chord) / s);

			// loechler pg 21 seems correct!
			for (i = 0; i < 3; i++)
			{
				v1dvl[i] = 1.0 / (lam * (1.0 + lam)) * sqrt(mu * (1.0 + x) / (2.0 * s * s * s * (L + x))) * ((r2[i] - r1[i]) + s * pow(1.0 + lam, 2) * (L + x) / (magr1 * (1.0 + x)) * r1[i]);
				// added v2
				v2dvl[i] = 1.0 / (lam * (1.0 + lam)) * sqrt(mu * (1.0 + x) / (2.0 * s * s * s * (L + x))) * ((r2[i] - r1[i]) - s * pow(1.0 + lam, 2) * (L + x) / (magr2 * (1.0 + x)) * r2[i]);
			}
			//fprintf(1, 'loe v1t %16.8f %16.8f %16.8f %16.8f\n', v1dvl, mag(v1dvl));
			//fprintf(1, 'loe v2t %16.8f %16.8f %16.8f %16.8f\n', v2dvl, mag(v2dvl));
		}  // if loops converged < 30

		//  if (fileout != null)
		//    fprintf(fileout, "%8.5f %3d\n", testamt, loops);

		//~ if (dtsec < 0.001)
			//~ fprintf(outfile, " \n");
		//~ else
			//~ fprintf(outfile, "x %3i %c %5i %12.5f %12.5f %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %12.9f \n",
			//~ nrev, dm, loops, dtsec, dtsec, y, f, v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);

		bigt = sqrt(8.0 / (s * s * s)) * dtsec;

		return uu.mat36(v1dvl, v2dvl);
	};  // lambertbattin
	Matrix61z TrajGenFlyby(Matrix61z Date1, Matrix61z Date2, int planet1, int planet2, gmpf mu, string sSV, string sCOE, char lambertSL, char lambertDR, gmpf odeStep, int Vflag, Matrix31z Vstart) {
		Uti u;
		Trajectory trj2;
		Orbital o;
	
		cout << "From planet " << planet1 << " to " << planet2 << endl << endl;
		gmpf JD1 = o.JD(Date1);
		gmpf JD2 = o.JD(Date2);
		u.DispNum(JD1, "JD1");
		u.DispNum(JD2, "JD2");
		Matrix61z SV1, SV2;
		if (!planet1 == 0)
		{
			SV1 = o.Eph(planet1, JD1);
		}
		else
		{
			SV1 << u.au2km() * 0.1610, u.au2km() * 0.7079, u.au2km() * 0.0004, 0, 0, 0;
		}
		if (!planet2 == 0)
		{
			SV2 = o.Eph(planet2, JD2);
		}
		else
		{
			SV2 << u.au2km() * 0.1610, u.au2km() * 0.7079, u.au2km() * 0.0004, 0, 0, 0;
		}
		Matrix31z R1, V1;
		R1 << SV1(0, 0), SV1(1, 0), SV1(2, 0);
		if (Vflag == 1) {
			V1 = Vstart;
		}
		else {
			V1 << SV1(3, 0), SV1(4, 0), SV1(5, 0);
		}
		Matrix31z R2, V2;
		R2 << SV2(0, 0), SV2(1, 0), SV2(2, 0);
		V2 << SV2(3, 0), SV2(4, 0), SV2(5, 0);
		Matrix61z VLam = lambertB(R1, R2, V1, lambertSL, lambertDR, 0, (JD2 - JD1) * 24.0 * 3600.0, mu);
		Matrix31z Vl11, Vl12;
		Vl11 << VLam(0, 0), VLam(1, 0), VLam(2, 0);
		Vl12 << VLam(3, 0), VLam(4, 0), VLam(5, 0);
		Matrix31z deltaV1 = Vl11 - V1;
		u.DispNum(u.nrm31(deltaV1), "DeltaV Start");
		Matrix61z intcon = SV1;
		intcon(3) = VLam(0);
		intcon(4) = VLam(1);
		intcon(5) = VLam(2);
		Matrix61z ode = trj2.RK4(intcon, 0.0, (JD2 - JD1) * 24.0 * 3600.0, odeStep, mu, sSV, sCOE);
		return VLam;
	};
	Matrix61z TrajGenArrival(Matrix61z Date1, Matrix61z Date2, int planet1, int planet2, gmpf mu, string sSV, string sCOE, char lambertSL, char lambertDR, gmpf odeStep, int Vflag, Matrix31z Vstart, gmpf Rpark, gmpf muD) {
		Uti u;
		Trajectory trj2;
		Orbital o;

		cout << "From planet " << planet1 << " to " << planet2 << endl << endl;
		gmpf JD1 = o.JD(Date1);
		gmpf JD2 = o.JD(Date2);
		u.DispNum(JD1, "JD1");
		u.DispNum(JD2, "JD2");
		Matrix61z SV1 = o.Eph(planet1, JD1);
		Matrix61z SV2 = o.Eph(planet2, JD2);
		Matrix31z R1, V1;
		R1 << SV1(0, 0), SV1(1, 0), SV1(2, 0);
		if (Vflag == 1) {
			V1 = Vstart;
		}
		else {
			V1 << SV1(3, 0), SV1(4, 0), SV1(5, 0);
		}
		Matrix31z R2, V2;
		R2 << SV2(0, 0), SV2(1, 0), SV2(2, 0);
		V2 << SV2(3, 0), SV2(4, 0), SV2(5, 0);
		Matrix61z VLam = lambertB(R1, R2, V1, lambertSL, lambertDR, 0, (JD2 - JD1) * 24.0 * 3600.0, mu);
		Matrix31z Vl11, Vl12;
		Vl11 << VLam(0, 0), VLam(1, 0), VLam(2, 0);
		Vl12 << VLam(3, 0), VLam(4, 0), VLam(5, 0);
		Matrix31z deltaV1 = Vl11 - V1;
		u.DispNum(u.nrm31(deltaV1), "DeltaV Start");
		u.DispNum(u.nrm31(Vl12 - V2), "Excess Velocity Arrival for Landing");

		gmpf Vp2 = sqrt(pow(u.nrm31(Vl12 - V2), 2) + ((2 * muD) / (Rpark)));
		gmpf Vpark2 = sqrt(muD / Rpark);
		gmpf DVatsat = Vp2 - Vpark2;
		u.DispNum(DVatsat, "Excess Velocity Arrival for Parking Orbit");

		Matrix61z intcon = SV1;
		intcon(3) = VLam(0);
		intcon(4) = VLam(1);
		intcon(5) = VLam(2);
		Matrix61z ode = trj2.RK4(intcon, 0.0, (JD2 - JD1) * 24.0 * 3600.0, odeStep, mu, sSV, sCOE);
		return VLam;
	};
};

int main(int argc, char** argv) {

	// Object Instantiation

	Uti u;
	Sun sun;
	Earth earth;
	Saturn saturn;
	Lambert l;
	Trajectory trj;
	Orbital o;

	// Defaults

	gmpf odeStep = 24.0 * 3600.0;
	gmpf Rpark = 377420;

	// Variables

	Matrix61z Tearth, Tjupiter, Tsaturn, Tvenus, Tmars, Tearth2, Tdsm, Tvenus2;;
	Matrix31z dum, Vstart;
	dum << 0, 0, 0;

	// Case 1 (EJSD)
	
		cout << "Case 1 (EJSD)" << endl << endl;

		Tearth << 2020, 4, 3, 0, 0, 0;
		Tjupiter << 2021, 8, 30, 0, 0, 0;
		Tsaturn << 2027, 4, 4, 0, 0, 0;

		Matrix61z vCase1Leg1 = l.TrajGenFlyby(Tearth, Tjupiter, 3, 5, sun.mu(), "SV_Case1_Leg1.dat", "COE_Case1_Leg1.dat", 's', 'd', odeStep, 0, dum);
		Vstart << vCase1Leg1(3), vCase1Leg1(4), vCase1Leg1(5);
		Matrix61z vCase1Leg2 = l.TrajGenArrival(Tjupiter, Tsaturn, 5, 6, sun.mu(), "SV_Case1_Leg2.dat", "COE_Case1_Leg2.dat", 's', 'd', odeStep, 1, Vstart, Rpark, saturn.mu());

		cout << "Case 1 Direct" << endl << endl;

		Matrix61z vCase1dir = l.TrajGenArrival(Tearth, Tsaturn, 3, 6, sun.mu(), "SV_Case1_dir.dat", "COE_Case1_dir.dat", 's', 'd', odeStep, 0, dum, Rpark, saturn.mu());
	
	// Case 2 (EJSD)
	
		cout << "Case 2 (EJSD)" << endl << endl;

		Tearth << 2021, 5, 4, 0, 0, 0;
		Tjupiter << 2022, 7, 24, 0, 0, 0;
		Tsaturn << 2028, 5, 4, 0, 0, 0;

		Matrix61z vCase2Leg1 = l.TrajGenFlyby(Tearth, Tjupiter, 3, 5, sun.mu(), "SV_Case2_Leg1.dat", "COE_Case2_Leg1.dat", 's', 'd', odeStep, 0, dum);
		Vstart << vCase2Leg1(3), vCase2Leg1(4), vCase2Leg1(5);
		Matrix61z vCase2Leg2 = l.TrajGenArrival(Tjupiter, Tsaturn, 5, 6, sun.mu(), "SV_Case2_Leg2.dat", "COE_Case2_Leg2.dat", 's', 'd', odeStep, 1, Vstart, Rpark, saturn.mu());

		cout << "Case 2 Direct" << endl << endl;

		Matrix61z vCase2Dir = l.TrajGenArrival(Tearth, Tsaturn, 3, 6, sun.mu(), "SV_Case2_dir.dat", "COE_Case2_dir.dat", 's', 'd', odeStep, 0, dum, Rpark, saturn.mu());
	
	// Case 3 (EMSD)
	
		cout << "Case 3 (EMSD)" << endl << endl;

		Tearth << 2033, 4, 20, 0, 0, 0;
		Tmars << 2033, 11, 5, 0, 0, 0;
		Tsaturn << 2038, 7, 20, 0, 0, 0;

		Matrix61z vCase3Leg1 = l.TrajGenFlyby(Tearth, Tmars, 3, 4, sun.mu(), "SV_Case3_Leg1.dat", "COE_Case3_Leg1.dat", 's', 'd', odeStep, 0, dum);
		Vstart << vCase3Leg1(3), vCase3Leg1(4), vCase3Leg1(5);
		Matrix61z vCase3Leg2 = l.TrajGenArrival(Tmars, Tsaturn, 4, 6, sun.mu(), "SV_Case3_Leg2.dat", "COE_Case3_Leg2.dat", 's', 'd', odeStep, 1, Vstart, Rpark, saturn.mu());

		cout << "Case 3 Direct" << endl << endl;

		Matrix61z vCase3Dir = l.TrajGenArrival(Tearth, Tsaturn, 3, 6, sun.mu(), "SV_Case3_dir.dat", "COE_Case3_dir.dat", 's', 'd', odeStep, 0, dum, Rpark, saturn.mu());
	
	// Case 4 (EMVSD)

		cout << "Case 4 (EMVSD)" << endl << endl;

		Tearth << 2020, 6, 20, 0, 0, 0;
		Tmars << 2020, 12, 31, 0, 0, 0;
		Tvenus << 2022, 3, 27, 0, 0, 0;
		Tsaturn << 2029, 9, 9, 0, 0, 0;

		Matrix61z vCase4Leg1 = l.TrajGenFlyby(Tearth, Tmars, 3, 4, sun.mu(), "SV_Case4_Leg1.dat", "COE_Case4_Leg1.dat", 's', 'd', odeStep, 0, dum);
		Vstart << vCase4Leg1(3), vCase4Leg1(4), vCase4Leg1(5);
		Matrix61z vCase4Leg2 = l.TrajGenFlyby(Tmars, Tvenus, 4, 2, sun.mu(), "SV_Case4_Leg2.dat", "COE_Case4_Leg2.dat", 's', 'd', odeStep, 1, Vstart);
		Vstart << vCase4Leg2(3), vCase4Leg2(4), vCase4Leg2(5);
		Matrix61z vCase4Leg3 = l.TrajGenArrival(Tvenus, Tsaturn, 2, 6, sun.mu(), "SV_Case4_Leg3.dat", "COE_Case4_Leg3.dat", 's', 'd', odeStep, 1, Vstart, Rpark, saturn.mu());

		cout << "Case 4 Direct" << endl << endl;

		Matrix61z vCase4Dir = l.TrajGenArrival(Tearth, Tsaturn, 3, 6, sun.mu(), "SV_Case4_dir.dat", "COE_Case4_dir.dat", 's', 'd', odeStep, 0, dum, Rpark, saturn.mu());
	
	// Case 5 (VEEGA)
		
		//cout << "Case 5 (VEEGA)" << endl << endl;
		//
		//Tearth << 1997, 10, 15, 9, 36, 0; // 2450736.9
		//Tvenus << 1998, 4, 26, 14, 24, 0; //2450930.1
		//Tdsm << 1998, 12, 3, 7, 12, 0; //2451150.8
		//Tvenus2 << 1999, 6, 24, 19, 12, 0; //2451354.3
		//Tearth2 << 1999, 8, 18, 4, 48, 0; //2451408.7
		//Tjupiter << 2000, 12, 30, 12, 0, 0; //2451909
		//Tsaturn << 2004, 7, 1, 9, 36, 0; //2453187.9

		//Matrix61z vCase5Leg1 = l.TrajGenFlyby(Tearth, Tvenus, 3, 2, sun.mu(), "SV_Case5_Leg1.dat", "COE_Case5_Leg1.dat", 'l', 'd', odeStep, 0, dum);
		//Vstart << vCase5Leg1(3), vCase5Leg1(4), vCase5Leg1(5);

		//Matrix61z vCase5Leg2 = l.TrajGenFlyby(Tvenus, Tdsm, 2, 0, sun.mu(), "SV_Case5_Leg2.dat", "COE_Case5_Leg2.dat", 's', 'd', odeStep, 1, Vstart);
		//Vstart << vCase5Leg2(3), vCase5Leg2(4), vCase5Leg2(5);

		//Matrix61z vCase5Leg3 = l.TrajGenFlyby(Tdsm, Tvenus2, 0, 2, sun.mu(), "SV_Case5_Leg3.dat", "COE_Case5_Leg3.dat", 's', 'd', odeStep, 1, Vstart);
		//Vstart << vCase5Leg3(3), vCase5Leg3(4), vCase5Leg3(5);

		//Matrix61z vCase5Leg4 = l.TrajGenFlyby(Tvenus2, Tearth2, 2, 3, sun.mu(), "SV_Case5_Leg4.dat", "COE_Case5_Leg4.dat", 's', 'd', odeStep, 1, Vstart);
		//Vstart << vCase5Leg4(3), vCase5Leg4(4), vCase5Leg4(5);

		//Matrix61z vCase5Leg5 = l.TrajGenFlyby(Tearth2, Tjupiter, 3, 5, sun.mu(), "SV_Case5_Leg5.dat", "COE_Case5_Leg5.dat", 's', 'd', odeStep, 1, Vstart);
		//Vstart << vCase5Leg5(3), vCase5Leg5(4), vCase5Leg5(5);

		//Matrix61z vCase5Leg6 = l.TrajGenArrival(Tjupiter, Tsaturn, 5, 6, sun.mu(), "SV_Case5_Leg6.dat", "COE_Case5_Leg6.dat", 's', 'd', odeStep, 1, Vstart, Rpark, saturn.mu());

		//cout << "Case 5 Direct" << endl << endl;

		//Matrix61z vCase5Dir = l.TrajGenArrival(Tearth, Tsaturn, 3, 6, sun.mu(), "SV_Case5_dir.dat", "COE_Case5_dir.dat", 's', 'd', odeStep, 0, dum, Rpark, saturn.mu());

		//// Case 6 (VEEGA No DSM)

		//cout << "Case 6 (VEEGA No DSM)" << endl << endl;

		//Tearth << 1997, 10, 15, 9, 36, 0; // 2450736.9
		//Tvenus << 1998, 4, 26, 14, 24, 0; //2450930.1
		//Tdsm << 1998, 12, 3, 7, 12, 0; //2451150.8
		//Tvenus2 << 1999, 6, 24, 19, 12, 0; //2451354.3
		//Tearth2 << 1999, 8, 18, 4, 48, 0; //2451408.7
		//Tjupiter << 2000, 12, 30, 12, 0, 0; //2451909
		//Tsaturn << 2004, 7, 1, 9, 36, 0; //2453187.9

		//Matrix61z vCase6Leg1 = l.TrajGenFlyby(Tearth, Tvenus2, 3, 2, sun.mu(), "SV_Case6_Leg1.dat", "COE_Case6_Leg1.dat", 'l', 'd', odeStep, 0, dum);
		//Vstart << vCase6Leg1(3), vCase6Leg1(4), vCase6Leg1(5);

		//Matrix61z vCase6Leg2 = l.TrajGenFlyby(Tvenus2, Tearth2, 2, 3, sun.mu(), "SV_Case6_Leg2.dat", "COE_Case6_Leg2.dat", 's', 'd', odeStep, 1, Vstart);
		//Vstart << vCase6Leg2(3), vCase6Leg2(4), vCase6Leg2(5);

		//Matrix61z vCase6Leg3 = l.TrajGenFlyby(Tearth2, Tjupiter, 3, 5, sun.mu(), "SV_Case6_Leg3.dat", "COE_Case6_Leg3.dat", 's', 'd', odeStep, 1, Vstart);
		//Vstart << vCase6Leg3(3), vCase6Leg3(4), vCase6Leg3(5);

		//Matrix61z vCase6Leg4 = l.TrajGenArrival(Tjupiter, Tsaturn, 5, 6, sun.mu(), "SV_Case6_Leg4.dat", "COE_Case6_Leg4.dat", 's', 'd', odeStep, 1, Vstart, Rpark, saturn.mu());

		//cout << "Case 6 Direct" << endl << endl;

		//Matrix61z vCase6Dir = l.TrajGenArrival(Tearth, Tsaturn, 3, 6, sun.mu(), "SV_Case6_dir.dat", "COE_Case6_dir.dat", 's', 'd', odeStep, 0, dum, Rpark, saturn.mu());

		cout << "END" << endl;

};
