#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <deque>
#include <string.h>
#include <sstream>

#if defined(__WIN32__)
#include <direct.h>
#include <windows.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif

#include <TROOT.h>
#include <TThread.h>
#include <TApplication.h>
#include <TRint.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <Rtypes.h>
#include <TVector.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Point2D.h>
#include <TRandom1.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TStyle.h"

double pulse_gauss_f(double *x, double *pars)
{
	return pars[0] * std::exp(-0.5* std::pow((x[0] - pars[1])/pars[2], 2));
}

double pulse_fitting_f(double *x, double *pars)
{
	//TODO:
}

TCanvas * canvas;
std::pair<double, double> time_window;
std::deque<std::deque<double>> raw_signals;

struct signal {
	std::deque<double> form;
	std::pair<double, double> xy_offset;
	std::string label;
	std::deque<double> linear_combination;
	long double Norm;
	std::pair<double, double> domain;
	std::pair<double, double> cursors;
	TF1 * fitter;
	TH1D* histogram;
	unsigned int n_bins;
	bool logscale;
	unsigned int pad;
};

std::deque<signal> signals;

std::deque<double> load_signal(std::string filename)
{
	std::deque<double> ret;
	std::ifstream str;
	str.open(filename, std::ios_base::binary);
	if (!str.is_open()) {
		std::cout << "Failed to open file \"" << filename << "\"" << std::endl;
		return ret;
	}
	double time;
	while (!str.eof()) {
		str.read((char*)&time, sizeof(double));
		ret.push_back(time);
	}
	str.close();
}

void init(void)
{
	//number of pads, raw signals and working signals may be different
	time_window = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	canvas = new TCanvas("Signal forms", "Signal forms", 1000, 900);
	raw_signals.resize(4);
	signals.resize(3);
	canvas->Divide(1, 3);
	raw_signals[0] = load_signal("test0.dat");
	raw_signals[1] = load_signal("test1.dat");
	raw_signals[2] = load_signal("test2.dat");
	raw_signals[3] = load_signal("test3.dat");
	
	signals[0].linear_combination = {1,0,0,0};
	signals[1].linear_combination = {0,1,0,0};
	signals[2].linear_combination = {0,0,1,1};
		
	signals[0].pad = 0;
	signals[1].pad = 1;
	signals[2].pad = 2;
	signals[0].n_bins = 300;
	signals[1].n_bins = 300;
	signals[2].n_bins = 300;
	signals[0].xy_offset = std::pair<double, double> (0,0);
	signals[1].xy_offset = std::pair<double, double> (0,0);
	signals[2].xy_offset = std::pair<double, double> (0,0);
	signals[0].label = "3PMT";
	signals[1].label = "PMT#1";
	signals[2].label = "SiPM";
	signals[0].domain = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	signals[1].domain = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	signals[2].domain = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	signals[0].cursors = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	signals[1].cursors = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	signals[2].cursors = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	
}


std::pair<double, double> x_limits(int channel)
{

}

std::pair<double, double> x_clipped_limits(int channel)
{

}

void set_domain(int ch, double from, double to) //-1==ch for all channels
{

}

void set_n_bins(int ch, double from, double to) //-1==ch for all channels
{

}

void add(int ch1, int ch2, double coeff) //set coef to 0 to remove sum
{

}

void set_offset(int ch, double dx, double dy) //-1 for all channels
{

}

void fit(int ch)
{

}

void set_fit_f(int ch, int type, std::pair<double, double> domain)
{

}

void set_fit_par_lims(int ch, int par, double from, double to)
{

}

double get_fit_par(int ch, int par)
{
	return 0;
}

TH1D *create_hist(int ch)
{
	return NULL;
}

void get_xy_axis(int pad, std::pair<double, double> &x_lims, std::pair<double, double> &y_lims)
{

}

TLegend* get_legend(int pad)
{
	return NULL;
}

void replot(void)
{

}

void script(void) 
{

}