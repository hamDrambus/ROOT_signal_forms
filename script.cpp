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
#include "TLegend.h"

double pulse_gauss_f(double *x, double *pars)
{
	return pars[0] * std::exp(-0.5* std::pow((x[0] - pars[1])/pars[2], 2));
}

double pulse_fitting_f(double *x, double *pars)
{
	//pars: {y0, A, x0, t_exp_rising, t_exp_falling}
	return pars[0] + pars[1]*(1- std::exp((pars[2]-x[0])/pars[3]))*std::exp((pars[2]-x[0])/pars[4]);
}

TCanvas * canvas;
int N_pads;
std::deque<Color_t> colors;

struct signal { //can be either raw or derived from raws
	std::deque<double> raw_signal;
	std::pair<double, double> xy_offset;
	std::string label;
	std::deque<double> linear_combination; //ignored if raw_signal is not empty
	long double Norm;
	std::pair<double, double> domain; //raw_signal (if present) is ignored outside the domain
	std::pair<double, double> cursors;
	TF1 * fitter; //used for signal operations when both raw signal and linear_combination are 0
	TH1D* histogram;
	unsigned int n_bins;
	bool logscale;
	unsigned int pad; //-1 for not displaying
};

std::deque<signal> signals;

std::deque<double> load_signal(std::string filename);
void init(void);
std::pair<double, double> x_raw_limits(int channel);
std::pair<double, double> x_clipped_limits(int channel); //does not account for resulting form being 0 after linear combination
void set_domain(int ch, double from, double to); //-1==ch for all channels
void set_n_bins(int ch, int N); //-1==ch for all channels
void add(int ch1, int ch2, double coeff); //set coef to 0 to remove sum. Interdependent linear combinations are not supported
void set_offset(int ch, double dx, double dy); //-1 for all channels
void fit(int ch);
void set_fit_f(int ch, int type, std::pair<double, double> domain);
void set_fit_par_lims(int ch, int par, double from, double to);
void set_fit_par(int ch, int par, double val);
double get_fit_par(int ch, int par);
void get_xy_axis(int pad, std::pair<double, double> &x_lims, std::pair<double, double> &y_lims);
TLegend* get_legend(int pad);
void FillHist(int ch, TH1D * temp_hist = NULL, std::pair<double, double> offset = std::pair<double, double>(0,0)); //temp_hist is necessary for arithmetic operation with signals
void replot(void);

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
	return ret;
}

void init(void)
{
	gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("");
	colors.push_back(kBlack);
	colors.push_back(kRed);
	colors.push_back(kBlue);
	colors.push_back(kGreen);
	colors.push_back(kMagenta);
	//number of pads, raw signals and working signals may be different
	canvas = new TCanvas("Signal forms", "Signal forms", 1000, 900);
	canvas->SetGrid();
	canvas->SetTicks();
	signals.resize(4);
	N_pads = 3;
	canvas->Divide(1, N_pads);
	signals[0].raw_signal = load_signal("test0.dat");
	signals[1].raw_signal = load_signal("test1.dat");
	signals[2].raw_signal = load_signal("test2.dat");
	
	signals[3].linear_combination = {1, 0, -1, 0};
		
	signals[0].pad = 0;
	signals[1].pad = 1;
	signals[2].pad = 2;
	signals[3].pad = 0;
	signals[0].n_bins = 300;
	signals[1].n_bins = 300;
	signals[2].n_bins = 300;
	signals[3].n_bins = 300;
	signals[0].xy_offset = std::pair<double, double> (0,0);
	signals[1].xy_offset = std::pair<double, double> (0,0);
	signals[2].xy_offset = std::pair<double, double> (0,0);
	signals[3].xy_offset = std::pair<double, double> (0,0);
	signals[0].label = "3PMT";
	signals[1].label = "PMT#1";
	signals[2].label = "SiPM";
	signals[3].label = "3PMT - SiPM";
	signals[0].domain = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	signals[1].domain = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	signals[2].domain = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	signals[3].domain = std::pair<double, double> (-DBL_MAX, DBL_MAX);
	signals[0].cursors = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	signals[1].cursors = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	signals[2].cursors = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	signals[3].cursors = std::pair<double, double>(-DBL_MAX, DBL_MAX);
}


std::pair<double, double> x_raw_limits(int channel)
{
	std::pair<double, double> ret(DBL_MAX, -DBL_MAX);
	if (channel>=signals.size() || channel<0) {
		std::cout<<"x_raw_limits: Error: channel is out of range"<<std::endl;
		return ret;
	}
	for (std::size_t i =0, i_end_ = signals[channel].raw_signal.size(); i!=i_end_; ++i) {
		if (signals[channel].raw_signal[i] < ret.first)
			ret.first = signals[channel].raw_signal[i];
		if (signals[channel].raw_signal[i] > ret.second)
			ret.second = signals[channel].raw_signal[i];
	}
	if (ret.first!=DBL_MAX) {
		ret.first += signals[channel].xy_offset.first;
		ret.second += signals[channel].xy_offset.first;
	}
	if (signals[channel].domain.first!=-DBL_MAX)
		ret.first = signals[channel].domain.first;
	if (signals[channel].domain.second!=DBL_MAX)
		ret.second = signals[channel].domain.second;
	return ret;
}

std::pair<double, double> x_clipped_limits(int channel) //does not account for resulting form being 0 after linear combination
{
	std::pair<double, double> ret(DBL_MAX, -DBL_MAX);
	if (channel>=signals.size() || channel<0) {
		std::cout<<"x_clipped_limits: Error: channel is out of range"<<std::endl;
		return ret;
	}
	ret = x_raw_limits(channel);
	if (ret.first!=DBL_MAX) {//there is raw signal for this channel (or fixed domain)
		return ret;
	}
	for (std::size_t i =0, i_end_ = signals[channel].linear_combination.size(); i!=i_end_; ++i) {
		std::pair<double, double> lims = ret;
		if (signals[channel].linear_combination[i]!=0 && i!=channel) {
			lims = x_clipped_limits(i);
		}
		if (lims.first!=DBL_MAX) {
			ret.first = std::min(ret.first, lims.first);
			ret.second = std::max(ret.second, lims.second);
		}
	}
	if (ret.first!=DBL_MAX) {
		ret.first += signals[channel].xy_offset.first;
		ret.second += signals[channel].xy_offset.first;
	}
	if (signals[channel].domain.first!=-DBL_MAX)
		ret.first = signals[channel].domain.first;
	if (signals[channel].domain.second!=DBL_MAX)
		ret.second = signals[channel].domain.second;
	return ret;
}

void set_domain(int ch, double from, double to) //-1==ch for all channels
{
	if (ch>=signals.size()) {
		std::cout<<"set_domain: Error: channel is out of range"<<std::endl;
		return;
	}
	if (ch<0) {
		for (std::size_t i =0, i_end_ = signals.size(); i!=i_end_; ++i)
			set_domain(i, from, to);
	}
	signals[ch].domain.first = std::min(from, to);
	signals[ch].domain.second = std::max(from, to);
}

void set_n_bins(int ch, int N) //-1==ch for all channels
{
	if (ch>=signals.size()) {
		std::cout<<"set_n_bins: Error: channel is out of range"<<std::endl;
		return;
	}
	if (ch<0) {
		for (std::size_t i =0, i_end_ = signals.size(); i!=i_end_; ++i)
			set_n_bins(i, N);
	}
	signals[ch].n_bins = std::max(N, 1);
}

void add(int ch1, int ch2, double coeff) //set coef to 0 to remove sum. Interdependent linear combinations are not supported
{
	if (ch1>=signals.size() || ch1<0) {
		std::cout<<"add: Error: channel1 is out of range"<<std::endl;
		return;
	}
	if (ch2>=signals.size() || ch2<0) {
		std::cout<<"add: Error: channel2 is out of range"<<std::endl;
		return;
	}
	if (ch1==ch2) {
		std::cout<<"add: Error: channel1 and channel2 cannot be the same"<<std::endl;
		return;
	}
	signals[ch1].linear_combination[ch2] = coeff;
}

void set_offset(int ch, double dx, double dy) //-1 for all channels
{
	if (ch>=signals.size()) {
		std::cout<<"set_offset: Error: channel is out of range"<<std::endl;
		return;
	}
	if (ch<0) {
		for (std::size_t i =0, i_end_ = signals.size(); i!=i_end_; ++i)
			set_offset(i, dx, dy);
	}
	signals[ch].xy_offset.first = dx;
	signals[ch].xy_offset.second = dy;
}

void fit(int ch)
{
	if (ch>=signals.size()||ch<0) {
		std::cout<<"fit: Error: channel is out of range"<<std::endl;
		return;
	}
	if ((signals[ch].fitter!=NULL)&&(NULL!=signals[ch].histogram))
		signals[ch].histogram->Fit(signals[ch].fitter, "RQ");
}

void set_fit_f(int ch, int type, std::pair<double, double> domain)
{
	if (ch>=signals.size()||ch<0) {
		std::cout<<"set_fit_f: Error: channel is out of range"<<std::endl;
		return;
	}
	switch (type)
	{
	case 0:	{
		if (NULL!=signals[ch].fitter)
			signals[ch].fitter->Delete();
		signals[ch].fitter = new TF1 ("gaus", pulse_gauss_f , domain.first, domain.second, 3);
		break;
	}
	case 1:	{
		if (NULL!=signals[ch].fitter)
			signals[ch].fitter->Delete();
		signals[ch].fitter = new TF1 ("pulse", pulse_fitting_f , domain.first, domain.second, 5);
		break;
	}
	default: {
		std::cout<<"set_fit_f: Error: fit function type is out of range"<<std::endl;
	}
	}
	return;
}

void set_fit_par_lims(int ch, int par, double from, double to)
{
	if (ch>=signals.size()||ch<0) {
		std::cout<<"set_fit_par_lims: Error: channel is out of range"<<std::endl;
		return;
	}
	if (NULL==signals[ch].fitter) {
		std::cout<<"set_fit_par_lims: Error: no fit function"<<std::endl;
		return;
	}
	signals[ch].fitter->SetParLimits(par, from, to);
	return;
}

void set_fit_par(int ch, int par, double val)
{
	if (ch>=signals.size()||ch<0) {
		std::cout<<"set_fit_par: Error: channel is out of range"<<std::endl;
		return;
	}
	if (NULL==signals[ch].fitter) {
		std::cout<<"set_fit_par: Error: no fit function"<<std::endl;
		return;
	}
	signals[ch].fitter->SetParameter(par, val);
	return;
}

double get_fit_par(int ch, int par)
{
	if (ch>=signals.size()||ch<0) {
		std::cout<<"set_fit_par: Error: channel is out of range"<<std::endl;
		return DBL_MAX;
	}
	if (NULL==signals[ch].fitter) {
		std::cout<<"set_fit_par: Error: no fit function"<<std::endl;
		return DBL_MAX;
	}
	return signals[ch].fitter->GetParameter(par);
}

void get_xy_axis(int pad, std::pair<double, double> &x_lims, std::pair<double, double> &y_lims)
{
	x_lims = std::pair<double, double> (DBL_MAX, -DBL_MAX);
	y_lims = std::pair<double, double> (DBL_MAX, -DBL_MAX);
	for (std::size_t ch, ch_end_=signals.size(); ch!=ch_end_; ++ch)	{
		if (pad==signals[ch].pad) {
			std::pair<double, double> x_lims_ch = x_clipped_limits(ch);
			x_lims.first = std::min(x_lims.first, x_lims_ch.first);
			x_lims.second = std::max(x_lims.second, x_lims_ch.second);
			if (NULL==signals[ch].histogram) {
				std::cout<<"get_xy_axis: Error: NULL histogram for ch"<<ch<<std::endl;
				continue;
			}
			y_lims.first = std::min(signals[ch].histogram->GetMinimum(-DBL_MAX), y_lims.first);
			y_lims.second = std::max(signals[ch].histogram->GetMaximum(DBL_MAX), y_lims.second);
		}
	}
}

TLegend* get_legend(int pad)
{
	return NULL;
}

void FillHist(int ch, TH1D * temp_hist, std::pair<double, double> offset) //temp_hist is necessary for arithmetic operation with signals
{
	if (ch>=signals.size()||ch<0) {
		std::cout<<"FillHist: Error: channel is out of range"<<std::endl;
		return;
	}
	if (NULL==signals[ch].histogram && temp_hist==NULL) {
		std::cout<<"FillHist: Error: NULL histogram"<<std::endl;
		return;
	}
	TH1D* hist = (NULL==temp_hist ? signals[ch].histogram : temp_hist);
	if (!signals[ch].raw_signal.empty()) {
		for (std::size_t i=0, i_end_ = signals[ch].raw_signal.size(); i!=i_end_; ++i)
			hist->Fill(signals[ch].raw_signal[i] + signals[ch].xy_offset.first + offset.first);
		for (int bin=1, bin_end_ = hist->GetNbinsX()+1; bin!=bin_end_; ++bin) { //ROOT's TH1 bin numbering is a bit tricky
			hist->AddBinContent(bin, offset.second);
		}
	} else {
		bool not_zero = false;
		for (std::size_t i=0, i_end_ = signals[ch].linear_combination.size(); i!=i_end_; ++i) {
			if (0!=signals[ch].linear_combination[i] && i!=ch) {
				not_zero = true;
				std::string name = "temp_" + std::to_string(ch) + "_" + std::to_string(i);
				if (gROOT->FindObject(name.c_str())!=NULL) {
					std::cout<<"FillHist: Error: histogram already exists - there is interdependency of signals"<<std::endl;
					continue;
				}
				TH1D * temp = new TH1D(name.c_str(), name.c_str(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmin());
				offset.first +=signals[ch].xy_offset.first;
				offset.second +=signals[ch].xy_offset.second;
				FillHist(i, temp, offset);
				if (false==hist->Add(temp, signals[ch].linear_combination[i])) {
					std::cout<<"FillHist: Error: failed to add histograms"<<std::endl;
				}
				temp->Delete();
			}
		}
		if (!not_zero) {
			if (signals[ch].fitter==NULL) {
				std::cout<<"FillHist: Error: empty signal ch "<<std::to_string(ch)<<": there is no raw signal, linear combination or fit function"<<std::endl;
				return;
			}
			for (int bin=1, bin_end_ = hist->GetNbinsX()+1; bin!=bin_end_; ++bin) { //ROOT's TH1 bin numbering is a bit tricky
				hist->SetBinContent(bin, offset.second + signals[ch].fitter->Eval(hist->GetBinCenter(bin) - offset.first));
			}
		}
	}
	return;
}

void replot(void)
{
	canvas->Clear();
	for (std::size_t ch, ch_end_=signals.size(); ch!=ch_end_; ++ch)	{
		std::pair<double, double> domain = x_clipped_limits(ch);
		if ((domain.first == -DBL_MAX)&&(domain.first == DBL_MAX)) {
			std::cout<<"replot: Error: no x axis limits for ch "<<ch<<std::endl;
			domain = std::pair<double, double> (0, 10);
		}
		if ((domain.first == -DBL_MAX)&&(domain.first != DBL_MAX)) {
			std::cout<<"replot: Error: no lower x axis limit for ch "<<ch<<std::endl;
			domain.first = domain.second-10;
		}
		if ((domain.first != -DBL_MAX)&&(domain.first == DBL_MAX)) {
			std::cout<<"replot: Error: no upper x axis limit for ch "<<ch<<std::endl;
			domain.second = domain.first+10;
		}
		if (NULL==signals[ch].histogram) {
			std::string name = "ch_" + std::to_string(ch);
			signals[ch].histogram = new TH1D (name.c_str(), name.c_str(), signals[ch].n_bins, domain.first, domain.second);
		} else {
			signals[ch].histogram->SetBins(signals[ch].n_bins, domain.first, domain.second);
		}
		FillHist(ch, NULL, signals[ch].xy_offset);
	}
	for (int pad = 0; pad<N_pads; ++pad) {
		std::pair<double, double> domain_x, domain_y;
		get_xy_axis(pad, domain_x, domain_y);
		canvas->cd(pad);
		std::string name = "frame_" + std::to_string(pad);
		TH2F* frame = new TH2F(name.c_str(), name.c_str(), 500, domain_x.first, domain_x.second, 500, domain_y.first*1.1, domain_y.second*1.1);
		frame->GetXaxis()->SetTitle("time [us]");
		frame->GetYaxis()->SetTitle("");
		frame->Draw();
		int counter = 0;
		for (std::size_t ch, ch_end_=signals.size(); ch!=ch_end_; ++ch)	{
			if (NULL!=signals[ch].histogram) {
				signals[ch].histogram->SetLineWidth(2);
				signals[ch].histogram->SetLineColor(colors[counter % colors.size()]);
				signals[ch].histogram->Draw("hist cpsame");
				++counter;
			}
			if (NULL!=signals[ch].fitter) {
				signals[ch].fitter->SetLineWidth(2);
				signals[ch].fitter->SetLineColor(colors[counter % colors.size()]);
				signals[ch].fitter->Draw("same");
				++counter;
			}
		}
		frame->Draw("sameaxis");
		TLegend* legend = get_legend(pad);
		if (NULL!=legend)
			legend->Draw("same");
	}
	canvas->Update();
}

void script(void) 
{
	init();
}
