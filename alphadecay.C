// Code creates synthetic data by filling a histogram with known Exponentiall modified Gaussians (EMG)
// Data is then fitted using a model which is the sum of EMGs using ROOT and displayed

#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include <vector>
#include <iostream>
using namespace std;

// EMG with no amplitude params: [0] = mean, [1] = sigma, [2] = tau 
// tau: tailing parameter
// Will use this function for random sampling to fill histogram
// Amplitude not defined yet, will be defined in fitting

double emg_shape(double *x, double *params) {

    double E = x[0];
    double mean = params[0];
    double sigma = params[1];
    double tau = params[2];

    if (sigma <= 0 || tau <= 0) return 0.0;
    double lambda = 1.0/tau;


    // Define the PDF of EMG 

    double exponent = lambda*(mean - E) + 0.5* pow(lambda,2)*pow(sigma,2);
    double erfc = TMath::Erfc(((mean + lambda*pow(sigma,2)-E)/(sqrt(2)*sigma)));

    double shape = lambda*0.5*TMath::Exp(exponent) * erfc;

    return shape;

}

// Model used for fitting 
int gNpeaks = 0; 
// params[0] = A1, params[1] = mean1, params[2] = sigma1, params[3] = tau1,
// params[4] = A2, params[5] = mean2, params[6] = sigma2, params[7] = tau2, etc
// params[-2] = background signal, params[-1] = linear background signal
double fitting_model(double *x, double *params) {
    double E = x[0];

    extern int gNpeaks;
    int N = gNpeaks;

    double y = 0;

    // Contributions of all the ith peaks into signal
    for (int i = 0; i < N; i++) {
        double Ai = params[4*i];
        double meani = params[4*i+1];
        double sigmai = params[4*i+2];
        double taui = params[4*i+3];

        double EMG_params[3];
        EMG_params[0] = meani;
        EMG_params[1] = sigmai;
        EMG_params[2] = taui;

        y += Ai * emg_shape(x,EMG_params);

    }

    y += params[4*N] + params[4*N+1] * E; // Adding the background signals

    return y;
}



void alphadecay() {

    const double E_min = 4000; // Boundaries of Energies in keV
    const double E_max = 8000;

    const int nbins = 2000; // number of bins in Histogram

    // Create Histogram
    TH1D *h = new TH1D("h","Simulated alpha decay; Energy(keV); Counts",nbins,E_min,E_max);

    TRandom3 rng(48922); // Fixed RNG seed so results can be reproduced

    struct PeakParams {
        double E;
        int counts;
        double sigma;
        double tau;
        const char* name;
    };

    vector<PeakParams> peaks;
    // {peak Energy,counts,sigma,tau, name} Peak energy in keV 
    // Took FWHM to be 30 keV
    double FWHM = 30;
    double u = pow(2*sqrt(2*log(2)),-1);
    cout << u;
    peaks.push_back({5685.37, 50000, FWHM*u, 10, "Ra-224 5685"});
    peaks.push_back({5448.6, 3000, FWHM*u, 10, "Ra-224 5448"});
    peaks.push_back({6288.08, 70000, FWHM*u, 10, "Rn-220 6288"});
    peaks.push_back({6778.3, 90000, FWHM*u, 10, "Po-216 6778"});

    // Some Values taken from paper "Novel device to study double-alpha decay at the FRS ion catcher" by Varga et al.

    int Npeaks = peaks.size();
    gNpeaks = Npeaks;

    vector<TF1*> peak_functions;

    for (int i = 0; i<Npeaks; i++){
        TF1* f = new TF1(Form("peak%d", i), emg_shape,E_min,E_max,3);
        f->SetParName(0,"mean"); f->SetParameter(0,peaks[i].E);
        f->SetParName(1,"sigma"); f->SetParameter(1,peaks[i].sigma);
        f->SetParName(2,"tau"); f->SetParameter(2,peaks[i].tau);

        peak_functions.push_back(f);
    };
    // Filling Histogram
    for (int i = 0; i<Npeaks; i++) {
        int Ni = peaks[i].counts;

        TF1* f = peak_functions[i];

        for (int e = 0;e<Ni;e++) {
            double E = f->GetRandom();
            h -> Fill(E);
        }


    }

    // Filling Histogram with background signals uniformly distributed like CMBR

    int bgcounts = 50000;

    for (int i = 0; i<bgcounts;i++) {
        
        double r = rng.Uniform();
        double E = E_min + r*(E_max - E_min);
        h -> Fill(E);
    }

    // Drawing 
    TCanvas* a = new TCanvas("a", "alpha simulation",1200,700);

    h -> Draw();

    // Fitting model to synthetic data
    int number_of_params = 4*Npeaks +2;
    TF1* model = new TF1("model",fitting_model,E_min,E_max,number_of_params);

    for (int i = 0;i<Npeaks;i++) {
        model -> SetParName(4*i, Form("A%d",i));
        model -> SetParName(4*i +1, Form("mean_%d",i));
        model -> SetParName(4*i+2, Form("sigma_%d",i));
        model -> SetParName(4*i+3, Form("tau_%d",i));

        // initial guesses based on known emission energies 

        model -> SetParameter(4*i,peaks[i].counts);
        model -> SetParameter(4*i+1,peaks[i].E);
        model -> SetParameter(4*i+2,peaks[i].sigma);
        model -> SetParameter(4*i+3,peaks[i].tau);

    }

    // Setting guesses for background radiation
    model -> SetParName(4*Npeaks, "Const_background");
    model -> SetParameter(4*Npeaks, 10);

    model -> SetParName(4*Npeaks + 1, "Slope_background");
    model -> SetParameter(4*Npeaks +1, 0);

    // Fit model

    h -> Fit("model","RL");

    // Draw fit over synthetic data

    model -> SetLineColor(kRed);
    model -> Draw("same");

    // Legend

    TLegend *l = new TLegend(0.6,1,0.6,0.6);
    l -> AddEntry(h,"Simulated Spectrum","l");
    l -> AddEntry(model,"Peaks fit","l");

    l -> Draw();

    // Save image

    a -> SaveAs("alpha_decay_simulation.png");

    // Printing fitted parameters

    cout << "Fitted Parameters:" << endl;
    for (int i = 0; i<number_of_params;i++) {
        cout << i << ":" << model->GetParName(i) << "=" <<model->GetParameter(i)
        << "±" << model->GetParError(i) << endl;
    }

    
}

