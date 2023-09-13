#include <math.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>

using namespace std;


const double lumi_mc16a = 36646.74;
const double lumi_mc16d = 44630.6;
const double lumi_mc16e = 58791.6;

void set_uo_flow(TH1D *hist) {
    int NBins = hist -> GetNbinsX();

    hist -> SetBinContent(1, hist -> GetBinContent(0) + hist -> GetBinContent(1));
    hist -> SetBinContent(NBins, hist -> GetBinContent(NBins) + hist -> GetBinContent(NBins + 1));

    hist -> SetBinError(1, sqrt(pow(hist -> GetBinError(0), 2) + pow(hist -> GetBinError(0), 2)));
    hist -> SetBinError(NBins, sqrt(pow(hist -> GetBinError(NBins), 2) + pow(hist -> GetBinError(NBins + 1), 2)));
}

void converter() {
    string outFName = "/workenv/converter/output/Zllgam.root";

    vector<string> inFilenames = {
        "/workenv/converter/source/base/user.akurova.MC16a.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/base/user.akurova.MC16a.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/base/user.akurova.MC16a.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/base/user.akurova.MC16d.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/base/user.akurova.MC16d.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/base/user.akurova.MC16d.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/base/user.akurova.MC16e.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/base/user.akurova.MC16e.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/base/user.akurova.MC16e.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.reproc-21-02-23_MxAOD.root/user.akurova.MxAOD.root",
    };

    int NBins = 20;
    double left_border = 0, right_border = 1000;

    double sum_of_weights_bk_xAOD, norm_koef, sum_weights, weight_coef;
    double weight;

    UInt_t n_jet, n_ph, n_mu, n_e_medium, ph_isem;
    double ph_pt, ph_phi, ph_eta, metTST_pt, jet_lead_pt, jet_sublead_pt, dphi_jj, dRj1gam, dphi_Zj, mT_Zg, soft_term_pt, metTSTsignif;
    double jet_lead_eta, jet_lead_phi, jet_lead_E, jet_sublead_eta, jet_sublead_phi, jet_sublead_E;
    double metTST_phi, ph_iso_et20, ph_iso_pt, ph_z_point;
    TLorentzVector met, ph, jet, jet2;

    TFile *outFile = new TFile(TString(outFName.data()), "RECREATE");

    TH1D *hist_SR = new TH1D ("SR", "SR", NBins, left_border, right_border);
    TH1D *hist_Wg = new TH1D ("Wg", "Wg", NBins, left_border, right_border);
    TH1D *hist_GammaJet = new TH1D ("GammaJet", "GammaJet", NBins, left_border, right_border);

    for (auto inFName : inFilenames) {
        TFile *inFile = new TFile(TString(inFName.data()), "READ");
        TTree* inTree = (TTree*)inFile->Get("output_tree");

        cout << inFName << endl;

        // расчёт коэффициента для веса =========================================
        TTree *inSWTree = (TTree*)inFile->Get("output_tree_sw");
        TTree *inNormTree = (TTree*)inFile->Get("norm_tree");

        inSWTree->SetBranchAddress("sum_of_weights_bk_xAOD", &sum_of_weights_bk_xAOD);
        inNormTree->SetBranchAddress("koef", &norm_koef);

        inNormTree->GetEntry(0);
        sum_weights = 0;
        for (int i=0; i<inSWTree->GetEntries(); i++) {
            inSWTree->GetEntry(i);
            sum_weights += sum_of_weights_bk_xAOD;
        }

        weight_coef = norm_koef/sum_weights;

        if (inFName.find("MC16a") != string::npos) weight_coef *= lumi_mc16a;
        else if (inFName.find("MC16d") != string::npos) weight_coef *= lumi_mc16d;
        else if (inFName.find("MC16e") != string::npos) weight_coef *= lumi_mc16e;
        else exit(1);
        // ======================================================================

        inTree -> SetBranchAddress("n_ph", &n_ph);
        inTree -> SetBranchAddress("n_jet", &n_jet);
        inTree -> SetBranchAddress("n_mu", &n_mu);
        inTree -> SetBranchAddress("n_e_looseBL", &n_e_medium);
        inTree -> SetBranchAddress("ph_isem", &ph_isem);

        inTree -> SetBranchAddress("ph_pt", &ph_pt);
        inTree -> SetBranchAddress("ph_phi", &ph_phi);
        inTree -> SetBranchAddress("ph_eta", &ph_eta);

        inTree -> SetBranchAddress("jet_lead_pt", &jet_lead_pt);  //leading jet p_x
        inTree -> SetBranchAddress("jet_lead_eta", &jet_lead_eta);  //p_y
        inTree -> SetBranchAddress("jet_lead_phi", &jet_lead_phi);  //p_z
        inTree -> SetBranchAddress("jet_lead_E", &jet_lead_E);    //E

        inTree -> SetBranchAddress("jet_sublead_pt", &jet_sublead_pt);  //leading jet p_x
        inTree -> SetBranchAddress("jet_sublead_eta", &jet_sublead_eta);  //p_y
        inTree -> SetBranchAddress("jet_sublead_phi", &jet_sublead_phi);  //p_z
        inTree -> SetBranchAddress("jet_sublead_E", &jet_sublead_E);    //E

        inTree -> SetBranchAddress("metTST_pt", &metTST_pt);  //MET p_x
        inTree -> SetBranchAddress("metTST_phi", &metTST_phi);  //p_y

        inTree -> SetBranchAddress("ph_iso_et20", &ph_iso_et20);
        inTree -> SetBranchAddress("ph_iso_pt", &ph_iso_pt);

        inTree -> SetBranchAddress("ph_z_point", &ph_z_point);
        inTree -> SetBranchAddress("metTSTsignif", &metTSTsignif);
        inTree -> SetBranchAddress("soft_term_pt", &soft_term_pt);

        inTree -> SetBranchAddress("weight", &weight);

        for (unsigned int i = 0; i < inTree->GetEntries(); i++) {
            inTree->GetEntry(i);

            weight *= weight_coef;

            int n_lep = n_e_medium + n_mu;
            double IsoVar = ph_iso_et20/ph_pt;

            jet.SetPtEtaPhiE(jet_lead_pt, jet_lead_eta, jet_lead_phi, jet_lead_E);
            jet2.SetPtEtaPhiE(jet_sublead_pt, jet_sublead_eta, jet_sublead_phi, jet_sublead_E);
            met.SetPtEtaPhiM(metTST_pt, 0, metTST_phi, 0);
            ph.SetPtEtaPhiE(ph_pt, ph_eta, ph_phi, ph_pt);

            // сюда нужно добавить отборы

            if(metTSTsignif >= 11 && n_lep == 0) hist_SR->Fill(ph_pt, weight);
            else if(metTSTsignif >= 11 && n_lep != 0) hist_Wg->Fill(ph_pt, weight);
            else if(metTSTsignif <= 11 && n_lep == 0) hist_GammaJet->Fill(ph_pt, weight);
        }

        inFile->Close();
    }
    set_uo_flow(hist_SR);
    set_uo_flow(hist_Wg);
    set_uo_flow(hist_GammaJet);

    outFile->Write("", TObject::kOverwrite);
    outFile->Close();
}