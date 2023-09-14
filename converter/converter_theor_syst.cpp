#include <math.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <TH1.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>

using namespace std;


const double lumi_mc16a = 36646.74;
const double lumi_mc16d = 44630.6;
const double lumi_mc16e = 58791.6;

unordered_map<std::string, int> theorSystMapping = {
    {"ABC", 5},
    {"DEF", 6},
    {"GHI", 7},
};

void set_uo_flow(TH1D *hist) {
    int NBins = hist -> GetNbinsX();

    hist -> SetBinContent(1, hist -> GetBinContent(0) + hist -> GetBinContent(1));
    hist -> SetBinContent(NBins, hist -> GetBinContent(NBins) + hist -> GetBinContent(NBins + 1));

    hist -> SetBinError(1, sqrt(pow(hist -> GetBinError(0), 2) + pow(hist -> GetBinError(1), 2)));
    hist -> SetBinError(NBins, sqrt(pow(hist -> GetBinError(NBins), 2) + pow(hist -> GetBinError(NBins + 1), 2)));
}

void converter_theor_syst() {
    string outFName = "/workenv/converter/output/Zllgam_theor_syst.root";

    vector<string> inFilenames = {
        "/workenv/converter/source/theor_syst/user.esoldato.MC16a.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/theor_syst/user.esoldato.MC16a.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/theor_syst/user.esoldato.MC16a.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/theor_syst/user.esoldato.MC16d.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/theor_syst/user.esoldato.MC16d.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/theor_syst/user.esoldato.MC16d.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/theor_syst/user.esoldato.MC16e.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/theor_syst/user.esoldato.MC16e.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/theor_syst/user.esoldato.MC16e.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.syst_TH_fR2_v0_MxAOD.root/user.akurova.MxAOD.root",
    };
    vector<string> systNames = {
        "ABC",
        "DEF",
        "GHI",
    };

    int NBins = 20;
    double left_border = 0, right_border = 1000;

    double sum_of_weights_bk_xAOD, norm_koef, sum_weights, weight_coef;
    double weight;

    UInt_t n_jet, n_ph, n_mu, n_e_medium, ph_isem;
    double ph_pt, ph_phi, ph_eta, metTST_pt, jet_lead_pt, jet_sublead_pt, dphi_jj, dRj1gam, dphi_Zj, mT_Zg, soft_term_pt, metTSTsignif;
    double jet_lead_eta, jet_lead_phi, jet_lead_E, jet_sublead_eta, jet_sublead_phi, jet_sublead_E;
    double metTST_phi, ph_iso_et20, ph_iso_pt, ph_z_point;
    vector<double>* weight_vec = new vector<double>;
    TLorentzVector met, ph, jet, jet2;

    Double_t xbins[11] = {150, 180, 210, 240, 270, 300, 340, 380, 430, 510, 600};

    for (auto systName : systNames) {
        int systWeightId = theorSystMapping[systName];

        TFile *outFile = new TFile(TString(outFName.data()), "UPDATE");

        cout << systName << endl;

        TH1D *hist_SR = new TH1D (TString(("SR_" + systName).data()), TString(("SR_" + systName).data()), NBins, left_border, right_border);
        TH1D *hist_Wg = new TH1D (TString(("Wg_" + systName).data()), TString(("Wg_" + systName).data()), NBins, left_border, right_border);
        TH1D *hist_GammaJet = new TH1D (TString(("GammaJet_" + systName).data()), TString(("GammaJet_" + systName).data()), NBins, left_border, right_border);

        hist_SR = dynamic_cast<TH1D*>(hist_SR->Rebin(NBins, "", xbins));
        hist_Wg = dynamic_cast<TH1D*>(hist_Wg->Rebin(NBins, "", xbins));
        hist_GammaJet = dynamic_cast<TH1D*>(hist_GammaJet->Rebin(NBins, "", xbins));

        for (auto inFName : inFilenames) {
            TFile* inFile = new TFile(TString(inFName.data()), "READ");
            TTree* inSysTree = (TTree*)inFile->Get("output_tree");

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

            inSysTree -> SetBranchAddress("n_ph", &n_ph);
            inSysTree -> SetBranchAddress("n_jet", &n_jet);
            inSysTree -> SetBranchAddress("n_mu", &n_mu);
            inSysTree -> SetBranchAddress("n_e_looseBL", &n_e_medium);
            inSysTree -> SetBranchAddress("ph_isem", &ph_isem);

            inSysTree -> SetBranchAddress("ph_pt", &ph_pt);
            inSysTree -> SetBranchAddress("ph_phi", &ph_phi);
            inSysTree -> SetBranchAddress("ph_eta", &ph_eta);

            inSysTree -> SetBranchAddress("jet_lead_pt", &jet_lead_pt);  //leading jet p_x
            inSysTree -> SetBranchAddress("jet_lead_eta", &jet_lead_eta);  //p_y
            inSysTree -> SetBranchAddress("jet_lead_phi", &jet_lead_phi);  //p_z
            inSysTree -> SetBranchAddress("jet_lead_E", &jet_lead_E);    //E

            inSysTree -> SetBranchAddress("jet_sublead_pt", &jet_sublead_pt);  //leading jet p_x
            inSysTree -> SetBranchAddress("jet_sublead_eta", &jet_sublead_eta);  //p_y
            inSysTree -> SetBranchAddress("jet_sublead_phi", &jet_sublead_phi);  //p_z
            inSysTree -> SetBranchAddress("jet_sublead_E", &jet_sublead_E);    //E

            inSysTree -> SetBranchAddress("metTST_pt", &metTST_pt);  //MET p_x
            inSysTree -> SetBranchAddress("metTST_phi", &metTST_phi);  //p_y

            inSysTree -> SetBranchAddress("ph_iso_et20", &ph_iso_et20);
            inSysTree -> SetBranchAddress("ph_iso_pt", &ph_iso_pt);

            inSysTree -> SetBranchAddress("ph_z_point", &ph_z_point);
            inSysTree -> SetBranchAddress("metTSTsignif", &metTSTsignif);
            inSysTree -> SetBranchAddress("soft_term_pt", &soft_term_pt);

            inSysTree -> SetBranchAddress("weight_vec", &weight_vec);


            for (unsigned int i = 0; i < inSysTree->GetEntries(); i++) {
                inSysTree -> GetEntry(i);

                weight = (*weight_vec)[systWeightId];

                int n_lep = n_e_medium + n_mu;
                double IsoVar = ph_iso_et20/ph_pt;

                jet.SetPtEtaPhiE(jet_lead_pt, jet_lead_eta, jet_lead_phi, jet_lead_E);
                jet2.SetPtEtaPhiE(jet_sublead_pt, jet_sublead_eta, jet_sublead_phi, jet_sublead_E);
                met.SetPtEtaPhiM(metTST_pt, 0, metTST_phi, 0);
                ph.SetPtEtaPhiE(ph_pt, ph_eta, ph_phi, ph_pt);

                // сюда нужно добавить отборы

                weight *= weight_coef;

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
        outFile -> Close();
    }
}
