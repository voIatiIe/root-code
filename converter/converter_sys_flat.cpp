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

void converter_sys_flat() {
    string outFName = "/workenv/converter/output/Zllgam_flat_syst.root";
    vector<string> inFilenames = {
        "/workenv/converter/source/flat_syst/user.esoldato.MC16a.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/flat_syst/user.esoldato.MC16a.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/flat_syst/user.esoldato.MC16a.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/flat_syst/user.esoldato.MC16d.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/flat_syst/user.esoldato.MC16d.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/flat_syst/user.esoldato.MC16d.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/flat_syst/user.esoldato.MC16e.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/flat_syst/user.esoldato.MC16e.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/flat_syst/user.esoldato.MC16e.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.FLAT_fR2_v4_MxAOD.root/user.akurova.MxAOD.root",
    };
    vector<string> treeNames = {
        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down",
        "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up",
        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1down",
        "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1up",
        "JET_JvtEfficiency__1down",
        "JET_JvtEfficiency__1up",
        "PRW_DATASF__1down",
        "PRW_DATASF__1up",
    };

    int NBins = 20;
    double left_border = 0, right_border = 1000;

    double sum_of_weights_bk_xAOD, norm_koef, sum_weights, weight_coef;
    double weight;

    UInt_t n_jet, n_lep;
    double ph_pt, ph_eta, metTST_pt, jet_lead_pt, jet_sublead_pt, dphi_jj, dRj1gam, dphi_Zj, mT_Zg, soft_term_pt, metTSTsignif;
    vector<double>* weight_vec_exp = new vector<double>;

    for (auto treeName : treeNames) {
        TFile *outFile = new TFile(TString(outFName.data()), "UPDATE");

        cout << treeName << endl;
        TH1D *hist_SR = new TH1D (TString(("SR_" + treeName).data()), TString(("SR_" + treeName).data()), NBins, left_border, right_border);
        TH1D *hist_Wg = new TH1D (TString(("Wg_" + treeName).data()), TString(("Wg_" + treeName).data()), NBins, left_border, right_border);
        TH1D *hist_GammaJet = new TH1D (TString(("GammaJet_" + treeName).data()), TString(("GammaJet_" + treeName).data()), NBins, left_border, right_border);

        for (auto inFName : inFilenames) {

            TFile* inFile = new TFile(TString(inFName.data()), "READ");
            TTree* inSysTree = (TTree*)inFile->Get(TString(("output_tree_sys_" + treeName).data()));
            TTree* inDefaultTree = (TTree*)inFile->Get("output_tree_sys_default");

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

            inDefaultTree -> SetBranchAddress("n_jet", &n_jet);
            inDefaultTree -> SetBranchAddress("n_lep", &n_lep);
            inDefaultTree -> SetBranchAddress("ph_et", &ph_pt);
            inDefaultTree -> SetBranchAddress("ph_eta", &ph_eta);
            inDefaultTree -> SetBranchAddress("metTST_pt", &metTST_pt);
            inDefaultTree -> SetBranchAddress("jet_lead_pt", &jet_lead_pt);
            inDefaultTree -> SetBranchAddress("jet_sublead_pt", &jet_sublead_pt);
            inDefaultTree -> SetBranchAddress("dphi_jj", &dphi_jj);
            inDefaultTree -> SetBranchAddress("dRj1gam", &dRj1gam);
            inDefaultTree -> SetBranchAddress("dphi_Zj", &dphi_Zj);
            inDefaultTree -> SetBranchAddress("mT_Zg", &mT_Zg);
            inDefaultTree -> SetBranchAddress("soft_term_pt", &soft_term_pt);
            inDefaultTree -> SetBranchAddress("metTSTsignif", &metTSTsignif);
            inDefaultTree -> SetBranchAddress("weight_vec_exp", &weight_vec_exp);

            inSysTree -> SetBranchAddress("weight", &weight);

            for (unsigned int i = 0; i < inSysTree->GetEntries(); i++) {
                inDefaultTree->GetEntry(i);
                inSysTree->GetEntry(i);

                weight *= weight_coef;

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
}
