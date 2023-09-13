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

void converter_non_flat_syst() {
    string outFName = "/workenv/converter/output/Zllgam_non_flat_syst.root";

    vector<string> inFilenames = {
        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16a.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16a.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16a.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16d.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16d.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16d.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",

        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16e.364504.Sherpa_222_NNPDF30NNLO_eegamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16e.364509.Sherpa_222_NNPDF30NNLO_mumugamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",
        "/workenv/converter/source/non_flat_syst/user.esoldato.MC16e.364514.Sherpa_222_NNPDF30NNLO_tautaugamma_pty_140_E_CMS.syst_nf_all_v2_MxAOD.root/user.akurova.MxAOD.root",
    };
    vector<string> tree_names = {
        "EG_RESOLUTION_ALL__1down",
        "EG_RESOLUTION_ALL__1up",
        "EG_SCALE_ALL__1down",
        "EG_SCALE_ALL__1up",
        // "JET_BJES_Response__1down",
        // "JET_BJES_Response__1up",
        // "JET_EffectiveNP_Detector1__1down",
        // "JET_EffectiveNP_Detector1__1up",
        // "JET_EffectiveNP_Detector2__1down",
        // "JET_EffectiveNP_Detector2__1up",
        // "JET_EffectiveNP_Mixed1__1down",
        // "JET_EffectiveNP_Mixed1__1up",
        // "JET_EffectiveNP_Mixed2__1down",
        // "JET_EffectiveNP_Mixed2__1up",
        // "JET_EffectiveNP_Mixed3__1down",
        // "JET_EffectiveNP_Mixed3__1up",
        // "JET_EffectiveNP_Modelling1__1down",
        // "JET_EffectiveNP_Modelling1__1up",
        // "JET_EffectiveNP_Modelling2__1down",
        // "JET_EffectiveNP_Modelling2__1up",
        // "JET_EffectiveNP_Modelling3__1down",
        // "JET_EffectiveNP_Modelling3__1up",
        // "JET_EffectiveNP_Modelling4__1down",
        // "JET_EffectiveNP_Modelling4__1up",
        // "JET_EffectiveNP_Statistical1__1down",
        // "JET_EffectiveNP_Statistical1__1up",
        // "JET_EffectiveNP_Statistical2__1down",
        // "JET_EffectiveNP_Statistical2__1up",
        // "JET_EffectiveNP_Statistical3__1down",
        // "JET_EffectiveNP_Statistical3__1up",
        // "JET_EffectiveNP_Statistical4__1down",
        // "JET_EffectiveNP_Statistical4__1up",
        // "JET_EffectiveNP_Statistical5__1down",
        // "JET_EffectiveNP_Statistical5__1up",
        // "JET_EffectiveNP_Statistical6__1down",
        // "JET_EffectiveNP_Statistical6__1up",
        // "JET_EtaIntercalibration_Modelling__1down",
        // "JET_EtaIntercalibration_Modelling__1up",
        // "JET_EtaIntercalibration_NonClosure_2018data__1down",
        // "JET_EtaIntercalibration_NonClosure_2018data__1up",
        // "JET_EtaIntercalibration_NonClosure_highE__1down",
        // "JET_EtaIntercalibration_NonClosure_highE__1up",
        // "JET_EtaIntercalibration_NonClosure_negEta__1down",
        // "JET_EtaIntercalibration_NonClosure_negEta__1up",
        // "JET_EtaIntercalibration_NonClosure_posEta__1down",
        // "JET_EtaIntercalibration_NonClosure_posEta__1up",
        // "JET_EtaIntercalibration_TotalStat__1down",
        // "JET_EtaIntercalibration_TotalStat__1up",
        // "JET_Flavor_Composition__1down",
        // "JET_Flavor_Composition__1up",
        // "JET_Flavor_Response__1down",
        // "JET_Flavor_Response__1up",
        // "JET_JER_DataVsMC_MC16__1down",
        // "JET_JER_DataVsMC_MC16__1up",
        // "JET_JER_EffectiveNP_1__1down",
        // "JET_JER_EffectiveNP_1__1up",
        // "JET_JER_EffectiveNP_2__1down",
        // "JET_JER_EffectiveNP_2__1up",
        // "JET_JER_EffectiveNP_3__1down",
        // "JET_JER_EffectiveNP_3__1up",
        // "JET_JER_EffectiveNP_4__1down",
        // "JET_JER_EffectiveNP_4__1up",
        // "JET_JER_EffectiveNP_5__1down",
        // "JET_JER_EffectiveNP_5__1up",
        // "JET_JER_EffectiveNP_6__1down",
        // "JET_JER_EffectiveNP_6__1up",
        // "JET_JER_EffectiveNP_7restTerm__1down",
        // "JET_JER_EffectiveNP_7restTerm__1up",
        // "JET_Pileup_OffsetMu__1down",
        // "JET_Pileup_OffsetMu__1up",
        // "JET_Pileup_OffsetNPV__1down",
        // "JET_Pileup_OffsetNPV__1up",
        // "JET_Pileup_PtTerm__1down",
        // "JET_Pileup_PtTerm__1up",
        // "JET_Pileup_RhoTopology__1down",
        // "JET_Pileup_RhoTopology__1up",
        // "JET_PunchThrough_MC16__1down",
        // "JET_PunchThrough_MC16__1up",
        // "JET_SingleParticle_HighPt__1down",
        // "JET_SingleParticle_HighPt__1up",
        // "MET_SoftTrk_ResoPara",
        // "MET_SoftTrk_ResoPerp",
        // "MET_SoftTrk_Scale__1down",
        // "MET_SoftTrk_Scale__1up",
        // "MUON_ID__1down",
        // "MUON_ID__1up",
        // "MUON_MS__1down",
        // "MUON_MS__1up",
        // "MUON_SAGITTA_RESBIAS__1down",
        // "MUON_SAGITTA_RESBIAS__1up",
        // "MUON_SAGITTA_RHO__1down",
        // "MUON_SAGITTA_RHO__1up",
        // "MUON_SCALE__1down",
        // "MUON_SCALE__1up"
    };

    int NBins = 20;
    double left_border = 0, right_border = 1000;

    double sum_of_weights_bk_xAOD, norm_koef, sum_weights, weight_coef;
    double weight;

    UInt_t n_jet, n_lep;
    double ph_pt, ph_eta, metTST_pt, jet_lead_pt, jet_sublead_pt, dphi_jj, dRj1gam, dphi_Zj, mT_Zg, soft_term_pt, metTSTsignif;

    for (auto treeName : tree_names) {
        TFile *outFile = new TFile(TString(outFName.data()), "UPDATE");

        cout << treeName << endl;
        TH1D *hist_SR = new TH1D (TString(("SR_" + treeName).data()), TString(("SR_" + treeName).data()), NBins, left_border, right_border);
        TH1D *hist_Wg = new TH1D (TString(("Wg_" + treeName).data()), TString(("Wg_" + treeName).data()), NBins, left_border, right_border);
        TH1D *hist_GammaJet = new TH1D (TString(("GammaJet_" + treeName).data()), TString(("GammaJet_" + treeName).data()), NBins, left_border, right_border);

        for (auto inFName : inFilenames) {
            TFile* inFile = new TFile(TString(inFName.data()), "READ");
            TTree* inSysTree = (TTree*)inFile->Get(TString(("output_tree_sys_" + treeName).data()));

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

            inSysTree -> SetBranchAddress("n_jet", &n_jet);
            inSysTree -> SetBranchAddress("ph_et", &ph_pt);
            inSysTree -> SetBranchAddress("ph_eta", &ph_eta);
            inSysTree -> SetBranchAddress("metTST_pt", &metTST_pt);
            inSysTree -> SetBranchAddress("n_lep", &n_lep);
            inSysTree -> SetBranchAddress("jet_lead_pt", &jet_lead_pt);
            inSysTree -> SetBranchAddress("jet_sublead_pt", &jet_sublead_pt);
            inSysTree -> SetBranchAddress("dphi_jj", &dphi_jj);
            inSysTree -> SetBranchAddress("dRj1gam", &dRj1gam);
            inSysTree -> SetBranchAddress("dphi_Zj", &dphi_Zj);
            inSysTree -> SetBranchAddress("mT_Zg", &mT_Zg);
            inSysTree -> SetBranchAddress("metTSTsignif", &metTSTsignif);
            inSysTree -> SetBranchAddress("soft_term_pt", &soft_term_pt);

            inSysTree -> SetBranchAddress("weight", &weight);

            for (unsigned int i = 0; i < inSysTree->GetEntries(); i++) {
                inSysTree -> GetEntry(i);

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
        outFile -> Close();
    }
}
