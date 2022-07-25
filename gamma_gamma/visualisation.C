#include "reconstructor.h" //<<-----ВСЁ ЗДЕСЬ!
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

TLorentzVector VectorInPolar(double rho, double theta, double phi, double mass)
{
    TVector3 p; 
    double Energy = Sqrt(rho*rho + mass*mass);
    p.SetMagThetaPhi(rho, theta, phi);
    return TLorentzVector(p, Energy);
}

double mass_by_id(int type) {
    switch(type){
        case  22: return 0;
        case 111: return 134.977;
        default: return 1E3;//proton
    }
}

//.....................................................
namespace gera_nm {

const unsigned int MAX_SIM = 100;
typedef struct {
    int nsim;
    int simtype[MAX_SIM];// particle ID from GEANT
    int simorig[MAX_SIM];// Origin of particle. 0 - from e+e- or particle ID from GEANT  
    Float_t simmom[MAX_SIM];
    Float_t simphi[MAX_SIM];
    Float_t simtheta[MAX_SIM];
    Float_t simvtx[MAX_SIM];// vertex X for origin
    Float_t simvty[MAX_SIM];//vertex Y
    Float_t simvtz[MAX_SIM];//vertex Z

} tree_data;
}


struct reco_data {
    int num_lines;
    double minSigma;
    std::vector<double> phi;
    std::vector<double> theta;
    std::vector<double> amp;
    double amp_chisq;
    std::vector<double> E; 
    
    TLorentzVector particle;
    double mass;
    double energy;

};

//*(*)




int get_sim_data(TTree *tree, tree_data& sim_orig)        //считываем оригинальные параметры МК-генератора
{  
    tree->SetBranchAddress("nsim", &sim_orig.nsim);
    tree->SetBranchAddress("simtype", &sim_orig.simtype);
    tree->SetBranchAddress("simorig", &sim_orig.simorig);
    tree->SetBranchAddress("simmom", &sim_orig.simmom);
    tree->SetBranchAddress("simphi", &sim_orig.simphi);
    tree->SetBranchAddress("simtheta", &sim_orig.simtheta);
    tree->SetBranchAddress("simvtx", &sim_orig.simvtx);
    tree->SetBranchAddress("simvty", &sim_orig.simvty);
    tree->SetBranchAddress("simvtz", &sim_orig.simvtz);
    return 0;
}

//*****



//перемещаем данные из реконструктора в структуру для хранения внутри дерева
void transport_to_reco(reco_data & reco, reconstructor& rcn, const tree_data &sim_orig);


//..............................................................................................

int visualisation() {

    const char* filename = "strips_run_gg";

    TFile *file = new TFile(std::string(filename).append(".root").c_str(), "read");
    if (!file->IsOpen()) {
        std::cerr << "FILE NOT FOUND\n";
        return 0;
    }
    TTree *tree = (TTree *)file->Get("tr_lxe");

    const strip_data *strd[] = {0, 0};
    tree_data sim_orig;


    tree->SetBranchAddress("strips_dir0", &strd[0]);
    tree->SetBranchAddress("strips_dir1", &strd[1]);
    get_sim_data(tree, sim_orig);
    
    int entr = 0;
    std::cout << "Total Num: "<< tree->GetEntries() << "\nEntry?\n";
    std::cin >> entr;
   
    tree->GetEntry( entr );

    reconstructor *reco = new reconstructor(strd);
    reco->find_first_approximation();
    //std::cout << sqrt(reco->appr[0]->minD) << " " << sqrt(reco->appr[1]->minD) << std::endl;
    
    auto h0 = get_hist("h0", strd[0]);
    auto h1 = get_hist("h1", strd[1]);
    h0->SetStats(false);
    h1->SetStats(false);

    TH2D *hic0 = nullptr, *hic1 = nullptr;
    {
        const char* dir = "Direction +";
        hic0 = new TH2D( dir, std::string(dir).append(";R;\\phi;a").c_str(), LAYERS, rho_begin, rho_end,  BINS_PER_LAYER, 0, 2*Pi() );
        auto clust = cluster_finder(strd[0]);
        for(int i=0; i<clust.size(); i++) {
            for(int j=0; j<clust[i].size(); j++) {
                hic0->Fill(rho[clust[i][j].layer], clust[i][j].start_phi, clust[i][j].amp);
            }
        }
    }
    {
        const char* dir = "Direction -";
        hic1 = new TH2D( dir, std::string(dir).append(";R;\\phi;a").c_str(), LAYERS, rho_begin, rho_end,  BINS_PER_LAYER, 0, 2*Pi() );
        auto clust = cluster_finder(strd[1]);
        for(int i=0; i<clust.size(); i++) {
            for(int j=0; j<clust[i].size(); j++) {
                hic1->Fill(rho[clust[i][j].layer], clust[i][j].start_phi, clust[i][j].amp);
            }
        }
    }
    hic0->SetStats(false);
    hic1->SetStats(false);

    TCanvas *canv__init = new TCanvas();
    canv__init->Divide(2, 1);

    canv__init->cd(1);
    //visualizate("direction +", strd[0]);
    h0->Draw("colz");

    canv__init->cd(2);
    //visualizate("direction -", strd[1]);
    h1->Draw("colz"); 

    
    TCanvas *canv0 = new TCanvas();
    canv0->Divide(2, 1);

    canv0->cd(1);
    h0->Draw("colz");
     for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[0]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }

    canv0->cd(2);
    h1->Draw("colz");
     for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[1]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }
    

    TCanvas *canv1 = new TCanvas();
    canv1->Divide(3, 1);
    canv1->cd(1);
    reco->approximate_N_lines(1);
    //std::cout << sqrt(reco->appr[0]->minD) << " " << sqrt(reco->appr[1]->minD) << std::endl;
    {
        hic0->Draw("colz");

        for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[0]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }
    }
    canv1->cd(2);
    reco->approximate_N_lines(2);
    //std::cout << sqrt(reco->appr[0]->minD) << " " << sqrt(reco->appr[1]->minD) << std::endl;
    {
        hic0->Draw("colz");

        for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[0]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }
    }

    canv1->cd(3);
    reco->approximate_N_lines(3);
    //std::cout << sqrt(reco->appr[0]->minD) << " " << sqrt(reco->appr[1]->minD) << std::endl;
    {
        hic0->Draw("colz");

        for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[0]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }
    }







    TCanvas *canv2 = new TCanvas();
    canv2->Divide(3, 1);
    canv2->cd(1);

    reco->approximate_N_lines(1);
    {
        hic1->Draw("colz");

        for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[1]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }
    }
    canv2->cd(2);
    reco->approximate_N_lines(2);
    {
        hic1->Draw("colz");

        for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[1]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }
    }

    canv2->cd(3);
    reco->approximate_N_lines(3);
    {
        hic1->Draw("colz");

        for(int k=0; k<reco->num_lines; k++) {
            TF1* f = &reco->appr[1]->func[k];
            f->SetLineWidth(2);
            f->DrawCopy("same");
        }
    }






    // {
    //     const char* dir = "direction -";
    //     TH2D *hist = new TH2D( dir, std::string(dir).append(";R;\\phi;a").c_str(), LAYERS, rho_begin, rho_end,  BINS_PER_LAYER, 0, 2*Pi() );
    //     auto clust = cluster_finder(strd[1]);
    //     for(int i=0; i<clust.size(); i++) {
    //         for(int j=0; j<clust[i].size(); j++) {
    //             hist->Fill(rho[clust[i][j].layer], clust[i][j].start_phi, clust[i][j].amp);
    //         }
    //     }
    //     hist->Draw("colz");
    //     for(int k=0; k<reco->num_lines; k++) {
    //         TF1* f = &reco->appr[1]->func[k];
    //         f->SetLineWidth(2);
    //         f->Draw("same");
    //     }
    // }
    

    return 0;
}





//............................................................









