
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
        default: return 0.511;//electron
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

void sort_orig(tree_data &simorig); //сортирует фотоны по энергиям

//..............................................................................................

int approximation__tree() {

    const char* filename = "strips_run22";

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
    
    std::cout << "Total Num: "<< tree->GetEntries() << "\nEntry?\n";
    TFile *fout = new TFile(std::string(filename).append("_reco.root").c_str(), "recreate");
    TTree *tout = new TTree("reco_lxe", "reco_lxe");

    reco_data reco[3];
    int clust_size0;
    int clust_size1;

    tout->Branch("nsim",    &sim_orig.nsim,   "nsim/I");
    tout->Branch("simtype", sim_orig.simtype, "simtype[nsim]/I");
    tout->Branch("simorig", sim_orig.simorig, "simorig[nsim]/I");
    tout->Branch("simmom",  sim_orig.simmom,  "simmom[nsim]/F");
    tout->Branch("simphi",  sim_orig.simphi,  "simphi[nsim]/F");
    tout->Branch("simtheta",sim_orig.simtheta,"simtheta[nsim]/F");

    //
    tout->Branch("clust_size0", &clust_size0, "clust_size0/I");
    tout->Branch("clust_size1", &clust_size1, "clust_size1/I");
    //
    tout->Branch("reco_one", &(reco[1]));
    tout->Branch("reco_two", &(reco[2]));
    tout->Branch("reco_optim", &(reco[0]));

    for(int i=0; i<tree->GetEntries(); i++) {

        tree->GetEntry(i);
        std::cout << std::endl << i << " ";

            // bool good_entry = true; //проверка что частицы не улетают из детектора
            //  for(int i=0; good_entry && i<sim_orig.nsim; i++) {
            //     if (sim_orig.simorig[i] == 0) continue;
                
            //     double my_theta = Pi()/2 - sim_orig.simtheta[i];
            //     double phi = sim_orig.simphi[i];
            //     double tg = Tan(my_theta);
                

            //     if (Abs(tg) > 50./38.) /*can't see*/ good_entry = false;
            //  }
            // if (!good_entry) continue;

        reconstructor rcn(strd);
        //if(rcn.clusters[0].empty() || rcn.clusters[1].empty()) continue;
        //sort_orig(sim_orig); //сортировка фотонов по убыванию энергии

        clust_size0 = rcn.clusters[0].size();
        clust_size1 = rcn.clusters[1].size();

        for(int nl=1; nl<=2; nl++) {
            rcn.approximate_N_lines(nl);
            transport_to_reco(reco[nl], rcn, sim_orig);
        }

        {
            rcn.find_first_approximation();
            transport_to_reco(reco[0], rcn, sim_orig);

            std::cout <<  " " << reco[0].num_lines << " | " << reco[0].energy <<  " " << reco[0].mass; 
        }

        tout->Fill();

    }
    std::cout << "\n\n";
    TCanvas *canv = new TCanvas("reco_numlines");
    tout->Draw("reco_optim.num_lines");
    canv->Write();

    fout->Write();
    fout->Close();
    return 0;
}





//............................................................


//перемещаем данные из реконструктора в структуру для хранения внутри дерева
void transport_to_reco(reco_data & reco, reconstructor& rcn, const tree_data & sim_orig) {
   
   reconstructor::coordinates coord = rcn.getCoordinates(); //ищет координаты, может переставлять прямые и амплитуды. Метод НЕ константный.
    reco.num_lines  = coord.phi.size();

    reco.phi.resize(reco.num_lines);
    reco.theta.resize(reco.num_lines);
    reco.amp.resize(reco.num_lines);

    reco.minSigma = Sqrt(rcn.minD());

    for(int j=0; j<reco.num_lines; j++) {
        reco.phi[j] = coord.phi[j];
        reco.theta[j] = Pi()/2. - ATan(coord.tg_theta[j]);
    }
    std::vector<double> amps = rcn.getAmplitudes();
    double ampSum = 0;
    for(int j=0; j<reco.num_lines; j++) { 
        reco.amp[j] = amps[j];
        ampSum += amps[j];
    }

    reco.amp_chisq = rcn.getAmplitudeChisq();
    
    reco.E.resize(reco.num_lines);
    reco.energy = Sqrt(  sim_orig.simmom[0]*sim_orig.simmom[0] + mass_by_id(sim_orig.simtype[0])*mass_by_id(sim_orig.simtype[0])  );

    reco.particle = TLorentzVector();
    for(int k=0; k< reco.num_lines; k++) {
        reco.E[k] = reco.energy /ampSum * reco.amp[k];
        reco.particle += VectorInPolar(reco.E[k], reco.theta[k], reco.phi[k],  0);
    }
    
    
    reco.mass   = reco.particle.M();
}




void sort_orig(tree_data &simorig) //сортирует фотоны по энергиям
{
    for(bool sorted = false; !sorted; ) {
        sorted = true;
        for(int i=1; i<simorig.nsim; i++) {
            if (simorig.simorig[i-1] != 0 && simorig.simorig[i] == 0 ) goto swap; 
            if ((simorig.simorig[i-1] != 0 && simorig.simorig[i] != 0) && simorig.simmom[i-1] < simorig.simmom[i]) goto swap;
            
            continue;
            swap: {
                sorted = false;
                std::swap(simorig.simorig[i-1], simorig.simorig[i]);
                std::swap(simorig.simtype[i-1], simorig.simtype[i]);
                std::swap( simorig.simmom[i-1],  simorig.simmom[i]);
                std::swap( simorig.simphi[i-1],  simorig.simphi[i]);
                std::swap(simorig.simtheta[i-1], simorig.simtheta[i]);
            }            
        }
    }
}
