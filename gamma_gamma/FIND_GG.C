#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <utility>

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"


using namespace TMath;



int FIND_GG(int ebeam)
{
    gSystem->Exec(Form("mkdir e%d", ebeam));
    gSystem->Exec(Form("rm -f e%d/*",ebeam)); 
    
    auto tr_ph = new TChain( "tr_ph" );
    tr_ph->Add( Form("root://cmd//scan2013_rho/scan2013_rho_tr_ph_fc_e%d_v8.root", ebeam) );
    
    std::cout << "BEGIN?\n";
    
    //see -> https://cmd.inp.nsk.su/docdb/files/0003/000381/003/ryzhenenkov_lum2021.pdf
    //select gamma gamma
    const char selection[] = "(nt == 0) && (nph == 2)   && (Abs(Abs(phphi[1] - phphi[0]) - Pi()) < 0.2) && (Abs(phth[1] + phth[0] - Pi()) < 0.6) &&  (  (Pi() - 1.0 > (phth[0] +Pi() - phth[1])/2) && ((phth[0] +Pi() - phth[1])/2 > 1.0)  )   &&   (0.5*emeas < phen[0]) && (0.5*emeas < phen[1]) && (phen[0] < 1.5*emeas) && (phen[1] < 1.5*emeas) ";
    
    std::cout << "BEFORE?\n";
    
    tr_ph->Draw( "evnum:runnum",  selection);
    
    std::cout << "DRAW?\n";

    auto selnum = tr_ph->GetSelectedRows();
    auto evnum = tr_ph->GetV1();
    auto runnum = tr_ph->GetV2();
    
    std::ofstream out(Form("scan2013_rho_tr_ph_fc_e%d_v8.txt", ebeam) );
    std::vector<std::pair<int, std::pair<std::string, int> > > selected;
    std::stringstream sout;
    
    
    int counter = 0;
    for(int i=0; i<selnum; i++, counter++) {
        if (i==0)
        {
             sout << '\"';
             out << runnum[i] << "    :   \"";
             sout << evnum[i];  out << evnum[i];
             continue;
        }   
        if (runnum[i] != runnum[i-1]) {
                 out << "\"     -n     " << counter << "\n\n\n";
                 sout << '\"';   // -n " << counter;          
                    
                 selected.push_back(std::make_pair(runnum[i-1], std::make_pair(sout.str(), counter)));
                 
                 counter = 0;
                 sout.str(std::string());
                 sout << '\"';
                 out << runnum[i] <<  "   :   \"";
                 sout << evnum[i]; out << evnum[i];
        } else { 
                 out <<  ';' <<  evnum[i]; 
                 sout <<  ';' <<  evnum[i]; 
        }
    }
        
    out << "\"     -n     " << counter << "\n\n\n";
    sout << '\"';   // -n " << counter;
        
    selected.push_back(std::make_pair(runnum[selnum-1], std::make_pair(sout.str(), counter)));

    out.close();


    std::sort(selected.begin(), selected.end(), [](const std::pair<int, std::pair<std::string, int>> & a, const std::pair<int,std::pair<std::string, int>> & b) -> bool {
            return a.second.second > b.second.second;
        });
    
    /*for(auto const &item : selected) {
           std::cout << item.first << ":    " << item.second.second << std::endl;
    }*/
    
    
    TPython::LoadMacro("raw2tree.py");
    
    int total_counter = 0;
    for(auto const &item : selected) {
        std::cout << item.first << ":    " << item.second.second << std::endl;
        TPython::Exec(Form("main(%d, %s, %d)", item.first, item.second.first.c_str(), item.second.second));
        gSystem->Exec(Form("mv tmp/online.raw.v1.%d.root e%d/", item.first, ebeam));
        total_counter += item.second.second;
        if (total_counter > 1000) break;
    }
    std::cout << "\n\ntotal count: "<< total_counter << "\n\n\n";
    gSystem->Exec(Form("hadd e%d/raw.online.e%d.root e%d/*", ebeam, ebeam, ebeam)); 
    gSystem->Exec(Form("mv scan2013_rho_tr_ph_fc_e%d_v8.txt e%d/", ebeam, ebeam) );
    
    return 0;
}



