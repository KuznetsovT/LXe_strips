#include <iostream>
#include <fstream>
#include <map>
#include <sstream>

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"


using namespace TMath;



int FIND_GG(int ebeam)
{
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
    std::map<int, std::string> selected;
    std::stringstream sout;
    
    
    int total_counter = 0;
    int counter = 0;
    for(int i=0; i<selnum; i++, counter++, total_counter++) {
        if (i==0)
        {
             sout << " -m \"";
             out << runnum[i] << "    :   \"";
             sout << evnum[i];  out << evnum[i];
             continue;
        }   
        if (runnum[i] != runnum[i-1]) {
                 out << "\"     -n     " << counter << "\n\n\n";
                 sout << "\" -n " << counter;          
                    
                 selected[runnum[i-1]] = sout.str();
                 if (total_counter > 1000) goto terminate;

                 counter = 0;
                 sout.str(std::string());
                 sout << " -m \"";
                 out << runnum[i] <<  "   :   \"";
                 sout << evnum[i]; out << evnum[i];
        } else { 
                 out <<  ';' <<  evnum[i]; 
                 sout <<  ';' <<  evnum[i]; 
        }
    }
        
    out << "\"     -n     " << counter << "\n\n\n";
    sout << "\" -n " << counter;
        
    selected[runnum[selnum-1]] = sout.str();

terminate:

    out.close();
    
    /*for(auto const &item : selected) {
           std::cout << item.first << ":    " << item.second << std::endl;
    }*/
    std::cout << "\ntotal: " <<  total_counter << "\n\n";

    TPython::LoadMacro("raw2tree.py");
    TPython::Exec(Form("main(%d, \'%s\')", selected.begin()->first, selected.begin()->second.c_str()));

    return 0;
}



