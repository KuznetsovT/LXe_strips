

#include "approximator.h"  //<<<----ФИТИРОВАНИЕ ЗДЕСЬ!

#include <algorithm>
#include <numeric>

#define THRESHOLD (3.25 *SIGMA)
#define MAX_TAN  55./38.



struct reconstructor {

    int num_lines = 0;
    const strip_data *strd[2];
    std::vector<std::vector<hist_bin>> clusters[2];
    approximator *appr[2] = {NULL, NULL};
      std::vector<approximator::amp_data> ampid[2]; // //отсортированные амплитуды и id 


    reconstructor(const strip_data *strd[2]) {
        this->strd[0] = strd[0];
        this->strd[1] = strd[1]; 
        clusters[0] = cluster_finder(strd[0]);
        clusters[1] = cluster_finder(strd[1]);
    }

    virtual ~reconstructor() { delete (appr[0]); delete (appr[1]); }


    double approximate_N_lines(int num_lines);
    double find_first_approximation();

    double minD() const {
        if (num_lines == 0) return 0;
        else return Max(appr[0]->minD, appr[1]->minD);
    }


    struct coordinates {
        std::vector<double> phi;
        std::vector<double> tg_theta;
    };
    struct coordinates getCoordinates();
    std::vector<double> getAmplitudes() const;
    double getAmplitudeChisq() const;
};







double reconstructor::approximate_N_lines(int num_lines) {
    if (clusters[0].empty() || clusters[1].empty()) {
        this->num_lines = 0;
        return 0;
    }
    this->num_lines = num_lines;
    delete (appr[0]); delete (appr[1]);
    
    appr[0] = new approximator(num_lines, strd[0], clusters[0], 0);
    appr[1] = new approximator(num_lines, strd[1], clusters[1], 1);
    
    double minD = Max(Sqrt(appr[0]->minimizeDisp()), Sqrt(appr[1]->minimizeDisp()));
      ampid[0] = appr[0]->calculate_amplitudes();
      ampid[1] = appr[1]->calculate_amplitudes(); //отсортированные амплитужы и id
    return minD;
}



double reconstructor::find_first_approximation() {
    if (clusters[0].empty() || clusters[1].empty()) {
        this->num_lines = 0;
        return 0;
    }
    double max_chi = 5*THRESHOLD;
    for(num_lines = 1;  true;  num_lines++) {
        
        delete (appr[0]); delete (appr[1]);
        
        appr[0] = new approximator(num_lines, strd[0], clusters[0], 0);
        appr[1] = new approximator(num_lines, strd[1], clusters[1], 1);
        
        double chi0 = Sqrt(appr[0]->minimizeDisp());
        double chi1 = Sqrt(appr[1]->minimizeDisp());
        max_chi = Max(chi0, chi1);
        
        if (max_chi < THRESHOLD ) break;
    }
      ampid[0] = appr[0]->calculate_amplitudes();
      ampid[1] = appr[1]->calculate_amplitudes(); //отсортированные амплитужы и id
    return max_chi;
}



std::vector<double> reconstructor::getAmplitudes() const {

    std::vector<double> sum(num_lines);
    for(int i=0; i<num_lines; i++) sum[i] = ( ampid[0][i].amp + ampid[1][i].amp ); 

    return sum;

}

double reconstructor::getAmplitudeChisq() const {

    double sum = 0;
    for(int i=0; i<num_lines; i++) {
        double diff = (ampid[0][i].amp - ampid[1][i].amp)/2;
        double mean = (ampid[0][i].amp + ampid[1][i].amp)/2;
        
        sum += diff*diff/mean;
    }
    return sum;
}



struct reconstructor::coordinates reconstructor::getCoordinates() {
    
    std::vector<double> params[2] = { std::vector<double>(num_lines),  std::vector<double>(num_lines) };  //вектора фитированных параметров phi +- tg
    for(int i=0; i<num_lines; i++) {
        params[0][i] = appr[0]->func[i].GetParameter(0);
        params[1][i] = appr[1]->func[i].GetParameter(0);
    }

    
    //так как все формулы выводятся с точностью до 2Pi, нужно проверять возможные ветви и менять на правильные
    bool GOOD = false;
    while(!GOOD) {
        //если тангенс получается больше нужного, мы можем попробовать изменить ветвь(добавить или отнять 2Pi)
        //В данном случае изменение на 4Pi нецелесообразно(хотя для очень длинных детекторов нужно смотреть когда n*2Pi станет больше MAX_TAN)
        for(int j=0; j<num_lines; j++) {
            
            
            //тут ПРОВЕРЯЕТСЯ возможность изменения ветви
            if (Abs( params[0][ampid[0][j].id] - params[1][ampid[1][j].id])/2 > MAX_TAN) {
                //check +2Pi leaf
                if(Abs( 2*Pi() + params[0][ampid[0][j].id] - params[1][ampid[1][j].id])/2 < MAX_TAN) {
                      //std::cout <<"can change_leaf ";
                    continue;
                }
                //check -Pi leaf
                if(Abs(-2*Pi() + params[0][ampid[0][j].id] - params[1][ampid[1][j].id])/2 < MAX_TAN) {
                      //std::cout <<"can change_leaf ";
                    continue;
                }
                //std::cout << "tg too big, next permutation\n";
                goto check_next_permutation; //если всё плохо, может мы перепутали амплитуды?
            }                    
        }                        
        //std::cout << "GOOD\n";
        for(int j=0; j<num_lines; j++) {
            //Если всё хорошо, МЕНЯЕМ ветви там где требуется
             if (Abs( params[0][ampid[0][j].id] - params[1][ampid[1][j].id])/2 > MAX_TAN) {
                //check +2Pi leaf
                if(Abs( 2*Pi() + params[0][ampid[0][j].id] - params[1][ampid[1][j].id])/2 < MAX_TAN) {
                      params[0][ampid[0][j].id] +=  2*Pi(); //не умаляя общности меняем параметр 0
                    continue;
                }
                //check -Pi leaf
                if(Abs(-2*Pi() + params[0][ampid[0][j].id] - params[1][ampid[1][j].id])/2 < MAX_TAN) {
                      params[0][ampid[0][j].id] += -2*Pi(); //не умаляя общности меняем параметр 0
                    continue;
                }
            } 
        }
        GOOD = true;
        break;           
        check_next_permutation:         
        if (std::next_permutation( ampid[1].begin(), ampid[1].end(), [](const approximator::amp_data & f, const approximator::amp_data & s){ return f.amp > s.amp; } )) continue; else break;
    }

    //!!!пока амплитуды на плохость не проверяем
    if (!GOOD) return coordinates(); //empty
    
    //находим координаты
    struct coordinates result = {std::vector<double>(num_lines), std::vector<double>(num_lines)};
    for(int j=0; j<num_lines; j++) {
        double phi = ( params[0][ampid[0][j].id] + params[1][ampid[1][j].id])/2;
        if (phi > 2*Pi()) phi -= 2*Pi();
        if (phi < 0)      phi += 2*Pi();

        result.phi[j] = phi;
        result.tg_theta[j] =  ( params[0][ampid[0][j].id] - params[1][ampid[1][j].id])/2;
    }
    return result;
} 

