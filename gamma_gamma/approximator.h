
#include "TH1D.h"
#include "TMath.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"

#include "TVirtualFitter.h"
#include "Math/PdfFunc.h"



#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>




/////   MAGIC  NUMBERS ?? //

#define LAYERS              7
#define BINS_PER_LAYER    312

#define TRIG_LEVEL   5.
#define SECOND_LEVEL 1.


#define SIGMA (0.035)


const double rho[LAYERS] = { 37.92, 39.96, 42.00, 44.04, 46.08, 48.12, 50.16 };  //радиусы цилиндров в сантиметрах
const double rho_begin = 37.0, rho_end = 51.0;

/////

namespace gera_nm {

struct tower_data {
    std::vector<int> ch;
    std::vector<double> amp;
    std::vector<double> rawamp;
};

struct tower_cluster {
    std::vector< std::vector<int>> towers;
    std::vector< int > index;
    std::vector< int > status;
    std::vector< int > type;
};

struct strip_data {
    std::vector<int> packedID;
    std::vector<int> innerID;
    std::vector<int> layer;

    std::vector<double> rho;
    std::vector<double> start_phi;

    // //std::vector<int> direction;
    std::vector<int> cluster_id;
    std::vector<double> amp;
    std::vector<double> sigma;
};

struct cross_data{
    std::vector<int> id1;
    std::vector<int> id2;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

}
using namespace gera_nm;
using namespace TMath;



////-- CLUSTER FINDER///  Нахождение кластеров и всё что с этим связано

//структура с информацией по элементам кластера
struct hist_bin {
    int layer; //номер уровня от 0 до 6 
    int phi;   //номер бина по углу
    double amp; //амплитуда
    double start_phi;  //для дальнейшего фитирования номера бина недостаточно - нужен стартовый угол
};

//алгоритм по поиску кластеров
std::vector<std::vector<hist_bin>> cluster_finder(const strip_data *strd);


//--------------------------------------------

//Аппроксимация прямыми

//структура, в которой будет храниться вся требуемая для фитирования информация
struct approximator {
public:
    const int num_lines;                               //количество прямых
    const strip_data *strd;                            //исходные данные из конвертора
    const std::vector<std::vector<hist_bin>> clusters; //фитирование будет проводиться по элементам из кластеров
    TF1 *func = NULL;
    const int direction;  //направление спиралей полосок
    TVirtualFitter *Minuit = NULL;
    double minD = -1;//посчитанный суммарный минимальный момент инерции

    approximator(int num_lines, const strip_data *strd, const std::vector<std::vector<hist_bin>> &clusters, int direction) 
        : num_lines(num_lines), strd(strd), clusters(clusters), direction(direction)
    {
        //инициализация функций
        func = new TF1[num_lines];
        /*
            alpha = phi0 + arcsin(x0/r) + z0/r + tg_theta*Sqrt( 1 - (x0*x0)/(r*r) )
            [1] - x0, [2] - +-z0, 
            в случае по умолчанию [0] - phi0 +- tg_theta, [3] = 0
            есть возможность задать [0] = phi0, [3] = tg_theta
        */
        for(int i=0; i<num_lines; i++) {
            
            func[i] = TF1("", "[0] + ASin([1]/x) + [2]/x + [3]*Sqrt( 1 - ([1]*[1])/(x*x) )", rho_begin, rho_end);
            func[i].SetParameters(2*Pi()/num_lines *i, 0., 0., 0.); //default parameters
        }

        Minuit = TVirtualFitter::Fitter(NULL, num_lines);
        static double arglist_for_quiet_mode[] = {-1}; 
        Minuit->ExecuteCommand("SET PRINTOUT", arglist_for_quiet_mode, 1);
        Minuit->ExecuteCommand("SET NOWARNINGS", NULL, 0);
    }

    approximator(int num_lines, const strip_data *strd, int direction) : 
        approximator(num_lines, strd, cluster_finder(strd), direction) {}



    virtual ~approximator() { delete[] func; }

    double getWeight() const;  //считает и возвращает амплитуду всех кластеров

    double minimizeDisp();  //выполнить подгонку, возвращает вероятность что данные описываются данным количеством прямых
    //double getProbability() const;   


    struct amp_data {double amp; int id;};
    std::vector<amp_data> calculate_amplitudes() const;   //рассчитать амплитуды сигналов по найденным прямым

    //ДАЛЕЕ ИДЁТ КОСТЫЛЬ. 
    //для фитирования используется статический метод fcn, который связывается с фитируемым объектом через статический указатель app
    static void fcn(Int_t &npar, double *gin, double &f, double *par, Int_t iflag);
    static approximator *app; //указатель на текущий фитируемый объект. нужен для связи fcn и экземпляром.
    
};

approximator* approximator::app = NULL;


int visualizate(const char dir[], const strip_data *strd);  //построить и вывести гистограмму
TH2D *get_hist(const char dir[], const strip_data *strd);   //построить и вернуть гистограмму
int visualisate_amp(const strip_data *strd[2]); //построить проекцию амплитуды на ось phi  


//.....................................................










std::vector<std::vector<hist_bin>> cluster_finder(const strip_data*strd) {

    std::vector<hist_bin> queue; 
    queue.reserve(2*LAYERS * BINS_PER_LAYER);   //очередь для поиска в ширину
    
    int status[LAYERS][BINS_PER_LAYER];         //статус каждой вершины
    struct {
        double amp;
        double start_phi;
    } info[LAYERS][BINS_PER_LAYER];              //информация о каждой вершине

    //инициализация info
    for(int l=0; l<LAYERS; l++) {
        for(int p=0; p<BINS_PER_LAYER; p++) {
            info[l][p] = {0, 2*Pi()/BINS_PER_LAYER *(p+0.5)}; //amp = 0, start_phi пусть примерно равен центру "бина"
            status[l][p] = -1;  //nо clusters
        }
    }
    for(unsigned i=0; i<strd->amp.size(); i++) {
        int l = strd->layer[i];
        int p = strd->start_phi[i]*BINS_PER_LAYER/(2*Pi());
        //если в один бин попадёт несколько точек, их start_phi будет определяться с учётом веса их амплитуд.
        info[l][p].start_phi = info[l][p].start_phi*info[l][p].amp + strd->start_phi[i]*strd->amp[i];
        info[l][p].amp       = info[l][p].amp + strd->amp[i];
        info[l][p].start_phi = info[l][p].start_phi / info[l][p].amp;
        

        queue.push_back({l, p, strd->amp[i], strd->start_phi[i]}); //инициализируем очередь для дальнейшего алгоритма
    }



    int total_clusters = -1;                                  //
    std::vector<std::vector<hist_bin>> clusters;              //выходной массив с кандидатами на кластеры
    
    //приготовили очередь для дальнейшего алгоритма
    //используем bfs? обход в ширину

    while( ! queue.empty() ) {
        hist_bin tag = queue.back(); 
        queue.pop_back();

        //смотрим зажигание текущей вершины
        if ( status[tag.layer][tag.phi] == -1 && tag.amp >= TRIG_LEVEL) {
            //у этой вершины нет зажённых соседей
            //создаём нового кандидата в кластер
            total_clusters++;
            status[tag.layer][tag.phi] = total_clusters;
            
            clusters.push_back({tag});       
        }

        if( status[tag.layer][tag.phi] != -1)
        {
            //check neighbours

            for(int new_layer = tag.layer - 2; new_layer <= tag.layer + 2; new_layer++)
            {
                if (new_layer < 0 || new_layer >= LAYERS ) continue; //check borders
                for(int new_phi = tag.phi - 3; new_phi <= tag.phi + 3; new_phi++) 
                {
                    /*check itself*/ if(new_layer == tag.layer && new_phi == tag.phi) continue;

                    int new_phi_absolute = (new_phi + BINS_PER_LAYER) % BINS_PER_LAYER;
                    hist_bin new_tag = {new_layer, new_phi_absolute, info[new_layer][new_phi_absolute].amp, info[new_layer][new_phi_absolute].start_phi};

                    if (status[new_layer][new_phi_absolute]  == -1  &&  new_tag.amp >=  SECOND_LEVEL )
                    {      
                        status[new_layer][new_phi_absolute] = status[tag.layer][tag.phi];
                        clusters[   status[tag.layer][tag.phi]    ].push_back(new_tag);

                        queue.push_back(new_tag);
                    }
                }
            }
        }
    }
    std::sort(clusters.begin(), clusters.end(), [](const std::vector<hist_bin> & a, const std::vector<hist_bin> & b){ return a.size() > b.size(); });
    

    //КОЛИЧЕСТВО ЭЛЕМЕНТОВ В КЛАСТЕРЕ
    //!!!!!!!!!!!!!!!
        while( !clusters.empty()  && clusters.back().size() < 4) clusters.pop_back();


    return clusters;
}




/////////////////////////////////////

double approximator::getWeight() const //возвращает суммарную амплитуду кластеров
{
    double w = 0;
    for(int i=0; i<clusters.size(); i++) {
        for(int j=0; j<clusters[i].size(); j++) w += clusters[i][j].amp;
    }
    return w;
}



double approximator::minimizeDisp() {
    

    approximator::app = this; //свяжем fcn и текущий объект
    

    Minuit->SetFCN( approximator::fcn );

    for(int i=0; i<num_lines; i++) {            //инициализируем параметры minuit
        auto c = std::string("c").append( std::to_string(i) );  //имена для параметров

        double default_phi = clusters[i % clusters.size()][i / clusters.size()].start_phi;  //стартовый угол берём из кластеров.
        //сначала берутся нулевые элементы, с которых начались кластеры, потом первые, вторые и т. д.   
        func[i].SetParameter(0, default_phi);
        
        Minuit->SetParameter( i, c.c_str(),  default_phi,  Pi(),  -3*Pi(),  3*Pi());
    }

    ///EXECUTE FITTING ALGORITHM
    Minuit->ExecuteCommand("SCAN", NULL, 0);
    Minuit->ExecuteCommand("MINIMIZE IMPROVE", NULL, 0);
    for(int i = 0; i < num_lines; i++) {
        double par = Minuit->GetParameter(i);
        while (par > 2*Pi()) par -= 2*Pi();
        while (par <  0) par += 2*Pi();
        func[i].SetParameter(0, par);
    }
    
    //КАЧЕСТВО ПОДГОНКИ
    
        double amin, edm, errdef;
        int nvpar, nparx;
        Minuit->GetStats(amin, edm, errdef, nvpar, nparx);
           
    approximator::app = NULL;
    return ( minD = ( amin/getWeight() ) *(SIGMA*SIGMA) );
}





//функция, вызываемая Minuit
void approximator::fcn(Int_t &npar, double *gin, double &f, double *par, Int_t iflag) {

    const int num_lines = approximator::app->num_lines;
    for(int i=0; i<num_lines; i++) {
        approximator::app->func[i].SetParameter(0, par[i]);
    }

    double chisq = 0;
    //проходим по всем элементам из кластеров
    for( int i=0; i< approximator::app->clusters.size();   i++ ) {
        for(int j=0; j<approximator::app->clusters[i].size(); j++) {

            double phi = approximator::app->clusters[i][j].start_phi;
            double amp = approximator::app->clusters[i][j].amp;
            int layer = approximator::app->clusters[i][j].layer;

            //находим минимальное расстояние от точки до одной из прямых

            double min_dist = 2*Pi(); //big number
            for(int k=0; k<num_lines; k++) {
                double dist = Abs(phi - approximator::app->func[k].Eval(rho[layer]));
                    while(dist > 2*Pi()) dist -= 2*Pi();
                    if (dist > Pi()) dist = 2*Pi() - dist;

                if (dist < min_dist) min_dist = dist;
            }
 
            chisq +=  (min_dist* min_dist)/(SIGMA*SIGMA)  * amp;
            //SIGMA в данном случае - магическое число
        }
    }
    f = chisq;
}





std::vector<approximator::amp_data> approximator::calculate_amplitudes() const   //рассчитать амплитуды сигналов по найденным прямым
{
    std::vector<approximator::amp_data> amp(num_lines); 
    
    for(int i=0; i<num_lines; i++) amp[i] = { 0, i };
    
    
    for(int i=0; i<clusters.size(); i++) {
        for(int j=0; j<clusters[i].size(); j++) {
            //для каждого сигнала смотрим, к какой прямой он пренадлежит - в такую и записываем амплитуду
            int id = -1;
            double min_dist = 2*Pi(); //big number
            for(int k=0; k<num_lines; k++) {
                double dist = Abs(clusters[i][j].start_phi - func[k].Eval( rho[clusters[i][j].layer] ));
                    while(dist > 2*Pi()) dist -= 2*Pi();
                    if (dist > Pi()) dist = 2*Pi() - dist;

                if (dist < min_dist) { min_dist = dist; id = k; }
            }

            amp[id].amp += clusters[i][j].amp;
        }
    }
    
    std::sort(amp.begin(), amp.end(), [](amp_data f, amp_data s){ return f.amp > s.amp; });

    return amp;
}



///////////////



//создаём двухмерную гистограмму для полосок в зависимости от направления
TH2D *get_hist(const char dir[], const strip_data *strd) {

    TH2D *hist = new TH2D( dir, std::string(dir).append(";R;\\phi;a").c_str(), LAYERS, rho_begin, rho_end,  BINS_PER_LAYER, 0, 2*Pi() );
    
    for(unsigned i=0; i<strd->amp.size(); i++) {
        //debug  // std::cout << strd->layer[i] << " | " << strd->start_phi[i] << " & " << strd->amp[i] << std::endl;
        hist->Fill(rho[ strd->layer[i] ], strd->start_phi[i], strd->amp[i]);
    }
    return hist;
}


//рисуем гистограмму
int visualizate(const char dir[], const strip_data *strd) {
    
    auto hist = get_hist(dir, strd);
    hist->Draw("colz");
    return 0;
}



int visualisate_amp(const strip_data *strd[2]) //построить проекцию амплитуды на ось phi  
{
    TCanvas *canv = new TCanvas();
        canv->Divide(2, 1);
        canv->cd(1);
        {
            TH1D *amp_phi = get_hist("project0", strd[0])->ProjectionY("dir0");

            amp_phi->Draw("HIST");
        }
        canv->cd(2);
        {
            TH1D *amp_phi = get_hist("project1", strd[1])->ProjectionY("dir1");

            amp_phi->Draw("HIST");
        }
    return 0;
}















