#ifndef GA
#define GA


#include <vector>
#include <iostream>
#include <set>
#include <map>

using namespace std;

namespace ga{

  struct Edge {
    Edge (int v1, int v2, int w12) : v1(v1), v2(v2), w12(w12) // how to  initialize struct
    {}
    int getNeighbor(int index) { // define function for Edge struct
      if (index == v1)
	return v2;
      else if (index == v2)
	return v1;
      else
	return -1;
    }
    int v1;
    int v2;
    int w12;
  };

  struct Chromosome{  // index corresponds to vertex 0 == S1, 1 == S2
    Chromosome (int gene_size) // necessary input
    : fitness(0.0f) // chromosome initialization
    { genes.clear(); genes.reserve(gene_size); } 
    void display() {
      for ( auto  &&gene : genes)
	cout << gene;
      cout <<endl;
  }
    float roulette;
    int index;  
    int fitness;
    std::vector<int> genes;
  };

  class FindMaxCut {
  public:
    FindMaxCut();
    ~FindMaxCut() {}

    void solve();//int popu_size, int crossover_method, float mut_prob,float crossover_prob, float converge_mut_prob);
    // public for now, for testing purposes
    void readinput(std::string inFile);
    bool isConverge();
    void readConfig(std::string config_file_name);
    void writeResult(std::string results_file_name);
    void writeResults2(std::string reults_file_namem, int pop_size, int crossover_method, float mut_prob, float converge_mut_prob, float crossover_prob);
    

  private:
    double getElapsedTime();
    void initialize_population(int size);
    void compute_fitness(shared_ptr <Chromosome>  chromosome );
    void update_fitness();
    pair<int, int> selection();
    vector<int> crossover_point(int mode, float prob);
    shared_ptr<Chromosome> crossover(int crossover_mode, float crossover_prob);
    void mutation (float mut_prob, shared_ptr<Chromosome> child );
    void replace(shared_ptr<Chromosome> child);
    vector<int> kpoint_crossover(int k);
    vector<int> uniform_crossover(float crossover_prob);

    // shared_variables
    vector<shared_ptr<Chromosome> > chromosomes_;
    vector<Edge> edges_;
    map<string, float> config_;
    int num_v_;
    int num_e_;
    int pop_size;
    chrono::time_point<chrono::system_clock> start_time_;
    int iteration;
    int num_pertubate;
    int best;
    int worst;
    float average;
    int worst_fit_i;
    int best_fit_i;
    int crossover_method;
    float mut_prob;
    float crossover_prob;
    float converge_mut_prob;

    
  };
    
      

}
#endif
