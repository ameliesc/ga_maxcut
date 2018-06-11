#ifndef GA
#define GA

#include <memory>
#include <vector>
#include <iostream>
#include <set>
#include <map>
#include <memory>
#include <ctime>
#include <chrono>
#include <random>
#include <queue>
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

  struct Vertex {
    Vertex (int ind) :ind(ind)
    {}
    void append(int id, int w){
      neighbours.push_back(make_pair(id,w));
    }
    int ind;
    vector<pair <int, int> > neighbours;
    bool seen;
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
    explicit FindMaxCut() :
    iteration(0),
    rd_(),
      gen_(rd_()),
      rand_real_(0,1),
      rand_int_(0,1),
      start_time_(chrono::system_clock::now())
	{
      srand(time(0));
    }
    ~FindMaxCut() { }

    void solve(std::string filename);//int popu_size, int crossover_method, float mut_prob,float crossover_prob, float converge_mut_prob);
    // public for now, for testing purposes
    void readinput(std::string inFile);
    bool isConverge();
    void writeResult(std::string results_file_name);
    void writeResults2(std::string results_file_name, int pop_size, int crossover_method, float mut_prob, float converge_mut_prob, float crossover_prob, std::string local_op);
    

  private:
    void reordering(int reordering_method);
    void BFSReordering();
    void DFSReordering();
    void localOptimizeAll();
    void localsearch(shared_ptr<Chromosome> chromosome);
    double getElapsedTime();
    void initialize_population(int size);
    void compute_fitness(shared_ptr<Chromosome> chromosome );
    void update_fitness();
    shared_ptr<Chromosome> sharing(double minDist, double share_degree,  shared_ptr<Chromosome> chromosome);
    int hamming_dist(shared_ptr<Chromosome> c1, shared_ptr<Chromosome> c2);
    void solveLocal();
    void selection();
    vector<int> crossover_point(int mode, float prob);
    shared_ptr<Chromosome> crossover(int crossover_mode, float crossover_prob);
    shared_ptr<Chromosome> best_solution_;
    void mutation (float mut_prob, shared_ptr<Chromosome> child );
    void replace(shared_ptr<Chromosome> child);
    vector<int> kpoint_crossover(int k);
    vector<int> uniform_crossover(float crossover_prob);


    // shared_variables
    vector<shared_ptr<Chromosome> > chromosomes_;
    vector<Edge> edges_;
    vector<Vertex> vertices_;
    vector<int> used_vertices_;
    pair<int, int> parents_index;
    //queue<int> used_vq_;
    // reordering related
    vector<int> origin_to_reorder_;
    vector<int> reorder_to_origin_;
    
    map<string, float> config_;
    int num_v_;
    int num_e_;
    int pop_size;
    chrono::time_point<chrono::system_clock> start_time_;
    int iteration;
    int num_pertubate;
    int reordering_method;
    float average;
    int worst_index_;
    int min_fitness_;
    int max_fitness_;
    //int worst_fit_i;
    //int best_fit_i;
    int crossover_method;
    vector<int> random_;
    float mut_prob;
    float crossover_prob;
    float converge_mut_prob;
    double sec;

    //random relatied 
    random_device rd_;
    mt19937 gen_;
    uniform_real_distribution<> rand_real_;
    uniform_int_distribution<> rand_int_;
  };
    
      

}
#endif
