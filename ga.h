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
  
  struct Vertex {
    Vertex (int v, int subgraph) : v(v), subgraph(subgraph) {}
    int v;
    int subgraph;
  };

  struct Chromosome{  // index corresponds to vertex 0 == S1, 1 == S2
    Chromosome (int gene_size) // necessary input
    : score(0), fitness(0.0f) // chromosome initialization
    { genes.clear(); genes.reserve(gene_size); } // ??????
    void display() {
      for ( auto  &&gene : genes)
	cout << gene;
      cout <<endl;
  }
    void display(set<int> cutting_points){
      auto itr = cutting_points.begin();
      for (size_t i = 0; i < genes.size(); ++i) {
	if (i == static_cast<size_t>(*itr)){
	  cout << " | ";
	  ++itr;
      }
	cout << genes[i];
      }
      cout << endl;
    }
    int score;
    int index;  
    float fitness;
    std::vector<int> genes;
};

 class FindMaxCut {
 public:
   FindMaxCut();
   ~FindMaxCut() {}
 // public for now, for testing purposes
   
   void readinput(std::string inFile);
   void solve(int popu_size, int crossover_method, float mut_prob,float crossover_prob, float converge_mut_prob);
   // solve
   // write Result

   // Log related methods??
   // helper functions
   void report_best(){
     shared_ptr<Chromosome> chrom;
     chrom = chromosomes_[best_fit_i];
     //cout<< "fitness:  " << chrom->fitness << endl;
     //cout << "chromosome";
     //chrom->display(); 
     //cout << "cutting_points" << chrom. 
   }
 private:
   // list all the GA functions used  - GA _ solver etc
   void initialize_population(int size);
   void compute_fitness(shared_ptr<Chromosome> chromosome);
   pair<int, int> selection();
   void replace(pair<shared_ptr<Chromosome>, shared_ptr<Chromosome> > children);
   pair<shared_ptr<Chromosome>, shared_ptr<Chromosome> > crossover(int crossover_mode, float crossover_prob);
   vector<int> uniform_crossover(float crossover_prob);
   vector<int> kpoint_crossover(int k);

   vector<int> crossover_point(int crossover_mode, float crossover_prob);
   pair<shared_ptr<Chromosome>, shared_ptr<Chromosome>> mutation(pair<shared_ptr<Chromosome>, shared_ptr<Chromosome>> children, float mut_prob);
   bool isConverge();
   void pertubate();
   // mutate
   //replace
   // fitness

   int num_v_;
   int iteration;
   int num_e_;
   int pop_size;
   int sum_fitness;
   int worst_fit_i;
   int second_worst_fit_i;
   int best_fit_i;
   int num_k_crossover;
   int crossover_mode;
   vector<float> roulette_fitness;
   int selection_pressure;
   vector<Edge> edges_;
   vector<Vertex> vertices_;
   vector<shared_ptr<Chromosome> > chromosomes_;
   // score keeping variables iteration, etc 
 };
  
}
#endif
