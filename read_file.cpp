#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>

#include "ga.h"

using namespace std;
using namespace ga;

namespace ga {
  FindMaxCut::FindMaxCut()
    : iteration(0){  srand(time(0));}
    
 

  void FindMaxCut::readinput(string inFilename){
    // return vertices and edges ?? save in graph sturcture?? 
    std::string line;
    ifstream inFile;
    inFile.open(inFilename);
    std::getline(inFile, line);
    std::istringstream iss(line);
    iss >> num_v_ >> num_e_;
    vertices_.reserve(num_v_);
    edges_.reserve(num_e_);
    
    inFile.close();
    inFile.open(inFilename);
    while (std::getline(inFile, line))
      {
	std::istringstream iss(line);
	int v_1, v_2, w_12;
	if (!(iss >> v_1 >> v_2 >> w_12)) { continue;}
	//cout  << v_1 << v_2 <<  w_12;
	Edge edge (v_1,v_2,w_12);
	edges_.emplace_back(edge);
      }
    // goes from 0 to num_edges-1 -> the indexing 
    inFile.close();
  }
  
  void FindMaxCut::initialize_population(int size){
  
    chromosomes_.clear();
    pop_size = size;
    chromosomes_.reserve(size);
    for (int i = 0; i < size; ++i){
      shared_ptr<Chromosome> chromosome(new Chromosome(num_v_));
      chromosome->index = i;
      
      for (int j = 0; j < num_v_; ++j) {
	chromosome->genes.push_back(rand() %2 );
      }
      chromosomes_.push_back(chromosome);
    }
  }

  void FindMaxCut::compute_fitness(shared_ptr<Chromosome> chromosome){
   float fitness = 0.0f;
   for (int i = 0; i < num_v_; ++i){ 
      if(chromosome->genes[i] == 0) {
	for (int j = 0;j < num_e_; ++j){
	  if (edges_[j].v1 == i+1){
	    if (chromosome->genes[edges_[j].v2 - 1] == 1){
	      fitness += edges_[j].w12;
	    }
	  }
	  else if (edges_[j].v2 == i+1){
	     if (chromosome->genes[edges_[j].v1 - 1] == 1){
	      fitness += edges_[j].w12;
	     }
	  }
	}
      }
   }
   chromosome->fitness = fitness;
  }

  pair<int, int>  FindMaxCut::selection(){ // for now only roulette wheel slection, returns single chromosome
    roulette_fitness.reserve(pop_size);
    int r = 0;
    int selected_index = 0;
    int best = chromosomes_[best_fit_i]->fitness;
    int worst = chromosomes_[worst_fit_i]->fitness;
    int selection_pressure = 3;
    float sum = 0.0f;
    float diff = best - worst;
    //cout << "diff" << diff << endl;
    //cout << "best : " << best;
    //cout << "worst : " << worst;
    for (int i = 0; i < pop_size; ++i){
      if (diff == 0.0f){
    	roulette_fitness.push_back(0.0f);
      }
      else{
    	float f_i = ((worst - chromosomes_[i]->fitness) + diff)/ (selection_pressure - 1 );
    	sum += f_i;
	//cout <<  f_i;
    	roulette_fitness.push_back(f_i);
	
      }
    }
    int index[] = {-1, -1};
    for(int i = 0; i < 2;) {
      int point = rand() % static_cast<int>(sum);
      float sum1 = 0.0f;
      for (size_t j = 0; j < chromosomes_.size(); ++j) {
	sum1 += chromosomes_[j]->fitness;
	if (point < static_cast<int>(sum1)) {
	  index[i] = j;
	  break;
	}
      }
      
      if (index[0] != index[1] || i == 0)
	++i;
    }
    return make_pair(index[0], index[1]);
  }
  
  vector<int>  FindMaxCut::crossover_point(int crossover_mode){
    vector<int>  crossover_point;
    int i,k;
    switch(crossover_mode){
    case 0: // on point crossover
      crossover_point.reserve(1);
      crossover_point.push_back(rand() % (num_v_ - 1) + 1);
    case 1: // uniform crossover but nothing nothing
	crossover_point = uniform_crossover(0.5f);// make vector contain 0 if parent 1 value, make vector = 1 if parent 2 value.
    case 2: // multi point cross over, don't know  how to implement
      k = 5;
      crossover_point = kpoint_crossover(k);
    
	//default:
	//throw "Not Implemented";
    }
    return crossover_point;
  }

  vector<int> FindMaxCut::kpoint_crossover(int k){
    vector<int> which_parent;
    int offset = num_v_/k;
    which_parent.reserve(k);
    int graph_n = 0;
    for( int i = 0; i < num_v_ ; ++i){
      if (i % offset == 0){
	which_parent.push_back(1-graph_n);
      }
      else{
	which_parent.push_back(graph_n);
      }
    }
    return which_parent;
  }
  
  vector<int> FindMaxCut::uniform_crossover(float crossover_prob){
    vector<int> which_parent;
    which_parent.reserve(num_v_);
    float mutation_probability = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    for (int i = 0; i < num_v_; ++i){
      if (mutation_probability >= crossover_prob )
	which_parent.push_back(0);
      else
        which_parent.push_back(1);
    }
    return which_parent;
  }

  pair<shared_ptr<Chromosome>, shared_ptr<Chromosome> > FindMaxCut::crossover(int crossover_mode){
    pair <int, int> indexes = selection();
    int i = indexes.first;
    int j = indexes.second;
    shared_ptr<Chromosome> child1(new Chromosome(num_v_));
    shared_ptr<Chromosome> child2(new Chromosome(num_v_));
    vector<int> c = crossover_point(crossover_mode);
    switch(crossover_mode){
    case 0: // one point crossover mode
      //int c = 2;
      // cout << "parents: " << endl;
      // chromosomes_[i]->display();
      // chromosomes_[j]->display();
      //cout << "crossover_point: " << c << endl;
      for (int t = 0; t< c[0]; ++t){
	child1->genes.push_back(chromosomes_[i]->genes[t]);
	child2->genes.push_back(chromosomes_[j]->genes[t]);
      }
      for (int t = c[0]; t < num_v_; ++t){
	child2->genes.push_back(chromosomes_[i]->genes[t]);
	child1->genes.push_back(chromosomes_[j]->genes[t]);
      }
      // cout << "children: " << endl;
      // child1->display();
      // child2->display(); 
      return make_pair(child1,child2);
    case 1: case 2: // kpoint crossover
      for (int t = 0; t < num_v_; ++t){
	if (c[t] == 0){
	  child1->genes.push_back(chromosomes_[i]->genes[t]);
	  child2->genes.push_back(chromosomes_[j]->genes[t]);
	}
	else if (c[t] == 1){
	  child1->genes.push_back(chromosomes_[j]->genes[t]);
	  child2->genes.push_back(chromosomes_[i]->genes[t]);
	}
      }
      return make_pair(child1,child2);
  }
  }



  pair<shared_ptr<Chromosome>, shared_ptr<Chromosome>> FindMaxCut::mutation(pair<shared_ptr<Chromosome>, shared_ptr<Chromosome>> children, float mut_prob){

    shared_ptr<Chromosome> child1 = children.first;
    shared_ptr<Chromosome> child2 = children.second;
    int r;
    mut_prob = mut_prob * 100;
    int p = (int)(mut_prob);
    for (int i; i < num_v_; ++i)
      {
	//srand(time(NULL));
	r = rand() % 101;
	 
	if (r <= p){
	  child1->genes[i] = 1 - child1->genes[i];
	  child2->genes[i] = 1 - child2->genes[i];
      }
      }
    //child1->display();
    //child2->display();  
    return make_pair(child1,child2);
  }

  // void FindMaxCut::pertubate() { /// ssame as mutation... 

  //   int size = pop_size;
  //   chromosomes_.clear();
  //   chromosomes_.reserve(size);
    
  //   for (int i = 0; i < size; ++i) {
  //     shared_ptr<Chromosome> chromosome(new Chromosome(num_v_));
  //     chromosome->index = i;
      
  //     for (int j = 0; j < num_v_; ++j) {
  // 	chromosome->genes.push_back(rand() % 2);
  //     }
      
  //     chromosomes_.push_back(chromosome);
  //   }

  // }
  
  

  void FindMaxCut::replace(pair<shared_ptr<Chromosome>, shared_ptr<Chromosome> > children){
    
    //children = crossover(); // onepoint chrosover
    shared_ptr<Chromosome> child1 = children.first;
    shared_ptr<Chromosome> child2 = children.second;
    child1->index = worst_fit_i;
    child2->index = second_worst_fit_i;
    chromosomes_[worst_fit_i] = child1;
    chromosomes_[second_worst_fit_i] = child2;
  }

  bool FindMaxCut::isConverge() {
    int num_converged_solution = 0;
    for (int i =0; i < pop_size; ++i)
      if (chromosomes_[i]->fitness == chromosomes_[best_fit_i]->fitness)
	++num_converged_solution;
    return num_converged_solution >= pop_size * 0.2;
  }
  

  void FindMaxCut::solve(){
    int POP_SIZE = 20;
    int num_pertubate = 0;
    
    pop_size = POP_SIZE;
    cout << "POPULATION_SIZE: " << pop_size << endl;
    cout << "CROSSOVER_METHOD:" << 0 << endl;
    cout << "MUTATION PROB" << 0.016 << endl;
    cout << "COVERGED MUTATION PROB" << 0.04 << endl;
    initialize_population(pop_size);
    int best = numeric_limits<int>::min();
    int worst = numeric_limits<int>::max();

    for( int it = 0; it < 1000000; it++){


      compute_fitness(chromosomes_[0]);
      sum_fitness = chromosomes_[0]->fitness;

      for (int i = 0; i < pop_size; ++i ){
	
	compute_fitness(chromosomes_[i]);
	if (chromosomes_[i]->fitness > best){
	  best = chromosomes_[i]->fitness;
	  best_fit_i = i;
	}
	
	if (chromosomes_[i]->fitness < best){
	  worst = chromosomes_[i]->fitness;
	  second_worst_fit_i = worst_fit_i; 
	  worst_fit_i = i;
	} 
	sum_fitness += chromosomes_[i]->fitness;
      }
      pair<shared_ptr<Chromosome>, shared_ptr<Chromosome> > children;
      children = crossover(0);
      
      if (isConverge()) {
	//break;
	
	//	cout << "CONVERGED" << endl;
	//cout << "MUTATING higher rate" << endl;
	children = mutation(children, 0.04f);
	// cout << "At iteration" << it << " ";
	// cout << "PERTUBATING.." << endl;
	// pertubate();
	// ++num_pertubate;
	// cout << "DONE." << endl;
      }
      else{
	children = mutation(children, 0.016f);
      }
      
      replace(children);
      // LOG
      if (it % 10000 == 0){
	
        cout << "ITERATION:" << it << endl;
      	cout << "################################################" << endl;
      	cout<< "fitness best:  " << chromosomes_[best_fit_i]->fitness << endl;
	chromosomes_[worst_fit_i]->display();
	cout<< "fitness worst:  " << chromosomes_[worst_fit_i]->fitness << endl;
      	chromosomes_[best_fit_i]->display();

      	cout<< "fitness0:  " << chromosomes_[0]->fitness << endl;
      	cout<< "fitness1:  " << chromosomes_[1]->fitness << endl;
      	cout<< "fitness2:  " << chromosomes_[2]->fitness << endl;
      	cout<< "fitness3:  " << chromosomes_[3]->fitness << endl;
	cout << "################################################" << endl;
      }
    }
  }
};
  
  


int main()
{
  ga::FindMaxCut find_maxcut;
  
  find_maxcut.readinput("unweighted_50.txt");
  find_maxcut.solve();
  
}