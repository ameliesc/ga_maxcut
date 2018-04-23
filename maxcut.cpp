#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>

#include "ga2.h"

using namespace std;
using namespace ga;

namespace ga {
  FindMaxCut::FindMaxCut()
    : iteration(0),
   start_time_(chrono::system_clock::now()){
    srand(time(0));
  }

  void FindMaxCut::readinput(string in_file_name) {
        string line;
        ifstream in_file(in_file_name);

        if (in_file.is_open()) {
            getline(in_file, line, ' ');
            num_v_ = stoi(line);
            getline(in_file, line, '\n');
            num_e_ = stoi(line);

            edges_.reserve(num_e_);

            int v1, v2, weight;
            while (getline(in_file, line, ' ')) {
                v1 = stoi(line);
                getline(in_file, line, ' ');
                v2 = stoi(line);
                getline(in_file, line, '\n');
                weight = stoi(line);

                ga::Edge edge(v1, v2, weight);
                edges_.emplace_back(edge);
            }

            in_file.close();
        } else
            cout << "Unable to open file" << endl;
    }

  // void FindMaxCut::readinput(string inFilename){
  //   // return vertices and edges ?? save in graph sturcture?? 
  //   std::string line;
  //   ifstream inFile;
  //   inFile.open(inFilename);
  //   std::getline(inFile, line);
  //   std::istringstream iss(line);
  //   iss >> num_v_ >> num_e_;
  //   edges_.reserve(num_e_);
    
  //   inFile.close();
  //   inFile.open(inFilename);
  //   while (std::getline(inFile, line))
  //     {
  // 	std::istringstream iss(line);
  // 	int v_1, v_2, w_12;
  // 	if (!(iss >> v_1 >> v_2 >> w_12)) { continue;}
  // 	//cout  << v_1 << v_2 <<  w_12;
  // 	Edge edge (v_1,v_2,w_12);
  // 	edges_.emplace_back(edge);
  //     }
  //     // goes from 0 to num_edges-1 -> the indexing 
  //   inFile.close();
  // }
  
  void FindMaxCut::initialize_population(int size){ //CLEAR
  
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

    void FindMaxCut::update_fitness(){
    worst_fit_i = 0;
    best_fit_i = 0;
    average = 0.0f;
    int sum = 0;
    best = numeric_limits<int>::min();
    worst = numeric_limits<int>::max();
    int worst_old = 0;

    for (int i = 0; i < pop_size; ++i ){
      compute_fitness(chromosomes_[i]);
      sum += chromosomes_[i]->fitness;
      if (chromosomes_[i]->fitness > best){
  	best = chromosomes_[i]->fitness;
  	best_fit_i = i;
      }
      
      if (chromosomes_[i]->fitness < worst){
  	worst_old  = worst_fit_i;
  	worst = chromosomes_[i]->fitness;
  	worst_fit_i = i;
  	//if (worst_old != worst_fit_i)
	//second_worst_fit_i = worst_old; 
	
      } 
    }
    average = sum/pop_size;
  }

    pair<int, int>  FindMaxCut::selection(){ // for now only roulette wheel slection, returns single chromosome // CLEAR
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
    for (auto &&chromosome : chromosomes_){
      if (diff == 0.0f){
    	chromosome->roulette = 0.0f;
      }
      else{
    	float f_i = ((worst - chromosome->fitness) + diff)/ (selection_pressure - 1 );
    	sum += f_i;
	chromosome->roulette = f_i;;
	
      }
    }
    int index[] = {-1, -1};
    for(int i = 0; i < 2;) {
      int point = rand() % static_cast<int>(sum);
      //cout << "SUM: " << sum << endl;
      float sum1 = 0.0f;
      for (size_t j = 0; j < chromosomes_.size(); ++j) {
	sum1 += chromosomes_[j]->roulette;
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

  vector<int>  FindMaxCut::crossover_point(int crossover_mode, float crossover_prob){
    vector<int>  crossover_point;
    crossover_point.clear();
    int i,k;
    if (crossover_mode ==  0) { // on point crossover
      crossover_point.reserve(1);
      crossover_point.push_back(rand() % (num_v_ - 1) + 1);
    }
    else if (crossover_mode == 1){ // uniform crossover but nothing nothing
      crossover_point = uniform_crossover(crossover_prob);// make vector contain 0 if parent 1 value, make vector = 1 if parent 2 value.
    }
    else if (crossover_mode == 2){// multi point cross over, don't know  how to implement
      k = 3;
      crossover_point = kpoint_crossover(k);
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

  shared_ptr<Chromosome> FindMaxCut::crossover(int crossover_mode, float crossover_prob){
    pair <int, int> indexes = selection();
    shared_ptr<Chromosome> child(new Chromosome(num_v_));

    int i = indexes.first;
    int j = indexes.second;
    vector<int> c = crossover_point(crossover_mode, crossover_prob);
    if (crossover_mode == 0 ){
      for (int t = 0; t< c[0]; ++t){
	child->genes.push_back(chromosomes_[i]->genes[t]);
      }
      for (int t = c[0]; t < num_v_; ++t){
	child->genes.push_back(chromosomes_[j]->genes[t]);
      }
    }
    else if (crossover_mode == 1 || crossover_mode ==2){ // kpoint crossover
	for (int t = 0; t < num_v_; ++t){
	  if (c[t] == 0){
	    child->genes.push_back(chromosomes_[i]->genes[t]);
	    // children.first->genes.push_back(chromosomes_[i]->genes[t]);
	    //children.second->genes.push_back(chromosomes_[j]->genes[t]);
	  }
	  else if (c[t] == 1){
	    child->genes.push_back(chromosomes_[j]->genes[t]);
	    //children.first->genes.push_back(chromosomes_[j]->genes[t]);
	    //children.second->genes.push_back(chromosomes_[i]->genes[t]);
	  }
	}
      }

    return child;
  }

  void FindMaxCut::mutation(float mut_prob,shared_ptr<Chromosome> child ){
    int r;
    float mutation_prob;
    mutation_prob = mut_prob * 100;
    int p = (int)(mutation_prob);
    for (auto &gene : child->genes)
      {
	r = rand() % 101;
	
	if (r <= p){
	  //cout << "mutationg" << endl;
	  gene = 1 - gene;
	}
      }
  }

  void FindMaxCut::writeResults2(string filename, int pop_size, int crossover_method, float mut_prob, float converge_mut_prob, float crossover_prob){
      ofstream outputFile;
      outputFile.open(filename);
      outputFile << "POP_SIZE" << pop_size << endl;
      outputFile << "CROSSOVER_METHOD:" << crossover_method << endl;
      outputFile << "MUTATION PROB: " << mut_prob << endl;
      outputFile  << "COVERGED MUTATION PROB: " << converge_mut_prob << endl;
      outputFile << "CROSSOVER_PROB: " << crossover_prob << endl;
      outputFile<< "fitness " << best<< endl;
      outputFile << "average " << average << endl;
      outputFile.close();
    }

  void FindMaxCut::writeResult(string result_file_name) {
        string line;
        ofstream result_file(result_file_name);

        if (result_file.is_open()) {
            bool first_flag = true;
            for (size_t i = 0; i < chromosomes_[best_fit_i]->genes.size(); ++i) {
                if (chromosomes_[best_fit_i]->genes[i] == 0) {
                    if (first_flag) {
                        result_file << (i+1);
                        first_flag = !first_flag;
                    } else
                        result_file << " " << (i+1);
                }
            }
            result_file.close();
        } else
            cout << "Unable to open file" << endl;

        if (config_["LOG"]) {
            for (size_t i = 0; i <chromosomes_[best_fit_i]->genes.size(); ++i) {
                if (chromosomes_[best_fit_i]->genes[i] == 0)
                    cout << i << " ";
            }
            cout << endl;
        }
    }

  void FindMaxCut::readConfig(string config_file_name) {
        string line;
        ifstream config_file(config_file_name);

        if (config_file.is_open()) {
            while (getline(config_file, line, ' ')) {
                string config_name = line;
                getline(config_file, line, '\n');
                float config_number = stof(line);
                config_[config_name] = config_number;
            }
            config_file.close();
        } else
            cout << "Unable to open file" << endl;

        if (config_["LOG"]) {
            for (auto itr = config_.begin(); itr != config_.end(); itr++)
                cout << itr->first << " : " << itr->second << endl;
            cout << endl;
        }
    }

  void FindMaxCut::replace(shared_ptr<Chromosome>  child){
    child->index = worst_fit_i;
    compute_fitness(child);
    chromosomes_[worst_fit_i] = child;
    //chromosomes_[second_worst_fit_i] = children.second; 
  }

  double FindMaxCut::getElapsedTime() {
    auto current_time = chrono::system_clock::now();
    chrono::duration<double> elapsed_time = current_time - start_time_;
    return elapsed_time.count();
  }

  bool FindMaxCut::isConverge() {
    int num_converged_solution = 0;
    for (int i =0; i < pop_size; ++i)
      if (chromosomes_[i]->fitness == chromosomes_[best_fit_i]->fitness)
	++num_converged_solution;
    return num_converged_solution >= pop_size * 0.2;
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

  
  void FindMaxCut::solve(){
    double sec;
    int num_pertubate;
    sec =  getElapsedTime();
    pop_size = 20;
    crossover_method = 0;
    mut_prob =  0.01; 
    crossover_prob = 0.5;
    converge_mut_prob = 0.5 ;
    initialize_population(pop_size);
    while (sec < 177.0){
      sec =  getElapsedTime();
      update_fitness();
      shared_ptr<Chromosome> child;
      child = crossover(0,0.5);
      if (isConverge()) {
	mutation(mut_prob,child);
	int j = 0;
	for(auto &chromosome :chromosomes_){
	  if (j % 2 == 0 && chromosome != chromosomes_[best_fit_i])
	    mutation(0.5,chromosome);
	}
      }
      else{
	mutation(mut_prob, child);
      }
      replace(child);
      ++iteration;
      if (sec > 170.0){
	    break; 
      }
    }

  }
};




int main(int argc, char **argv){
  if (argc != 3) {
    std::cerr << "Usage: ./maxcut {data_in} {data_out}" << std::endl;
    return -1;
  }
   ga::FindMaxCut find_maxcut;
   find_maxcut.readinput(argv[1]);
   find_maxcut.solve();
   find_maxcut.writeResult(argv[2]);
   

   return 0;}
