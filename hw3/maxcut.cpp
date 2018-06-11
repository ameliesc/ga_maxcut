#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <queue>
#include <stack>
#include <set>
#include <stdlib.h>
#include "ga.h"

using namespace std;
using namespace ga;

// #include <fstream>
// #include <algorithm>
// #include <ctime>
// #include <queue>
// #include <stack>
// #include <set>
// #include <algorithm>

// #include <stdlib.h>

// #include "ga.h"

// using namespace std;
// using namespace ga;



namespace ga {
  void FindMaxCut::readinput(string in_file_name) {
        string line;
        ifstream in_file(in_file_name);

        if (in_file.is_open()) {
            getline(in_file, line, ' ');
            num_v_ = stoi(line);
            getline(in_file, line, '\n');
            num_e_ = stoi(line);

            edges_.reserve(num_e_);
	    vertices_.reserve(num_v_);

	    for (int i = 0; i < num_v_; ++i){ // for somme reason size needs to be 101 otherwise segmentationfault
	      Vertex vertex(i);
	      vertices_.push_back(vertex);
	    }
	    
            int v1, v2, weight;
            while (getline(in_file, line, ' ')) {
	      v1 = stoi(line) - 1;
	      getline(in_file, line, ' ');
	      v2 = stoi(line) - 1;
	      getline(in_file, line, '\n');
	      weight = stoi(line);

              ga::Edge edge(v1, v2, weight);
              edges_.emplace_back(edge);
	      vertices_[v1].append(v2, weight);
	      vertices_[v2].append(v1, weight);
	      
		
            }

            in_file.close();
        } else
            cout << "Unable to open file" << endl;
    }

  void FindMaxCut::reordering(int reordering_method) {
    switch (reordering_method) {
    case 0:
      BFSReordering();
      return;
    case 1:
      DFSReordering();
      return;
    default:
      throw "Not implemented";
    }
  }

  void FindMaxCut::solveLocal() {
    random_.clear();
    random_.reserve(num_v_);
    for (int i = 0; i < num_v_; ++i){
      if (sec >= 490.0){
	
	break; 
      }
    random_.push_back(i);
    }
    initialize_population(pop_size);
    localOptimizeAll();
    
    cout << "Total iter : " << iteration << endl;
  }

  void FindMaxCut::localOptimizeAll() {
    for (int a = 0; a < pop_size; ++a) {
      random_shuffle(random_.begin(), random_.end());
      bool improved = true;
      
      while (improved) {
	improved = false;
	for (int i = 0; i < num_v_; ++i) {
	  int sum_weights = 0;
	  int max_j = vertices_[random_[i]].neighbours.size();
	  for (int j = 0; j < max_j; ++j) {
	    sec =  getElapsedTime();
	    int index_neigh = vertices_[random_[i]].neighbours[j].first;
	    int weight_neigh = vertices_[random_[i]].neighbours[j].second;
	    
	    // check whether in same group or not
	    if (chromosomes_[a]->genes[random_[i]] == chromosomes_[a]->genes[index_neigh]) {
	      sum_weights += weight_neigh;
	    } else {
	      sum_weights -= weight_neigh;
	    }

	    
	  }
	  
	  if (sum_weights > 0) {
	    chromosomes_[a]->genes[random_[i]] = 1 - chromosomes_[a]->genes[random_[i]];
	    improved = true;

	  }
	  //if (sec >= 490.0){
	  // goto endofloop;
	      
	  //}
	}
      }
      //endofloop:
      //compute_fitness(chromosomes_[a]);
      
    }
  }
  
  void FindMaxCut::BFSReordering() {
    // initialize variables for reordering
    for (int i = 0; i < num_v_; ++i) {
      origin_to_reorder_.push_back(-1);
      reorder_to_origin_.push_back(-1);
    }
    
    // initialize
    for (int i = 0; i < num_v_; ++i)
      vertices_[i].seen = false;
    
    queue<int> q;
    
    int root = rand() % num_v_;
    vertices_[root].seen = true;
    q.push(root);
    
    // do BFS
    int reorder_counter = 0;
    while (!q.empty()) {
      int current = q.front();
      q.pop();
      
      origin_to_reorder_[current] = reorder_counter;
      reorder_to_origin_[reorder_counter] = current;
      
      for (pair<int, int> neighbor : vertices_[current].neighbours) {
	int n = neighbor.first;
	
	if (!vertices_[n].seen) {
	  vertices_[n].seen = true;
	  q.push(n);
	}
      }
      
      ++reorder_counter;
    }
    
    
    // reordering
    vector<Vertex> vertices_reord;
    vector<Edge> edges_reord;
    
    vertices_reord.reserve(num_v_);
    edges_reord.reserve(num_e_);
    
    for (int i = 0; i < num_v_; ++i) {
      Vertex vertex(i);
      vertices_reord.push_back(vertex);
    }
    
    for (Edge edge : edges_) {
      int r_v1 = origin_to_reorder_[edge.v1];
      int r_v2 = origin_to_reorder_[edge.v2];
      int w12 = edge.w12;
      
      // update edge
      Edge e(r_v1, r_v2, w12);
      edges_reord.push_back(e);
      
      // update vertex
      vertices_reord[r_v1].append(r_v2, w12);
      vertices_reord[r_v2].append(r_v1, w12);
    }
    
    vertices_ = vertices_reord;
    edges_ = edges_reord;
  }
  
  void FindMaxCut::DFSReordering() {
    // initialize variables for reordering
    for (int i = 0; i < num_v_; ++i) {
      origin_to_reorder_.push_back(-1);
      reorder_to_origin_.push_back(-1);
    }
    
    // initialize
    for (int i = 0; i < num_v_; ++i)
      vertices_[i].seen = false;
    
    stack<int> s;
    
    int root = rand() % num_v_;
    s.push(root);
    
    // do DFS
    int reorder_counter = 0;
    while (!s.empty()) {
      int current = s.top();
      s.pop();
      
      
      if (!vertices_[current].seen) {
	origin_to_reorder_[current] = reorder_counter;
	reorder_to_origin_[reorder_counter] = current;
	
	vertices_[current].seen = true;
	
	for (pair<int, int> neighbour : vertices_[current].neighbours) {
	  int n = neighbour.first;
	  s.push(n);
	}
	++reorder_counter;
      }
    }
    // reordering
    vector<Vertex> vertices_reord;
    vector<Edge> edges_reord;
    
    vertices_reord.reserve(num_v_);
    edges_reord.reserve(num_e_);
    
    for (int i = 0; i < num_v_; ++i) {
      Vertex vertex(i);
      vertices_reord.push_back(vertex);
    }
    
    for (Edge edge : edges_) {
      int r_v1 = origin_to_reorder_[edge.v1];
      int r_v2 = origin_to_reorder_[edge.v2];
      int w12 = edge.w12;
      
      // update edge
      Edge e(r_v1, r_v2, w12);
      edges_reord.push_back(e);
      
      // update vertex
      vertices_reord[r_v1].append(r_v2, w12);
      vertices_reord[r_v2].append(r_v1, w12);
    }
    
    vertices_ = vertices_reord;
    edges_ = edges_reord;
  }

  void FindMaxCut::localsearch(shared_ptr<Chromosome> chromosome) {
    //      int score_difference = 0;

        if (iteration  % 1000 == 0)
            random_shuffle(random_.begin(), random_.end());

        bool improved = true;

        while (improved) {
	  improved = false;
            for (int i = 0; i < num_v_; ++i) {
                int score_difference = 0;
                size_t max_j = vertices_[random_[i]].neighbours.size();
                for (size_t j = 0; j < max_j; ++j) {
                    int neighbor_index = vertices_[random_[i]].neighbours[j].first;
                    int weight = vertices_[random_[i]].neighbours[j].second;
                    if (chromosome->genes[random_[i]] == chromosome->genes[neighbor_index]) { // f_x(v)

                        score_difference += weight;
                    } else {
                        score_difference -= weight;
                    }
                }
                if (score_difference > 0) {
                    chromosome->genes[random_[i]] = 1 - chromosome->genes[random_[i]];
                    improved = true;
                }
            }

	    //  for (int i = 0; i < num_v_; ++i) {
            //     int score_difference2 = 0;
            //     size_t max_j = vertices_[random_[i]].neighbours.size();
            //     for (size_t j = 0; j < max_j; ++j) {
            //         int neighbor_index = vertices_[random_[i]].neighbours[j].first;
            //         int weight = vertices_[random_[i]].neighbours[j].second;
            //         if (chromosome->genes[random_[i]] == chromosome->genes[neighbor_index]) { // f_x(v)

            //             score_difference2 += weight;
            //         } else {
            //             score_difference2 -= weight;
            //         }
            //     }
            //     if (score_difference2 > 0) {
            //         chromosome->genes[random_[i]] = 1 - chromosome->genes[random_[i]];
            //         improved = true;
            //     }
            // }
	    //cout << score_difference << "\n";

        }
	
	//compute_fitness(chromosome);
    }
  
  // void  FindMaxCut::localsearch(shared_ptr<Chromosome> chromosome){
  //    // need to change this to random index generator somehow
  //   //double start;
  //   //double current;
  //   int sum_weights = 0;
  //   if (iteration  % 1000 == 0)
  //     random_shuffle(random_.begin(), random_.end());
  //   // start = getElapsedTime();
  //   bool improved = true;
    
  //   while(improved){
  //     improved = false;
  //     for(int i = 0; i < num_v_; ++i){
  // 	sum_weights = 0; 
  // 	size_t num_neigh = vertices_[random_[i]].neighbours.size();
  // 	  for(size_t j = 0; j < num_neigh; ++j){
  // 	    int index_neigh = vertices_[random_[i]].neighbours[j].first;
  // 	    int weight_neigh = vertices_[random_[i]].neighbours[j].second;
	   
  // 	    if (chromosome->genes[random_[i]] == chromosome->genes[index_neigh]){
  // 	      sum_weights += weight_neigh;
  // 	      }
  // 	    else {	      
  // 	      sum_weights -=weight_neigh;
  // 	    }
  // 	  }
  // 	    if (sum_weights > 0) {
  // 	      //chromosome2  = chromosome;
  // 	      chromosome->genes[random_[i]] = 1 - chromosome->genes[random_[i]];
  // 	      improved = true;
  // 	    }

  // 	    //if (i == 0 && sum_weights <= 0){
  // 	    //  improved = false;
  // 	    // goto endofloop;
  // 	    //}


  // 	    //cout << "sum_weights:" << sum_weights << "\n";
  //     }
  //     cout << "sum_weights:" << sum_weights << "\n";
  //     cout << "imporved" << improved << "\n";
  //     //cout << "improve:" << improved << "\n";
  //     //cout << sum_weights > 0 << "\n";

  //   }
  // //endofloop:
  // //int a;
  //   //compute_fitness(chromosome);
  // }

  // void FindMaxCut::loalsearchwrapper(){
  //   for (int i = 0; i < pop_size; ++i){
  //     chromosomes_[i] = localsearch(chromosomes_[i]);
  //   }
  // }

      
  void FindMaxCut::initialize_population(int size){
  
    chromosomes_.clear();
    pop_size = size;
    chromosomes_.reserve(size);
    for (int i = 0; i < size; ++i){
      shared_ptr<Chromosome> chromosome(new Chromosome(num_v_));
      chromosome->index = i;
      
      for (int j = 0; j < num_v_; ++j) {
  	chromosome->genes.push_back(rand_int_(gen_));
      }
      chromosomes_.push_back(chromosome);
    }
  }



  shared_ptr<Chromosome> FindMaxCut::sharing(double minDist, double share_degree, shared_ptr<Chromosome> chromosome){
    double denominator = 1;
    double dist;
    for(int j = 0; j < pop_size; ++j){
      dist = hamming_dist(chromosome, chromosomes_[j]);
	if (dist < minDist){
	  denominator += (1 - (dist/share_degree));
	}
    }
    chromosome->fitness = chromosome->fitness/denominator;
    return chromosome;
  }

  int FindMaxCut::hamming_dist(shared_ptr<Chromosome> c1, shared_ptr<Chromosome> c2){
    int hamming_dist = 0;
    for(int i = 0; i < num_v_; ++i){
      if(c1->genes[i] != c2->genes[i]){
	hamming_dist +=1;
      }
    }
    return hamming_dist;
  }
  void FindMaxCut::compute_fitness(shared_ptr<Chromosome> chromosome){
    //float fitness = 0.0f;ss
   chromosome->fitness = 0;
   for (int i = 0; i < num_e_; ++i)
     if (chromosome->genes[edges_[i].v1] != chromosome->genes[edges_[i].v2])
       chromosome->fitness += edges_[i].w12;

   
   // for (int i = 0; i < num_v_; ++i){ 
   //    if(chromosome->genes[i] == 0) {
   // 	for (int j = 0;j < num_e_; ++j){
   // 	  if (edges_[j].v1 == i+1){
   // 	    if (chromosome->genes[edges_[j].v2 - 1] == 1){
   // 	      fitness += edges_[j].w12;
   // 	    }
   // 	  }
   // 	  else if (edges_[j].v2 == i+1){
   // 	     if (chromosome->genes[edges_[j].v1 - 1] == 1){
   // 	      fitness += edges_[j].w12;
   // 	     }
   // 	  }
   // 	}
   //    }
   // }
   // chromosome->fitness = fitness;
  }
  
  void FindMaxCut::update_fitness(){
    int worst_fit_i = 0;
    int best_fit_i = 0;
    average = 0.0f;
    int sum = 0;
    int best = 0; numeric_limits<int>::min();
    int worst =0; numeric_limits<int>::max();
    //int worst_old = 0;



      for (int i = 0; i < pop_size; ++i) {
	//	compute_fitness(chromosomes_[i]);
	best_fit_i = best > chromosomes_[i]->fitness ? best_fit_i : chromosomes_[i]->index;
	worst_fit_i = worst < chromosomes_[i]->fitness ? worst_fit_i : chromosomes_[i]->index;
        worst = worst < chromosomes_[i]->fitness ? worst : chromosomes_[i]->fitness;
	best = best  > chromosomes_[i]->fitness ? best : chromosomes_[i]->fitness;
      }
    // for (int i = 0; i < pop_size; ++i ){
    //   compute_fitness(chromosomes_[i]);
    //   //chromosomes_[i] = sharing(10.0,3000.00,chromosomes_[i]);
    //   sum += chromosomes_[i]->fitness;
    //   if (chromosomes_[i]->fitness > best){
    // 	best = chromosomes_[i]->fitness;
    // 	best_fit_i = i;
    //   }
      
    //   if (chromosomes_[i]->fitness < worst){
    // 	//worst_old  = worst_fit_i;
    // 	worst = chromosomes_[i]->fitness;
    // 	worst_fit_i = i;
    // 	//if (worst_old != worst_fit_i)
    // 	//second_worst_fit_i = worst_old; 
	
    //   } 
    // }
    average = sum/pop_size;
    worst_index_ = worst_fit_i;
    min_fitness_ = worst;
    max_fitness_ = best;

    shared_ptr<Chromosome> best_solution(new Chromosome(num_v_));
    best_solution->fitness = chromosomes_[best_fit_i]->fitness;
    best_solution->genes = chromosomes_[best_fit_i]->genes;
    best_solution_ = best_solution;

    // if (chromosomes_[best_fit_i]->fitness > best_solution_->fitness) {
    // best_solution_->fitness = chromosomes_[best_fit_i]->fitness;
    // best_solution_->genes = chromosomes_[best_fit_i]->genes;
    //}
  }

  //   pair<int, int>  FindMaxCut::selection(){ // for now only roulette wheel slection, returns single chromosome // CLEAR
  //     //int r = 0;
  //     //int best = chromosomes_[best_fit_i]->fitness;
  //     //int worst = chromosomes_[worst_fit_i]->fitness;
  //   int selection_pressure = 3;
  //   double sum = 0.0;
  //   //float diff = max_fitness_ - min_fitness_;
  //   //cout << "diff" << diff << endl;
  //   //cout << "best : " << best;
  //   //cout << "worst : " << worst;
  //   for (int i = 0; i < pop_size; ++i){
  //     double f_i =  (chromosomes_[i]->fitness - min_fitness_) + (max_fitness_ - min_fitness_ + 0.001f) / (selection_pressure - 1);
  // 	if (chromosomes_[i]->fitness < 1.0f)
  // 	  chromosomes_[i]->fitness = 1.0f;
  //   	sum += f_i;
  // 	chromosomes_[i]->roulette = f_i;
	
  //     }
  //   int index[] = {-1, -1};
  //   for(int i = 0; i < 2;) {
  //     uniform_real_distribution<> rand_real(0,sum);
  //     //double point =
  // 	//int point = rand() % static_cast<int>(sum);
  //     //cout << "SUM: " << sum << endl;
  //     //cout << point << endl;
  //     double sum1 = 0.0;
  //     for (size_t j = 0; j < num_v_; ++j) {
  //     	sum1 += chromosomes_[j]->roulette;
  //     	if  ( rand_real(gen_)<  sum1) {
  // 	  index[i] = j;
  // 	  break;
  // 	}
  //     }
      
  //     if (index[0] != index[1] || i == 0)
  //     	++i;
  //   }
  //   return make_pair(index[0], index[1]);
  //   //return make_pair(1,2);
  // }

      pair<int, int > FindMaxCut::selection() {

            double sum_of_roulette = 0.0;
            for (int i = 0; i < pop_size; ++i) {
	      chromosomes_[i]->roulette = (chromosomes_[i]->fitness - min_fitness_) + (max_fitness_ - min_fitness_ + 0.001f) / (3 - 1);
                if (chromosomes_[i]->roulette < 1.0f)
                    chromosomes_[i]->roulette = 1.0f;
                sum_of_roulette += chromosomes_[i]->roulette;
            }

            // selection first parent
            int index[] = {pop_size-1 , pop_size-1};
            int point;
            double sum;

            point = rand() % static_cast<int>(sum_of_roulette);
            sum = 0.0f;
            for (size_t i = 0; i < chromosomes_.size(); ++i) {
                sum += chromosomes_[i]->roulette;
                if (point <= static_cast<int>(sum)) {
                    index[0] = i;
                    break;
                }
            }

            // select second parent
            do {
                point = rand() % static_cast<int>(sum_of_roulette);
                sum = 0.0f;
                for (size_t i = 0; i < chromosomes_.size(); ++i) {
                    sum += chromosomes_[i]->roulette;
                    if (point <= static_cast<int>(sum)) {
                        index[1] = i;
                        break;
                    }
                }
            } while (index[0] == index[1]);

	    return make_pair(index[0], index[1]);
            //return make_pair(solutions_[index[0]], solutions_[index[1]]);
      //else
  //  throw "Not implemented";
    }

  vector<int>  FindMaxCut::crossover_point(int crossover_mode, float crossover_prob){
    vector<int>  crossover_point;
    crossover_point.clear();
    int k;
    if (crossover_mode ==  0) { // one point crossover
      crossover_point.reserve(1);
      crossover_point.push_back(rand() % (num_v_ - 1) + 1);
    }
    else if (crossover_mode == 1){ // uniform crossover 
      crossover_point = uniform_crossover(crossover_prob);// make vector contain 0 if parent 1 value, make vector = 1 if parent 2 value.
    }
    else if (crossover_mode == 2){
      k = 6;
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
    //cout << "parent 1 : " << i << endl;
    int j = indexes.second;
    //cout << "parent 2 : " << j << endl;
    vector<int> c = crossover_point(crossover_mode, crossover_prob);
    if (crossover_mode == 0 ){
      for (int t = 0; t< c[0]; ++t){
	child->genes.push_back(chromosomes_[i]->genes[t]);
      }
      for (int t = c[0]; t < num_v_; ++t){
	child->genes.push_back(chromosomes_[j]->genes[t]);
      }
    }
    else if (crossover_mode ==2){ // kpoint crossover
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
    compute_fitness(child);
    return child;
  }

  void FindMaxCut::mutation(float mut_prob,shared_ptr<Chromosome> child ){
    //float mutation_prob;
    //mutation_prob = mut_prob * 100;
    //int p = (int)(mutation_prob);
    for (int i = 0; i < num_v_; ++i)
      {	
	if (rand_real_(gen_) <= mut_prob){
	  //cout << "mutationg" << endl;
	  child->genes[i] = 1 - child->genes[i];
	  //gene = 1 - gene;
	}
      }
  }

  void FindMaxCut::writeResults2(string filename, int pop_size, int crossover_method, float mut_prob, float converge_mut_prob, float crossover_prob, string local_op ){
    ofstream outputFile(filename, ios::app);

    //outputFile.open(filename, ios::app);
    if (!outputFile) {
      cerr << "can't open output file" << endl;
    }
    std::cout << pop_size << endl;
    std::cout << max_fitness_ << endl;
    std::cout << average << endl;
    //outputFile.open(filename);
    outputFile << "POP_SIZE" << pop_size << endl;
    outputFile << "CROSSOVER_METHOD:" << crossover_method << endl;
    outputFile << "MUTATION PROB: " << mut_prob << endl;
    outputFile  << "COVERGED MUTATION PROB: " << converge_mut_prob << endl;
    outputFile << "CROSSOVER_PROB: " << crossover_prob << endl;
    outputFile << "LOCAL OPTIMIZATION: " << local_op << endl;
    outputFile<< "fitness " << max_fitness_<< endl;
    outputFile << "average " << average << endl;
    //outputFile.close();
    }

  void FindMaxCut::writeResult(string result_file_name) {
        string line;
        ofstream result_file(result_file_name);

        if (result_file.is_open()) {
            bool first_flag = true;
            for (size_t i = 0; i < best_solution_->genes.size(); ++i) {
                if (best_solution_->genes[i] == 0) {
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
    }


  void FindMaxCut::replace(shared_ptr<Chromosome>  child){
    child->index = worst_index_;
    //compute_fitness(child);
    chromosomes_[worst_index_] = child;
    //update_fitness();
    
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
      if (chromosomes_[i]->fitness == max_fitness_)
	++num_converged_solution;
    return num_converged_solution >= pop_size * 0.5;
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

  
  void FindMaxCut::solve(string filestring){
    random_.clear();
    random_.reserve(num_v_);
    for (int i = 0; i < num_v_; ++i)
      random_.push_back(i);
    string filename = filestring + "_no_reorder_localsearch_crossover2_pop100.txt";
    sec =  getElapsedTime();
    pop_size = 50;
    crossover_method = 2; //
    mut_prob =  0.015;
    reordering_method = 0;
    crossover_prob = 0.5;
    string local_search = "true";
    //reordering(reordering_method);
    initialize_population(pop_size);
    int num_restart = 0;
    while (sec <= 490.0){   
      sec =  getElapsedTime();
      update_fitness();
      shared_ptr<Chromosome> child;
     
      if (isConverge() && num_restart == 5 ) {
      	for(auto &chromosome :chromosomes_){
      	  if (chromosome != best_solution_){
      	  mutation(0.10,chromosome);
      	  }
      	}
      	update_fitness();
      }
      else if(isConverge()){
      	for(int i = 0; i < pop_size; ++i){
      	  if(chromosomes_[i] != best_solution_)
      	    {
      	      for(int j = 0; j < num_v_; ++j){
      		chromosomes_[i]->genes[j] = rand_int_(gen_);
      	      }
      	    }
      	  update_fitness();
      	}
      	num_restart += 1; 
      }

      child = crossover(crossover_method,0.5);
      mutation(mut_prob, child);
      //cout << "fitness before: "<<  child->fitness << "\n";
      localsearch(child);
      //cout << "best solution after: "<< best_solution_->fitness << "\n";
      replace(child);
      ++iteration;

      if (sec >= 490.0){
      	break; 
      }
    }
    cout << "writing results:";
    writeResults2( filename,  pop_size, crossover_method,  mut_prob,  converge_mut_prob, crossover_prob, local_search);
  }
};




int main(int argc, char **argv){
  if (argc != 3) {
    std::cerr << "Usage: ./maxcut {data_in} {data_out}" << std::endl;
    return -1;
  }
   ga::FindMaxCut find_maxcut;
   find_maxcut.readinput(argv[1]);
   find_maxcut.solve(argv[1]);
   
   find_maxcut.writeResult(argv[2]);
   

   return 0;}
