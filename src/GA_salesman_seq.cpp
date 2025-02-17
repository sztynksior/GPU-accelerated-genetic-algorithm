#include <fstream>
#include <iostream>
#include <list>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <random>
#include <algorithm>
#include <math.h>
#include <unordered_map>
#include <sstream>
#include "parameters.h"

struct Solution {
    int chromosome[CHROMOSOME_LENGTH];
    int cost;

    static bool compare_solutions_descending(const Solution &solution_a, const Solution &solutiion_b) {
        return solution_a.cost > solutiion_b.cost;
    }

    void print_solution() {
        std::cout << "Cost:" << this->cost << " Path: " << STARTING_VERTICE;
        for (int vertice : this->chromosome) std::cout << "-" << vertice;
        std::cout << "-" << STARTING_VERTICE << std::endl;
    }
};

// for simplicity we assume that optimal path goes through vertices in order
int calculate_optimal_cost(int graph[VERTICES][VERTICES]) {
    int optimal_cost = 0;
    for (int i = 0; i < VERTICES; i++) {
        optimal_cost += graph[i][(i + 1) % VERTICES];
    }
    return optimal_cost;
}

void get_graph_from_file(std::string file_path, int graph[VERTICES][VERTICES]) {
    std::ifstream graph_file(file_path);
    std::string line;
    std::string path_length;
    if (graph_file.is_open()) {
        for (int n = 0; std::getline(graph_file, line); n++) {
            std::stringstream line_stream(line);
            for (int m = 0; std::getline(line_stream, path_length, ';'); m++) {
                graph[n][m] = std::stoi(path_length);
            }
        }
    } else {
        std::cout << "Problem opening graph file!";
    }
    graph_file.close();
}

void init_population(std::vector<Solution> &population) {
    // we assume that salesman have to start and end at STARTING_VERTICE, 
    // so places that he can visit in different order are vertices in {0,...,VERTICES - 1} / STARTING_VERTICE
    std::vector<int> places_to_visit(CHROMOSOME_LENGTH);
    int skip = 0;
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        if (i == STARTING_VERTICE) { // when i == STARTING_VERTICE skip it
            skip = 1;
        }
        places_to_visit[i] = i + skip;   
    }
    // each solution is a permutation of vertices in {0,...,VERTICES - 1} / STARTING_VERTICE
    std::random_device rd;
    std::mt19937 generator(rd());
    for (int i = 0; i < POPUPLATION_SIZE; i++) {
        std::shuffle(places_to_visit.begin(), places_to_visit.end(), generator);
        for (int j = 0; j < CHROMOSOME_LENGTH ; j++) {
            population[i].chromosome[j] = places_to_visit[j];
        }
    }
}

void calculate_cost(std::vector<Solution> &population, const int graph[VERTICES][VERTICES]) {
    // cost of each solution is total path length, lower is better
    for (int i = 0; i < POPUPLATION_SIZE; i++) {
        // salesman have to start and end at vartice STARTING_VERTICE
        population[i].cost = graph[STARTING_VERTICE][population[i].chromosome[0]];
        population[i].cost += graph[population[i].chromosome[CHROMOSOME_LENGTH - 1]][STARTING_VERTICE];
        for (int j = 0; j < CHROMOSOME_LENGTH - 1; j++) {
            int from = population[i].chromosome[j];
            int to = population[i].chromosome[j + 1];
            population[i].cost += graph[from][to];
        }
    }
}

void rank_selection(std::vector<Solution> &population, Solution const* parents[POPUPLATION_SIZE]) {
    int ranks_sum = (1 + POPUPLATION_SIZE) * POPUPLATION_SIZE / 2 ;
    for (int i = 0; i < POPUPLATION_SIZE; i++) {
        // genreatete random value from range 0 - 1 which represents choice of a parent
        double solution_choice = ((double)rand()) / (RAND_MAX + 1);
        // calculate rank of a parent given cumulative probability formula, eg: p_0 in [0; 1/21), p_1 in [1/21; 3/21), ... , p_6 in [15/21; 1)
        int parent_rank = std::floor((-1 + std::sqrt(8 * solution_choice * ranks_sum + 1)) / 2);
        parents[i] = &population[parent_rank];
    }
}

void perform_crossover(Solution const &parent_1, Solution const &parent_2, Solution &ofspring) {
    int num_of_genes_from_par1 = CHROMOSOME_LENGTH / 2;
    int num_of_genes_from_par2 = CHROMOSOME_LENGTH - num_of_genes_from_par1;
    int first_parent_substring_offset = num_of_genes_from_par2 / 2;
    std::unordered_map<int, int> copied_alleles; // store already visited veritices, we do not want any duplicates
    
    // copy substring of first parent into ofspring, starting at position 'first_parent_sunbstring_offset' with length 'num_of_genes_from_par1'
    for (int j = 0; j < num_of_genes_from_par1; j++) {
        int allele_to_copy = parent_1.chromosome[first_parent_substring_offset + j];
        ofspring.chromosome[first_parent_substring_offset + j] = allele_to_copy;
        copied_alleles.insert(std::make_pair(allele_to_copy, allele_to_copy));
    }
    
    // copy remaning alleles from second parent preserving their order and not allowing for any duplicates
    int ofspring_copy_index = 0;
    for (int j = 0; j < CHROMOSOME_LENGTH; j++) {
        int allele_to_copy = parent_2.chromosome[j];
        if (copied_alleles.find(allele_to_copy) == copied_alleles.end()) {
            ofspring.chromosome[ofspring_copy_index] = allele_to_copy;
            ofspring_copy_index++;
            if (ofspring_copy_index == first_parent_substring_offset) {
                ofspring_copy_index += num_of_genes_from_par1;
            }
        }
    }
}

void linear_order_crossover(Solution const* parents[POPUPLATION_SIZE], std::vector<Solution> &new_population) {
    for (int i = 0; i < POPUPLATION_SIZE; i += 2) {
        perform_crossover(*(parents[i]), *(parents[i + 1]), new_population[i]);
        perform_crossover(*(parents[i + 1]), *(parents[i]), new_population[i + 1]);
    }

    if (POPUPLATION_SIZE % 2 == 1) {
        perform_crossover(*(parents[POPUPLATION_SIZE - 1]), *(parents[0]), new_population[POPUPLATION_SIZE - 1]);
    }
}
struct Triple {
    int pos1;
    int pos2;
    int sol;
    Triple(int pos1, int pos2, int sol) {
        this->pos1 = pos1;
        this->pos2 = pos2;
        this->sol = sol;
    }
};

void inversion_mutation(std::vector<Solution> &new_population, float mutation_probability) {
    std::vector<Triple> mutated;
    int swap_pos_1 = -1;
    int swap_pos_2 = -1;
    for (int i = 0; i < POPUPLATION_SIZE; i++) {
        float random = ((float)rand()) / RAND_MAX;
        if (random < mutation_probability) {
            while (swap_pos_1 == swap_pos_2) {
                swap_pos_1 = rand() % CHROMOSOME_LENGTH;
                swap_pos_2 = rand() % CHROMOSOME_LENGTH;
            }
            int buff = new_population[i].chromosome[swap_pos_1];
            new_population[i].chromosome[swap_pos_1] = new_population[i].chromosome[swap_pos_2];
            new_population[i].chromosome[swap_pos_2] = buff;
            Triple trip(swap_pos_1, swap_pos_2, i);
            mutated.push_back(trip);
            swap_pos_1 = -1;
            swap_pos_2 = -1;
        }
    }
    int x;
}

int main() {
    srand(time(0));
    int graph[VERTICES][VERTICES];
    get_graph_from_file("..\\data\\graph.txt", graph);
    int optimal_cost = calculate_optimal_cost(graph);
    int max_deviation = optimal_cost * (1 + DEVIATION);

    // use pointers to allow shallow copying
    std::vector<Solution> *population = new std::vector<Solution>(POPUPLATION_SIZE);
    std::vector<Solution> *new_population = new std::vector<Solution>(POPUPLATION_SIZE);
    Solution const* parents[POPUPLATION_SIZE];

    init_population(*population);
    calculate_cost(*population, graph);

    // sort population in ascending order to have the worst solution at the beginning (rank = 1) and the best solution at the end (rank = POPULATION_SIZE)
    std::sort((*population).begin(), (*population).end(), Solution::compare_solutions_descending);
    for (int i = 0; i < NUMBER_OF_POPULATIONS && (*population)[POPUPLATION_SIZE - 1].cost > max_deviation; i++) {
        rank_selection(*population, parents);
        linear_order_crossover(parents, *new_population);
        inversion_mutation(*new_population, MUTATION_PROBABILITY);
        std::swap(population, new_population); // use pointers to avoid unecessary memory operations (shallow copy)
        calculate_cost(*population, graph);
        std::sort((*population).begin(), (*population).end(), Solution::compare_solutions_descending);

        if (i % POPULATION_INFO == 0) {
            std::cout << "Population: " << i << " ";
            (*population)[POPUPLATION_SIZE - 1].print_solution();
        }
    }
    
    std::cout << "Optimal cost: " << optimal_cost << " Max deviation: " << max_deviation << " Best solution:" << std::endl;
    (*population)[POPUPLATION_SIZE - 1].print_solution();

    delete population;
    delete new_population;
}