#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <fstream>
#include "parameters.h"

int main() {
    int graph[VERTICES][VERTICES]; 

   srand(time(0));

   for (int n = 0; n < VERTICES; n++) {
        for (int m = 0; m < VERTICES; m++) {
            if (m == n) {
                graph[n][m] = 0;
            } else if (n - m > 0) {
                graph[n][m] = graph[m][n];
            } else if (n - m == -1) {
                graph[n][m] = rand() % 100 + 50;
            } else {
                graph[n][m] = rand() % 100 + 150;
            }
        }
    }

    graph[0][VERTICES - 1] = rand() % 100 + 50;
    graph[VERTICES - 1][0] = graph[0][VERTICES - 1];

    std::ofstream graph_file("..\\data\\graph.txt");

    if (!graph_file) {
        std::cout << "Error creating graph text file!";
    }

    for (int n = 0; n < VERTICES; n++) {
        for (int m = 0; m < VERTICES; m++) {
            graph_file << graph[n][m] << ";";
        }
        graph_file << "\n";
    }

    graph_file.close();
}
