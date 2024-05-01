#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
using namespace std;
struct Bee {
    vector<int> route;
    double fitness; //represents the distance travelled by a bee in its journey.
};

class BeeColonyOptimization {
private:
    vector<vector<double>> distances;
    int num_cities; //size of the distances matrix
    int num_bees; // number of bees we have.
    int max_trials; // max trials or the iterations of the bco algorithm.
    Bee best_solution; // as the name it stores the best_solution as a bee.
    vector<Bee> population; //stores each bee in the hive.

public:
    BeeColonyOptimization(const vector<vector<double>>& dist, int n_bees, int trials) 
        : distances(dist), num_cities(dist.size()), num_bees(n_bees),
         max_trials(trials) {}

    void initialize_population() {
        population.resize(num_bees);
        for (int i = 0; i < num_bees; ++i) {
            population[i].route.resize(num_cities);
            for (int j = 0; j < num_cities; ++j) {
                population[i].route[j] = j;
            }
            random_shuffle(population[i].route.begin() + 1, population[i].route.end()); //randomly shuffling the entire row of positions of bees.
        }
    }

    /*
    This Function is calculating the total distance travelled in a tour[i] for some i
    where tour[i] represents the order of travelling of cities.
    */
    double tour_distance(const vector<int>& tour) {
        double dist = 0.0;
        for (int i = 0; i < num_cities; ++i) {
            int from = tour[i];
            int to = tour[(i + 1) % num_cities];
            dist += distances[from][to];
        }
        return dist;
    }

    // This founction
    void evaluate_population() {
        for (int i = 0; i < num_bees; ++i) {
            population[i].fitness = tour_distance(population[i].route);
        }
    }

    void send_employee_bees() {
        for (int i = 0; i < num_bees; ++i) {
            int curr_bee = i;
            int neighbor_idx = rand() % num_bees; //choosing a random neighbouring bee for the current bee.
            while (neighbor_idx == curr_bee) {
                neighbor_idx = rand() % num_bees; // basically making sure that the neighbouring bee is not the current bee.
            }

            Bee new_solution = population[curr_bee]; // takes the current bee.
            int idx1 = rand() % num_cities;
            int idx2 = rand() % num_cities;
            swap(new_solution.route[idx1], new_solution.route[idx2]); //randomly swapping two cities in the path of the current bee.
            new_solution.fitness = tour_distance(new_solution.route); // now calculating the fitness of the new path or the new bee.

            if (new_solution.fitness < population[curr_bee].fitness) {
                population[curr_bee] = new_solution; //now if the distance is less than updating the old one otherwise it'll remain the same.
            }
        }
    }

    void send_onlooker_bees() {
        vector<int> tour_lengths(num_bees);
        double total_fitness = 0.0; // total distance travelled by all the bees.

        //Calculating the tour length of each and every bee.
        for (int i = 0; i < num_bees; ++i) {
            tour_lengths[i] = population[i].fitness;
            total_fitness += tour_lengths[i];
        }

        //Determining an estimate of the how much each bee can travel in comparison to others.
        vector<double> probabilities(num_bees);
        for (int i = 0; i < num_bees; ++i) {
            probabilities[i] = tour_lengths[i] / total_fitness;
        }

        for (int i = 0; i < num_bees; ++i) {
            double r = (double)rand() / RAND_MAX;
            double accum_prob = 0.0;
            int chosen_bee = 0;
            //This ensures that bees who have travelled less will likely continue further.
            for (int j = 0; j < num_bees; ++j) {
                accum_prob += probabilities[j];
                if (r <= accum_prob) {
                    chosen_bee = j;
                    break;
                }
            }

            int neighbor_idx = rand() % num_bees;
            while (neighbor_idx == chosen_bee) {
                neighbor_idx = rand() % num_bees;
            }
            
            //This part is similar to the employeer bees part.
            Bee new_solution = population[chosen_bee];
            int idx1 = rand() % num_cities;
            int idx2 = rand() % num_cities;
            swap(new_solution.route[idx1], new_solution.route[idx2]);
            new_solution.fitness = tour_distance(new_solution.route);

            if (new_solution.fitness < population[chosen_bee].fitness) {
                population[chosen_bee] = new_solution;
            }
        }
    }

    void send_scout_bees() {
        for (int i = 0; i < num_bees; ++i) {
            if (rand() % 100 < 5) { // 5% probability
                Bee new_solution;
                //This part is same as in initialize_population.
                new_solution.route.resize(num_cities);
                for (int j = 0; j < num_cities; ++j) {
                    new_solution.route[j] = j;
                }
                random_shuffle(new_solution.route.begin() + 1, new_solution.route.end());
                new_solution.fitness = tour_distance(new_solution.route);
                population[i] = new_solution;
            }
        }
    }

//Here the main optimization occurs by interating throught the whole process max_trials times.
    Bee optimize() {
        srand(time(0));
        initialize_population();
        evaluate_population();

        best_solution = population[0];
        for (int trial = 0; trial < max_trials; ++trial) {
            if((trial+1)%1000==0){
                cout<<"Processing the "<<trial+1<<"th iteration of bee optimization."<<endl;
            }
            send_employee_bees();
            send_onlooker_bees();
            send_scout_bees();
            evaluate_population();
            //Here we are finalizing the best solution found by the bees in the max_trials.
                for (int i = 0; i < num_bees; ++i) {
                    if (population[i].fitness < best_solution.fitness) {
                        best_solution = population[i];
                    }
                }
        }
        return best_solution;
    }
};
//reading from the .txt file and extracting the distance matrix.
vector<vector<double>> read_distance_matrix(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    int num_cities;
    file >> num_cities;
    vector<vector<double>> distances(num_cities, vector<double>(num_cities));
    for (int i = 0; i < num_cities; ++i) {
        for (int j = 0; j < num_cities; ++j) {
            file >> distances[i][j];
        }
    }

    file.close();
    return distances;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return EXIT_FAILURE;
    }
    
    int numberOfBees=2000;
    int maxTrials=10000;

    vector<vector<double>> distances = read_distance_matrix(argv[1]);
    int num_cities = distances.size();
    BeeColonyOptimization bco(distances, numberOfBees, maxTrials);
    time_t start,end;
    time(&start);
    Bee best_solution = bco.optimize();
    time(&end);

    cout << "Best solution : ";
    int start_city=best_solution.route[0];
    for (int city : best_solution.route) {
        cout << city +1 << " --> ";
    }
    cout<<start_city+1;
    cout << "\nTotal distance : " << best_solution.fitness << endl;
    double time_taken=double(end-start);
    cout<<"Total time taken : "<<time_taken<<" sec"<<endl;
    return 0;
}
