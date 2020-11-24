#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>

int main(int argc, const char *argv[]) {
    if (argc < 4) {
      cout << "./MaxService graph.txt param.txt result.txt MVE SAE algorithm timeLimit" << endl;
        return 0;
    } else {
        auto *graph = new Graph(argv[1], argv[2], argv[3]);
        string mve = "1", sae = "1";

	cout << "Begin the preprocessing" << endl;
	int preprocessingTime;
	auto start = chrono::steady_clock::now();
	if (mve.compare(argv[4]) == 0) graph->MVE(argv[3], argv[8]);
	if (sae.compare(argv[5]) == 0) graph->SAE(argv[3], argv[8]);
	graph->finishPreprocessing(argv[3], mve.compare(argv[4]) == 0, sae.compare(argv[5]) == 0);
	auto end = chrono::steady_clock::now();
	preprocessingTime = chrono::duration_cast<chrono::seconds>(end - start).count();

	string rnfOne = "rnf-one", rnfAll = "rnf-all", mod = "model", lr = "lr";	
	auto *model = new Model(graph);
	if(mod.compare(argv[6]) == 0) {
	  model->initialize();
	  model->initModel();
	  model->solve(argv[7]);
	  model->writeSolution(argv[3], preprocessingTime);
	} else if (lr.compare(argv[6]) == 0) {
	  model->initializeRnf();
	  model->initModel();
	  model->solveLinear(argv[7]);
	  model->writeSolutionLinear(argv[3], preprocessingTime);
	} else {
	  int timeLimit;
	  stringstream tl(argv[7]);
	  tl >> timeLimit;
	  
	  model->initializeRnf();
	  model->initModel();
	  int fo = -1;
	  if (rnfAll.compare(argv[6]) == 0) fo = model->relaxAndFix(timeLimit, true);
	  else fo = model->relaxAndFix(timeLimit, false);
	  model->writeSolutionRnf(argv[3], preprocessingTime, fo);
	}
    }

    return 0;
}
