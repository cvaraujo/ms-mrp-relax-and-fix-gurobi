//
// Created by carlos on 06/03/19.
//

#ifndef MRP_MODEL_H
#define MRP_MODEL_H

#include <iostream>
#include <vector>
#include "string"
#include <iomanip>
#include <bits/ios_base.h>
#include <algorithm>
#include <fstream>
#include <gurobi_c++.h>
#include "Graph.h"
#include <boost/algorithm/string.hpp>

using namespace std;

class Model {
  Graph *graph;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  vector<vector<vector<GRBVar>>> f;
  vector<vector<GRBVar>> y;
  vector<GRBVar> z;
  int rnfTime = 0;

  void preprocessing();

  void objectiveFunction();

  void rootFlow();

  void flowConservation();

  void terminalsFlow();

  void relXandY();

  void onlyWay();

  void maxArcs();

  void limDelayAndJitter();

  void limVariation();

  void primeToTerminals();

  void nonTerminalsLeafs();

public:
  Model(Graph *graph);

  void initialize();

  void initializeRnf();

  void initModel();
    
  int relaxAndFix(int timeLimit, bool rnfAll);

  void solve(string timeLimit);

  void solveLinear(string timeLimit);

  void writeSolution(string instance, int preprocessingTime);

  void writeSolutionLinear(string instance, int preprocessingTime);

  void writeSolutionRnf(string instance, int preprocessingTime, int obj);

};


#endif //MRP_MODEL_H
