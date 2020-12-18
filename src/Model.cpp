//
// Created by carlos on 06/03/19.
//

#include <chrono>
#include "../headers/Model.h"

Model::Model(Graph *graph) {
  if (graph != nullptr) {
    this->graph = graph;
  } else exit(EXIT_FAILURE);
}

void Model::initialize() {
  int o, d, n = graph->getN(), m = graph->getM();
  try {

    env.set("LogFile", "MS_mip.log");
    env.start();

    f = vector<vector<vector<GRBVar>>>(n, vector<vector<GRBVar>>(n, vector<GRBVar>(n)));
    y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
    z = vector<GRBVar>(n);

    char name[40];
    for (o = 0; o < n; o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	sprintf(name, "y_%d_%d", o, d);
	y[o][d] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
	for (int k: graph->DuS) {
	  sprintf(name, "f_%d_%d_%d", o, d, k);
	  if (!graph->removedF[o][d][k]) this->f[o][d][k] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
	  else this->f[o][d][k] = model.addVar(0.0, 0.0, 0, GRB_BINARY, name);
	}
      }
    }

   
    for (auto i : graph->terminals) {
      sprintf(name, "z_%d", i);
      z[i] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
    }
    
    // ifstream sol;
    // sol.open("solution.sol");
    
    // int i, j;
    // while(!sol.eof()) {
    //     sol >> i >> j;
    //     model.addConstr(y[i][j] == 1);
    // }

    model.update();
    cout << "Create variables" << endl;
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void Model::initializeRnf() {
  int o, d, n = graph->getN(), m = graph->getM();
  try {
    env.set("LogFile", "MS_mip.log");
    env.start();

    f = vector<vector<vector<GRBVar>>>(n, vector<vector<GRBVar>>(n, vector<GRBVar>(n)));
    y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
    z = vector<GRBVar>(n);

    char name[20];
    for (o = 0; o < n; o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	sprintf(name, "y_%d_%d", o, d);
	y[o][d] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
	for (int k: graph->DuS) {
	  sprintf(name, "f_%d_%d_%d", o, d, k);
	  if (!graph->removedF[o][d][k]) this->f[o][d][k] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
	  else this->f[o][d][k] = model.addVar(0.0, 0.0, 0, GRB_CONTINUOUS, name);
	}
      }
    }

    for (auto i : graph->terminals) {
      sprintf(name, "z_%d", i);
      z[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
    }

    model.update();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void Model::initModel() {
  // ifstream file;
  // file.open("solution.sol");
  // string line;
  // vector<string> token;
  // int i , j;
  cout << "Begin the model creation" << endl;
  objectiveFunction();
  rootFlow(), flowConservation(), terminalsFlow();
  relXandY(), maxArcs();
  limDelayAndJitter();
  limVariation();
  primeToTerminals(); // Opcional, mas efetua cortes
  nonTerminalsLeafs();
  cout << "All done!" << endl;
}

int Model::relaxAndFix (int time, bool rnfAll) {
  auto start = chrono::steady_clock::now();
  cout << rnfAll << endl;

  for (auto k : graph->DuS) {
    if (graph->removed[k] || graph->noPath[k]) {
      f[graph->getRoot()][0][k].set(GRB_DoubleAttr_LB, 1.0);
      f[0][k][k].set(GRB_DoubleAttr_LB, 1.0);   
    }
  }
 
  graph->orderingPaths(true);
  solve(to_string(time));
  
  int i, j, k, n = graph->getN();
  double fo = model.get(GRB_DoubleAttr_ObjVal);

  if (rnfAll) {
    for (auto q : graph->terminals)
      z[q].set(GRB_CharAttr_VType, GRB_BINARY);
  }
  model.update();
    
  // for (auto k : graph->terminals)
  //   cout << k << " - " << z[k].get(GRB_CharAttr_VType) << endl;
  // getchar();
  
  for (auto p : graph->pathValues) {
    k = p.first;
    if (graph->noPath[k]) continue;
    for (i = 0; i < n; i++) {
      for (auto arc : graph->arcs[i]) {
	j = arc->getD();
	if (!graph->removedF[i][j][k])
	  f[i][j][k].set(GRB_CharAttr_VType, GRB_BINARY);
      }
    }

    if (!rnfAll) z[k].set(GRB_CharAttr_VType, GRB_BINARY);
    
    model.update();
    solve(to_string(time));

    if (model.get(GRB_IntAttr_Status) == 3) {
      cout << "Infeasible" << endl;
      break;
    }
    
    fo = model.get(GRB_DoubleAttr_ObjVal);
    auto end = chrono::steady_clock::now();
    rnfTime = chrono::duration_cast<chrono::seconds>(end - start).count();
    if (rnfTime >= time) return fo;

    for (i = 0; i < n; i++) {
      for (auto arc : graph->arcs[i]) {
	j = arc->getD();
	if (!graph->removedF[i][j][k]) {
	  if (f[i][j][k].get(GRB_DoubleAttr_X) > 0.9) {
	    f[i][j][k].set(GRB_DoubleAttr_LB, 1.0);
	  } else {
	    f[i][j][k].set(GRB_DoubleAttr_UB, 0.0);
	  }
	}
      }
    }
    model.update();
  }
  solve(to_string(time));
  fo = model.get(GRB_DoubleAttr_ObjVal);
  cout << "Final result: " << fo << endl;
  return fo;
}

void Model::objectiveFunction() {
  GRBLinExpr objective;
  for (auto k : graph->terminals) objective += z[k];
  model.setObjective(objective, GRB_MINIMIZE);
  cout << "Objective Function was added successfully!" << endl;
}

void Model::rootFlow() {
  int o, d, root = graph->getRoot();
  for (auto k : graph->terminals) {
    GRBLinExpr flowExpr, rootExpr;
    for (o = 0; o < graph->getN(); o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	if (o == root) flowExpr += f[root][d][k];
	else if (d == root) rootExpr += f[o][root][k];
      }
    }
    model.addConstr((flowExpr - rootExpr) == 1, "root_flow_all_" + to_string(k));
  }
  model.update();
  cout << "Flow on root node" << endl;
}

void Model::flowConservation() {
  int o, d, root = graph->getRoot();
  for (auto k : graph->DuS) {
    for (int j = 0; j < graph->getN(); j++) {
      if (j != root && j != k) {
	GRBLinExpr flowIn, flowOut;
	for (o = 0; o < graph->getN(); o++) {
	  for (auto *arc : graph->arcs[o]) {
	    d = arc->getD();
	    if (o == j) flowOut += f[j][d][k];
	    if (d == j) flowIn += f[o][j][k];
	  }
	}
	model.addConstr((flowIn - flowOut) == 0, "flow_conservation_" + to_string(j) + "_" + to_string(k));
      }
    }
  }
  model.update();
  cout << "Flow conservation" << endl;
}

void Model::terminalsFlow() {
  int o, d;
  for (auto k : graph->DuS) {
    GRBLinExpr flowIn, flowOut;
    for (o = 0; o < graph->getN(); o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	if (o == k) flowOut += f[k][d][k];
	if (d == k) flowIn += f[o][k][k];
      }
    }
    model.addConstr((flowOut - flowIn) == -1, "flow_on_terminals_" + to_string(k));
  }
  model.update();
  cout << "Flow on terminals" << endl;
}

void Model::relXandY() {
  int o, d;
  for (o = 0; o < graph->getN(); o++) {
    for (auto *arc : graph->arcs[o]) {
      d = arc->getD();
      for (auto k : graph->DuS) {
	model.addConstr(f[o][d][k] <= y[o][d], "f_and_y_relation_" + to_string(o) + "_" + to_string(d) + "_" + to_string(k));
      }
    }
  }
  model.update();
  cout << "f and Y relation" << endl;
}

void Model::maxArcs() {
  GRBLinExpr totalArcs;
  for (int o = 0; o < graph->getN(); o++) {
    for (auto *arc : graph->arcs[o]) {
      totalArcs += y[arc->getO()][arc->getD()];
    }
  }
  model.addConstr(totalArcs == (graph->getN() - 1), "maximum_of_arcs");

  model.update();
  cout << "maximum of arcs in the tree" << endl;
}

void Model::limDelayAndJitter() {
  int o, d, paramDelay, paramJitter;
  for (auto k : graph->terminals) {
    GRBLinExpr limDelay, limJitter;
    for (o = 0; o < graph->getN(); o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	limDelay += arc->getDelay() * f[o][d][k];
	limJitter += arc->getJitter() * f[o][d][k];
      }
    }
    paramDelay = graph->getParamDelay(), paramJitter = graph->getParamJitter();
    model.addConstr(limDelay <= (paramDelay + (graph->getBigMDelay() - paramDelay) * z[k]),
		    "delay_limit_" + to_string(k));
    model.addConstr(limJitter <= (paramJitter + (graph->getBigMJitter() - paramJitter) * z[k]),
		    "jitter_limit_" + to_string(k));
  }
  model.update();
  cout << "Delay and Jitter limits" << endl;
}

void Model::limVariation() {
  int o, d, bigMK, bigML;
  for (auto k : graph->terminals) {
    for (auto l : graph->terminals) {
      if (k != l) {
	GRBLinExpr delayVariation;
	for (o = 0; o < graph->getN(); o++) {
	  for (auto *arc : graph->arcs[o]) {
	    d = arc->getD();
	    delayVariation += arc->getDelay() * (f[o][d][k] - f[o][d][l]);
	  }
	}
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
	model.addConstr(delayVariation <= graph->getParamVariation() + bigMK * z[k] + bigML * z[l],
			"limit_of_variation_between_pairs_" + to_string(k) + "_" + to_string(l));
      }
    }
  }
  model.update();
  cout << "Delay variation limits" << endl;
}

void Model::primeToTerminals() {
   for (auto k : graph->terminals)
    model.addConstr(z[k] >= f[0][k][k], "prime_to_terminals_" + to_string(k));
  model.update();
  cout << "S' to terminals" << endl;
}

void Model::nonTerminalsLeafs() {
  model.addConstr(y[graph->getRoot()][0] == 1);
//   for (auto q : graph->DuS) {
//     for (auto e : graph->DuS) {
//       if (e != q) {
// 	model.addConstr(f[0][q][e] == 0, "non_terminals_leafs_" + to_string(q) + "_" + to_string(e));
//       }
//     }
//   }
//  model.update();
  cout << "Non terminals are leafs" << endl;
}

void Model::solve(string timeLimit) {
  try {
    model.set("TimeLimit", timeLimit);
    // model.computeIIS();
    // model.set("OutputFlag", "0");
    model.update();
    model.write("model.lp");
    model.optimize();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }
}

void Model::solveLinear(string timeLimit) {
  try {
    model.set("TimeLimit", timeLimit);
    //model.set("OutputFlag", "0");
    model.relax();
    model.update();
    model.write("model.lp");
    model.optimize();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }
}

void Model::writeSolution(string instance, int preprocessingTime) {
  try {
    ofstream output;
    output.open(instance, ofstream::app);
    output << "Prep. Time: " << preprocessingTime << endl;
    double ub = model.get(GRB_DoubleAttr_ObjVal), lb = model.get(GRB_DoubleAttr_ObjBound);
    output << "UB: " << ub << endl;
    output << "LB: " << lb << endl;
    if (ub != 0) output << "gap: " << (ub - lb) / ub << endl;
    
        
    output << "N. Nodes: " << model.get(GRB_DoubleAttr_NodeCount) << endl;
    output << "Runtime: " << model.get(GRB_DoubleAttr_Runtime) << endl;

    output << "----- Solution -----" << endl;
    for (auto i : graph->terminals)
      if (z[i].get(GRB_DoubleAttr_X) > 0.5)
	output << i << endl;
    output.close();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }

}

void Model::writeSolutionLinear(string instance, int preprocessingTime) {
  try {
    ofstream output;
    output.open(instance, ofstream::app);
    output << "Prep. Time: " << preprocessingTime << endl;
    double lb = model.get(GRB_DoubleAttr_ObjVal);
    output << "LB: " << lb << endl;
    output << "Runtime: " << model.get(GRB_DoubleAttr_Runtime)  << endl;
    output.close();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }

}

void Model::writeSolutionRnf(string instance, int preprocessingTime, int obj) {
  ofstream output;
  output.open(instance, ofstream::app);
  output << "Prep. Time: " << preprocessingTime << endl;
  output << "UB: " << obj << endl;
  output << "Runtime: " << rnfTime  << endl;
  output.close();
}
