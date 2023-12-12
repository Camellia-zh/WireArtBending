#include "DataProcessor.h"
#include "WireBending.h"
#include <iostream>
#include <igl/readPLY.h>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <utility>
#include "GCoptimization.h" 
#include "GeneAlg.h"
#include "BeamSearch.h"
#include <chrono>
#include "build/fitting.h"

//std::chrono::duration<double> timeTry = std::chrono::duration<double>::zero();
//std::chrono::duration<double> timeCC1 = std::chrono::duration<double>::zero();
//std::chrono::duration<double> timeCC2 = std::chrono::duration<double>::zero();////croÅö×²¼ì²âÊ±¼ä

//void lisansanwei()
//{
//
//}