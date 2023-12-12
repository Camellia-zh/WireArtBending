#pragma once
#ifndef BEAMSEARCH_H
#define BEAMSEARCH_H

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <string>
#include <tuple>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include "WireBending.h"
#include "build/fitting.h"

//
////beam search without weight
//class State {
//public:
//    std::vector<double> sequence;
//    std::vector<double> type;
//    double score;
//
//    State(const std::vector<double>& seq, double sc, const std::vector<double>& ty) : sequence(seq), score(sc), type(ty) {}
//
//    bool operator<(const State& other) const {
//        return score > other.score;
//    }
//
//    bool operator==(const State& other) const {
//        return sequence == other.sequence && type == other.type && score == other.score;
//    }
//
//};
//class BeamSearch {
//public:
//    BeamSearch(const std::vector<double>& distances, const std::vector<double>& angles);
//    std::vector<std::vector<double>> runSearch(int sequence_beam_width);
//
//    std::vector<double> getBendingSequence();
//    std::vector<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<int, int>>> getBendingSpace(int angle_start, int angle_end);
//
//
//    int beam_width;
//    int max_length;
//    std::vector<double> distances;
//    std::vector<double> angles;
//    //std::vector<std::pair<int, int>> bendingrange;
//    std::vector<double> bendingsequence;
//    std::vector<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<int, int>>> bendingrange;
//
//
//
//private:
//    //int beamWidth_;
//};


////beam search with weight
//class State {
//public:
//    std::vector<double> sequence;
//    std::vector<double> type;
//    double score;
//    double score2;  // 新添加的 score2 成员变量
//
//    State(const std::vector<double>& seq, double sc, const std::vector<double>& ty, double sc2) : sequence(seq), score(sc), type(ty), score2(sc2) {}
//
//    bool operator<(const State& other) const {
//        if (score == other.score) {
//            return score2 < other.score2;
//        }
//        return score > other.score;
//    }
//
//    bool operator==(const State& other) const {
//        return sequence == other.sequence && type == other.type && score == other.score && score2 == other.score2;
//    }
//};
//class BeamSearch {
//public:
//    BeamSearch(const std::vector<double>& distances, const std::vector<double>& angles, const std::vector<double>& collisionAnglesSeq);
//    std::vector<std::vector<double>> runSearch(int sequence_beam_width);
//
//    std::vector<double> getBendingSequence();
//    std::vector<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<int, int>>> getBendingSpace(int angle_start, int angle_end);
//    int countMatchingOnes(const std::vector<double>& collisionAnglesSeq, const std::vector<double>& type);
//
//    int beam_width;
//    int max_length;
//    std::vector<double> distances;
//    std::vector<double> angles;
//    std::vector<double> collisionAnglesSeq;
//    //std::vector<std::pair<int, int>> bendingrange;
//    std::vector<double> bendingsequence;
//    std::vector<std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<int, int>>> bendingrange;
//    
//
//
//private:
//    //int beamWidth_;
//};


//beam search without weight  freeform
class State {
public:
    std::vector<double> sequence;
    std::vector<double> type;
    double score;

    State(const std::vector<double>& seq, double sc, const std::vector<double>& ty) : sequence(seq), score(sc), type(ty) {}

    bool operator<(const State& other) const {
        return score > other.score;
    }

    bool operator==(const State& other) const {
        return sequence == other.sequence && type == other.type && score == other.score;
    }

};
class BeamSearch {
public:
    BeamSearch(const std::vector<ShapeElement>& ShapeElements, const std::vector<double>& angles);
    std::vector<std::vector<double>> runSearch(int sequence_beam_width);


    int beam_width;
    int max_length;
    std::vector<double> distances;
    std::vector<double> angles;
    std::vector<double> curvature;
    std::vector<bool> segArc;
    std::vector<ShapeElement> ShapeElements;
    //std::vector<std::pair<int, int>> bendingrange;

private:
    //int beamWidth_;
};

#endif 
