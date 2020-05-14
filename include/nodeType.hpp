#pragma once
enum class NodesType : unsigned int { OWN, HAlO, ALL };

/*
 *  Struct represents info about any nodes vector:
 * - nodesType - type of nodes to operate with
 * - nodesNumber - number of nodes to operate with
 */

struct NodesInfo {

    NodesInfo(){};

    NodesInfo(int number, NodesType type) : nodesNumber(number), nodesType(type){};

    int nodesNumber;
    NodesType nodesType;
};
