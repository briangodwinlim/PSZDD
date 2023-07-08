#include <set>
#include <map>
#include <cmath>
#include <iomanip>

// TdZdd
#include <tdzdd/DdEval.hpp>
#include <tdzdd/DdSpecOp.hpp>
#include <tdzdd/DdStructure.hpp>
#include <tdzdd/util/Graph.hpp>
#include <tdzdd/util/IntSubset.hpp>
#include <tdzdd/spec/DegreeConstraint.hpp>
#include <tdzdd/spec/FrontierBasedSearch.hpp>

using namespace tdzdd;

const int INF = INT_MAX / 2;    // Infinity placeholder

// Stochastic graph class with edge failure probabilities and source vertices
class StochGraph: public Graph {
private:
    std::set<std::string> vertices, sources;
    std::map<std::pair<std::string,std::string>,double> probabilities, distances;

public:
    void addEdge(std::string vertexName1, std::string vertexName2, double prob = 1.0, double dist = 1.0) {
        Graph::addEdge(vertexName1, vertexName2);
        
        assert(0.0 <= prob && prob <= 1.0);
        probabilities[std::make_pair(vertexName1, vertexName2)] = prob;
        probabilities[std::make_pair(vertexName2, vertexName1)] = prob;

        assert(dist > 0.0);
        distances[std::make_pair(vertexName1, vertexName2)] = dist;
        distances[std::make_pair(vertexName2, vertexName1)] = dist;
    }

    double getProb(EdgeNumber e) const {
        return probabilities.at(edgeName(e));
    }

    double getDist(EdgeNumber e) const {
        return distances.at(edgeName(e));
    }

    void addVertex(std::string vertex) {
        vertices.insert(vertex);
    }

    std::set<std::string> getVertices() const {
        return vertices;
    }

    void addSource(std::string source) {
        sources.insert(source);
    }

    std::set<std::string> getSources() const {
        return sources;
    }
};

// DdEval for finding the shortest distance
class MinDist: public DdEval<MinDist,double> {
private:
    int const numEdges;
    StochGraph const& graph;

public:
    MinDist(StochGraph const& graph) 
        : graph(graph), numEdges(graph.edgeSize()) {
    }

    void evalTerminal(double &v, int id) const {
        v = id ? 0.0 : INF;
    }

    void evalNode(double &v, int level, DdValues<double,2> const& values) const {
        v = std::min(values.get(0), values.get(1) + graph.getDist(numEdges - level));
    }
};

// DdEval for calculating the per-vertex path survival probabilities
class PathSurvival: public DdEval<PathSurvival,double> {
private:
    int const numEdges;
    StochGraph const& graph;

public:
    PathSurvival(StochGraph const& graph) 
        : graph(graph), numEdges(graph.edgeSize()) {
    }

    void evalTerminal(double &v, int id) const {
        v = id ? 1.0 : 0.0;
    }

    void evalNode(double &v, int level, DdValues<double,2> const& values) const {
        v = values.get(0) + graph.getProb(numEdges - level) * values.get(1);
    }
};


std::map<std::string,bool> opt;
std::string options[][2] = { //
        {"graph", "Dump input graph to STDOUT in DOT format"}, //
        {"adjusted", "Output adjusted per-vertex path survival probabilities"}, //
    };

void usage(char const* cmd) {
    std::cerr << "usage: " << cmd
              << " <edgelist> <probabilities> <sources> [ <options>... ]\n";
    std::cerr << "options\n";
    for (unsigned i = 0; i < sizeof(options) / sizeof(options[0]); ++i) {
        std::cerr << "  -" << options[i][0];
        for (unsigned j = options[i][0].length(); j < 10; ++j) {
            std::cerr << " ";
        }
        std::cerr << ": " << options[i][1] << "\n";
    }
}

// Global variables
StochGraph graph;                                                       // Stochastic graph object
double density;                                                         // Graph density
int numEdges, numVertices;                                              // Number of edges and vertices
std::string edgelistFilename, probabilitiesFilename, sourcesFilename;   // Input filenames

int main(int argc, char *argv[]) {
    for (unsigned i = 0; i < sizeof(options) / sizeof(options[0]); ++i) {
        opt[options[i][0]] = false;
    }

    // Process inputs
    try {
        for (int i = 1; i < argc; ++i) {
            std::string s = argv[i];
            if (s[0] == '-') {
                s = s.substr(1);
                if (opt.count(s)) {
                    opt[s] = true;
                }
                else {
                    throw std::exception();
                }
            }
            else if (edgelistFilename.empty()) {
                edgelistFilename = s;
            }
            else if (probabilitiesFilename.empty()) {
                probabilitiesFilename = s;
            }
            else if (sourcesFilename.empty()) {
                sourcesFilename = s;
            }
            else {
                throw std::exception();
            }
        }
    }
    catch (std::exception& e) {
        usage(argv[0]);
        return 1;
    }

    MessageHandler::showMessages();
    MessageHandler mh;
    mh.begin("Started\n");

    // Read inputs
    try {
        std::ifstream edgelistFile, probabilitiesFile, sourcesFile;
        edgelistFile.open(edgelistFilename);
        probabilitiesFile.open(probabilitiesFilename);
        sourcesFile.open(sourcesFilename);
        std::string edge, probability, source;
        
        // Get edges and probabilities
        getline(probabilitiesFile, probability);
        std::istringstream probabilitySS(probability);
        while (getline(edgelistFile, edge)) {
            std::istringstream edgeSS(edge);
            std::string u, v, p;
            edgeSS >> u >> v;
            probabilitySS >> p;

            graph.addVertex(u);
            graph.addVertex(v);
            graph.addEdge(u, v, std::stod(p), 1.0);  // Can be physical distances
        }
        
        // Get sources
        getline(sourcesFile, source);
        std::istringstream sourceSS(source);
        std::string s;
        while (getline(sourceSS, s, ' ')) {
            graph.addSource(s);
        }
        
        // Close all files
        edgelistFile.close();
        probabilitiesFile.close();
        sourcesFile.close();

        // Update graph
        graph.update();

        // Output graph information
        numVertices = graph.vertexSize();
        numEdges = graph.edgeSize();
        density = (2.0 * numEdges) / (numVertices * (numVertices - 1.0));
        mh << "#vertex = " << numVertices << ", #edge = " << numEdges << ", density = " << density << "\n";

        // Check graph size
        if (numEdges == 0) 
            throw std::runtime_error("ERROR: The graph is empty!");

        // Output graph as dot format
        if (opt["graph"]) graph.dump(std::cout);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    // Path survival probabilities
    try {
        std::map<std::string,double> min_path;     // Shortest distance to any source
        std::map<std::string,double> per_vertex;   // Per-vertex path survival probabilities
        for (std::string v: graph.getVertices()) {
            min_path[v] = INF;
            per_vertex[v] = 0.0;
        }

        // Initialize Frontier-based search and Degree constraint
        FrontierBasedSearch fbs(graph, 1);
        IntRange ZeroOrTwo(0, 2, 2);
        IntRange OnlyOne(1, 1);
        IntRange OnlyZero(0, 0);
        DegreeConstraint dc(graph, &ZeroOrTwo);
        for (std::string s: graph.getSources()) dc.setConstraint(s, &OnlyZero);
        
        // Construct ZDDs
        MessageHandler::showMessages(false);    // Turn off ZDD construction messages
        for (std::string v: graph.getVertices()) {
            if (graph.getSources().count(v)) continue;
            
            // Set degree constaint for non-source vertex
            long long numPaths = 0;
            dc.setConstraint(v, &OnlyOne);

            for (std::string s: graph.getSources()) {
                // Set degree constraint for source vertex
                dc.setConstraint(s, &OnlyOne);
                
                // Construct ZDD
                DdStructure<2> dd(zddIntersection(fbs, dc));
                dd.zddReduce();

                // Add number of single-source paths to source
                numPaths += dd.evaluate(ZddCardinality<long long,2>());

                // Compute path survival probabilities
                per_vertex[v] += dd.evaluate(PathSurvival(graph));

                // Update shortest path to source
                min_path[v] = std::min(min_path[v], dd.evaluate(MinDist(graph)));

                // Reset degree constriant
                dc.setConstraint(s, &OnlyZero);
            }
            // Reset degree constraint
            dc.setConstraint(v, &ZeroOrTwo);

            // Average path survival probabilties
            per_vertex[v] /= numPaths;
        }
        MessageHandler::showMessages();

        // Per-vertex path survival probabilities
        mh << "Per-Vertex Path Survival Probability:\n";
        if (!opt["adjusted"]) {
            mh << std::left << std::setw(5) << "v" << std::setw(10) << "PV" << "\n";
            for (std::string v: graph.getVertices()) {
                mh << std::left 
                   << std::setw(5) << v 
                   << std::setw(10) << per_vertex[v] << "\n"; 
            }
        }
        else {
            mh << std::left << std::setw(5) << "v" << std::setw(10) << "PV" << std::setw(10) << "Adj. PV" << "\n";
            for (std::string v: graph.getVertices()) {
                mh << std::left 
                   << std::setw(5) << v 
                   << std::setw(10) << per_vertex[v]
                   << std::setw(10) << pow(per_vertex[v], density) << "\n"; 
            }
        }

        // Total path survival probability
        double numerator = 0.0;
        double denominator = 0.0;
        for (std::string v: graph.getVertices()) {
            if (graph.getSources().count(v)) continue;
            numerator += min_path[v] * pow(per_vertex[v], density);
            denominator += min_path[v];
        }
        mh << "Total Path Survival Probability: " << numerator / denominator << "\n";
    }
    catch (std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    mh.end("finished");
    return 0;
}