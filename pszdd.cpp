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

// Stochastic graph class
// Features: edge survival probability, edge length
// Features: source vertices
class StochGraph: public Graph {
private:
    std::set<std::string> vertices, sources;
    std::map<std::pair<std::string,std::string>,double> probabilities, lengths;

public:
    void addEdge(std::string vertexName1, std::string vertexName2, double prob = 1.0, double len = 1.0) {
        Graph::addEdge(vertexName1, vertexName2);
        addVertex(vertexName1);
        addVertex(vertexName2);
        
        assert(0.0 <= prob && prob <= 1.0);
        probabilities[std::make_pair(vertexName1, vertexName2)] = prob;

        assert(len > 0.0);
        lengths[std::make_pair(vertexName1, vertexName2)] = len;
    }

    double getProb(EdgeNumber e) const {
        return probabilities.at(edgeName(e));
    }

    double getLen(EdgeNumber e) const {
        return lengths.at(edgeName(e));
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

    void readEdges(std::string const& filename, bool len) {
        std::ifstream is(filename.c_str(), std::ios::in);
        readEdges(is, len);
        Graph::update();
    }

    void readSources(std::string const& filename) {
        std::ifstream is(filename.c_str(), std::ios::in);
        std::string s;
        while (getline(is, s, ' ')) {
            addSource(s);
        }
    }

private:
    void readEdges(std::istream& is, bool len) {
        std::string s;
        std::map<std::string,std::string> edge;

        while (is) {
            char ss = is.get();

            if (isspace(ss)) {
                if (!s.empty()) {
                    if (!edge.count("v1")) {
                        edge["v1"] = s;
                    }
                    else if (!edge.count("v2")) {
                        edge["v2"] = s;
                    }
                    else if (!edge.count("p")) {
                        edge["p"] = s;
                    }
                    else if (!edge.count("l") and len) {
                        edge["l"] = s;
                    }
                    else {
                        throw std::runtime_error("ERROR: Exceeded expected number of tokens in a line");
                    }
                    s.clear();
                }

                if (ss == '\n') {
                    // Default values
                    edge["l"] = edge.count("l") ? edge["l"] : "1.0";
                    if (edge.count("v1") && edge.count("v2") && edge.count("p")) {
                        addEdge(edge["v1"], edge["v2"], std::stod(edge["p"]), std::stod(edge["l"])); 
                        edge.clear();
                    }
                    else {
                        throw std::runtime_error("ERROR: Less than minimum number of tokens in a line");
                    }
                }
            }
            else {
                s += ss;
            }
        }
    }
};

// DdEval for finding the shortest distance
class MinDist: public DdEval<MinDist,long double> {
private:
    int const numEdges;
    StochGraph const& graph;

public:
    MinDist(StochGraph const& graph) 
        : graph(graph), numEdges(graph.edgeSize()) {
    }

    void evalTerminal(long double &v, int id) const {
        v = id ? 0.0 : INF;
    }

    void evalNode(long double &v, int level, DdValues<long double,2> const& values) const {
        v = std::min(values.get(0), values.get(1) + graph.getLen(numEdges - level));
    }
};

// DdEval for calculating the per-vertex path survival probabilities
class PathSurvival: public DdEval<PathSurvival,double> {
private:
    bool const prob;
    int const numEdges;
    StochGraph const& graph;

public:
    PathSurvival(StochGraph const& graph, bool prob = true) 
        : graph(graph), numEdges(graph.edgeSize()), prob(prob) {
    }

    void evalTerminal(double &v, int id) const {
        v = id ? 1.0 : 0.0;
    }

    void evalNode(double &v, int level, DdValues<double,2> const& values) const {
        double v0 = values.get(0);
        double v1 = values.get(1);
        if (prob) v1 *= graph.getProb(numEdges - level);
        v = v0 + v1;
    }
};


std::map<std::string,bool> opt;
std::string options[][2] = { //
        {"len", "Flag including edge length in edgelist"}, //
        {"graph", "Dump input graph to STDOUT in DOT format"}, //
        {"adjusted", "Output adjusted per-vertex path survival probabilities"}, //
        {"nodes", "Output total ZDD nodes per vertex"}, //
        {"paths", "Output total single-source paths per vertex"}, //
    };

void usage(char const* cmd) {
    std::cerr << "usage: " << cmd
              << " <edgelist> <sources> [ <options>... ]\n";
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
StochGraph graph;                               // Stochastic graph object
double density;                                 // Graph density
int numEdges, numVertices;                      // Number of edges and vertices
std::string edgelistFilename, sourcesFilename;  // Input filenames

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
        // Read files
        graph.readEdges(edgelistFilename, opt["len"]);
        graph.readSources(sourcesFilename);

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
        std::map<std::string,long long> num_path;       // Number of single-source paths
        std::map<std::string,long double> sum_path;     // Sum of path survival probabilities
        std::map<std::string,long double> min_path;     // Shortest distance to any source
        std::map<std::string,long long> num_node;       // Number of ZDD nodes
        std::map<std::string,long double> per_vertex;   // Per-vertex path survival probabilities
        for (std::string v: graph.getVertices()) {
            num_path[v] = 0;
            sum_path[v] = 0;
            min_path[v] = INF;
            num_node[v] = 0;
            per_vertex[v] = 0;
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
            dc.setConstraint(v, &OnlyOne);

            for (std::string s: graph.getSources()) {
                // Set degree constraint for source vertex
                dc.setConstraint(s, &OnlyOne);
                
                // Construct ZDD
                DdStructure<2> dd(zddIntersection(fbs, dc));
                dd.zddReduce();

                // Compute number of single-source paths
                num_path[v] += dd.evaluate(PathSurvival(graph, false));

                // Compute sum of path survival probabilities
                sum_path[v] += dd.evaluate(PathSurvival(graph, true));

                // Update shortest distance to any source
                min_path[v] = std::min(min_path[v], dd.evaluate(MinDist(graph)));

                // Record number of ZDD nodes
                num_node[v] += dd.size();

                // Reset degree constriant
                dc.setConstraint(s, &OnlyZero);
            }
            // Reset degree constraint
            dc.setConstraint(v, &ZeroOrTwo);

            // Per-vertex path survival probability
            per_vertex[v] = sum_path[v] / num_path[v];
        }
        MessageHandler::showMessages(true);

        // Per-vertex path survival probabilities
        mh << "Per-Vertex Path Survival Probability (PV PSP):\n";
        mh << std::left << std::setw(10) << "v" << std::setw(15) << "PV PSP";
        if (opt["adjusted"]) mh << std::left << std::setw(15) << "Adj. PV PSP";
        if (opt["nodes"]) mh << std::left << std::setw(15) << "ZDD Nodes";
        if (opt["paths"]) mh << std::left << std::setw(15) << "SS Paths";
        mh << "\n";
        for (std::string v: graph.getVertices()) {
            mh << std::left << std::setw(10) << v << std::setw(15) << per_vertex[v];
            if (opt["adjusted"]) mh << std::left << std::setw(15) << pow(per_vertex[v], density);
            if (opt["nodes"]) mh << std::left << std::setw(15) << num_node[v];
            if (opt["paths"]) mh << std::left << std::setw(15) << num_path[v];
            mh << "\n"; 
        }

        // Total path survival probability
        long double numerator = 0.0;
        long double denominator = 0.0;
        for (std::string v: graph.getVertices()) {
            if (graph.getSources().count(v)) continue;
            numerator += min_path[v] * pow(per_vertex[v], density);
            denominator += min_path[v];
        }
        long double total = numerator / denominator;
        mh << "Total Path Survival Probability (Total PSP): " << total << "\n";
    }
    catch (std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    mh.end("finished");
    return 0;
}