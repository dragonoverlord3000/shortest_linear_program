#include <iostream>
using namespace std;

#include "slp/algorithm.hpp"

int main() {
    cout << "Z2 Greedy Potential" << endl;
    slp::Z2Matrix G = {{1,1,0,0,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},{0,1,0,0,0,0,1,0,0},{1,0,0,0,0,0,1,0,0},{1,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0},{1,0,0,0,0,1,1,0,0},{0,0,0,1,1,1,0,0,0},{0,0,0,0,1,0,1,0,0},{0,1,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0,0},{0,1,0,0,0,0,0,1,0},{0,0,0,0,0,0,1,1,0},{0,1,0,0,1,0,0,0,1},{0,0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,1,1},{0,1,1,0,1,1,0,0,0},{0,0,0,0,0,1,0,0,0},{1,0,1,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,1},{0,0,1,0,0,1,0,0,1},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0,1}}; //{{1, 1, 1, 0}, {0, 1, 1, 1}};

    slp::Options options = {true, 0.2, slp::SearchStrategy::GreedyPotential, 12};
    slp::Result result = slp::gf2::run(G, options);

    cout << "Additions before: " << result.additions_before << endl;
    cout << "Additions after: " << result.additions_after << endl;
    cout << "Transformation:" << endl;
    for(std::pair<std::size_t,std::size_t>& t: result.method.additions) 
        cout << "{" << t.first << ", " << t.second << "}" << endl;


    cout << "-----------------------------------------" << endl;
    cout << "Z2 Backtrack Potential" << endl;
    options.strategy = slp::SearchStrategy::BacktrackingPotential;
    result = slp::gf2::run(G, options);
    
    cout << "Additions before: " << result.additions_before << endl;
    cout << "Additions after: " << result.additions_after << endl;
    cout << "Transformation:" << endl;
    for(std::pair<std::size_t,std::size_t>& t: result.method.additions) 
        cout << "{" << t.first << ", " << t.second << "}" << endl;

    return 0;
}
