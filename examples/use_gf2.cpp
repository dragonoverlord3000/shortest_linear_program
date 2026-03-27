#include "slp/types.hpp"
#include "slp/algorithm.hpp"

#include <iostream>
#include <chrono>
using namespace std;


int main() {
    cout << "Z2 Greedy Potential" << endl;
    slp::Z2Matrix G = {{1,1,0,0,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},{0,1,0,0,0,0,1,0,0},{1,0,0,0,0,0,1,0,0},{1,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0},{1,0,0,0,0,1,1,0,0},{0,0,0,1,1,1,0,0,0},{0,0,0,0,1,0,1,0,0},{0,1,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0,0},{0,1,0,0,0,0,0,1,0},{0,0,0,0,0,0,1,1,0},{0,1,0,0,1,0,0,0,1},{0,0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,1,1},{0,1,1,0,1,1,0,0,0},{0,0,0,0,0,1,0,0,0},{1,0,1,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,1},{0,0,1,0,0,1,0,0,1},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0,1}}; 

    slp::Options options = {true, slp::SearchStrategy::GreedyPotential, 0.2, slp::ReachableStrategy::BacktracingSparseAware, 12};
    auto t0 = std::chrono::steady_clock::now();
    slp::Result result = slp::gf2::run(G, options);
    auto t1 = std::chrono::steady_clock::now();

    cout << "Additions before: " << result.additions_before << endl;
    cout << "Additions after: " << result.additions_after << endl;
    cout << "Transformation:" << endl;
    for(std::pair<std::size_t,std::size_t>& t: result.method.additions) 
        cout << "{" << t.first << ", " << t.second << "}" << endl;
    cout << "Duration: " << static_cast<std::chrono::nanoseconds>(t1 - t0) << " ns" << endl;


    cout << "-----------------------------------------" << endl;
    cout << "Z2 Backtrack Potential" << endl;
    options.strategy = slp::SearchStrategy::BacktrackingPotential;

    t0 = std::chrono::steady_clock::now();
    result = slp::gf2::run(G, options);
    t1 = std::chrono::steady_clock::now();
    
    cout << "Additions before: " << result.additions_before << endl;
    cout << "Additions after: " << result.additions_after << endl;
    cout << "Transformation:" << endl;
    for(std::pair<std::size_t,std::size_t>& t: result.method.additions) 
        cout << "{" << t.first << ", " << t.second << "}" << endl;
    cout << "Duration: " << static_cast<std::chrono::nanoseconds>(t1 - t0) << " ns" << endl;

    cout << "-----------------------------------------" << endl;
    cout << "Z2 Boyar Peralta" << endl;
    options.strategy = slp::SearchStrategy::BoyarPeralta;

    t0 = std::chrono::steady_clock::now();
    result = slp::gf2::run(G, options);
    t1 = std::chrono::steady_clock::now();
    
    cout << "Additions before: " << result.additions_before << endl;
    cout << "Additions after: " << result.additions_after << endl;
    cout << "Transformation:" << endl;
    for(std::pair<std::size_t,std::size_t>& t: result.method.additions) 
        cout << "{" << t.first << ", " << t.second << "}" << endl;
    cout << "Duration: " << static_cast<std::chrono::nanoseconds>(t1 - t0) << " ns" << endl;


    return 0;
}
