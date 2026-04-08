#include "slp/types.hpp"
#include "slp/algorithm.hpp"

#include <iostream>
#include <chrono>
using namespace std;


int main() {
    cout << "Greedy Potential" << endl;
    slp::Z2Matrix G = {{1,0,1,0,1}, {0,1,1,0,0}, {1,0,1,1,1}, {1,1,0,1,1}}; // {{1,1,1,0}, {0,1,1,1}}; // {{1,1,0,0,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},{0,1,0,0,0,0,1,0,0},{1,0,0,0,0,0,1,0,0},{1,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0},{1,0,0,0,0,1,1,0,0},{0,0,0,1,1,1,0,0,0},{0,0,0,0,1,0,1,0,0},{0,1,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0,0},{0,1,0,0,0,0,0,1,0},{0,0,0,0,0,0,1,1,0},{0,1,0,0,1,0,0,0,1},{0,0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,1,1},{0,1,1,0,1,1,0,0,0},{0,0,0,0,0,1,0,0,0},{1,0,1,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,1},{0,0,1,0,0,1,0,0,1},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0,1}}; 

    slp::Options options;
    options.verbose = true;
    options.strategy = slp::SearchStrategy::GreedyPotential;
    options.alpha = 0.2;
    options.reachable_strategy = slp::ReachableStrategy::BacktracingSparseAware;
    options.nearest = 0;
    options.max_level = 12;


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
    cout << "Backtrack Potential" << endl;
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
    cout << "Boyar Peralta" << endl;
    options.strategy = slp::SearchStrategy::BP;

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
    cout << "Random Normal Boyar Peralta" << endl;
    options.strategy = slp::SearchStrategy::RNBP;

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

    cout << "-----------------------------------------" << endl;
    cout << "A1 Boyar Peralta" << endl;
    options.strategy = slp::SearchStrategy::A1;

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
    cout << "A2 Boyar Peralta" << endl;
    options.strategy = slp::SearchStrategy::A2;

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
    cout << "Paar1" << endl;
    options.strategy = slp::SearchStrategy::Paar1;

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
