#include <Rcpp.h>
using namespace Rcpp;

/*
Apply viterbi algorithm to label a sequence of symbols with states.

input[pos]
symbol_cost[state, symbol]
transition_cost[prev_state, this_state]
initial_cost[state]
final_cost[state]

We find a sequence of states that minimizes the cost, where the cost is the sum of:

- cost of initial state
- cost of transition from initial state to first position
- cost at position 1
- cost of state transition to next position
- cost at position 2
...
- cost of state transition to final state
- cost of final state
*/
// [[Rcpp::export]]
IntegerVector viterbi(
        IntegerVector input, 
        NumericMatrix symbol_cost, 
        NumericMatrix transition_cost,
        NumericVector initial_cost,
        NumericVector final_cost) {
    int n_pos = input.length();
    int n_state = symbol_cost.nrow();
    if (transition_cost.nrow() != n_state ||
        transition_cost.ncol() != n_state ||
        initial_cost.length() != n_state ||
        final_cost.length() != n_state)
        stop("Wrong sized arguments.");
    
    IntegerMatrix came_from(n_pos+1, n_state);
    NumericVector score(n_state);
    NumericVector new_score(n_state);
    
    for(int state=0;state<n_state;state++)
        score[state] = initial_cost[state];
    
    for(int pos=0;;pos++) {
        // Transition to new state
        for(int new_state=0;new_state<n_state;new_state++) {
            int best_old_state = 0;
            double best_score = score[0] + transition_cost(0,new_state);
            for(int old_state=1;old_state<n_state;old_state++) {
                double this_score = score[old_state] + transition_cost(old_state,new_state);
                if (this_score < best_score) {
                    best_score = this_score;
                    best_old_state = old_state;
                }
            }
            
            new_score[new_state] = best_score;
            came_from(pos, new_state) = best_old_state;
        }
        
        if (pos >= n_pos) 
            break;
        
        // Cost of symbol in new state
        for(int state=0;state<n_state;state++) {
            score[state] = new_score[state] + symbol_cost(state,input[pos]);
        }
    }
    
    for(int state=0;state<n_state;state++)
        score[state] = new_score[state] + final_cost[state];

    // Start from highest score
    int current_state = 0;
    for(int state=1;state<n_state;state++)
        if (score[state] < score[current_state])
            current_state = state;
    
    // Rewind down best path
    IntegerVector assignment(n_pos);
    for(int pos=n_pos;pos>0;) {
        current_state = came_from(pos, current_state);
        pos--;
        assignment[pos] = current_state;
    }
    
    return assignment;
}
