#include <Rcpp.h>
using namespace Rcpp;

/*
Given read sequence and quality strings.
Return a good point to clip off the end of the read, i.e. the number of characters to retain.
Maximizes the number of high-quality non-G bases minus the number of low quality non-G bases times penalty.
*/
// [[Rcpp::export]]
double quality_clip(std::string seq, std::string qual, char c, double penalty) {
    double best_score = 0;
    int best_end = 0;
    double score = 0;
    int n = seq.size();
    for(int i = 0;i<n;i++) {
        if (seq[i] != 'G') {
            if (qual[i] < c)
                score -= penalty;
            else
                score += 1;
        }
        if (score >= best_score) {
            best_end = i;
            best_score = score;
        }
    }
    
    return best_end+1;
}

/*
In a string, up to the nth character,
Find the best span containing mostly the c character.
score = matches - penalty * mixmatches
Returns 1-based range of the span, and the score.
*/
// [[Rcpp::export]]
NumericVector scan(std::string s, char c, int n, double penalty) {
    int start = 0;
    double score = 0.0;
    double best_score = 0.0;
    int best_start = 0, best_end = 0;
    for(int i=0;i<n;i++) {
        if (score <= 0) {
            start = i;
            score = 0;
        } 
        
        if (s[i] == c)
            score += 1;
        else
            score -= penalty;
        
        if (score >= best_score) {
            best_score = score;
            best_start = start;
            best_end = i+1;
        }
    }
    
    return NumericVector::create(best_start+1, best_end, best_score);
}

/*
In a string, up to the nth character,
Find the best span containing mostly the c character, then possibly continuing into a suffix sequence.
The suffix sequence may extend beyond the nth character.
score = matches - penalty * mixmatches
Returns 1-based range of the span, the suffix match length, and the score.
*/
// [[Rcpp::export]]
NumericVector scan_suffix(std::string s, char c, std::string suffix, int n, double penalty, double suffix_penalty) {
    int start = 0;
    double score = 0.0;
    int n_suffix = suffix.size();
    int n_s = s.size();
    
    double best_score = 0.0;
    int best_start = 0, best_end = 0, best_suffix = 0;
    for(int i=0;i<n_s;i++) {
        if (score <= 0) {
            start = i;
            score = 0;
        } 
        
        // Look for the suffix
        // Could add early abort...
        double suffix_score = score;
        for(int j=0;j<n_suffix && i+j<n_s;j++) {
            if (s[i+j] == suffix[j])
                suffix_score += 1;
            else
                suffix_score -= suffix_penalty;
            
            if (suffix_score >= best_score) {
                best_score = suffix_score;
                best_start = start;
                best_end = i;
                best_suffix = j+1;
            }
        }
        
        // Stop at end of good quality
        if (i >= n) break;
        
        if (s[i] == c)
            score += 1;
        else
            score -= penalty;
        
        if (score >= best_score) {
            best_score = score;
            best_start = start;
            best_end = i+1;
            best_suffix = 0;
        }
        
    }
    
    return NumericVector::create(best_start+1, best_end, best_suffix, best_score);
}


/*
In a string, strictly from the specified start, up to the nth character,
Find the best span containing mostly the c character.
score = matches - penalty * mismatches
Returns 1-based range of the span, and the score.
*/
// [[Rcpp::export]]
NumericVector scan_from(std::string s, char c, int start, int n, double penalty) {
    // Convert arguments to 0-based
    start -= 1;
    
    double score = 0.0;
    double best_score = 0.0;
    int best_end = start;
    
    for(int i=start;i<n;i++) {
        if (s[i] == c)
            score += 1;
        else
            score -= penalty;
        
        if (score >= best_score) {
            best_score = score;
            best_end = i+1;
        }
    }
    
    return NumericVector::create(start+1, best_end, best_score);
}