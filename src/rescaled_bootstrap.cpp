// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

NumericMatrix resample_stratum(CharacterVector clusters, int B);

void resample_with_zeros(std::map<std::string, int> &counts,
                         CharacterVector cluster_names,
                         int m_h);
void count_clusters(std::map<std::string, int> &counts,
                    CharacterVector these_clusters,
                    CharacterVector cluster_names);
void print_cluster_counts(std::map<std::string, int> &counts);

// CharacterVector ids - the id of each row; we return this along with the
//                       results purely for convenience
// CharacterVector clusters - for each obs, the cluster it is in
// NumericVector weights - for each obs, its analysis weight
// int B - the number of bootstrap reps
//NumericMatrix resample_stratum(CharacterVector clusters,
//
// [[Rcpp::export]]
NumericMatrix resample_stratum(CharacterVector clusters, int B) {

    // get the set of clusters and the count of obs in each one
    IntegerVector cluster_set = table(clusters);
    CharacterVector cluster_names = cluster_set.names();

    int n_h = cluster_set.size();
    int m_h = n_h - 1;

    int num_rows = clusters.size();

    // this matrix will hold the rescaling factors for the weights
    // for each obs (row) and boot strap rep (column)
    NumericMatrix result_scales = NumericMatrix(num_rows, B);

    // uncomment to print useful debugging messages
    /*
    Rprintf("\n cluster table: \n\n");
    Rf_PrintValue (cluster_set);

    Rprintf("\n number of clusters: %d\n", n_h);

    Rprintf("\n names attribute: \n\n");
    Rf_PrintValue(cluster_names);
    */

    std::map<std::string, int> resampled_clusters;

    double n_h_factor = (double) n_h / (double) m_h;

    // resample clusters and compute rescaling factors once for each bootstrap rep
    for (int b = 0; b < B; b++) {

        //Rprintf("starting iteration %d...\n", b+1);

        // resamplea  set of clusters
        resample_with_zeros(resampled_clusters, cluster_names, m_h);

        // go through each row and
        // compute resampled weights for this cluster resample
        for (int r=0; r < num_rows; r++) {

            double this_scale = n_h_factor *
                                ((double) resampled_clusters[as<std::string>(clusters[r])]);

            result_scales(r, b) = this_scale;

        }

        // uncomment for debugging
        /*
        Rprintf("\n cluster resamples: \n");
        print_cluster_counts(resampled_clusters);
        */

        // clear out the resampled clusters
        resampled_clusters.clear();

    }

    return(result_scales);

}

// resample a set of clusters
//
// cluster_names - a CharacterVector of cluster names to be resampled
// m_h - the number of resamples to take
void resample_with_zeros(std::map<std::string, int> &counts,
                         CharacterVector cluster_names,
                         int m_h) {

    RNGScope scope;

    // resample the clusters
    CharacterVector these_clusters = RcppArmadillo::sample(cluster_names,
                                                           m_h,
                                                           TRUE);

    // tabulate the results, including zero counts for those clusters that were not
    // drawn in the resample
    count_clusters(counts, these_clusters, cluster_names);

    return;

}

void count_clusters(std::map<std::string, int> &counts,
                    CharacterVector these_clusters,
                    CharacterVector cluster_names) {

    int n = these_clusters.size();

    for(int i=0; i < n; i++) {
        counts[as<std::string>(these_clusters[i])]++;
    }

    //print_cluster_counts(counts);

    // add clusters with no counts
    CharacterVector::iterator it;
    for(it = cluster_names.begin(); it != cluster_names.end(); ++it) {

        // if this name not in map, make an entry and set it to zero
        if(counts.find(as<std::string>(*it)) == counts.end()) {
            counts[as<std::string>(*it)] = 0;
        }

    }

    //print_cluster_counts(counts);

    return;
}

// print out a table of resampled cluster counts;
// useful for debugging
void print_cluster_counts(std::map<std::string, int> &counts) {

    Rprintf("\nSampled cluster counts:\n");

    std::map<std::string, int>::iterator it;
    int total = 0;

    for(it = counts.begin(); it != counts.end(); ++it) {
        Rprintf("%s => %d\n", it->first.c_str(), it->second);
        total = total + it->second;
    }

    Rprintf("\n==========\ntotal: %d\n\n", total);

}


