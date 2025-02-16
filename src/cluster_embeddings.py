import os
import argparse
import flexynesis
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Cluster flexynesis embeddings", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--embeddings", help="Path to train/test embeddings", type=str, required = True)
    parser.add_argument("--clin", help="Path to clin table", type=str, required = True)
    args = parser.parse_args()
    
    E = pd.read_csv(args.embeddings, index_col=0)
    cluster_labels, G, partition = flexynesis.get_optimal_clusters(E, min_k = 18, max_k = 24)
    
    labels = pd.read_csv(args.clin, index_col=0)
    known_labels = [labels.loc[x]['cohort'] for x in E.index]
    
    print(flexynesis.compute_ami_ari(known_labels, cluster_labels))

    clusters = pd.DataFrame({'cluster': cluster_labels, 'sample': E.index, 'label': known_labels})

    clusters.to_csv("clusters.csv", header=True)

if __name__ == "__main__":
    main()
