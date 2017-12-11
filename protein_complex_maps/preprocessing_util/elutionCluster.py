import os.path
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

def main():
    
    parser = argparse.ArgumentParser(description="Cluster a single or multiple joined elution files")
    parser.add_argument("--input_elutions", action='store', dest='input_elutions', required=True, nargs='+', help="Input csv elution files")
    parser.add_argument("--outfile", action="store", dest="outfile", required=True, help="Name of outfile. Must end with .xlsx")
    parser.add_argument("--hclust_method", action="store",dest="hclust_method",required=False,default='average',choices=["single","complete","average","weighted","centroid","median","ward"],help="Method for hierarchical clustering of rows")
    parser.add_argument("--hclust_metric", action="store",dest="hclust_metric",required=False,default="euclidean",choices=["euclidean","canberra","braycurtis","pearson","spearman"],help="Distance or correlation metric used for clustering")
    parser.add_argument("--path_to_annotations", action='store', dest='path_to_annotations', required=False, default=None, help="Path to annotations file. A tab-delimited file with protein ids in first column and annotations in second")
    parser.add_argument("--elution_separator", action='store', dest='sep', required=False, default="\t", help="Separator used for input elution profiles, default='\t'")
    
    args = parser.parse_args()
    
    if args.path_to_annotations != None:
        assert os.path.isfile(args.path_to_annotations), "File {} not found".format(args.path_to_annotations)
        
    assert args.outfile.split(".")[-1] == "xlsx", "Outfile must be .xlsx - sorry"
    
    ## Read in and join files
    
    joined_elutions = pd.DataFrame()
    for f in args.input_elutions:
        #df = pd.read_table(f,index_col=0)
        df = pd.read_csv(f,index_col=0, sep=args.sep)
        if "TotalCount" in df.columns:
            df.drop("TotalCount",inplace=True,axis=1)
        joined_elutions = joined_elutions.join(df, how='outer')

    # is this kosher? Put it in because I got and error clustering: ValueError: The condensed distance matrix must contain only finite values.
    #kdrew: should be, likely due to proteins missing in one or the other joined elution profiles from the outer join above, zero makes sense if not seen
    joined_elutions.fillna(0,inplace=True) 
    
    print "Clustering method={}, metric={}".format(args.hclust_method,args.hclust_metric)
        
    labels_index = pd.Series(joined_elutions.index,index=range(len(joined_elutions))) # Create indexer to manage scipy output
    
    if args.hclust_metric == "pearson":
        corr_mat = np.around( np.corrcoef(joined_elutions), decimals = 4)
        dist_mat = (1. - corr_mat) / 2. # BJL: Invert correlations and move to 0-1 axis, i.e. convert them to distances
        distArray = dist_mat[ np.triu_indices(dist_mat.shape[0],1) ]
        link = hclust.linkage(distArray, method=args.hclust_method)
    elif args.hclust_metric == "spearman":
        corr_mat = np.around( stats.spearmanr(joined_elutions.T)[0], decimals = 4 )
        dist_mat = (1. - corr_mat) / 2.
        distArray = dist_mat[ np.triu_indices(dist_mat.shape[0],1) ]
        link = hclust.linkage(distArray, method=args.hclust_method)
    elif args.hclust_metric in ("euclidean","canberra","braycurtis"):
        link = hclust.linkage(joined_elutions, method=args.hclust_method, metric=args.hclust_metric)
    else:
        raise Exception("Clustering metric unknown: {}".format(args.hclust_metric))# shouldn't be necessary b/c arparse should catch it with "choices", but whatever
    
    leaves_list = hclust.leaves_list(link)
    ordered_index = labels_index.iloc[leaves_list] # reorder index from hier clustering
    ordered_index.to_csv("ordered_index.csv")

    column_order = list(joined_elutions.columns)

    print "Summing"
    colSums = pd.Series(joined_elutions.sum(),name="colSum") # add to df later
    joined_elutions["rowSum"] = joined_elutions.sum(axis=1)
    clustered_index_name = "_".join(["order",args.hclust_metric,args.hclust_method])
    joined_elutions[clustered_index_name] = range(len(ordered_index))
    column_order = ['rowSum', clustered_index_name] + column_order
    
    #kdrew: create empty annotation dataframe in case no annotations are read in
    if args.path_to_annotations != None:
        print "Adding annotations"
        annots = pd.read_table(args.path_to_annotations,index_col=0)
        #print "length of annotation df: {}".format(len(annots)) # Debugging
        joined_elutions = annots.join(joined_elutions,how='right')
        column_order = list(annots.columns) + column_order
    
    # reorder by hier clustering and ordering column
    clustered_df = joined_elutions.reindex( ordered_index )
    clustered_df[clustered_index_name] = range(len(ordered_index))
    
    # reorder columns
    #column_order = annots.columns + ["rowSum",clustered_index_name] + [i for i in clustered_df.columns if i not in annots.columns + ["rowSum",clustered_index_name]] # haaacky bleerggh (BJL)
    reordered_df = clustered_df[column_order]
    
    # add column sum row
    outdf = reordered_df.append(colSums) # BJL: Must add annots, then row sum, the cluster order, then col sums -- in that order
    
    print "Writing excel file"
    outdf.to_excel(args.outfile)
    
if __name__ == '__main__':
    main()
