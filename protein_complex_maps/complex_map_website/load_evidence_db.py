
import argparse
import numpy as np
import itertools as it

import pandas as pd
import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    #kdrew: should really generalize the bait prey files and make them optional, also need to reorganize the readin through pandas
    #kdrew: but just trying to get the db loaded

    parser = argparse.ArgumentParser(description="Loads edge sql tables from input files")
    parser.add_argument("--edge_file", action="store", dest="edge_file", required=True, 
                                    help="Filename edge table")
    #parser.add_argument("--feature_matrix", action="store", dest="feature_matrix", required=False, default=None,
    #                                help="Filename of feature matrix (format: acc1,acc2,mann_correlation,cong_RFraw,hart_wmm_score,in_train_set,in_test_set)")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    #Michaelis_pairs = set()
    #Hart_pairs = set()
    #Humphreys_pairs = set()
    #if args.feature_matrix != None:
    #    featmat = pd.read_csv(args.feature_matrix)


    edge_table_file = open(args.edge_file,"rb")
    for line in edge_table_file.readlines():

        #kdrew: if header do not parse
        #id1     score   fractions       bioplex hein    bioplex_prey    hein_prey
        if 'score' in line:
            continue

        print line
        split_line = line.split(',')
        print split_line

        #kdrew: example: 0_P25359,0_P38792,P25359,P38792,0_P25359 (pp) 0_P38792,P25359,P38792,0.7471802830696106,pp,mann_correlation cong_RFraw hart_wmm_score

        prot1 = split_line[2]
        prot2 = split_line[3]

        try:
            score = float(split_line[7])
        except ValueError:
            score = np.nan

        evidence_str = split_line[9].strip()
        evidence_dict = dict()

        p1 = db.session.query(cdb.Protein).filter_by(uniprot_acc=prot1).first()
        p2 = db.session.query(cdb.Protein).filter_by(uniprot_acc=prot2).first()

        if p1 and p2:
            #kdrew: enforce order on protein ids
            if p2.id < p1.id:
                p2, p1 = p1, p2
            print "protein id1: %s" % p1.id
            print "protein id2: %s" % p2.id

            edge = db.session.query(cdb.Edge).filter_by(protein_key=p1.id).filter_by(protein_key2=p2.id).first()

            #kdrew: temp comment out for testing/debug
            #edge = cdb.get_or_create(db, cdb.Edge, 
            #                            protein_key = p1.id,
            #                            protein_key2 = p2.id,
            #                            score = score,
            #                            )
            #db.session.add(edge)
            #db.session.commit()

            if edge:
                evidence = cdb.get_or_create(db, cdb.Evidence,
                                            edge_key = edge.id,
                                            evidence_type = evidence_str
                                            )
                db.session.add(evidence)
                db.session.commit()
            else:
                print("edge not found: %s , %s" % (prot1, prot2))

        else:
            print "Cannot find proteins %s" % (id1)


if __name__ == "__main__":
    main()


