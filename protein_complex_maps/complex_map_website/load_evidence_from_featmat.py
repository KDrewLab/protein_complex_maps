
import argparse
import gzip
import csv

import numpy as np

import protein_complex_maps.complex_map_website.complex_db as cdb

#kdrew: feature matrix column -> evidence name, in display order
EVIDENCE_COLUMNS = [('mann_correlation', 'Michaelis_correlation'),
                    ('cong_RFraw', 'Humphreys_RF'),
                    ('hart_wmm_score', 'Hart_wmm')]

def confidence(value, low_cutoff, high_cutoff):
    if value >= high_cutoff:
        return 'high'
    elif value < low_cutoff:
        return 'low'
    return 'medium'

def main():

    parser = argparse.ArgumentParser(description="Loads edge evidence (with confidence bins) from the feature matrix used for predictions")
    parser.add_argument("--feature_matrix", action="store", dest="feature_matrix", required=True,
                                    help="Filename of feature matrix, csv, optionally gzipped (format: acc1,acc2,mann_correlation,cong_RFraw,hart_wmm_score); a blank value means no data for that pair")
    parser.add_argument("--humphreys_high", action="store", dest="humphreys_high", type=float, default=0.99,
                                    help="Humphreys_RF values >= this are high confidence (default: 0.99)")
    parser.add_argument("--humphreys_low", action="store", dest="humphreys_low", type=float, default=0.3,
                                    help="Humphreys_RF values < this are low confidence (default: 0.3)")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    #kdrew: evidence is fully replaced on each load (also picks up schema changes to the evidence table)
    cdb.Evidence.__table__.drop(db.engine, checkfirst=True)
    db.create_all()

    acc_dict = dict()
    for p in db.session.query(cdb.Protein).all():
        acc_dict[p.id] = p.uniprot_acc

    #kdrew: a pair can be stored as more than one edge (both protein orders), so map pair -> list of edge ids
    pair_edges = dict()
    for e in db.session.query(cdb.Edge).all():
        pair = frozenset([acc_dict[e.protein_key], acc_dict[e.protein_key2]])
        pair_edges.setdefault(pair, []).append(e.id)

    if args.feature_matrix.endswith('.gz'):
        featmat_file = gzip.open(args.feature_matrix, 'rb')
    else:
        featmat_file = open(args.feature_matrix, 'rb')
    reader = csv.reader(featmat_file)
    header = next(reader)
    column_index = [header.index(col) for col, _ in EVIDENCE_COLUMNS]

    pair_values = dict()
    for row in reader:
        pair = frozenset(row[:2])
        if pair in pair_edges:
            pair_values[pair] = [float(row[i]) if row[i] != '' else None for i in column_index]
    featmat_file.close()

    print("edges: %s, distinct pairs: %s, pairs found in feature matrix: %s" % (sum([len(v) for v in pair_edges.values()]), len(pair_edges), len(pair_values)))

    #kdrew: Michaelis_correlation and Hart_wmm: top/bottom 25% of values over the map's edges
    #kdrew: Humphreys_RF: fixed cutoffs, its values below 0.3 are no more likely than background to be true interactions (benchmark check, 2026-11)
    cutoffs = dict()
    for i, (_, name) in enumerate(EVIDENCE_COLUMNS):
        if name == 'Humphreys_RF':
            cutoffs[name] = (args.humphreys_low, args.humphreys_high)
        else:
            values = [v[i] for v in pair_values.values() if v[i] != None]
            cutoffs[name] = (np.percentile(values, 25), np.percentile(values, 75))
        print("%s: low < %.4g, high >= %.4g" % (name, cutoffs[name][0], cutoffs[name][1]))

    evidence_rows = []
    for pair, values in pair_values.items():
        for edge_id in sorted(pair_edges[pair]):
            for (_, name), value in zip(EVIDENCE_COLUMNS, values):
                if value == None:
                    continue
                evidence_rows.append({'edge_key': edge_id, 'evidence_type': name, 'value': value,
                                      'confidence': confidence(value, cutoffs[name][0], cutoffs[name][1])})

    #kdrew: insert in edge order then type order, so Edge.evidences (ordered by id) displays types consistently
    type_order = dict([(name, i) for i, (_, name) in enumerate(EVIDENCE_COLUMNS)])
    evidence_rows.sort(key=lambda r: (r['edge_key'], type_order[r['evidence_type']]))
    db.session.execute(cdb.Evidence.__table__.insert(), evidence_rows)
    db.session.commit()

    print("loaded %s evidence rows" % len(evidence_rows))
    for _, name in EVIDENCE_COLUMNS:
        counts = dict()
        for r in evidence_rows:
            if r['evidence_type'] == name:
                counts[r['confidence']] = counts.get(r['confidence'], 0) + 1
        print("%s: high %s, medium %s, low %s" % (name, counts.get('high', 0), counts.get('medium', 0), counts.get('low', 0)))

if __name__ == "__main__":
    main()
