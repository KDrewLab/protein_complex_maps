
from flask import Flask
#from flask.ext.sqlalchemy import SQLAlchemy
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import func, or_, and_
from sqlalchemy.orm import subqueryload

import itertools as it

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db/test.db'

app.config['SECRET_KEY'] = 'please, tell nobody'

db = SQLAlchemy(app)

def get_db():
    return db

def get_app():
    return app

def get_or_create(db, model, **kwargs):
    session = db.session

    table_columns = [x.key for x in model.__table__.columns]
    missing_args = set(kwargs.keys()) - set(table_columns)
    if len(missing_args) > 0:
        print "Warning: ignoring keywords because not in Model: %s" % ','.join(missing_args)

        #kdrew: remove from kwargs
        for ma in missing_args:
            del kwargs[ma]

    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        session.commit()
    return instance


class Complex(db.Model):
    """A single complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_id = db.Column(db.Integer, unique=True, index=True)
    #kdrew: uses table name for ProteinComplexMapping class (annoying sqlalchemy magic)
    proteins = db.relationship('Protein', secondary='protein_complex_mapping', back_populates='complexes')
    enrichments = db.relationship('ComplexEnrichment')
    rnp_stats = db.relationship('ComplexRNPStats')
    corum_id = db.Column(db.Integer)
    humap_id = db.Column(db.Integer)
    rnp_label = db.Column(db.String(255))

    def complex_link(self,):
        retstr = "<a href=displayComplexes?complex_key=%s>%s</a>" % (self.complex_id, self.complex_id)
        return retstr

    def complex_link_button(self,):
        retstr = "<a href=displayComplexes?complex_key=%s class='btn btn-default btn-sm' role='button'>%s</a>" % (self.complex_id, self.complex_id)
        return retstr

    def corum_complex_link(self,):
        retstr = ""
        if self.corum_id != -1:
            retstr = retstr + "<a href=http://mips.helmholtz-muenchen.de/corum/#?complexID=%s target=\"_blank\">CORUM</a> " % (self.corum_id)
        return retstr

    def humap_complex_link(self,):
        retstr = ""
        if self.humap_id != -1:
            retstr = retstr + "<a href=http://hu.proteincomplexes.org/displayComplexes?complex_key=%s target=\"_blank\">hu.MAP</a> " % (self.humap_id)
        return retstr

    #kdrew: this used to generate every pairwise combination of the
    #complex's proteins and query the edge table once per pair (O(n^2)
    #queries just to find which pairs happen to have an edge) -- for a
    #148-protein complex that's ~10,900 individual queries before even
    #getting to get_proteins()/evidences per edge. This schema has no
    #direct complex->edge mapping table (unlike humap2/humap3v1, which
    #already moved to one), so instead this fetches every edge touching
    #ANY of the complex's proteins in one query (indexed on protein_key/
    #protein_key2), then filters in Python for edges where BOTH ends are
    #complex members -- one query total instead of one per candidate
    #pair, plus subqueryload to batch-load evidences instead of one
    #query per edge in complex.html's "for evidence in edge.evidences".
    def edges(self,):
        protein_ids = [p.id for p in self.proteins]
        id_set = set(protein_ids)
        candidates = db.session.query(Edge).filter(
            or_(Edge.protein_key.in_(protein_ids), Edge.protein_key2.in_(protein_ids))
        ).options(subqueryload(Edge.evidences)).all()
        es = [e for e in candidates if e.protein_key in id_set and e.protein_key2 in id_set]

        return sorted(set(es), key=lambda es: es.score, reverse=True)

    #kdrew: rnp select is a class of rnps which we have data for in DIFFRAC and most of the subunits show signal
    def is_rnp_select(self,):
        return (self.is_rnp() and self.rnp_stats[0].sec2_cntl_median_correlation >= 0.75 and self.rnp_stats[0].sec2_fraction_fdr5 > 0.5)

    def is_rnp(self,):
        return self.rnp_label != ''

    def proteinsOrderedByEvidence(self,):
        return sorted(self.proteins, key=lambda x: len(x.rbp_evidences), reverse=True)
        
class Gene(db.Model):
    """A gene"""
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.String(63), index=True)
    genename = db.Column(db.String(255), index=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'))


class Protein(db.Model):
    """A single protein"""
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.String(63), index=True)
    uniprot_acc = db.Column(db.String(63), index=True)
    #genename = db.Column(db.String(255))
    proteinname = db.Column(db.String(255))
    uniprot_url = db.Column(db.String(255))
    #kdrew: uses table name for ProteinComplexMapping class (annoying sqlalchemy magic)
    complexes = db.relationship('Complex', secondary='protein_complex_mapping',  back_populates='proteins')
    genenames = db.relationship('Gene')

    rbp_evidences = db.relationship('RBPEvidence')
    rbp_stats = db.relationship('RBPStats')

    def genename(self,):
        gnames = [g for g in self.genenames]
        if len(gnames) > 0:
            return gnames[0].genename
        else:
            return self.gene_id

    def uniprot_link(self,):
        if self.uniprot_url != "":
            retstr = "<a href=https://www.uniprot.org/uniprot/%s target=\"_blank\">%s</a>" % (self.uniprot_acc, 'UniProt')
        else:
            retstr = ""
        return retstr

    def ncbi_link(self,):
        retstr = "<a href=https://www.ncbi.nlm.nih.gov/gene/%s target=\"_blank\">%s</a>" % (self.gene_id, 'NCBI')
        return retstr

class Edge(db.Model):
    """A protein protein edge"""
    id = db.Column(db.Integer, primary_key=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'), index=True )
    protein_key2 = db.Column(db.Integer, db.ForeignKey('protein.id'), index=True )
    score = db.Column(db.Float)

    evidences = db.relationship('EdgeEvidence')

    def get_proteins(self,):
        #kdrew: complex.html calls this twice per edge (once per protein
        #column), and the old Protein.id.in_([...]) query always hit the
        #DB fresh regardless -- .get() is identity-map-aware, so on pages
        #where the complex's proteins were already loaded earlier in the
        #same request (comp.sorted_proteins(), rendered above this table),
        #these become free, in-memory lookups instead of new queries.
        prot1 = db.session.query(Protein).get(self.protein_key)
        prot2 = db.session.query(Protein).get(self.protein_key2)
        return [prot1, prot2]

class EdgeEvidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'))
    evidence_type = db.Column(db.String(255))
    
#kdrew: this table holds evidence of the protein being an rna binder
class RBPEvidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'))
    evidence = db.Column(db.String(255))
    evidence_type = db.Column(db.String(255))

    #low_throughput                                    
    #high_throughput                                   
    #computational                                     
    #rnabinding                                        
    #rnabinding_uniprot                                
    #trendel                                           
    #queiroz                                           
    #huang                                             
    #bao                                               
    #HEK293-RIC_Hs_Baltz2012                           
    #HuH7-IC_Hs_Beckmann2015                           
    #HeLa-RNPxl_Hs_Kramer2014                          
    #HeLa-IC_Hs_Castello2012                           
    #HeLa-RBDmap_Hs_Castello2016                       
    #K562-serIC-chr_Hs_Conrad2016_chr                  
    #K562-serIC_Hs_Conrad2016                          
    #annotation                                        

    #annotated                                         

class RBPStats(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'))
    #,diffrac,diffrac_percent,diffrac_normalized,emd,pearsonr,mean_abundance,annotated,zscore,sliding_zscore,pvalues,pvalues_fdrcor,sliding_pvalues,sliding_pvalues_fdrcor,low_throughput,high_throughput,computational,rnabinding,rnabinding_uniprot,trendel,queiroz,huang,bao,HEK293-RIC_Hs_Baltz2012,HuH7-IC_Hs_Beckmann2015,HeLa-RNPxl_Hs_Kramer2014,HeLa-IC_Hs_Castello2012,HeLa-RBDmap_Hs_Castello2016,K562-serIC-chr_Hs_Conrad2016_chr,K562-serIC_Hs_Conrad2016,Enzyme,Metabolism,"""Metabolic.Enzyme""",annotation,Entry name,Status,Protein names,Gene names,Organism,Cross-reference (GeneID),Gene names  (primary ),Keywords,Gene ontology (biological process),Gene ontology (cellular component),Gene ontology (molecular function),ht_go_uniprot_sec2_rbp

    diffrac = db.Column(db.Float) 
    diffrac_percent = db.Column(db.Float) 
    diffrac_normalized = db.Column(db.Float) 
    emd = db.Column(db.Float) 
    pearsonr = db.Column(db.Float) 
    mean_abundance = db.Column(db.Float) 
    zscore = db.Column(db.Float) 
    sliding_zscore = db.Column(db.Float) 
    pvalues = db.Column(db.Float) 
    pvalues_fdrcor = db.Column(db.Float) 
    sliding_pvalues = db.Column(db.Float) 
    sliding_pvalues_fdrcor = db.Column(db.Float) 

    #kdrew: tests if protein is not significant but below the 50% FDR
    def is_diffrac_50(self,):
        if self.sliding_pvalues_fdrcor > 0.05 and self.sliding_pvalues_fdrcor < 0.5:
            return True


class ProteinComplexMapping(db.Model):
    """A mapping between proteins and complexes"""
    __tablename__ = 'protein_complex_mapping'
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'), primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'), primary_key=True)

class ComplexRNPStats(db.Model):
    """RNP stats for a Complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'))
    #,sec2_cntl_median_correlation,sec2_rnase_median_correlation,sec2_cntl_num_ids,sec2_cntl_fraction_ids,sec2_rnase_num_ids,sec2_rnase_fraction_ids,sec2_num_fdr05,sec2_fraction_fdr05,sec2_num_fdr2,sec2_fraction_fdr2,sec2_num_annotated,sec2_fraction_annotated,sec2_num_rnabinding,sec2_fraction_rnabinding,sec2_num_rnabinding_uniprot,sec2_fraction_rnabinding_uniprot,sec2_num_lt,sec2_fraction_lt,sec2_num_ht,sec2_fraction_ht,sec2_num_unannotated,sec2_fraction_unannotated,sec2_num_ht_go_uniprot_sec2_rbp,sec2_fraction_ht_go_uniprot_sec2_rbp,Significant_RNP,Significant_20FDR_08Corr_RNP,High_Throughput_RNP,GO_Uniprot_RNP,RNP_label,complex_id,sec2_logratio_psm,sec2_bw_exps_median_correlation,sec2_num_fdr5,sec2_fraction_fdr5
    

    Significant_RNP = db.Column(db.Boolean)
    Significant_20FDR_08Corr_RNP = db.Column(db.Boolean)
    High_Throughput_RNP = db.Column(db.Boolean)
    GO_Uniprot_RNP = db.Column(db.Boolean)
    #RNP_label = db.Column(db.String(255))

    sec2_cntl_median_correlation = db.Column(db.Float) 
    sec2_rnase_median_correlation = db.Column(db.Float) 
    sec2_bw_exps_median_correlation = db.Column(db.Float) 

    sec2_cntl_num_ids = db.Column(db.Integer)
    sec2_rnase_num_ids = db.Column(db.Integer)
    sec2_num_fdr05 = db.Column(db.Integer)
    sec2_num_fdr2 = db.Column(db.Integer)
    sec2_num_annotated = db.Column(db.Integer)
    sec2_num_rnabinding = db.Column(db.Integer)
    sec2_num_rnabinding_uniprot = db.Column(db.Integer)
    sec2_num_lt = db.Column(db.Integer)
    sec2_num_ht = db.Column(db.Integer)
    sec2_num_unannotated = db.Column(db.Integer)
    sec2_num_ht_go_uniprot_sec2_rbp = db.Column(db.Integer)
    sec2_num_fdr5 = db.Column(db.Integer) 

    sec2_cntl_fraction_ids= db.Column(db.Float) 
    sec2_rnase_fraction_ids= db.Column(db.Float) 
    sec2_fraction_fdr05 = db.Column(db.Float) 
    sec2_fraction_fdr2 = db.Column(db.Float) 
    sec2_fraction_fdr5 = db.Column(db.Float) 
    sec2_fraction_annotated = db.Column(db.Float) 
    sec2_fraction_rnabinding = db.Column(db.Float) 
    sec2_fraction_rnabinding_uniprot= db.Column(db.Float) 
    sec2_fraction_lt = db.Column(db.Float) 
    sec2_fraction_ht = db.Column(db.Float) 
    sec2_fraction_unannotated = db.Column(db.Float) 
    sec2_fraction_ht_go_uniprot_sec2_rbp= db.Column(db.Float) 

    sec2_logratio_psm = db.Column(db.Float) 

class ComplexEnrichment(db.Model):
    """Annotation Enrichment for a Complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'))
    #  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list
    corr_pval = db.Column(db.Float)
    t_count = db.Column(db.Integer)
    q_count = db.Column(db.Integer)
    qandt_count = db.Column(db.Integer)
    qandt_by_q = db.Column(db.Float)
    qandt_by_t = db.Column(db.Float)
    term_id = db.Column(db.String(255))
    t_type = db.Column(db.String(63))
    t_group = db.Column(db.Integer)
    t_name = db.Column(db.String(255))
    depth_in_group = db.Column(db.Integer)
    qandt_list = db.Column(db.String(255))

    def get_proteins(self,):
        proteins = db.session.query(Protein).filter(Protein.uniprot_acc.in_(self.qandt_list.split(','))).all()
        return proteins





