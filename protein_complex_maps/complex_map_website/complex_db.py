
from flask import Flask
#from flask.ext.sqlalchemy import SQLAlchemy
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import func, or_, and_
from sqlalchemy.orm import subqueryload

import itertools as it

app = Flask(__name__)
#app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db/test.db'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db/humap1.db'

app.config['SECRET_KEY'] = 'please, tell nobody'

db = SQLAlchemy(app)

def get_db():
    return db

def get_app():
    return app

def get_or_create(db, model, **kwargs):
    session = db.session
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


    def complex_link(self,):
        retstr = "<a href=displayComplexes?complex_key=%s>%s</a>" % (self.complex_id, self.complex_id)
        return retstr

    #kdrew: this used to generate every pairwise combination of the
    #complex's proteins and query the edge table once per pair (O(n^2)
    #queries just to find which pairs happen to have an edge) -- for a
    #150-protein complex that's ~11,000 individual queries before even
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
    annotation_score = db.Column(db.Integer)

    def genename(self,):
        gnames = [g for g in self.genenames]
        if len(gnames) > 0:
            return gnames[0].genename
        else:
            return self.gene_id

    def uniprot_link(self,):
        if self.uniprot_url != "":
            retstr = "<a href=%s target=\"_blank\">%s</a>" % (self.uniprot_url, 'UniProt')
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

    evidences = db.relationship('Evidence')

    def get_proteins(self,):
        #kdrew: complex.html calls this twice per edge (once per protein
        #column), and the old Protein.id.in_([...]) query always hit the
        #DB fresh regardless -- .get() is identity-map-aware, so on pages
        #where the complex's proteins were already loaded earlier in the
        #same request (comp.proteins, rendered above this table), these
        #become free, in-memory lookups instead of new queries.
        prot1 = db.session.query(Protein).get(self.protein_key)
        prot2 = db.session.query(Protein).get(self.protein_key2)
        return [prot1, prot2]

class Evidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'))
    evidence_type = db.Column(db.String(255))
    

class ProteinComplexMapping(db.Model):
    """A mapping between proteins and complexes"""
    __tablename__ = 'protein_complex_mapping'
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'), primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'), primary_key=True)

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





