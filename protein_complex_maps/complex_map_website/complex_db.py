
from flask import Flask
#from flask.ext.sqlalchemy import SQLAlchemy
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import func, or_, and_

import datetime as dt

import itertools as it

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db/humap3v1.db'

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
    humap_id = db.Column(db.String(63), unique=True, index=True)
    #kdrew: uses table name for ProteinComplexMapping class (annoying sqlalchemy magic)
    proteins = db.relationship('Protein', secondary='protein_complex_mapping', back_populates='complexes')
    #kdrew: uses table name for EdgeComplexMapping class (annoying sqlalchemy magic)
    edges = db.relationship('Edge', secondary='edge_complex_mapping', back_populates='complexes')
    enrichments = db.relationship('ComplexEnrichment')
    top_rank = db.Column(db.Integer)
    average_precision = db.Column(db.Float)
    average_precision_rank = db.Column(db.Integer)
    average_window_precision = db.Column(db.Float)
    average_window_precision_rank = db.Column(db.Integer)

    def complex_link(self,):
        #retstr = "<a href=displayComplexes?complex_key=%s>%s</a>" % (self.complex_id, self.complex_id)
        retstr = "<a href=displayComplexes?complex_key=%s>%s</a>" % (self.humap_id, self.humap_id)
        return retstr

    def sorted_proteins(self,):
        return sorted(self.proteins, key=lambda p: len([g for g in p.genenames]), reverse=True)
            

    #kdrew: a bit of a bottleneck when serving pages, has to generate all combinations of proteins in complex and then search for combination in edge table
    #kdrew: seems like there should be a better way of doing this, either 1) combining protein keys as a single index or 2) storing mapping between complex and edges directly
    def all_edges(self,):
        #es = []
        #for prot1, prot2 in it.combinations(self.proteins,2):
        #    #kdrew: edge table enforces order
        #    if prot2.id < prot1.id:
        #        prot2, prot1 = prot1, prot2
        #    edge = db.session.query(Edge).filter( and_(Edge.protein_key == prot1.id, Edge.protein_key2 == prot2.id) ).first()
        #    if edge != None:
        #        es.append(edge)

        #edges = [db.session.query(Edge).filter((and_(Edge.protein_key == prot1.id, Edge.protein_key2 == prot2.id) | and_(Edge.protein_key == prot2.id,Edge.protein_key2 == prot1.id))).first() for prot1, prot2 in it.combinations(self.proteins,2)]
        #es = [e for e in edges if e != None]

        return sorted(list(set(self.edges)), key=lambda es: es.score, reverse=True)

        
class Gene(db.Model):
    """A gene"""
    id = db.Column(db.Integer, primary_key=True)
    #gene_id = db.Column(db.String(63), index=True)
    genename = db.Column(db.String(255), index=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'))


class ProteinDiseaseAnnotation(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'))
    disease_database = db.Column(db.String(255))
    disease_annotation = db.Column(db.Text)


class Protein(db.Model):
    """A single protein"""
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.String(63), index=True)
    uniprot_acc = db.Column(db.String(63), index=True)
    #genename = db.Column(db.String(255))
    proteinname = db.Column(db.String(255))
    uniprot_url = db.Column(db.String(255))
    annotation_score = db.Column(db.Integer)
    #kdrew: uses table name for ProteinComplexMapping class (annoying sqlalchemy magic)
    complexes = db.relationship('Complex', secondary='protein_complex_mapping',  back_populates='proteins')
    genenames = db.relationship('Gene')
    disease_annotations = db.relationship('ProteinDiseaseAnnotation')

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
    #kdrew: uses table name for EdgeComplexMapping class (annoying sqlalchemy magic)
    complexes = db.relationship('Complex', secondary='edge_complex_mapping',  back_populates='edges')
    score = db.Column(db.Float)
    #kdrew: precision should be calculated from sklearn precision_recall_curve
    precision = db.Column(db.Float)
    #kdrew: window_precision should be calculated for a window surrounding the given edge's svm confidence score
    window_precision = db.Column(db.Float)

    evidences = db.relationship('Evidence')
    prothd = db.relationship('ProtHD')

    def get_proteins(self,):
        #prot1 = db.session.query(Protein).filter(Protein.id==self.protein_key).first()
        #prot2 = db.session.query(Protein).filter(Protein.id==self.protein_key2).first()
        prots = db.session.query(Protein).filter(Protein.id.in_([self.protein_key,self.protein_key2])).all()
        return prots

    def get_prothd_score(self,):
        return self.prothd[0].prothd_score

class Evidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'))
    evidence_type = db.Column(db.String(255))
    
class ProtHD(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'))
    prothd_score = db.Column(db.Float)

class ProteinComplexMapping(db.Model):
    """A mapping between proteins and complexes"""
    __tablename__ = 'protein_complex_mapping'
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'), primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'), primary_key=True)

class EdgeComplexMapping(db.Model):
    """A mapping between edges and complexes"""
    __tablename__ = 'edge_complex_mapping'
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'), primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'), primary_key=True)


class ComplexEnrichment(db.Model):
    """Annotation Enrichment for a Complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'))
    #  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list
    #index,complex_id,description,effective_domain_size,goshv,group_id,intersection_genes,intersection_size,name,native,p_value,parents,precision,query,query_size,recall,significant,source,source_order,term_size

    description = db.Column(db.BLOB)
    effective_domain_size = db.Column(db.Integer)
    goshv = db.Column(db.Integer)
    group_id = db.Column(db.Integer)
    intersection_genes = db.Column(db.BLOB)
    intersection_size = db.Column(db.Integer) 
    name = db.Column(db.String(255)) 
    native = db.Column(db.String(255))  
    p_value = db.Column(db.Float)
    parents = db.Column(db.String(255)) 
    precision = db.Column(db.Float)
    query = db.Column(db.String(255))  
    query_size = db.Column(db.Integer) 
    recall = db.Column(db.Float) 
    significant = db.Column(db.String(255))  
    source = db.Column(db.String(255))  
    source_order = db.Column(db.Integer) 
    term_size = db.Column(db.Integer) 

    #corr_pval = db.Column(db.Float)
    #t_count = db.Column(db.Integer)
    #q_count = db.Column(db.Integer)
    #qandt_count = db.Column(db.Integer)
    #qandt_by_q = db.Column(db.Float)
    #qandt_by_t = db.Column(db.Float)
    #term_id = db.Column(db.String(255))
    #t_type = db.Column(db.String(63))
    #t_group = db.Column(db.Integer)
    #t_name = db.Column(db.String(255))
    #depth_in_group = db.Column(db.Integer)
    #qandt_list = db.Column(db.String(255))

    def get_proteins(self,):
        proteins = db.session.query(Protein).filter(Protein.uniprot_acc.in_(self.intersection_genes.split())).all()
        c_proteins = db.session.query(Complex).filter(Complex.id == self.complex_key).first().proteins
        return set(proteins).intersection(c_proteins)

class SessionComplexCSV(db.Model):
    """Temporary storage for session data"""
    id = db.Column(db.Integer, primary_key=True)
    remote_addr = db.Column(db.String(255))
    query_id = db.Column(db.String(255))
    #store_date = db.Column(db.TimeStamp)
    store_date = db.Column(db.DateTime(timezone=True), server_default=func.now())
    complex_csv = db.Column(db.Text)


    def delete_expired_entries(self):
        yesterdays_datetime = dt.datetime.now() - dt.timedelta(days=1)
        #print dt.datetime() - 1
        db.session.query(SessionComplexCSV).filter(SessionComplexCSV.store_date < yesterdays_datetime ).delete(synchronize_session=False)


